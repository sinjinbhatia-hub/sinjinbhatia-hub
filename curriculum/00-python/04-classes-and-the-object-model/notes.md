# PY-4: Classes and the Object Model — Dunder Methods, Dataclasses, Protocols, and a Dual-Number Differentiator

**Time:** ~10 hours · **Prerequisites:** PY-3 (LA-1 helps for the `Vec` example) · **Unlocks:** PY-5, PY-6, and, through the dual numbers of §3.7, LA-14 and ML-2

> **Conventions in this module.** Python 3.11+; every output below came from a real run under Python 3.13. `>>>` blocks are interactive transcripts; other blocks are scripts followed by their output. A block that starts `class Vec:` and shows only some methods means *add these inside the class you already have*; `project.py` holds the complete classes. A **dunder** is a method named like `__add__` that Python calls for you. `__repr__` returns text you could type back in; `__eq__` and `__hash__` are designed together. A **dual number** is $a + b\varepsilon$ with $\varepsilon^2 = 0$; $\varepsilon$ is never machine epsilon (`eps`, $2.2\times10^{-16}$). The classes are pure Python; NumPy and SymPy only check them. Ids change between runs.

## 0. Where this sits

PY-3 gave you functions as values. A class is the other way to bundle state with behavior, and in Python it is also how your own types obey the language's operators. You need PY-1's names-and-objects model and PY-2's data-model table (§3.7). The payoffs: PY-5's exceptions are classes and its tests check the contracts of §3.2; PY-6's `ndarray` implements the dunders of §3.3–§3.4; and the dual numbers of §3.7 are forward-mode automatic differentiation, which LA-14 reads as a Jacobian and ML-2 and ML-7 run backwards as backpropagation.

## 1. The problem that forced this idea

You want to add two 3-vectors, scale one, and take a length. The obvious representation is a list.

```python
v = [1.0, 2.0, 3.0]
w = [4.0, 5.0, 6.0]
print(v + w)
print(v * 2)
```
```text
[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
[1.0, 2.0, 3.0, 1.0, 2.0, 3.0]
```

The list answers `+` and `*` with *its* meanings (concatenate, repeat) and offers no length. Writing `add(v, w)` and `scale(2, v)` turns `2*v + w` into `add(scale(2, v), w)`. Operators belong to types, not to the language: `v + w` asks `v`'s type what `+` means. So define a type whose answers are yours. The forcing problem is **operator overloading for a type the language does not ship**.

The second problem looks unrelated. You need $f'(x)$ for $f(x) = x\sin x/(1 + x^2)$ at $x = 1.3$ without differentiating by hand (a Newton solver or optimizer needs it for thousands of different $f$). The standard shortcut is a difference quotient. Its error against the exact derivative (SymPy, 30 digits):

```python
import math
import sympy as sp

x = sp.symbols("x")
exact = float(sp.diff(x * sp.sin(x) / (1 + x**2), x).subs(x, 1.3).evalf(30))
f = lambda t: t * math.sin(t) / (1 + t**2)
print(f"exact f'(1.3) = {exact:.16f}")
print(f"{'h':>7} {'forward':>10} {'central':>10}     |error| of the difference quotient")
for k in (2, 5, 8, 11):
    h = 10.0 ** -k
    fwd = (f(1.3 + h) - f(1.3)) / h
    cen = (f(1.3 + h) - f(1.3 - h)) / (2 * h)
    print(f"{h:7.0e} {abs(fwd - exact):10.2e} {abs(cen - exact):10.2e}")
```
```text
exact f'(1.3) = 0.0373942117663269
      h    forward    central     |error| of the difference quotient
  1e-02   3.41e-03   1.17e-05
  1e-05   3.43e-06   5.19e-12
  1e-08   9.29e-09   6.51e-09
  1e-11   3.65e-06   1.90e-06
```

The error falls, bottoms out, and *rises*: a V. Truncation error shrinks with $h$ (like $h$ forward, $h^2$ central), but the rounding error of $f(x+h) - f(x)$, about `eps` times $|f|$, is divided by $h$ and grows as $h$ shrinks (PY-1 §3.5). The best you can reliably do is about $10^{-8}$ (forward) or $10^{-11}$ (central), at a step size you must guess (exercise B7 derives it). §3.7 builds a number that returns the derivative with error near $10^{-16}$ and no step size. Both problems have one fix: a new kind of object that behaves under the operators Python already has.

## 2. Intuition

**A class is a blueprint; an instance is a dict plus a pointer.**

```
v = Vec(1, 2, 3)

 v ──► ┌ instance ──────────────────────┐      ┌ class object Vec ───────────────┐
       │ __class__ ───────────────────────────►│ __init__ __repr__ __add__ __eq__ │
       │ __dict__ = {'_c': (1.0,2.0,3.0)}│      │ __mul__  __rmul__ dot  cross ... │
       └────────────────────────────────┘      └─────────────────────────────────┘
```

Data lives in the instance; behavior lives in the class. `v.dot(w)` looks in `v.__dict__`, finds no `dot`, then finds the function in `Vec`.

**Operators are method calls with special names.** `a + b` is `type(a).__add__(a, b)`, so any object can be taught to add, compare, index, iterate or be called:

| You write | Python evaluates | Section |
|---|---|---|
| `a + b`, `-a`, `abs(a)` | `type(a).__add__(a, b)`, `__neg__`, `__abs__`; if `a` cannot, `type(b).__radd__(b, a)` | §3.3 |
| `a == b`, `a < b`, `hash(a)`, `repr(a)` | `__eq__`, `__lt__`, `__hash__`, `__repr__` | §3.2 |
| `len(a)`, `a[i]`, `x in a`, `for x in a` | `__len__`, `__getitem__`, `__contains__`, `__iter__` | §3.4 |
| `a(x)`, `a @ b`, `a += b` | `__call__`, `__matmul__`, `__iadd__` (else `a = a + b`) | §3.3–3.4 |

**Duck typing.** Python never asks "is this a `Duck`?" when you define a function; it asks "can it quack?" when you call it. `x * x + 1` works for ints, floats, `Fraction`s and arrays because each answers `*` and `+`. Once `Dual` answers `+ - * /`, every function written with those operators accepts it, unedited.

**Guess before §3.7.** Taylor: $f(x + h) = f(x) + f'(x)h + \tfrac12 f''(x)h^2 + \cdots$. Suppose a "number" $\varepsilon$ has $\varepsilon^2 = 0$. The series stops after two terms: $f(x + \varepsilon) = f(x) + f'(x)\,\varepsilon$, value and derivative together. Take $u + u'\varepsilon$ and $v + v'\varepsilon$, with $u', v'$ the derivatives of two functions, and guess their product as $\text{value} + \text{coefficient}\cdot\varepsilon$. The rule for $(uv)'$ should be in the coefficient.

## 3. Formal development

### 3.1 Defining a class; what `self` is; where attributes live

Four rules fix the semantics.

1. A `class` statement runs its body once in a fresh namespace. The names it binds (mostly functions) become attributes of a **class object**, bound to the class name.
2. `Vec(1, 2)` makes an empty instance (`__new__`, rarely overridden) and calls `Vec.__init__(instance, 1, 2)`, which fills it in and returns `None`.
3. **Reading** `obj.name` searches: (a) a *data descriptor* on the class, such as a property (§3.5); (b) the instance's `__dict__`; (c) the class, then its parents in method-resolution order, where a function found becomes a *bound method*; (d) `__getattr__`, if defined. Otherwise `AttributeError`.
4. **Writing** `obj.name = value` stores into `obj.__dict__` (unless a data descriptor intercepts). It never writes to the class.

▶ **Try it first.** Predict every printed value, then explain `count` versus `log`.

```python
class Tally:
    count = 0
    log = []

    def __init__(self, name):
        self.name = name

    def hit(self):
        self.count += 1
        self.log.append(self.name)


a, b = Tally("a"), Tally("b")
a.hit(); a.hit(); b.hit()
print(a.count, b.count, Tally.count)
print(a.log, b.log, Tally.log, a.log is b.log)
print(a.__dict__, b.__dict__)
```

■ **What just happened.**

```text
2 1 0
['a', 'a', 'b'] ['a', 'a', 'b'] ['a', 'a', 'b'] True
{'name': 'a', 'count': 2} {'name': 'b', 'count': 1}
```

`self.count += 1` is `self.count = self.count + 1`: the read finds the class attribute `0` (rule 3c), the write creates an *instance* attribute (rule 4), so `Tally.count` never moves. `self.log.append(...)` is a read followed by a mutation of whatever was found, and the only `log` is the class's single list. Class attributes are for constants and shared immutable defaults; anything an instance mutates is created in `__init__`.

▶ **Build it yourself (0): objects from dicts.** An instance is a dict holding its class and attributes; a class is a dict of functions; a "method call" looks the function up in the class and passes the instance first.

```python
def new(cls, *args):
    obj = {"cls": cls, "attrs": {}}
    cls["__init__"](obj, *args)
    return obj


def send(obj, name, *args):
    return obj["cls"][name](obj, *args)


def point_init(self, x, y):
    self["attrs"].update(x=x, y=y)


def point_norm(self, k):
    return k * (self["attrs"]["x"] ** 2 + self["attrs"]["y"] ** 2) ** 0.5


Point = {"__init__": point_init, "norm": point_norm}
p = new(Point, 3, 4)
print(send(p, "norm", 2), p["cls"] is Point)
```
```text
10.0 True
```

`self` is whatever `send` passes first. Python does the same: `p.norm(2)` finds `norm` through `type(p)` and calls it with `p` in front, so **`obj.method(args)` is `Class.method(obj, args)`**. A *bound method* object does the insertion (Example 4.1 takes one apart); every `self` confusion is this one fact.

### 3.2 Showing, comparing, hashing: `__repr__`, `__eq__`, `__hash__`

▶ **Build it yourself (1a): `Vec`.** Store components in a tuple, which is immutable, so the vector cannot change under a set or dict that holds it.

```python
class Vec:
    def __init__(self, *components):
        self._c = tuple(float(c) for c in components)

    def __repr__(self):
        return f"Vec({', '.join(map(repr, self._c))})"


v = Vec(1, 2, 3)
print(v)
print(v is Vec(1, 2, 3), v == Vec(1, 2, 3))
```
```text
Vec(1.0, 2.0, 3.0)
False False
```

The last `False`: without `__eq__`, `==` is identity. First, what are `repr` and `str` for?

▶ **Try it first.** Predict the four lines.

```python
class A:
    def __repr__(self):
        return "A()"


class B(A):
    def __str__(self):
        return "a B"


a, b = A(), B()
print(a, b)
print([a, b])
print(f"{a} {b} {b!r}")
print(str([b]), repr(b))
```

■ **What just happened.**

```text
A() a B
[A(), A()]
A() a B A()
[A()] A()
```

`print` and f-strings call `str(x)`, which uses `__str__` and falls back to `__repr__`, so `a` prints as `A()`. Containers and the REPL use `repr`; `B` inherited `A.__repr__`. `repr` is for programmers: unambiguous and, where possible, text you could type back in (`eval(repr(x)) == x`). `str` is for end users. Write `__repr__` first; add `__str__` only if the friendly text differs.

**Equality and the hash contract.**

1. `a == b` calls `type(a).__eq__(a, b)`. If that returns the singleton `NotImplemented` (a returned value, not the exception `NotImplementedError`), Python tries `type(b).__eq__(b, a)`; if that also declines, `==` falls back to identity. Return `NotImplemented` for operand types you do not understand, never `False`. `!=` is derived from `==`.
2. **The contract:** `a == b` implies `hash(a) == hash(b)`. Sets and dicts find an item by its hash first (PY-2 §3.4); equal keys with different hashes sit in different slots and never meet.
3. Defining `__eq__` in a class body *without* `__hash__` sets `__hash__` to `None`: instances become unhashable. A by-value `__eq__` with the inherited by-address hash would break statement 2, so Python refuses.
4. Define `__hash__` only for objects immutable in the fields `__eq__` compares, and hash exactly those fields. Mutable objects (list, dict, `Matrix` in §3.8) stay unhashable.

```python
class Vec:
    def __str__(self):
        return "(" + ", ".join(f"{c:g}" for c in self._c) + ")"

    def __eq__(self, other):
        if not isinstance(other, Vec):
            return NotImplemented
        return self._c == other._c

    def __hash__(self):
        return hash(self._c)


v = Vec(1, 2, 3)
print(v, repr(v))
print(v == Vec(1.0, 2.0, 3.0), v == [1, 2, 3], v != Vec(1, 2, 4))
print(eval(repr(v)) == v, len({v, Vec(1, 2, 3), Vec(3, 2, 1)}))
```
```text
(1, 2, 3) Vec(1.0, 2.0, 3.0)
True False True
True 2
```

`v == [1, 2, 3]` is `False` because both sides declined and identity was the fallback. Example 4.2 shows the error when `__hash__` is forgotten.

**Ordering.** `a < b` calls `type(a).__lt__(a, b)`, then the mirror `type(b).__gt__(b, a)` if that declines. Define `__eq__` and one ordering method and `functools.total_ordering` supplies the rest. Define an order only when one is natural: vectors have none (sorting `Vec`s raises `TypeError`, correctly), rationals do (§3.6).

### 3.3 Arithmetic: operators are method calls

The rules for a binary operator, written for `+` (data model reference §3.3.8 has the full list):

1. Python calls `type(a).__add__(a, b)`. If that is missing or returns `NotImplemented`, and the types differ, it calls the **reflected** method `type(b).__radd__(b, a)`. If that also fails: `TypeError: unsupported operand type(s)`.
2. Exception: if `type(b)` is a *subclass* of `type(a)` and overrides `__radd__`, `b.__radd__(a)` is tried **first**.
3. `a += b` calls `type(a).__iadd__(a, b)` if it exists (it should mutate `a` and `return self`); if it is missing or declines, Python evaluates `a = a + b` by rules 1–2. Immutable types never define `__iadd__`.
4. The family: `__sub__`, `__mul__`, `__truediv__`, `__pow__`, `__matmul__` (for `@`) and so on, each with `r` and `i` twins, plus `__neg__`, `__abs__`. Lookup is on the *type*, so assigning `v.__add__ = f` does not change `v + w`.

▶ **Build it yourself (1b): arithmetic for `Vec`.** Add `+`, `-`, unary `-` and scalar `*`. Test `v * 2`, and before running it predict whether `2 * v` works.

```python
import math
import numbers


class Vec:
    def __add__(self, other):
        if not isinstance(other, Vec):
            return NotImplemented
        return Vec(*(a + b for a, b in zip(self._c, other._c, strict=True)))

    def __sub__(self, other):
        return self + (-other) if isinstance(other, Vec) else NotImplemented

    def __neg__(self):
        return Vec(*(-a for a in self._c))

    def __mul__(self, k):
        if not isinstance(k, numbers.Real):
            return NotImplemented
        return Vec(*(k * a for a in self._c))


v, w = Vec(1, 2, 3), Vec(4, 5, 6)
print(v + w, v - w, -v, v * 2)
```
```text
(5, 7, 9) (-3, -3, -3) (-1, -2, -3) (2, 4, 6)
```

```pycon
>>> 2 * v
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: unsupported operand type(s) for *: 'int' and 'Vec'
```

`numbers.Real` accepts `int`, `float`, `Fraction` and NumPy scalars: duck typing against the numeric tower, not a list of concrete types. `zip(..., strict=True)` makes a dimension mismatch a `ValueError`.

The failure is rule 1: `2 * v` asks `int.__mul__(2, v)`, which declines; Python looks for `Vec.__rmul__(v, 2)`, finds none, and gives up. Scalar multiplication commutes, so the reflected method is the same function. Step 4 adds it, with division, length, dot and cross products.

```python
class Vec:
    __rmul__ = __mul__

    def __truediv__(self, k):
        if not isinstance(k, numbers.Real):
            return NotImplemented
        return Vec(*(a / k for a in self._c))

    def __abs__(self):
        return math.hypot(*self._c)

    def dot(self, other):
        return sum(a * b for a, b in zip(self._c, other._c, strict=True))

    def __matmul__(self, other):
        return self.dot(other) if isinstance(other, Vec) else NotImplemented

    def cross(self, other):
        if len(self._c) != 3 or len(other._c) != 3:
            raise ValueError("cross product needs two 3-vectors")
        (a1, a2, a3), (b1, b2, b3) = self._c, other._c
        return Vec(a2 * b3 - a3 * b2, a3 * b1 - a1 * b3, a1 * b2 - a2 * b1)


u, v = Vec(1, 2, 3), Vec(4, 5, 6)
print(2 * u + v / 2, abs(Vec(3, 4)), u @ v, u.cross(v))
print(u.cross(v) @ u, u.cross(v) @ v)
```
```text
(4, 6.5, 9) 5.0 32.0 (-3, 6, -3)
0.0 0.0
```

■ **What just happened.**

These are LA-1's definitions as code: componentwise sum and scalar multiple (LA-1 §3.1), $\mathbf u\cdot\mathbf v = \sum_i u_iv_i$ and $\|\mathbf v\| = \sqrt{\mathbf v\cdot\mathbf v}$ (LA-1 §3.5); the cross product is orthogonal to both factors, as the last line shows. `project.py` Part 1 checks Cauchy–Schwarz (LA-1 Theorem 3.11) and the parallelogram law on thousands of random pairs. Every method returns a new `Vec` or `NotImplemented`; none mutates.

### 3.4 Containers, iteration and `__call__`

Brackets, `len`, `in`, `for` and `bool` are protocols (PY-2 §3.7). What arrives in `__getitem__` is **whatever was written between the brackets, as one object.**

```python
class Probe:
    def __getitem__(self, key):
        return key


p = Probe()
print(p[3], p[1:5:2], p[::-1], p[1, 2], p[:, 0], p[..., None])
```
```text
3 slice(1, 5, 2) slice(None, None, -1) (1, 2) (slice(None, None, None), 0) (Ellipsis, None)
```

A slice is the object `slice(start, stop, step)`; `a[i, j]` passes the tuple `(i, j)`. NumPy's `a[:, 0]` is the fifth case: PY-6 indexing is this protocol with more axes. A sequence implements `__len__` and `__getitem__` (an `int` gives an item, a `slice` the same kind of sequence); `__setitem__` and `__delitem__` make it mutable.

```python
class Squares:
    def __len__(self):
        return 4

    def __getitem__(self, i):
        if not 0 <= i < 4:
            raise IndexError(i)
        return i * i


s = Squares()
print(s[3], list(s), 9 in s, list(reversed(s)), bool(s))
```
```text
9 [0, 1, 4, 9] True [9, 4, 1, 0] True
```

We wrote no `__iter__`, `__contains__`, `__reversed__` or `__bool__`, yet all work. With no `__iter__`, Python calls `__getitem__(0)`, `__getitem__(1)`, ... until `IndexError`; `in` falls back to iteration and `==`; `reversed` uses `__len__` and `__getitem__`; `bool` falls back to `len(x) != 0`. An explicit `__iter__` is clearer and faster: it returns an iterator (an object with `__next__`; `iter(...)` or a generator function supplies one, PY-2 §3.7–3.8). `Vec` needs one so that `x, y, z = v` works, and a `__bool__`: with neither `__len__` nor `__bool__`, every instance is truthy, zero vector included.

▶ **Build it yourself (1c): iteration for `Vec`.**

```python
class Vec:
    def __iter__(self):
        return iter(self._c)

    def __len__(self):
        return len(self._c)

    def __bool__(self):
        return any(self._c)


x, y, z = Vec(1, 2, 3)
print(x, y, z, len(Vec(1, 2)), bool(Vec(0, 0)), bool(Vec(0, 1)), max(Vec(3, 9, 4)))
```
```text
1.0 2.0 3.0 2 False True 9.0
```

**Callable objects.** `obj(args)` calls `type(obj).__call__(obj, args)`. A class with `__call__` is a function that remembers state, the object twin of PY-3's closure, and `callable(obj)` is `True`. `Polynomial` in §3.9 uses it so that a polynomial can be passed wherever a function is expected.

### 3.5 Properties, `@classmethod`, `@staticmethod`

**Properties** make attribute syntax run code: validation on assignment, or a value computed on demand. Callers write `r.ohms = 5`, never `r.set_ohms(5)`, so you can begin with a plain attribute and add checks later without touching a caller.

```python
class Resistor:
    def __init__(self, ohms):
        self.ohms = ohms                      # goes through the setter below

    @property
    def ohms(self):
        return self._ohms

    @ohms.setter
    def ohms(self, value):
        if value <= 0:
            raise ValueError(f"resistance must be positive, got {value}")
        self._ohms = float(value)

    @property
    def conductance(self):                    # computed, no setter: read-only
        return 1 / self._ohms

    def __repr__(self):
        return f"Resistor({self._ohms!r})"


r = Resistor(50)
print(r.ohms, r.conductance)
for action in (lambda: setattr(r, "ohms", -5), lambda: setattr(r, "conductance", 3)):
    try:
        action()
    except (ValueError, AttributeError) as err:
        print(type(err).__name__ + ":", err)
```
```text
50.0 0.02
ValueError: resistance must be positive, got -5
AttributeError: property 'conductance' of 'Resistor' object has no setter
```

The underscore in `_ohms` is a convention, not enforcement. A property is a **descriptor**: an object stored in the class with a `__get__` method and, if it also has `__set__`, a *data descriptor*, which is why rule 3a of §3.1 lets it outrank the instance `__dict__`. A plain function has `__get__` but no `__set__`; that `__get__` is what turns `p.norm` into a bound method (Example 4.1).

**`@classmethod`** receives the class instead of the instance. Its use is the *alternative constructor*; because it calls `cls(...)`, a subclass inherits it and gets its own type back. **`@staticmethod`** receives nothing: an ordinary function that lives in the class namespace because it belongs there conceptually (§3.6's `_lowest_terms`). If unsure, use a module-level function.

▶ **Build it yourself (1d): finishing `Vec`.**

```python
class Vec:
    @property
    def x(self):
        return self._c[0]

    @property
    def y(self):
        return self._c[1]

    @property
    def z(self):
        return self._c[2]

    @classmethod
    def from_polar(cls, r, theta):
        return cls(r * math.cos(theta), r * math.sin(theta))


p = Vec.from_polar(2, math.pi / 3)
print(p, p.x, p.y, abs(p))
```
```text
(1, 1.73205) 1.0000000000000002 1.7320508075688772 2.0
```

### 3.6 `Rational`: equality, hashing and order for a number type

▶ **Build it yourself (2): exact fractions.** A rational number is a class of pairs ($\tfrac12 = \tfrac24$). To make `==` a comparison of two integers, store one representative: lowest terms, positive denominator. Predict first: `Rational(3) == 3` should be `True`. What does the contract of §3.2 then demand of `hash(Rational(3))`?

1. **Normalize** by $g = \gcd(n, d)$, using $-g$ when $d < 0$ so the sign moves to the numerator. `operator.index` accepts ints (and NumPy ints) but rejects floats, so `Rational(0.5)` fails loudly.
2. **Equality and hash.** Compare the normalized pairs. Since `Rational(3) == 3`, the contract forces `hash(Rational(3)) == hash(3)`: hash the bare numerator when $d = 1$.
3. **Order.** For $d, d' > 0$: $\tfrac nd < \tfrac{n'}{d'}$ exactly when $nd' < n'd$. Define `__lt__`; `total_ordering` does the rest.
4. **Arithmetic.** $\tfrac nd + \tfrac{n'}{d'} = \tfrac{nd' + n'd}{dd'}$, $\tfrac nd\cdot\tfrac{n'}{d'} = \tfrac{nn'}{dd'}$, and division multiplies by the reciprocal. One decorator promotes an `int` operand to a `Rational` and declines everything else; reflected versions are then one line each.

```python
import functools
import operator
from math import gcd


def _coerced(method):
    """Accept a Rational or an int (promoted to Rational); decline anything else."""
    @functools.wraps(method)
    def wrapper(self, other):
        if isinstance(other, numbers.Integral):
            other = Rational(other)
        elif not isinstance(other, Rational):
            return NotImplemented
        return method(self, other)
    return wrapper


@functools.total_ordering
class Rational:
    def __init__(self, n, d=1):
        n, d = operator.index(n), operator.index(d)
        if d == 0:
            raise ZeroDivisionError("Rational with zero denominator")
        self._n, self._d = Rational._lowest_terms(n, d)

    @staticmethod
    def _lowest_terms(n, d):
        g = gcd(n, d) if d > 0 else -gcd(n, d)
        return n // g, d // g

    numerator = property(lambda self: self._n)
    denominator = property(lambda self: self._d)

    def __repr__(self):
        return f"Rational({self._n}, {self._d})"

    def __str__(self):
        return str(self._n) if self._d == 1 else f"{self._n}/{self._d}"

    @_coerced
    def __eq__(self, other):
        return (self._n, self._d) == (other._n, other._d)

    def __hash__(self):
        return hash(self._n) if self._d == 1 else hash((self._n, self._d))

    @_coerced
    def __lt__(self, other):
        return self._n * other._d < other._n * self._d

    def __bool__(self):
        return self._n != 0
```

```python
class Rational:
    def __neg__(self):
        return Rational(-self._n, self._d)

    @_coerced
    def __add__(self, other):
        return Rational(self._n * other._d + other._n * self._d, self._d * other._d)

    __radd__ = __add__

    @_coerced
    def __sub__(self, other):
        return self + (-other)

    @_coerced
    def __rsub__(self, other):
        return other + (-self)

    @_coerced
    def __mul__(self, other):
        return Rational(self._n * other._n, self._d * other._d)

    __rmul__ = __mul__

    @_coerced
    def __truediv__(self, other):
        return Rational(self._n * other._d, self._d * other._n)

    @_coerced
    def __rtruediv__(self, other):
        return other / self

    def __pow__(self, k):
        if not isinstance(k, int):
            return NotImplemented
        return Rational(self._n ** k, self._d ** k) if k >= 0 else Rational(self._d ** -k, self._n ** -k)


a = Rational(6, -8)
print(a, repr(a), Rational(0, 5), Rational(2, 4) == Rational(1, 2))
print(sorted([Rational(1, 2), Rational(1, 3), Rational(-1, 2)]))
print(Rational(3) == 3, hash(Rational(3)) == hash(3), {3: "x"}[Rational(3)])
print(Rational(1, 3) + Rational(1, 6), 1 - Rational(1, 3), 2 * Rational(3, 4), 3 / Rational(1, 2), Rational(2, 3) ** -2)
```
```text
-3/4 Rational(-3, 4) 0 True
[Rational(-1, 2), Rational(1, 3), Rational(1, 2)]
True True x
1/2 2/3 3/2 6 9/4
```

■ **What just happened.**

Every operator is one expression on normalized pairs; the decorator carries the type checking. The library's `fractions.Fraction` also equals *floats* (`Fraction(1, 2) == 0.5`) and so hashes like them: Python's numeric tower gives equal numbers of any type equal hashes. Ours equals ints only, so we owe only the int hash. `project.py` Part 2 checks all four operations, order and equality against `Fraction` on 1000 random pairs.

### 3.7 `Dual`: a number that carries its derivative

Let $\mathbb{D}$ be the expressions $a + b\varepsilon$, $a, b \in \mathbb{R}$, added and multiplied like polynomials in $\varepsilon$ with one extra rule, $\varepsilon^2 = 0$. (The matrix $E = \begin{pmatrix}0&1\\0&0\end{pmatrix}$ has $E^2 = 0$ and $a + b\varepsilon \leftrightarrow aI + bE$ is a faithful copy: $E$ is the nilpotent part of a Jordan block, LA-11. `project.py` checks that matrix products reproduce (3.1).)

▶ **Try it first.** On paper, expand $(a + b\varepsilon)(c + d\varepsilon)$. Then find $p + q\varepsilon$ with $(c + d\varepsilon)(p + q\varepsilon) = a + b\varepsilon$ (division). Compare both with rules from calculus.

**Derivation of the algebra.**

1. *Sum:* $(a + b\varepsilon) \pm (c + d\varepsilon) = (a \pm c) + (b \pm d)\varepsilon$.
2. *Product:* $(a + b\varepsilon)(c + d\varepsilon) = ac + (ad + bc)\,\varepsilon + bd\,\varepsilon^2 = ac + (ad + bc)\,\varepsilon.$ $\tag{3.1}$
3. *Quotient, $c \ne 0$:* $(c + d\varepsilon)(p + q\varepsilon) = cp + (cq + dp)\varepsilon$ must equal $a + b\varepsilon$, so $p = a/c$ and $q = (b - dp)/c = (bc - ad)/c^2$:
   $$\frac{a + b\varepsilon}{c + d\varepsilon} = \frac ac + \frac{bc - ad}{c^2}\,\varepsilon. \tag{3.2}$$
4. *Reading the coefficients.* Put $a = u(x)$, $b = u'(x)$, $c = v(x)$, $d = v'(x)$. Then (3.1) is $(uv)' = u'v + uv'$ and (3.2) is $(u/v)' = (u'v - uv')/v^2$. If every value travels as the pair (value, derivative), arithmetic maintains both.
5. *Functions.* If $f$ has a convergent Taylor series at $a$, then $f(a + b\varepsilon) = \sum_k f^{(k)}(a)(b\varepsilon)^k/k!$, and every term with $k \ge 2$ contains $\varepsilon^2 = 0$:
   $$f(a + b\varepsilon) = f(a) + f'(a)\,b\,\varepsilon. \tag{3.3}$$
   This *is* the chain rule: if $a + b\varepsilon = g(x) + g'(x)\varepsilon$, (3.3) gives $f(g(x)) + f'(g(x))g'(x)\,\varepsilon$. For integer $n$, $(a + b\varepsilon)^n = a^n + na^{n-1}b\,\varepsilon$ (binomial theorem; negative $n$ by (3.2)).

**Theorem 3.1 (forward-mode differentiation).** Let $F(x)$ be built from $x$, real constants, $+,-,\times,\div$, integer powers and functions $f$ satisfying (3.3). Evaluating $F$ at $x + 1\cdot\varepsilon$ returns $F(x) + F'(x)\,\varepsilon$.

*Proof.* Induction on the expression tree. The leaf $x$ becomes $x + 1\varepsilon$ and a constant $c$ becomes $c + 0\varepsilon$. If $U, V$ evaluate to $u + u'\varepsilon$ and $v + v'\varepsilon$ with $u', v'$ their true derivatives, then by step 1, (3.1), (3.2) and (3.3) the expressions $U \pm V$, $UV$, $U/V$, $U^n$, $f(U)$ evaluate to a value plus the sum, product, quotient or chain rule applied to $(u, u', v, v')$, which is the true derivative of each. ∎

■ **What just happened.**

Automatic differentiation is the product rule stored in a tuple. There is no symbolic algebra and no step size, hence no truncation error; each floating-point operation still rounds, so Example 4.4 reports an error near $10^{-16}$, not $0$. Where $F$ is not differentiable (`abs` at 0, a branch), you get the derivative of the branch taken.

▶ **Build it yourself (3): the `Dual` class.** Predict, using (3.1) and (3.2), the results of `x * x * x`, `1 / x` and `(x + 1) / (x - 1)` for `x = Dual(3.0, 1.0)`. Operators must accept plain numbers (`2 * x`, `x + 1`): the lifting decorator of `Rational` does it, promoting `n` to `Dual(n, 0)`. Equality takes only another `Dual`, so unlike `Rational` it owes no int hash.

```python
def _dual_operand(method):
    """Promote a real number to Dual(number, 0); decline anything else."""
    @functools.wraps(method)
    def wrapper(self, other):
        if isinstance(other, numbers.Real):
            other = Dual(other, 0.0)
        elif not isinstance(other, Dual):
            return NotImplemented
        return method(self, other)
    return wrapper


class Dual:
    """a + b*eps with eps**2 = 0: a value a and its derivative b."""
    __slots__ = ("a", "b")

    def __init__(self, a, b=0.0):
        self.a, self.b = a, b

    def __repr__(self):
        return f"Dual({self.a!r}, {self.b!r})"

    def __eq__(self, other):
        if not isinstance(other, Dual):
            return NotImplemented
        return (self.a, self.b) == (other.a, other.b)

    def __hash__(self):
        return hash((self.a, self.b))

    def __neg__(self):
        return Dual(-self.a, -self.b)

    @_dual_operand
    def __add__(self, o):
        return Dual(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    @_dual_operand
    def __sub__(self, o):
        return self + (-o)

    @_dual_operand
    def __rsub__(self, o):
        return o + (-self)

    @_dual_operand
    def __mul__(self, o):                                         # (3.1)
        return Dual(self.a * o.a, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    @_dual_operand
    def __truediv__(self, o):                                     # (3.2)
        return Dual(self.a / o.a, (self.b * o.a - self.a * o.b) / (o.a * o.a))

    @_dual_operand
    def __rtruediv__(self, o):
        return o / self

    def __pow__(self, n):                                         # integer powers, step 5
        if not isinstance(n, int):
            return NotImplemented
        if n == 0:
            return Dual(1.0, 0.0)
        return Dual(self.a ** n, n * self.a ** (n - 1) * self.b)


x = Dual(3.0, 1.0)
print(x * x * x, 1 / x, (x + 1) / (x - 1), 2 * x - 5, x ** 3)
```
```text
Dual(27.0, 27.0) Dual(0.3333333333333333, -0.1111111111111111) Dual(2.0, -0.5) Dual(1.0, 2.0) Dual(27.0, 27.0)
```

Now the functions. By (3.3) each needs only its own derivative. A factory builds the Dual-aware version of `f` from `f` and `df`; because `df` is itself Dual-aware, the chain rule keeps working however deeply `Dual`s nest, which second derivatives need ("Your turn").

```python
from types import SimpleNamespace


def _elementary(f, df):
    """Dual-aware f: f(a + b eps) = f(a) + df(a) * b eps  (3.3)."""
    def wrapped(x):
        if isinstance(x, Dual):
            return Dual(wrapped(x.a), df(x.a) * x.b)
        return f(x)
    return wrapped


dmath = SimpleNamespace()                       # stands in for a module dmath.py
dmath.sin = _elementary(math.sin, lambda a: dmath.cos(a))
dmath.cos = _elementary(math.cos, lambda a: -dmath.sin(a))
dmath.exp = _elementary(math.exp, lambda a: dmath.exp(a))
dmath.log = _elementary(math.log, lambda a: 1 / a)
dmath.sqrt = _elementary(math.sqrt, lambda a: 0.5 / dmath.sqrt(a))

print(dmath.sin(Dual(0.5, 1.0)), (math.sin(0.5), math.cos(0.5)))
print(dmath.sin(Dual(1.0, 1.0) ** 2), (math.sin(1.0), 2 * math.cos(1.0)))
```
```text
Dual(0.479425538604203, 0.8775825618903728) (0.479425538604203, 0.8775825618903728)
Dual(0.8414709848078965, 1.0806046117362795) (0.8414709848078965, 1.0806046117362795)
```

The second line is the chain rule, $\frac{d}{dx}\sin(x^2) = 2x\cos(x^2)$ at $x = 1$. A formula written `x * dmath.sin(x) / (1 + x**2)` works for floats (`dmath.sin` falls through to `math.sin`) and for `Dual`s, because `*`, `/`, `+`, `**` are duck-typed: nobody edited the formula.

**Many inputs.** For $F: \mathbb{R}^n \to \mathbb{R}^m$, seed input $j$ with $1\varepsilon$ and the others with $0\varepsilon$; the $\varepsilon$ parts of the outputs are column $j$ of the Jacobian $\partial F_i/\partial x_j$ (LA-14): $n$ forward passes (Example 4.7). Reverse mode (ML-2, ML-7) runs the chain rule from the outputs backward and gets the gradient of a scalar loss ($m = 1$) in one pass whatever $n$ is: backpropagation.

### 3.8 `Matrix`: a mutable container

▶ **Build it yourself (4): a matrix class.** Unlike `Vec`, a `Matrix` supports `A[i, j] = x`, so it is *mutable*: unhashable, and `A += B` should update in place (rule 3 of §3.3). `solve` is LA-2's algorithm: elimination with partial pivoting (the largest entry in the column is the pivot) stores the multipliers in $L$ and gives $PA = LU$ (LA-2 Theorem 3.13); then forward substitution $L\mathbf y = P\mathbf b$ (LA-2 (3.4)) and back substitution $U\mathbf x = \mathbf y$ (LA-2 (3.3)).

```python
class Matrix:
    def __init__(self, rows):
        self._rows = [[float(x) for x in row] for row in rows]
        if len({len(r) for r in self._rows}) > 1:
            raise ValueError("rows have different lengths")

    shape = property(lambda self: (len(self._rows), len(self._rows[0])))
    T = property(lambda self: Matrix(zip(*self._rows)))

    def __repr__(self):
        return f"Matrix({self._rows!r})"

    def __eq__(self, other):
        if not isinstance(other, Matrix):
            return NotImplemented
        return self._rows == other._rows

    __hash__ = None                                  # mutable, so unhashable

    def __getitem__(self, key):
        if isinstance(key, tuple):                   # A[i, j]
            i, j = key
            return self._rows[i][j]
        return Vec(*self._rows[key])                 # A[i] is row i

    def __setitem__(self, key, value):
        i, j = key
        self._rows[i][j] = float(value)

    def __iadd__(self, other):
        if not isinstance(other, Matrix):
            return NotImplemented
        if self.shape != other.shape:
            raise ValueError(f"shapes {self.shape} and {other.shape} differ")
        for r, s in zip(self._rows, other._rows):
            r[:] = [a + b for a, b in zip(r, s)]
        return self                                  # forgetting this binds the name to None

    def __matmul__(self, other):
        if isinstance(other, Vec):
            return Vec(*(sum(a * b for a, b in zip(row, other)) for row in self._rows))
        if not isinstance(other, Matrix):
            return NotImplemented
        if self.shape[1] != other.shape[0]:
            raise ValueError(f"shapes {self.shape} and {other.shape} do not align")
        cols = list(zip(*other._rows))
        return Matrix([[sum(a * b for a, b in zip(row, col)) for col in cols] for row in self._rows])
```

The `lu` method is the part to type yourself:

```python
class Matrix:
    def lu(self):
        """Return (perm, L, U) with P A = L U; row i of P A is row perm[i] of A."""
        n = self.shape[0]
        if self.shape[1] != n:
            raise ValueError("LU needs a square matrix")
        U = [row[:] for row in self._rows]
        L = [[0.0] * n for _ in range(n)]
        perm = list(range(n))
        for k in range(n):
            p = max(range(k, n), key=lambda i: abs(U[i][k]))        # partial pivoting
            if U[p][k] == 0.0:
                raise ValueError("matrix is singular")
            if p != k:
                U[k], U[p] = U[p], U[k]
                L[k], L[p] = L[p], L[k]                              # swap the multipliers found so far
                perm[k], perm[p] = perm[p], perm[k]
            for i in range(k + 1, n):
                L[i][k] = U[i][k] / U[k][k]                          # the multiplier
                for j in range(k, n):
                    U[i][j] -= L[i][k] * U[k][j]
        for i in range(n):
            L[i][i] = 1.0
        return perm, Matrix(L), Matrix(U)

    def solve(self, b):
        """Solve A x = b: L y = P b forward, then U x = y backward."""
        perm, L, U = self.lu()
        n, b = len(perm), list(b)
        y = [0.0] * n
        for i in range(n):
            y[i] = float(b[perm[i]]) - sum(L._rows[i][j] * y[j] for j in range(i))
        x = [0.0] * n
        for i in reversed(range(n)):
            x[i] = (y[i] - sum(U._rows[i][j] * x[j] for j in range(i + 1, n))) / U._rows[i][i]
        return Vec(*x)


A = Matrix([[1, 2], [3, 4]])
alias = A
A += Matrix([[10, 0], [0, 10]])                  # __iadd__ mutates and returns self
print(alias is A, alias, A[1, 0], A[1], A.T[0, 1])
v = Vec(1, 2)
w = v
v += Vec(1, 1)                                   # Vec has no __iadd__: v = v + Vec(1, 1), a new object
print(w, v, w is v)
try:
    {A: 1}
except TypeError as err:
    print("TypeError:", err)
```
```text
True Matrix([[11.0, 2.0], [3.0, 14.0]]) 3.0 (3, 14) 3.0
(1, 2) (2, 3) False
TypeError: unhashable type: 'Matrix'
```

Example 4.5 checks `solve` against NumPy; `project.py` Part 4 does 200 random $6\times6$ systems.

### 3.9 `@dataclass`: the boilerplate you should not write

Most classes hold data, and for those you keep writing `__init__`, `__repr__`, `__eq__`, then `__hash__` and `__lt__`. `@dataclass` generates them from annotated class-level names.

▶ **Build it yourself (5): what `@dataclass` writes for you.** Write by hand what a record `Measurement(value, unit, time)` needs: constructor, round-tripping repr, equality, a consistent hash, and `<` comparing the fields as a tuple. About twenty lines: what you are about to be spared. Then the rules:

1. **Fields** are annotated names in the class body; declaration order is the `__init__` parameter order. Annotations are not enforced.
2. **Defaults** may be simple values. A list, dict or set default raises `ValueError` (the shared-object trap of §3.1); write `field(default_factory=list)` for a fresh one per instance.
3. **Options.** `frozen=True`: assignment raises `FrozenInstanceError`, and with the default `eq=True` a consistent `__hash__` is generated. `eq=True` without `frozen` sets `__hash__` to `None` (§3.2). `order=True` generates `<`, `<=`, `>`, `>=` comparing the fields as a tuple in declaration order, so **field order is sort order**. `slots=True` and `kw_only=True` also exist.
4. **`__post_init__`** runs after the generated `__init__` (validation, derived fields); in a frozen class assign with `object.__setattr__`. `replace(obj, f=new)` returns a modified copy.

```python
from dataclasses import dataclass, field, replace, FrozenInstanceError


@dataclass(frozen=True, order=True)
class Measurement:
    value: float
    unit: str
    time: float


m = Measurement(21.5, "degC", 3.0)
twin = Measurement(21.5, "degC", 3.0)
print(m, m == twin, hash(m) == hash(twin), replace(m, value=22.0))
try:
    m.value = 0.0
except FrozenInstanceError as err:
    print("FrozenInstanceError:", err)
try:
    @dataclass
    class Run:
        name: str
        readings: list = []
except ValueError as err:
    print("ValueError:", err)


@dataclass
class Run:
    name: str
    readings: list[Measurement] = field(default_factory=list)


r1, r2 = Run("a"), Run("b")
r1.readings.append(m)
print(r1, r2.readings, r1 == Run("a", [m]))
```
```text
Measurement(value=21.5, unit='degC', time=3.0) True True Measurement(value=22.0, unit='degC', time=3.0)
FrozenInstanceError: cannot assign to field 'value'
ValueError: mutable default <class 'list'> for field readings is not allowed: use default_factory
Run(name='a', readings=[Measurement(value=21.5, unit='degC', time=3.0)]) [] True
```

A `Polynomial` shows `__post_init__`, `__call__` and the payoff of duck typing. Coefficients are ascending, $p(x) = c_0 + c_1x + c_2x^2 + \cdots$, evaluated by Horner's rule $c_0 + x(c_1 + x(c_2 + \cdots))$, which uses only `*` and `+`:

```python
@dataclass(frozen=True)
class Polynomial:
    coeffs: tuple

    def __post_init__(self):
        c = [float(a) for a in self.coeffs]
        while len(c) > 1 and c[-1] == 0.0:                    # canonical form: no trailing zeros
            c.pop()
        object.__setattr__(self, "coeffs", tuple(c))

    def __call__(self, x):
        acc = 0.0
        for c in reversed(self.coeffs):
            acc = acc * x + c
        return acc

    def deriv(self):
        return Polynomial(tuple(k * c for k, c in enumerate(self.coeffs))[1:] or (0.0,))


p = Polynomial((1, -3, 0, 2, 0))                              # 1 - 3x + 2x^3
print(p, p(2.0), p.deriv(), p.deriv()(2.0))
print(p(Dual(2.0, 1.0)))
```
```text
Polynomial(coeffs=(1.0, -3.0, 0.0, 2.0)) 11.0 Polynomial(coeffs=(-3.0, 0.0, 6.0)) 21.0
Dual(11.0, 21.0)
```

`Polynomial` was written without knowing `Dual` exists, yet `p(Dual(2.0, 1.0))` returns $p(2)$ and $p'(2)$, agreeing with `p.deriv()(2.0)`. `project.py` Part 5 compares it with `numpy.polynomial.Polynomial`; exercise B6 adds `+` and `*`.

**When not to use one:** when identity or behavior defines the object (a file handle, a `Matrix` with invariants to maintain, an `__init__` that does real work), when you need control over `==`, or when the class is mostly methods. For records it is the default.

### 3.10 Inheritance, abstract base classes, protocols, composition

**Inheritance and `super()`.** `class Child(Parent)` gives `Child` the parent's attributes by rule 3c of §3.1; a method defined in `Child` *overrides* the parent's, and `isinstance(x, Parent)` is true for a `Child`. `super().method()` calls the next implementation along the **method resolution order**, `Child.__mro__`, which is not "the parent" once there are several bases: in a diamond `Both(Left, Right)` over a common `Base`, the `super()` call inside `Left` reaches `Right`, a class `Left` has never heard of. That makes cooperative multiple inheritance work, and it is also why inheritance couples classes tightly.

**The fragile base class.** A subclass depends on how its parent calls its own methods, which is rarely documented. Count the writes into a dict:

```python
from collections import UserDict


class CountingDict(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.writes = 0

    def __setitem__(self, key, value):
        self.writes += 1
        super().__setitem__(key, value)


class CountingUserDict(UserDict):
    def __init__(self, *args, **kwargs):
        self.writes = 0                           # first: UserDict.__init__ may already call __setitem__
        super().__init__(*args, **kwargs)

    def __setitem__(self, key, value):
        self.writes += 1
        super().__setitem__(key, value)


for cls in (CountingDict, CountingUserDict):
    d = cls()
    d["a"] = 1
    d.update(b=2, c=3)
    d.setdefault("e", 5)
    print(cls.__name__, "stored", len(d), "keys but counted", d.writes, "writes")
```
```text
CountingDict stored 4 keys but counted 1 writes
CountingUserDict stored 4 keys but counted 4 writes
```

CPython's `dict.update` and `setdefault` do not call your `__setitem__`; `UserDict`, built on the abstract `MutableMapping`, routes everything through it. Do not subclass `dict`, `list` or `set` to change behavior. Subclass the `collections.UserX` or abstract classes, or **compose**: hold a `dict` inside your object and expose only the methods you want.

**Abstract base classes.** `abc.ABC` with `@abstractmethod` declares an interface: a subclass that leaves an abstract method unimplemented cannot be instantiated. The classes in `collections.abc` supply *mixin methods for free*: implement `__getitem__` and `__len__` and `Sequence` gives you `__contains__`, `__iter__`, `__reversed__`, `index` and `count`.

```python
from collections.abc import Iterable, Sequence


class Samples(Sequence):
    def __init__(self, s):
        self._s = list(s)

    def __len__(self):
        return len(self._s)

    def __getitem__(self, i):
        return self._s[i]

    def __repr__(self):
        return f"Samples({self._s!r})"


x = Samples([5, 7, 7, 9])
print(x.index(7), x.count(7), 9 in x, list(reversed(x)))
try:
    Sequence()
except TypeError as err:
    print("TypeError:", err)
print(isinstance(x, Iterable), isinstance(Vec(1, 2), Iterable), isinstance(3, Iterable))
```
```text
1 2 True [9, 7, 7, 5]
TypeError: Can't instantiate abstract class Sequence without an implementation for abstract methods '__getitem__', '__len__'
True True False
```

The last line bridges to protocols: `Vec` does not inherit from `Iterable`, yet `isinstance(Vec(1, 2), Iterable)` is `True`, because `Iterable` only checks for an `__iter__`. Inherit from an ABC when you plug into a framework that supplies the rest (`Sequence`, `MutableMapping`, `Exception`).

**Protocols: structural typing.** `typing.Protocol` declares "anything with these methods", with no inheritance. Static checkers (PY-5) verify call sites; at run time nothing is checked unless you add `@runtime_checkable`, and then `isinstance` tests only that the names exist, not their signatures.

```python
from typing import Protocol, runtime_checkable
import numpy as np


@runtime_checkable
class Evaluable(Protocol):
    def __call__(self, x: float) -> float: ...


class Reading:                                    # NumPy's protocol: "I can turn myself into an array"
    def __array__(self, dtype=None, copy=None):
        return np.array([1.0, 2.0, 3.0])


print(isinstance(Polynomial((1, 2)), Evaluable), isinstance(math.sin, Evaluable), isinstance(3, Evaluable))
print(np.sum(Reading()), np.asarray(Vec(1, 2, 3)).shape)
```
```text
True True False
6.0 ()
```

A `Polynomial` (a class) and `math.sin` (a built-in function) share no ancestor, and both satisfy `Evaluable`. NumPy accepts `Reading` because `np.asarray` looks for `__array__`; an object with `__len__` and `__getitem__` it reads as a sequence, and anything else as an opaque scalar. `Vec` has no `__getitem__` (it is "Your turn" 1), so NumPy sees a 0-dimensional object array, shape `()`. PY-6 builds on this.

▶ **Build it yourself (6): inheritance versus composition.** Integration rules approximate $\int_a^b f$ from samples $y_i = f(a + ih)$, $h = (b-a)/n$: the trapezoid rule $h\big(\tfrac12 y_0 + y_1 + \cdots + y_{n-1} + \tfrac12 y_n\big)$ and Simpson's rule ($n$ even) $\tfrac h3\big(y_0 + 4y_1 + 2y_2 + \cdots + 4y_{n-1} + y_n\big)$. Design 1: an abstract `Integrator` builds the grid and each subclass supplies a rule. Design 2: one function takes the rule as an argument. Write both first; predict which is shorter and what happens to the error when $n$ doubles.

```python
from abc import ABC, abstractmethod


class Integrator(ABC):
    def __init__(self, n):
        self.n = n

    def integrate(self, f, a, b):
        h = (b - a) / self.n
        return self._combine([f(a + i * h) for i in range(self.n + 1)], h)

    @abstractmethod
    def _combine(self, ys, h): ...


class Trapezoid(Integrator):
    def _combine(self, ys, h):
        return h * (sum(ys) - 0.5 * (ys[0] + ys[-1]))


class Simpson(Integrator):
    def _combine(self, ys, h):
        return h / 3 * (ys[0] + ys[-1] + 4 * sum(ys[1:-1:2]) + 2 * sum(ys[2:-1:2]))


def trapezoid(ys, h):
    return h * (sum(ys) - 0.5 * (ys[0] + ys[-1]))


def simpson(ys, h):
    return h / 3 * (ys[0] + ys[-1] + 4 * sum(ys[1:-1:2]) + 2 * sum(ys[2:-1:2]))


def integrate(f, a, b, n, rule):
    h = (b - a) / n
    return rule([f(a + i * h) for i in range(n + 1)], h)


for name, cls, rule in (("trapezoid", Trapezoid, trapezoid), ("simpson", Simpson, simpson)):
    ns = (8, 16, 32)
    errs = [abs(cls(n).integrate(math.sin, 0.0, math.pi) - 2.0) for n in ns]
    same = all(cls(n).integrate(math.sin, 0.0, math.pi) == integrate(math.sin, 0.0, math.pi, n, rule) for n in ns)
    print(f"{name:9s} errors {[f'{e:.2e}' for e in errs]}  ratios {[round(a / b, 1) for a, b in zip(errs, errs[1:])]}  same: {same}")
```
```text
trapezoid errors ['2.58e-02', '6.43e-03', '1.61e-03']  ratios [4.0, 4.0]  same: True
simpson   errors ['2.69e-04', '1.66e-05', '1.03e-06']  ratios [16.2, 16.1]  same: True
```

■ **What just happened.**

Both designs return identical floats, and the error ratio when $n$ doubles is $\approx 4$ for the trapezoid rule (error $\propto h^2$) and $\approx 16$ for Simpson ($\propto h^4$): each rule's order, read off a table. As design, the function version is half the length, a rule can be tested alone on a list of samples, a new rule is one function (even a `lambda`), and nothing needs a constructor. The class version buys an `isinstance` check and a home for state, but its base class fixes the skeleton (a uniform grid) for every subclass; a rule with non-uniform nodes, Gauss–Legendre say, means rewriting the base. Use classes when integrators carry state (a tolerance, an evaluation counter), and even then let an `Integrator` *have* a rule rather than *be* one. **Inherit for "is-a" with substitutability and shared behavior, as `Samples` is a `Sequence`; otherwise compose.**

### 3.11 `__slots__`, and rules of thumb

Each instance normally carries a `__dict__` for its attributes. `__slots__ = ("a", "b")`, as in `Dual`, replaces it with fixed storage for exactly those names: smaller instances, and a typo creates no attribute. It pays off for millions of small objects and is premature elsewhere; `project.py` Part 5 measures the saving.

```python
d = Dual(1.0, 2.0)
print(hasattr(d, "__dict__"), hasattr(Vec(1, 2), "__dict__"))
try:
    d.c = 3
except AttributeError as err:
    print("AttributeError:", err)
```
```text
False True
AttributeError: 'Dual' object has no attribute 'c' and no __dict__ for setting new attributes
```

**Rules of thumb.** One responsibility and few methods. Prefer immutable state (a frozen dataclass, a tuple inside): hashable and safe to share. Write `__repr__` first. Create per-instance state in `__init__`, never in the class body. Reach for a dataclass before hand-writing `__init__`. Do not subclass built-in containers. Start with plain attributes and add `@property` when you need validation. Test the contracts: `eval(repr(x)) == x`, equal objects hash equal, `NotImplemented` for foreign operands.

## 4. Worked examples

### Example 4.1: `obj.method()` is `Class.method(obj)`

```python
class Point:
    def __init__(self, x, y):
        self.x, self.y = x, y

    def norm(self):
        return (self.x ** 2 + self.y ** 2) ** 0.5


p = Point(3, 4)
bound = p.norm                                       # attribute lookup builds a bound method
print(bound.__self__ is p, bound.__func__ is Point.__dict__["norm"], "norm" in p.__dict__)
print(p.norm(), Point.norm(p), Point.__dict__["norm"].__get__(p, Point)())
print(f"id(p) = {id(p)}, id(bound.__self__) = {id(bound.__self__)}, id(Point.norm) = {id(Point.norm)}")
try:
    p.norm(5)
except TypeError as err:
    print("TypeError:", err)
```
```text
True True False
5.0 5.0 5.0
id(p) = 140261945239472, id(bound.__self__) = 140261945239472, id(Point.norm) = 140261902338464
TypeError: Point.norm() takes 1 positional argument but 2 were given
```

`p.norm` finds the function in the *class* (not in `p.__dict__`: one function serves every instance). A function is a non-data descriptor, so lookup calls `function.__get__(p, Point)`, which returns a bound method holding `__self__` (here `p`) and `__func__`; calling it puts `p` in front. The three spellings on the second line are one call, and `p` and `bound.__self__` have the same id. The error counts `self`: two arguments, `p` and `5`, for a function that takes one.

### Example 4.2: `__eq__` without `__hash__`

```python
class Pt:
    def __init__(self, x, y):
        self.x, self.y = x, y

    def __eq__(self, other):
        return isinstance(other, Pt) and (self.x, self.y) == (other.x, other.y)


class BadPt(Pt):
    __hash__ = object.__hash__                       # restores the by-address hash: legal, wrong


class GoodPt(Pt):
    def __hash__(self):
        return hash((self.x, self.y))


print(Pt(1, 2) == Pt(1, 2), Pt.__hash__)
try:
    {Pt(1, 2)}
except TypeError as err:
    print("TypeError:", err)
print(BadPt(1, 2) in {BadPt(1, 2)}, GoodPt(1, 2) in {GoodPt(1, 2)})
```
```text
True None
TypeError: unhashable type: 'Pt'
False True
```

Defining `__eq__` set `Pt.__hash__` to `None`, so `Pt` objects cannot go in a set or be dict keys. The message ("unhashable type") never mentions `__eq__`, which is why it confuses. `BadPt` silences the error by restoring the address hash; now two equal points land in different slots and `BadPt(1, 2) in {BadPt(1, 2)}` is `False`: the contract is broken. `GoodPt` hashes the fields `__eq__` compares.

### Example 4.3: Operator dispatch, traced

```python
class Left:
    def __add__(self, o):
        print("  Left.__add__"); return NotImplemented
    def __radd__(self, o):
        print("  Left.__radd__"); return "from Left.__radd__"


class Right:
    def __radd__(self, o):
        print("  Right.__radd__"); return "from Right.__radd__"


class SubLeft(Left):
    def __radd__(self, o):
        print("  SubLeft.__radd__"); return "from SubLeft.__radd__"


print(Left() + Right())
print(1 + Right())
print(Left() + SubLeft())
try:
    Left() + Left()
except TypeError as err:
    print("TypeError:", err)
```
```text
  Left.__add__
  Right.__radd__
from Right.__radd__
  Right.__radd__
from Right.__radd__
  SubLeft.__radd__
from SubLeft.__radd__
  Left.__add__
TypeError: unsupported operand type(s) for +: 'Left' and 'Left'
```

Line 1: `Left.__add__` declines, so `Right.__radd__` answers. Line 2: `int.__add__` declines and `Right.__radd__` answers; this is `2 * v`. Line 3: `SubLeft` subclasses the left operand's type and overrides `__radd__`, so it is asked *first*. Line 4: same types, so no reflected call, and `TypeError`.

### Example 4.4: One pass of `Dual` against SymPy

Differentiate $f(x) = x\sin x/(1 + x^2)$ at $x = 1.3$. Seed $x = 1.3 + 1\varepsilon$. By the product rule $u = x\sin x$ has slope $\sin x + x\cos x$, and $v = 1 + x^2$ has slope $2x = 2.6$.

```python
import sympy as sp

x = Dual(1.3, 1.0)
u = x * dmath.sin(x)
v = 1 + x ** 2
y = u / v
print("u =", u)
print("v =", v)
print("y =", y)
t = sp.symbols("t")
F = t * sp.sin(t) / (1 + t ** 2)
value, slope = float(F.subs(t, 1.3)), float(sp.diff(F, t).subs(t, 1.3))
print(f"SymPy: f = {value!r}, f' = {slope!r}")
print(f"|Dual - SymPy| = {abs(y.a - value):.1e}, {abs(y.b - slope):.1e}")
assert abs(y.a - value) < 1e-14 and abs(y.b - slope) < 1e-14
```
```text
u = Dual(1.2526256410423509, 1.3113066626291565)
v = Dual(2.6900000000000004, 2.6)
y = Dual(0.4656600896068218, 0.03739421176632707)
SymPy: f = 0.46566008960682176, f' = 0.0373942117663269
|Dual - SymPy| = 5.6e-17, 1.7e-16
```

By (3.2), $y' = (u'v - uv')/v^2 = (1.31131 \cdot 2.69 - 1.25263 \cdot 2.6)/2.69^2 = 0.27059/7.2361 = 0.037394$, matching the output. One evaluation, no step size, an error of one or two roundings (near $10^{-16}$) against $10^{-8}$ at best for the forward difference of §1.

### Example 4.5: A `Matrix` LU solve against NumPy

Solve $A\mathbf x = \mathbf b$ with $A = \begin{pmatrix}2&1&1\\4&-6&0\\-2&7&2\end{pmatrix}$, $\mathbf b = (5, -2, 9)^\top$. Partial pivoting takes the entry $4$ in row 2 as the first pivot, so `perm` starts `[1, 0, ...]`.

```python
import numpy as np

A = Matrix([[2, 1, 1], [4, -6, 0], [-2, 7, 2]])
b = [5, -2, 9]
perm, L, U = A.lu()
PA = Matrix([A._rows[i] for i in perm])
print("perm", perm)
print("L", L)
print("U", U)
print("P A == L U:", PA == L @ U)
y = L.solve([b[i] for i in perm])                 # L y = P b  (L is unit lower triangular)
x = U.solve(y)                                    # U x = y
print("y", y, " x", x, " numpy", np.linalg.solve(np.array(A._rows), np.array(b)))
print("residual", A @ A.solve(b) - Vec(*b))
```
```text
perm [1, 0, 2]
L Matrix([[1.0, 0.0, 0.0], [0.5, 1.0, 0.0], [-0.5, 1.0, 1.0]])
U Matrix([[4.0, -6.0, 0.0], [0.0, 4.0, 1.0], [0.0, 0.0, 1.0]])
P A == L U: True
y (-2, 6, 2)  x (1, 1, 2)  numpy [1. 1. 2.]
residual (0, 0, 0)
```

By hand: $P\mathbf b = (-2, 5, 9)$, forward substitution gives $\mathbf y = (-2, 6, 2)$, back substitution gives $x_3 = 2$, $x_2 = (6 - 2)/4 = 1$, $x_1 = (-2 + 6)/4 = 1$. The second column offers $4$ and $4$ as pivot, a tie resolved to the first. All entries are exact binary fractions, so $PA = LU$ holds with `==`; in general use a tolerance.

### Example 4.6: `order=True` sorts records by field order

```python
from operator import attrgetter

readings = [Measurement(21.5, "degC", 3.0), Measurement(19.0, "degC", 1.0),
            Measurement(21.5, "degC", 0.5), Measurement(20.2, "degC", 2.0)]
for m in sorted(readings):
    print(m)
print("by time:", [m.time for m in sorted(readings, key=attrgetter("time"))])
```
```text
Measurement(value=19.0, unit='degC', time=1.0)
Measurement(value=20.2, unit='degC', time=2.0)
Measurement(value=21.5, unit='degC', time=0.5)
Measurement(value=21.5, unit='degC', time=3.0)
by time: [0.5, 1.0, 2.0, 3.0]
```

The generated `<` compares `(value, unit, time)` as a tuple. Two readings tie on value $21.5$ and on unit, so time decides: $0.5$ before $3.0$. To sort by time the field order is wrong, so pass `key=attrgetter("time")` or declare `time` first.

### Example 4.7: Newton's method with no hand-written derivative; a Jacobian (from mechanics)

PY-1's `newton(f, df, x0)` needed `df` typed in. With `Dual` the derivative comes from `f` itself. Apply it to $\cos x = x$ (PY-1) and to Kepler's equation $E - e\sin E = M$ for the eccentric anomaly $E$ of an orbit with eccentricity $e = 0.6$ and mean anomaly $M = 1$ rad.

```python
def newton_ad(f, x0, tol=1e-14, max_iter=50):
    """Newton's method where f'(x) is the Dual part of f(Dual(x, 1))."""
    x = x0
    for k in range(1, max_iter + 1):
        y = f(Dual(x, 1.0))
        step = y.a / y.b
        x -= step
        if abs(step) <= tol:
            return x, k
    raise RuntimeError("Newton did not converge")


root, k1 = newton_ad(lambda x: dmath.cos(x) - x, 1.0)
e, M = 0.6, 1.0
E, k2 = newton_ad(lambda E: E - e * dmath.sin(E) - M, M)
print(root, k1)
print(E, k2, E - e * math.sin(E) - M)


def jacobian(f, x):
    """Rows of the Jacobian of f: R^n -> R^m, one forward pass per input."""
    cols = []
    for j in range(len(x)):
        out = f([Dual(xi, 1.0 if i == j else 0.0) for i, xi in enumerate(x)])
        cols.append([y.b for y in out])
    return [list(row) for row in zip(*cols)]


polar = lambda p: [p[0] * dmath.cos(p[1]), p[0] * dmath.sin(p[1])]
J = jacobian(polar, [2.0, 0.5])
print(J, J[0][0] * J[1][1] - J[0][1] * J[1][0])
```
```text
0.7390851332151607 5
1.5997485482275295 6 0.0
[[0.8775825618903728, -0.958851077208406], [0.479425538604203, 1.7551651237807455]] 2.0
```

The first root is PY-1's $0.7390851332151607$; the Kepler residual prints as $0$. The map $(r, \theta) \mapsto (r\cos\theta, r\sin\theta)$ has Jacobian $\begin{pmatrix}\cos\theta & -r\sin\theta\\ \sin\theta & r\cos\theta\end{pmatrix}$ (LA-14), with determinant $r$, the area factor of double integrals in polar coordinates. Two forward passes produced it; the determinant prints as $2.0 = r$.

## 5. Why it matters and where it connects

- **PY-5:** exceptions are classes (§3.10 is their inheritance story); its tests check the contracts of §3.2–§3.3; its type hints are what `Protocol` needs.
- **PY-6:** `numpy.ndarray` implements `__getitem__` with the tuple-of-slices keys of §3.4, `__matmul__` and in-place `__iadd__` (§3.8); `Matrix` is its slow ancestor, and `project.py` Part 4 measures the gap. **PY-9:** a `DataFrame` answers `[]`, `+` and `==` through the same dunders.
- **LA-14 and ML-2:** the Jacobian of Example 4.7 is forward-mode AD, and one pass per input is why reverse mode wins for a scalar loss with many parameters. **ML-7** derives backpropagation as that computation reversed; exercise C1 builds a small reverse-mode class.
- **LA-2 and LA-11:** `Matrix.lu` is LA-2's algorithm; $\varepsilon$ is LA-11's nilpotent Jordan block. **CM-1:** its `euler_step`, `rk4_step` and `symplectic_euler_step` take the acceleration as an argument: the composition design of §3.10.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 are predict-the-output, bug-fix and short-class drills; do all. B1–B7 are the module; B1 (`Vec` indexing), B3 (more `dmath` functions against SymPy), B4 (second derivatives by nesting `Dual`) and B7 (extend `project.py`) are essential. C1–C3: a reverse-mode differentiator, an integrator design with measurements, and the cost of pure-Python matrix multiplication.

## 7. Coding project

`project.py` builds `Vec`, `Rational`, `Dual` with `dmath`, `Matrix`, `Measurement`, `Polynomial` and both integrator designs, and tests them:

1. **`Vec`** against LA-1: Cauchy–Schwarz, the parallelogram law and Lagrange's identity on 2000 random pairs; dispatch.
2. **`Rational`** against `fractions.Fraction` on 1000 random pairs; the hash contract.
3. **`Dual`** on ten functions, each against SymPy (`lambdify`) and central differences; the error-versus-$h$ plot (a V for finite differences, a flat line near $10^{-16}$ for `Dual`) saved to `figures/`; a Jacobian; Newton with `Dual`; the matrix form of $\varepsilon$.
4. **`Matrix`** LU against `numpy.linalg.solve` on 200 random $6\times6$ systems, and a measured pure-Python versus NumPy multiply.
5. Every number quoted in §4 asserted.

**Build it yourself:** before opening the file, write `Dual` (arithmetic plus `dmath.sin`, `dmath.exp`) and test it on one function against SymPy. Then "Your turn": `Vec.__getitem__` with slices, a `Complex` class against the built-in, `Dual.__pow__` for non-integer exponents, and second derivatives by nesting.

## 8. Where you'll struggle

- **`self`.** It confuses for about a week, then it is obvious. The cure is §3.1: `obj.method(a)` is `Class.method(obj, a)`. "Takes 1 positional argument but 2 were given" after a one-argument call means you forgot `self` in the `def`.
- **A half-working operator type.** Without `__rmul__`, `v * 2` works and `2 * v` does not. Worse: if `__add__` *raises* `TypeError` (or `__eq__` returns `False`) for operands it does not understand instead of returning `NotImplemented`, the other operand never gets its turn. NumPy adds a quiet version: `array * v` treats `v` as a scalar and returns an array of `Vec`s with no error; `__array_ufunc__ = None` on your class makes NumPy decline.
- **`__eq__` without `__hash__`.** The error ("unhashable type") appears where you build a set, far from the class. Define both together or neither (Example 4.2).
- **Inheritance by reflex.** People subclass because that is what they were taught, then fight the parent. Composition is usually right (§3.10).
- **Hand-written boilerplate.** A dataclass removes about 80% of what simple classes contain, and people still type `__init__`. A mutable dataclass default raises because the alternative is one list shared by every instance.
- **Small traps.** `__iadd__` without `return self` rebinds `A` to `None` on `A += B`; a mutable class attribute is shared by all instances (§3.1).
- **Dual numbers feel like a trick.** They are the product rule stored in a tuple, and seeing that is the point. If Theorem 3.1 still feels like magic, differentiate $x\cdot x\cdot x$ by hand beside the trace in §3.7 and match the numbers.

## 9. Check yourself

1. Write `p.norm(5)` as a call of the underlying function. Where is the function stored?
2. In what order does Python search for `obj.name`? Why does a property beat an entry in `obj.__dict__`?
3. What are `repr` and `str` for, and which does the REPL use?
4. State the `__eq__`/`__hash__` contract. What happens to `__hash__` if you define only `__eq__`?
5. Trace `2 * v` for a `Vec` with and without `__rmul__`. When is the right operand tried first?
6. Why return `NotImplemented` instead of raising `TypeError` or returning `False`?
7. What object does `m[:, 2]` pass to `__getitem__`?
8. Expand $(a + b\varepsilon)(c + d\varepsilon)$ and name the rule. What is $\sin(a + b\varepsilon)$?
9. Why does `Polynomial.__call__` work on a `Dual` unedited?
10. For $f:\mathbb{R}^n \to \mathbb{R}$, how many forward passes give the gradient, and how many reverse passes?
11. Give one case for inheritance, one for composition, and the failing `dict` example.

## 10. Pair with

- **Ramalho, *Fluent Python* 2e, ch. 1, 11–13, 16.** Ch. 1 is the data model in one chapter; 11–12 build Pythonic objects and sequences (our `Vec`, `Samples`); 13 covers protocols and ABCs; 16 covers operator overloading and its dispatch rules. The best single companion.
- **Python Data Model reference, §3.3 (docs.python.org/3/reference/datamodel.html).** The authority for §3.1–§3.4; check it when a dunder surprises you.
- **Downey, *Think Python* 3e, ch. 15–18.** The gentlest introduction to classes and inheritance, with exercises.
- **Python Tutorial §9; the `dataclasses` documentation.** Class syntax in brief; every dataclass option.
- **Slatkin, *Effective Python* 2e, ch. 5.** Short numbered rules on classes and interfaces; a good review checklist.
- **Baydin et al., "Automatic differentiation in machine learning: a survey," §3.1 (reference only).** The dual-number construction of §3.7 and the forward/reverse comparison you need for ML-7.
