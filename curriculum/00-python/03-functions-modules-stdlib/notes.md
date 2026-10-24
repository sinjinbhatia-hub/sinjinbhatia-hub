# PY-3: Functions as Values, Closures and Decorators, Recursion, Modules and the Standard Library

**Time:** ~8 hours · **Prerequisites:** PY-2 · **Unlocks:** PY-4, PY-5, PY-10

> **Conventions in this module.** Python 3.11+; every output below is the real output of a run under Python 3.13. `>>>` blocks are interactive transcripts; other blocks are scripts followed by their output; `$` blocks are shell sessions (in a scratch directory, or in this module's folder when they mention `numtools`). Timings come from one run on one machine: trust ratios, not microseconds. $F_n$ is the Fibonacci sequence ($F_0 = 0$, $F_1 = 1$); $\varphi = (1+\sqrt5)/2$. Standard library only in the notes; `project.py` adds NumPy and Matplotlib. The package `numtools/` in this folder is the code deliverable (§3.8).

## 0. Where this sits

PY-1 gave you scope (LEGB), the mutable-default trap and a first closure; PY-2 ended with `sorted(key=...)` and generators. This module treats the function as an object. PY-4 builds on it (`@property` and `@dataclass` are decorators, §3.5), PY-5 (`@pytest.fixture` is a decorator) and PY-10 (caching, profiling and packaging all start here). Later courses use it without comment: `scipy.integrate.quad(f, a, b)` takes a function, and `import` is how any library reaches you.

## 1. The problem that forced this idea

Numerical code is full of "do this to that function": integrate $f$, find a root of $g$, time $h$, cache the expensive $k$. Given the PY-1 root finders and three test functions, you can paste a near-identical block for each pair, or write a table of functions and a double loop:

```python
import math

from numtools import bisect, secant

cases = {"cos x - x": (lambda x: math.cos(x) - x, 0.0, 1.0),
         "x^2 - 2": (lambda x: x * x - 2, 1.0, 2.0),
         "x^3 - 2x - 5": (lambda x: x**3 - 2 * x - 5, 2.0, 3.0)}
solvers = {"bisect": bisect, "secant": secant}      # same call signature: (f, a, b, tol=...)
for cname, (f, a, b) in cases.items():
    for sname, solve in solvers.items():
        root, iterations = solve(f, a, b, tol=1e-10)
        print(f"{cname:>13} {sname:>7}  root {root:.10f}  iterations {iterations:>2}")
```
```text
    cos x - x  bisect  root 0.7390851332  iterations 34
    cos x - x  secant  root 0.7390851332  iterations  6
      x^2 - 2  bisect  root 1.4142135623  iterations 34
      x^2 - 2  secant  root 1.4142135624  iterations  7
 x^3 - 2x - 5  bisect  root 2.0945514816  iterations 34
 x^3 - 2x - 5  secant  root 2.0945514815  iterations  7
```

A tenth function or a sixth solver costs one dictionary entry. In the loop, `solve` is a function taken out of a dict and called. That is the forcing problem. Three more pressures of the same kind drive the module: **doing something around every call** (time, count, retry, cache) without editing each function, which needs a function that takes a function and returns a wrapped one (§3.5, §3.6); **configuring a function once and calling it later** ("bisection with `tol=1e-6`" should itself be a function of $(f, a, b)$), which needs a function that remembers values (§3.4); and **splitting code across files** that can be found by name and imported without side effects (§3.8), with the standard library (§3.9) as code you never have to write. Recursion (§3.7) is here because it is where "what is the call stack?" first has to be answered.

## 2. Intuition

**A function is an object that a name points to.** `def` builds a function object and binds a name to it. Another name can point to it, a dict can hold it, a call can return it.

```
  def square(x): ...        square ───►  ⟨ function object: code x*x, __defaults__, __closure__, __dict__ ⟩
  g = square                g      ───►  (the same object)
  ops = {"sq": square}      ops["sq"] ─►  (the same object)
```

**A closure is a function with a backpack.** A function built inside another and using one of its variables keeps that variable in a small box, a *cell*. The outer call ends and its frame disappears; the cell lives as long as the inner function.

**A decorator takes a function and returns a function** (usually a wrapper): `@timer` above a `def` means "build the function, pass it through `timer`, bind the name to the result".

**Calling a function pushes a frame; returning pops it.** A frame holds the call's own parameters and locals and where to resume in the caller. Recursion stacks frames of one function with different `n`:

```
        top ►  fib(1)    n = 1
               fib(2)    n = 2   (waiting for fib(1))
               fib(3)    n = 3   (waiting for fib(2))
               fib(4)    n = 4   (waiting for fib(3))
               <module>
```

**A module is a file that is run once on import and becomes an object** whose attributes are the file's top-level names. Python keeps every loaded module in a dict, `sys.modules`; a second `import` finds it there. A package is a folder of modules with an `__init__.py`.

Guess before reading on. (1) What do the functions in `[lambda: i for i in range(3)]` return when called? (§3.4) (2) Naive recursive `fib(30)` makes roughly how many calls: 30, 900, or millions? (§3.7)

## 3. Formal development

### 3.1 Functions are objects

`def` is a statement that runs: it builds a function object and binds a name. The object has a type, a name, a docstring, defaults, compiled code, and a `__dict__` to which you can add attributes.

```pycon
>>> def area(r, pi=3.14159):
...     """Area of a circle of radius r."""
...     return pi * r * r
...
>>> type(area), area.__name__, area.__doc__
(<class 'function'>, 'area', 'Area of a circle of radius r.')
>>> area.__defaults__, area.__code__.co_varnames
((3.14159,), ('r', 'pi'))
>>> area.units = "m^2"; area.__dict__
{'units': 'm^2'}
>>> g = area
>>> g(2), g is area, callable(area), callable(3)
(12.56636, True, True, False)
```

▶ **Try it first.** Predict the printed values.

```python
def f():
    return 1

g = f

def f():
    return 2

print(g(), f())
del f
print(g())
```
```text
1 2
1
```

■ **What just happened.** A name is a label. The second `def` made a new object and moved the label `f` to it; `g` still labels the first. `del f` removes a label, not an object. That is why a function can sit in a list or dict and be called later: `table[key](args)` is a lookup followed by a call. The most useful consequence is the **dispatch table**, a dict from keys to functions that replaces a long `if/elif` chain; Example 4.1 builds one.

### 3.2 Higher-order functions

A **higher-order function** takes a function as an argument, returns one, or both. You have used one since PY-2: `sorted(xs, key=f)` calls `f` once per element and sorts by the results (`min` and `max` take `key` too). `map(f, xs)` and `filter(pred, xs)` are lazy iterators. `functools.reduce(f, xs)` folds a sequence: `reduce(f, [a, b, c])` is `f(f(a, b), c)`.

▶ **Try it first.** Predict all four lines of output.

```python
from functools import reduce

nums = [1, 2, 3, 4]
squares = map(lambda x: x * x, nums)
print(list(squares), list(squares))
print(reduce(lambda acc, x: acc * 10 + x, nums))
words = ["pear", "Fig", "apple"]
print(sorted(words), sorted(words, key=str.lower), sorted(words, key=len))
```
```text
[1, 4, 9, 16] []
1234
['Fig', 'apple', 'pear'] ['apple', 'Fig', 'pear'] ['Fig', 'pear', 'apple']
```

■ **What just happened.** The second `list(squares)` is empty: the `map` was used up. `reduce` computed $((1\cdot10+2)\cdot10+3)\cdot10+4 = 1234$. Plain `sorted` puts uppercase first (code point order); `key=str.lower` and `key=len` sort by derived values used only for comparison. `str.lower` is a method looked up on the class: a function value you can pass.

**`lambda`** builds a function from one expression (no statements, no `return`; its name is `'<lambda>'` in tracebacks). Use it for short `key=` functions and callbacks; if you are about to write `name = lambda ...`, write a `def` (PEP 8). The `operator` module has named versions (`itemgetter(1)`, `attrgetter("x")`, `add`). **`map`/`filter` against comprehensions:** `[f(x) for x in xs]` and `[x for x in xs if p(x)]` are usually clearer and give a reusable list; use `map(float, tokens)` when the function already exists. `reduce` is rarely clearest (`sum`, `math.prod`, `max`, `any`, `all` and `str.join` cover most folds); Horner's rule from PY-1, `reduce(lambda acc, c: acc * x + c, coeffs)`, is a genuine use.

▶ **Build it yourself (4): `compose` and `pipeline`.** Write `compose(f, g)` returning the function $h = f \circ g$, $h(x) = f(g(x))$, then `pipeline(*funcs)`, which applies functions left to right. Predict first: does `compose` return a value or a function?

1. A **function**: $x$ does not exist yet. Define an inner `h(*args, **kwargs)` that returns `f(g(*args, **kwargs))` and return `h`. It uses `f` and `g` from the enclosing call: a closure (§3.4).
2. `pipeline` keeps the tuple `funcs` in the closure and feeds the value through them in a loop.
3. The same fold is a `reduce`: `reduce(lambda acc, f: f(acc), funcs, x)`.

```python
import math
from functools import reduce

def compose(f, g):
    def h(*args, **kwargs):
        return f(g(*args, **kwargs))
    return h

def pipeline(*funcs):
    return lambda x: reduce(lambda acc, f: f(acc), funcs, x)

to_kelvin = pipeline(str.strip, float, lambda c: c + 273.15)
print(compose(math.sin, math.sqrt)(2.0), to_kelvin("  21.5 "), compose(len, str.split)("a b c d"))
```
```text
0.9877659459927356 294.65 4
```

■ **What just happened.** Both take functions and return a new one, which can go to `sorted(key=...)`, into a table, or into `compose` again. `numtools/functional.py` writes `pipeline` with a loop and lets the first function take several arguments.

### 3.3 Arguments: `*args`, `**kwargs`, keyword-only, positional-only

A `def` header can contain five kinds of parameter, in this order:

```
def f(pos_only, /, normal, *args, kw_only, **kwargs):
```

- Before `/`: **positional-only**; the caller cannot use the name, so you may rename it later (`inspect.signature(divmod)` shows `(x, y, /)`).
- Between `/` and `*`: **normal** (positional or keyword).
- `*args` collects extra positional arguments into a **tuple**; a bare `*` only ends the positional list.
- After `*`: **keyword-only**, so call sites read well: `bisect(f, 0, 1, tol=1e-9)`.
- `**kwargs` collects extra keyword arguments into a **dict**.

At the call site, `f(*seq)` spreads a sequence into positional arguments and `f(**mapping)` a dict into keyword arguments. Wrappers use both to forward what they receive (§3.5).

▶ **Try it first.** Predict each result; two lines raise an error, so predict the messages.

```pycon
>>> def f(a, b=2, /, c=3, *args, d, e=5, **kw):
...     return a, b, c, args, d, e, kw
...
>>> f(1, d=4)
(1, 2, 3, (), 4, 5, {})
>>> f(1, 2, 3, 4, 5, d=6, z=7)
(1, 2, 3, (4, 5), 6, 5, {'z': 7})
>>> f(*[1, 2], **{"d": 0})
(1, 2, 3, (), 0, 5, {})
>>> f(a=1, d=2)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: f() missing 1 required positional argument: 'a'
>>> f(1)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: f() missing 1 required keyword-only argument: 'd'
```

■ **What just happened.** `c` is the only parameter that accepts either style; surplus positionals land in `args`, unknown keywords in `kw`. `a` is positional-only, so the keyword `a=1` does not bind it: with `**kw` present it is swallowed by `kw`, and `a` is reported missing. `d` has no default, so omitting it is an error; Python allows a keyword-only parameter without a default after one with a default.

**Defaults are evaluated once**, when the `def` runs, and stored in `f.__defaults__` (see `area` above and PY-1 Example 4.6). So `def f(x, lst=[])` is still wrong: every call that omits `lst` shares one list. The repair is a sentinel: `def f(x, lst=None)` with `if lst is None: lst = []` in the body (or a private `_MISSING = object()` when `None` is a legitimate argument, compared with `is`).

### 3.4 Closures

A **free variable** of a function is a name it reads but does not bind itself, bound in an enclosing *function* scope (the E of LEGB, PY-1). A **closure** is a function object plus the cells holding its free variables: `f.__code__.co_freevars` lists the names, `f.__closure__` holds the cells, and `cell.cell_contents` reads one.

```pycon
>>> def make_adder(n):
...     def add(x):
...         return x + n
...     return add
...
>>> add5, add7 = make_adder(5), make_adder(7)
>>> add5(10), add7(10)
(15, 17)
>>> add5.__code__.co_freevars, add5.__closure__[0].cell_contents
(('n',), 5)
>>> add5.__closure__[0] is add7.__closure__[0]
False
```

Each call of `make_adder` made its own cell. Two rules (Language Reference, Execution model, "Resolution of names") matter in practice. **(i)** A free variable is looked up in the enclosing scope *when the inner function runs*: a closure captures the *variable*, not its value at creation. **(ii)** Assigning to a name inside a function makes it local unless declared `nonlocal` (or `global`); reading a free variable needs no declaration, rebinding does (Example 4.2).

▶ **Try it first.** This is guess (1) of §2. Predict the list, and whether the lambdas share a cell.

```pycon
>>> fs = [lambda: i for i in range(3)]
>>> [f() for f in fs]
[2, 2, 2]
>>> fs[0].__closure__[0] is fs[2].__closure__[0]
True
```

■ **What just happened.** The answer is `[2, 2, 2]`. There is **one** variable `i`, all three lambdas hold its cell, and when they are called `i` is 2. This is the **late-binding trap**; it bites in any loop that builds callbacks.

▶ **Build it yourself (3): predict, run, fix.** Build $x^1, x^2, x^3$ as functions so that evaluating at $x = 2$ gives `[2, 4, 8]`. Predict what the plain `lambda` gives, then fix it three ways.

1. **Plain closure** (broken): `lambda x: x ** k` in a loop over `k`.
2. **Default argument:** `lambda x, k=k: x ** k`. The default is evaluated when the `lambda` expression runs, once per iteration, and stored in that function's `__defaults__`. Weakness: a caller can overwrite it.
3. **`functools.partial`:** `partial(power, k=k)` stores the *value* of `k` now, in a callable that remembers the function and its fixed arguments.
4. **Factory function:** `make_power(k)`; each call has a fresh scope, hence a fresh cell (`make_adder` again).

```python
from functools import partial

def power(x, k):
    return x ** k

def make_power(k):
    return lambda x: x ** k

variants = {"late": [lambda x: x ** k for k in (1, 2, 3)],
            "default": [lambda x, k=k: x ** k for k in (1, 2, 3)],
            "partial": [partial(power, k=k) for k in (1, 2, 3)],
            "factory": [make_power(k) for k in (1, 2, 3)]}
for name, fs in variants.items():
    print(f"{name:>8}", [f(2) for f in fs])
```
```text
    late [8, 8, 8]
 default [2, 4, 8]
 partial [2, 4, 8]
 factory [2, 4, 8]
```

■ **What just happened.** Only the plain closure fails. Prefer `partial` or a factory in new code; you will meet the default-argument form in old code. `partial` also configures solvers: `fast = partial(bisect, tol=1e-6)` is a function of $(f, a, b)$ with the tolerance baked in.

### 3.5 Decorators

▶ **Build it yourself (1): from a pasted timing block to `@timer`.** The first thing everyone writes:

```python
import time

def slow_sum(n):
    """Sum of squares below n."""
    return sum(i * i for i in range(n))

start = time.perf_counter()
result = slow_sum(10**6)
print(f"slow_sum took {time.perf_counter() - start:.3f} s, result {result}")
```
```text
slow_sum took 0.088 s, result 333332833333500000
```

Five steps turn that into a tool.

1. **The pain.** Three timing lines pasted around every call of every function.
2. **Take the function as an argument.** `timed(func, *args)` holds the one copy of the logic, but every call site must change to `timed(slow_sum, n)`.
3. **Return a function with the same call signature.** `timer(func)` returns a `wrapper(*args, **kwargs)` that times and calls `func`. Rebind the name and call sites stay as they are: `slow_sum = timer(slow_sum)`.
4. **`@` is sugar for that rebinding.** `@timer` above `def slow_sum` means `slow_sum = timer(slow_sum)` after the `def`.
5. **Repair what the wrapper broke.** The name now points to `wrapper`. Predict `slow_sum.__name__` and `slow_sum.__doc__` before running.

```python
import time

def timer(func):
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        print(f"{func.__name__} took {time.perf_counter() - start:.3f} s")
        return result
    return wrapper

@timer
def slow_sum(n):
    """Sum of squares below n."""
    return sum(i * i for i in range(n))

print(slow_sum(10**6))
print(slow_sum.__name__, repr(slow_sum.__doc__))
```
```text
slow_sum took 0.075 s
333332833333500000
wrapper None
```

`functools.wraps(func)` copies `__name__`, `__qualname__`, `__doc__`, `__module__` and `__dict__` from `func` onto the wrapper and sets `wrapper.__wrapped__ = func`, so `help` and `inspect.signature` see through it. The final `timer` (as in `numtools/decorators.py`) stores durations in a list on the wrapper (functions have a `__dict__`) and uses `try/finally` so that a call that raises is still recorded:

```python
import functools
import time

def timer(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        try:
            return func(*args, **kwargs)
        finally:
            wrapper.timings.append(time.perf_counter() - start)
    wrapper.timings = []
    return wrapper

@timer
def slow_sum(n):
    """Sum of squares below n."""
    return sum(i * i for i in range(n))

for _ in range(3):
    slow_sum(10**5)
print(slow_sum.__name__, repr(slow_sum.__doc__), len(slow_sum.timings), slow_sum.__wrapped__.__name__)
```
```text
slow_sum 'Sum of squares below n.' 3 slow_sum
```

■ **What just happened.** The wrapper pattern is always: take `*args, **kwargs`, do something before, call the original, do something after, return its result. New decorators forget one of three things: the `return`, the argument forwarding, or `wraps`. **`wraps` is not optional:** without it tracebacks, `help()` and test reports all say `wrapper`. Stacked decorators apply bottom-up (`@a @b def f` is `f = a(b(f))`), and the outermost wrapper runs first.

**Decorators with arguments** have three levels, because `@retry(3)` is evaluated first and must *return* a decorator. Outer level: configuration. Middle: the function. Inner: each call.

```python
import functools

def retry(n):
    """Call the function up to n times in total; re-raise the last ConnectionError."""
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            for attempt in range(1, n + 1):
                try:
                    return func(*args, **kwargs)
                except ConnectionError:
                    if attempt == n:
                        raise
        return wrapper
    return decorator

attempts = []

@retry(3)
def flaky():
    attempts.append(1)
    if len(attempts) < 3:
        raise ConnectionError(f"attempt {len(attempts)} failed")
    return "ok"

print(flaky(), len(attempts))
attempts.clear()
try:
    retry(2)(flaky.__wrapped__)()           # the decorator syntax written out, on the original
except ConnectionError as err:
    print("gave up:", err)
```
```text
ok 3
gave up: attempt 2 failed
```

`numtools.retry(n, exceptions=(Exception,), delay=0.0)` adds the exception types to catch (never a bare `except:`) and an optional pause. A decorator that only *registers* a function and returns it unchanged, as in Example 4.1, needs no wrapper.

### 3.6 Caching: `memoize`, `lru_cache`, and hashability

A **pure** function (same arguments, same result, no side effects) never needs to compute an answer twice.

▶ **Build it yourself (2): `@memoize`.** Predict first: to get `fib(25)` by the textbook recursion, how many times does the body run, and how many with a cache?

1. Keep a dict `cache` in the closure. The key is the tuple `args` plus the sorted `kwargs` items; all of it must be **hashable** (PY-2 §3.4).
2. If the key is missing, call the function and store the result; return `cache[key]`.
3. Apply it to a recursive function. The body calls `fib_memo` by its global name, which now points to the wrapper, so every recursive call goes through the cache.
4. Replace it with `functools.cache`, the C-implemented standard version.

```python
import functools

calls = {"plain": 0, "memo": 0}

def memoize(func):
    cache = {}
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        key = (args, tuple(sorted(kwargs.items())))
        if key not in cache:
            cache[key] = func(*args, **kwargs)
        return cache[key]
    wrapper.cache = cache
    return wrapper

def fib_plain(n):
    calls["plain"] += 1
    return n if n < 2 else fib_plain(n - 1) + fib_plain(n - 2)

@memoize
def fib_memo(n):
    calls["memo"] += 1
    return n if n < 2 else fib_memo(n - 1) + fib_memo(n - 2)

print(fib_plain(25), fib_memo(25), calls, len(fib_memo.cache))
```
```text
75025 75025 {'plain': 242785, 'memo': 26} 26
```

■ **What just happened.** The plain version ran its body 242,785 times to produce a number that depends on 26 distinct values ($F_0, \dots, F_{25}$); the memoized one ran it 26 times, once per distinct argument. §3.7 derives the first count.

`functools.cache` is `lru_cache(maxsize=None)`; `lru_cache(maxsize=N)` keeps at most $N$ entries and discards the **l**east **r**ecently **u**sed. Both add `cache_info()` and `cache_clear()`. Two requirements and a trap:

```pycon
>>> import functools
>>> @functools.lru_cache(maxsize=None)
... def total(xs):
...     return sum(xs)
...
>>> total((1, 2, 3)), total((1, 2, 3)), total.cache_info()
(6, 6, CacheInfo(hits=1, misses=1, maxsize=None, currsize=1))
>>> total([1, 2, 3])
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: unhashable type: 'list'
>>> @functools.cache
... def zeros(n):
...     return [0] * n
...
>>> a = zeros(3)
>>> a.append(9)
>>> zeros(3)
[0, 0, 0, 9]
```

**Arguments must be hashable** (a list fails; pass a tuple). **The function must be pure,** or the cache serves stale answers. **The trap:** the cache returns the *stored object*, not a copy, so mutating a returned list changes every later answer (PY-1's aliasing). Keys compare by `==`, so `f(1, 2)` and `f(1.0, 2.0)` share an entry unless `typed=True`.

### 3.7 Recursion

A **recursive** function calls itself. It is correct when it has (1) a **base case** answered directly, (2) a **recursive case** that reduces the problem to *smaller* instances, and (3) a measure of size that strictly decreases. Drop (1) or (3) and you get a `RecursionError`.

▶ **Build it yourself (5): the call stack for `fib(4)`.** Before running anything, draw the tree of calls: `fib(4)` calls `fib(3)` and then `fib(2)`, down to `fib(1)` and `fib(0)`. How many calls? Then check with a decorator that prints each call, indented by depth.

1. `trace(func)` prints `name(arg)` at the current depth, increments the depth, calls the function, decrements it, and prints the result. The depth is rebound, so it needs `nonlocal` (§3.4).
2. Apply it to `fib`. The left operand `fib(n - 1)` is evaluated *completely* before `fib(n - 2)` starts; every `->` line is a frame being popped.
3. Apply it again with a cache underneath (`@trace` above `@functools.cache`).

```python
import functools

def trace(func):
    depth = 0
    @functools.wraps(func)
    def wrapper(n):
        nonlocal depth
        print(f"{'  ' * depth}{func.__name__}({n})")
        depth += 1
        try:
            result = func(n)
        finally:
            depth -= 1
        print(f"{'  ' * depth}-> {result}")
        return result
    return wrapper

@trace
def fib(n):
    return n if n < 2 else fib(n - 1) + fib(n - 2)

fib(4)
print("--- with a cache underneath ---")

@trace
@functools.cache
def fib(n):
    return n if n < 2 else fib(n - 1) + fib(n - 2)

fib(4)
```
```text
fib(4)
  fib(3)
    fib(2)
      fib(1)
      -> 1
      fib(0)
      -> 0
    -> 1
    fib(1)
    -> 1
  -> 2
  fib(2)
    fib(1)
    -> 1
    fib(0)
    -> 0
  -> 1
-> 3
--- with a cache underneath ---
fib(4)
  fib(3)
    fib(2)
      fib(1)
      -> 1
      fib(0)
      -> 0
    -> 1
    fib(1)
    -> 1
  -> 2
  fib(2)
  -> 1
-> 3
```

■ **What just happened.** Naive `fib(4)` made 9 calls (`fib(2)` twice, `fib(1)` three times). At the first `-> 1` the stack holds four frames, the picture of §2: depth is $n$ although the calls are far more numerous. With the cache, the repeated `fib(2)` and `fib(1)` are single lines with no children: 7 calls, 5 of which ran the body (the distinct arguments $0, \dots, 4$).

**Counting calls exactly.** Let $C(n)$ be the number of calls (including the first) of naive `fib(n)`: $C(0) = C(1) = 1$ and $C(n) = C(n-1) + C(n-2) + 1$.

*Claim:*
$$C(n) = 2F_{n+1} - 1. \tag{3.1}$$
*Proof by induction.* $C(0) = 1 = 2F_1 - 1$ and $C(1) = 1 = 2F_2 - 1$. If the claim holds for $n - 1$ and $n - 2$, then $C(n) = (2F_n - 1) + (2F_{n-1} - 1) + 1 = 2F_{n+1} - 1$. ∎

By (3.1) and $F_n \sim \varphi^n / \sqrt5$, calls grow by $\varphi \approx 1.618$ per unit of $n$, so $C(n+5)/C(n) \to \varphi^5 \approx 11.09$. With a cache the body runs once per distinct argument ($n + 1$ times) and the calls total $2n - 1$: exponential becomes linear.

```python
def fib_calls(n):
    return 1 if n < 2 else fib_calls(n - 1) + fib_calls(n - 2) + 1

counts = [fib_calls(n) for n in (20, 25, 30)]
print(counts, [2 * f - 1 for f in (10946, 121393, 1346269)], round(counts[1] / counts[0], 2), "vs phi^5 =", round(((1 + 5**0.5) / 2) ** 5, 2))
```
```text
[21891, 242785, 2692537] [21891, 242785, 2692537] 11.09 vs phi^5 = 11.09
```

(The middle list is $2F_{21} - 1$, $2F_{26} - 1$, $2F_{31} - 1$.) About 2.7 million calls for `fib(30)` answers guess (2) of §2. `project.py` Part 4 times the recursion and asserts that the time ratio between $n = 25$ and $n = 20$ is within a factor of 2 of $\varphi^5$; it prints a value between 10 and 12 on the machine used for these notes.

**The recursion limit.** Every call holds a frame in memory, so CPython caps the depth, by default at 1000.

```python
import sys

def depth(n):
    return 0 if n == 0 else 1 + depth(n - 1)

print(sys.getrecursionlimit(), depth(900))
try:
    depth(5000)
except RecursionError as err:
    print(type(err).__name__ + ":", err)
```
```text
1000 900
RecursionError: maximum recursion depth exceeded
```

CPython does not optimize tail calls: every level costs a frame. `sys.setrecursionlimit` raises the cap but trades a clean error for memory use and, on some paths through C code, a crashed interpreter. **Python is not the language for deep recursion**; a depth proportional to the input size is a bug waiting for a large input. A cache does not change the depth (Example 4.4).

▶ **Build it yourself (6): binary search, invariant first.** Given a sorted list `a` and a value `x`, find the smallest index `i` with `a[i] >= x` (or `len(a)`). Before writing code, write what stays true: **everything left of `lo` is `< x`; everything from `hi` on is `>= x`.**

1. *Start.* `lo, hi = 0, len(a)`: both parts are empty, so the statement holds.
2. *Step.* While `lo < hi`, let `mid = (lo + hi) // 2`, so `lo <= mid < hi`. If `a[mid] < x`, sortedness gives `a[j] < x` for all `j <= mid`: set `lo = mid + 1`. Otherwise `a[j] >= x` for all `j >= mid`: set `hi = mid`. Both branches preserve the invariant.
3. *Progress.* The size $s$ = `hi - lo` becomes at most $\lfloor s/2 \rfloor$, so the loop ends after at most $\lfloor \log_2 n \rfloor + 1$ steps.
4. *Done.* At `lo == hi` everything before `lo` is `< x` and everything from `lo` on is `>= x`: `lo` is the answer.
5. *Recursive version.* The invariant is the contract of `search(lo, hi)`: base case `lo == hi`, recursive case a call on the remaining half. Depth $O(\log n)$.

```python
import math
from bisect import bisect_left

def lower_bound(a, x):
    lo, hi, steps = 0, len(a), 0
    while lo < hi:                      # invariant: a[:lo] < x <= a[hi:]
        steps += 1
        mid = (lo + hi) // 2
        if a[mid] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo, steps

def lower_bound_rec(a, x, lo=0, hi=None):
    hi = len(a) if hi is None else hi
    if lo == hi:
        return lo
    mid = (lo + hi) // 2
    if a[mid] < x:
        return lower_bound_rec(a, x, mid + 1, hi)
    return lower_bound_rec(a, x, lo, mid)

a = [1, 3, 3, 5, 8, 8, 8, 13]
for x in (0, 3, 4, 14):
    i, steps = lower_bound(a, x)
    assert i == lower_bound_rec(a, x) == bisect_left(a, x)
    print(f"x = {x:>2}: index {i}, {steps} steps")
print("bound floor(log2 n) + 1 =", math.floor(math.log2(len(a))) + 1)
```
```text
x =  0: index 0, 4 steps
x =  3: index 1, 4 steps
x =  4: index 3, 3 steps
x = 14: index 8, 3 steps
bound floor(log2 n) + 1 = 4
```

■ **What just happened.** The loop and the recursion are one algorithm; the invariant makes either correct. `project.py` checks both against `bisect_left` on thousands of random lists. In real code use `bisect.bisect_left`.

**Divide and conquer: merge sort.** Sort the two halves, then merge. A merge costs at most one comparison per output item, there are $\lceil \log_2 n \rceil$ levels of halving, and each level merges $n$ items in total: $T(n) = 2T(n/2) + n$, at most $n \lceil \log_2 n \rceil$ comparisons. Merging with `<=` keeps equal items in input order (a **stable** sort).

```python
def merge_sort(xs):
    if len(xs) <= 1:
        return list(xs)
    mid = len(xs) // 2
    left, right = merge_sort(xs[:mid]), merge_sort(xs[mid:])
    out, i, j = [], 0, 0
    while i < len(left) and j < len(right):
        if left[i] <= right[j]:
            out.append(left[i])
            i += 1
        else:
            out.append(right[j])
            j += 1
    return out + left[i:] + right[j:]

print(merge_sort([5, 2, 9, 1, 5, 6, 0]), merge_sort([]))
```
```text
[0, 1, 2, 5, 5, 6, 9] []
```

`project.py` checks it against `sorted` on $10^4$ items; the built-in (Timsort, in C) wins by a large constant, but the point is the shape of the recursion.

**The cofactor determinant** is the standard correct-but-hopeless recursion: LA-6 eq. (3.12) expands $\det A$ into $n$ determinants of size $n - 1$, costing $T(n) = nT(n-1) + n$ multiplications, which LA-6 §3.5 solves as $n!\sum_{k=1}^{n-1} 1/k!$ (Example 4.7). **Mutual recursion** (`is_even` calling `is_odd` calling `is_even`) is rare and has all the limits above. Recursion fits recursive data (nested lists, trees, parsers) and halving algorithms, where the depth is small; for a "loop in disguise" such as `xs[0] + sum_rec(xs[1:])`, use a loop. Any recursion becomes a loop over an explicit stack (a list you push to and pop from) when depth is the problem; "Your turn" asks for one.

**Generators, revisited as resumable functions.** A generator function's *frame survives between calls*: at `yield`, the frame (locals and position) is kept on the heap instead of being popped, and `next` resumes it. `yield from` delegates to another iterable, which makes recursive traversal short:

```pycon
>>> def count_to(n):
...     i = 0
...     while i < n:
...         i += 1
...         yield i
...
>>> g = count_to(3)
>>> next(g), next(g)
(1, 2)
>>> g.gi_frame.f_locals
{'n': 3, 'i': 2}
>>> def flatten(x):
...     for item in x:
...         if isinstance(item, list):
...             yield from flatten(item)
...         else:
...             yield item
...
>>> list(flatten([1, [2, [3, [4, 5]], 6], [], 7]))
[1, 2, 3, 4, 5, 6, 7]
```

The frozen frame shows `i` at 2, waiting for the next `next`. `flatten` is recursion on nested data, yielding lazily with no intermediate lists.

### 3.8 Modules, packages, and the import system

**What `import m` does**, in order:

1. Look for `"m"` in `sys.modules`, the dict of every module loaded so far. If found, bind the name `m` and stop.
2. Otherwise find the file: the finders walk the list `sys.path` (directories, in order) for `m.py` or a package folder `m/`.
3. Create an empty module object and register it in `sys.modules["m"]` **before** running anything.
4. Execute the file's top-level code, top to bottom, in the module's namespace. Each top-level `def`, `class` and assignment becomes an attribute of the module.
5. Bind the name `m` in the importer.

If step 4 raises, the half-built entry is removed. Registering first is why a circular import sees a partly built module instead of looping. Because of step 1, **a module's top-level code runs once per process, however many files import it.**

▶ **Try it first.** Predict every line `python3 main.py` prints.

```python
# file: chatty.py
print("chatty.py is running; __name__ =", __name__)
VALUE = 42
```

```python
# file: main.py
import sys
import chatty
import chatty
from chatty import VALUE
import chatty as c

print(VALUE, c is chatty, "chatty" in sys.modules, sys.modules["chatty"] is chatty)
print(sorted(name for name in vars(chatty) if not name.startswith("__")))
```
```text
chatty.py is running; __name__ = chatty
42 True True True
['VALUE']
```

■ **What just happened.** The file ran once, on the first `import`; the later forms found the module in `sys.modules`. The module is an ordinary object whose attributes are its top-level names.

**`from m import *` is banned in this course.** It dumps every public name of `m` into your namespace, so a reader cannot tell where a name came from, and it silently replaces names you already had, built-ins included:

```pycon
>>> from os import *
>>> open("notes.txt")
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: open() missing required argument 'flags' (pos 2)
```

The built-in `open` has been replaced by `os.open`, and the error never mentions the import. A module's `__all__` lists what `import *` exports; its better job is documenting the public API.

**`if __name__ == "__main__":`** Python sets `__name__` to `"__main__"` in the file it was asked to run and to the module's own name when it is imported. The guard lets one file be a library and a program:

```python
# file: tool.py
def double(x):
    return 2 * x

if __name__ == "__main__":
    print("running as a script:", double(21))
```

```console
$ python3 tool.py
running as a script: 42
$ python3 -c "import tool; print(tool.double(4), tool.__name__)"
8 tool
$ python3 -m tool
running as a script: 42
```

`python -m tool` runs a module *by name*, searching `sys.path` like an import; prefer it for anything inside a package.

**`sys.path`: "it works in the REPL but not from the script".** `sys.path[0]` depends on how Python was started: for `python script.py` it is **the directory containing the script**; for `python -m pkg.mod`, `python -c ...` and the REPL it is **the current working directory**. "Module not found when the file is right there" is almost always about where Python was run from. (The folders below have no `__init__.py`; such folders still import, as *namespace packages*, which is why `-m scripts.run` works.)

```python
# file: proj/tools/helpers.py
def hello():
    return "hello from tools.helpers"
```

```python
# file: proj/scripts/run.py
from tools.helpers import hello

print(hello())
```

```console
$ cd proj && python3 scripts/run.py
Traceback (most recent call last):
  File "proj/scripts/run.py", line 2, in <module>
    from tools.helpers import hello
ModuleNotFoundError: No module named 'tools'
$ cd proj && python3 -m scripts.run
hello from tools.helpers
```

The first command puts `proj/scripts/` first on `sys.path`, and `tools` is not in it; the second puts `proj/` first. Fixes, in order of preference: run modules with `python -m` from the project root; install your package so it is importable everywhere (`pip install -e .`, PY-10); set `PYTHONPATH` for a one-off. Do not edit `sys.path` inside library code.

**Packages.** A folder with an `__init__.py` is a *package*: importing it runs `__init__.py`, and the names defined there are the package's attributes. Modules inside a package import each other with **relative imports**: `from .roots import bisect` means "the module `roots` next to me", resolved against the package name. A package may contain subpackages (folders with their own `__init__.py`), imported as `numtools.sub.mod`; inside one, `from ..roots import bisect` climbs one level.

▶ **Build it yourself (7): the `numtools` package.** The finished package is in this folder. To build your own: move PY-1's `bisect`, `newton`, `secant` to `roots.py`, PY-2's sparse polynomial to `polynomials.py`, your decorators to `decorators.py`, then:

1. Add `__init__.py` that re-exports the public names (`from .roots import bisect, newton, secant`, and so on) and lists them in `__all__`; users write `from numtools import bisect`.
2. Use relative imports between the package's own modules.
3. Add `__main__.py`, which `python -m numtools` runs with the package as context (below).
4. Run from the folder that **contains** `numtools/`.

`numtools/__main__.py`:

```python
"""``python -m numtools`` runs this file with numtools as its package."""
import sys

from .cli import main

sys.exit(main())
```

```console
$ ls numtools
__init__.py
__main__.py
cli.py
decorators.py
functional.py
polynomials.py
roots.py
$ python3 -m numtools --version
numtools 0.1
$ python3 -c "import numtools; print(numtools.__all__)"
['bisect', 'newton', 'secant', 'poly_add', 'poly_mul', 'poly_eval', 'poly_deriv', 'timer', 'memoize', 'retry', 'compose', 'pipeline']
$ python3 numtools/cli.py
Traceback (most recent call last):
  File "numtools/cli.py", line 9, in <module>
    from . import __version__
ImportError: attempted relative import with no known parent package
```

The last command shows the other half of relative imports: run directly, `cli.py` is a top-level script with no package, so `from . import __version__` has nothing to be relative to. `project.py` is a script *next to* `numtools/`, so Python puts this folder first on `sys.path` when you run it; in an interactive session started elsewhere the import can fail, so `project.py` inserts its own folder explicitly.

### 3.9 The standard library tour

One example each; the documentation has the rest.

**Numbers.** `math`, `cmath` (complex), `fractions` (exact rationals), `decimal` (exact decimal arithmetic at a chosen precision), `statistics` (`stdev` is the sample version with $n-1$, `pstdev` the population one).

```pycon
>>> import math, cmath, statistics
>>> from fractions import Fraction
>>> from decimal import Decimal
>>> math.hypot(3, 4), math.comb(5, 2), math.isclose(0.1 + 0.2, 0.3), cmath.sqrt(-4)
(5.0, 10, True, 2j)
>>> Fraction(1, 3) + Fraction(1, 6), Decimal("0.1") + Decimal("0.2") == Decimal("0.3")
(Fraction(1, 2), True)
>>> data = [2, 4, 4, 4, 5, 5, 7, 9]
>>> statistics.mean(data), statistics.median(data), statistics.pstdev(data), round(statistics.stdev(data), 4)
(5, 4.5, 2.0, 2.1381)
```

**`random`.** The module-level functions share one hidden global generator; a `random.Random(seed)` instance is an independent generator with the same methods (`randint`, `sample`, `gauss`), and two with the same seed give the same stream. Pass an instance (or a NumPy `Generator`, PY-6) as an argument. The footgun is a library that seeds the *global* generator at import, because step 4 of the import runs top-level code in the importer's process:

```python
# file: noisy.py
import random

random.seed(0)          # runs at import time: seeds the shared global generator
```

```python
import random

random.seed(123)        # you carefully seed your experiment...
import noisy            # ...and importing this module re-seeds the global generator with 0
print(random.random() == random.Random(0).random())
```
```text
True
```

`True`: your seed was overwritten. Seed in `main`, or create `rng = random.Random(seed)` and pass it down.

**`time`, `datetime`, `os`, `sys`, `textwrap`, `pprint`.** `time.perf_counter()` measures intervals (it is what `timer` uses); `time.time()` is wall-clock time, which can jump. `datetime` does calendar dates, ISO strings and `timedelta` arithmetic. `sys.argv` holds the command-line words, `sys.exit(status)` ends the program, `os.environ` is a dict of environment variables; for paths use `pathlib` (PY-2 §3.10), not `os.path`.

```pycon
>>> import os, textwrap, pprint
>>> from datetime import datetime, timedelta, timezone
>>> d = datetime.fromisoformat("2024-03-01T12:00:00")
>>> d + timedelta(days=30), (datetime(2024, 3, 1) - datetime(2024, 2, 1)).days
(datetime.datetime(2024, 3, 31, 12, 0), 29)
>>> d.replace(tzinfo=timezone.utc).timestamp(), os.environ.get("NO_SUCH_VARIABLE", "unset")
(1709294400.0, 'unset')
>>> print(textwrap.fill("Python is a language that lets you work quickly and integrate systems.", width=30))
Python is a language that lets
you work quickly and integrate
systems.
>>> pprint.pprint({"gamma": "a string", "alpha": list(range(10)), "beta": {"x": 1, "y": [1, 2, 3]}}, width=50)
{'alpha': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
 'beta': {'x': 1, 'y': [1, 2, 3]},
 'gamma': 'a string'}
```

▶ **Build it yourself (8): a command-line tool with `argparse`, in four steps.** A program receives its arguments as the list `sys.argv` of strings; `argparse` converts and validates them and writes the usage message.

1. **Bare `sys.argv`.** `a, b = float(sys.argv[1]), float(sys.argv[2])` works for good input. Wrong calls end in an `IndexError` or `ValueError` traceback pointing into your code, and there is no `--help`.
2. **Typed positionals and options.** `parser.add_argument("a", type=float)` converts and validates, and `-h/--help` appears for free. A name starting with `--` is optional and has a default; `action="store_true"` makes a flag; `choices=[...]` restricts values. (`%(default)s` in a help string is argparse's own placeholder for the default.)
3. **Make it testable.** Put the work in `main(argv=None)`, which calls `parser.parse_args(argv)` (`None` means `sys.argv[1:]`) and returns an exit status; end the file with `if __name__ == "__main__": sys.exit(main())`.
4. **Subcommands** come from `add_subparsers`; `set_defaults(func=...)` on each stores the function to call, a dispatch table (§3.1) that argparse builds for you.

`numtools/cli.py` has exactly this shape, and Example 4.6 runs it: predict the output of its error cases before you look.

**`re`.** Regular expressions describe text patterns. Write them as **raw strings** (`r"\d+"`) so Python does not interpret the backslashes first. `match` anchors at the start, `search` finds the first match anywhere, `fullmatch` requires the whole string, `findall` returns all matches, `sub` replaces. A match object is truthy and `None` falsy: test before calling `.group`. Parentheses make *groups*.

```pycon
>>> import re
>>> re.match(r"\d+", "42 apples"), re.match(r"\d+", "apples 42")
(<re.Match object; span=(0, 2), match='42'>, None)
>>> re.search(r"\d+", "apples 42").group(), re.findall(r"\d+", "a1 b22 c333")
('42', ['1', '22', '333'])
>>> re.sub(r"\s+", " ", "too   many    spaces"), re.match(r"(\w+)=(\d+)", "alpha=7").groups()
('too many spaces', ('alpha', '7'))
```

### 3.10 Environments: `venv` and `pip`

Different projects need different versions of the same package, and the system Python belongs to the operating system. A **virtual environment** is a private folder with its own interpreter link and its own `site-packages`; what you install while it is active goes there only.

```console
$ python3 -m venv .venv
$ .venv/bin/python -c "import sys; print(sys.prefix != sys.base_prefix)"
True
$ .venv/bin/python -m pip freeze
```

`sys.prefix != sys.base_prefix` tests "am I in a virtual environment". `pip freeze` lists what you installed, with exact versions, and prints nothing in a fresh one. The daily routine needs a network, so it is not run here:

```text
python3 -m venv .venv                       # once per project
source .venv/bin/activate                   # each new shell (Windows: .venv\Scripts\activate)
python -m pip install numpy scipy           # installs into .venv only
python -m pip freeze > requirements.txt     # record exact versions; commit this file
python -m pip install -r requirements.txt   # reproduce the environment elsewhere
```

Use `python -m pip` so that the `pip` belongs to the interpreter you mean. Do not `pip install` into the system Python (many systems refuse: PEP 668). Keep `.venv/` out of version control; commit `requirements.txt`.

## 4. Worked examples

### Example 4.1: A calculator as a dispatch table

The `if/elif` version has a branch per operator, and a new operator edits the function. A table separates data (which names exist) from logic (look up, convert, call). A decorator factory registers functions and returns them unchanged: a decorator with no wrapper.

```python
import math
import operator

OPS = {"+": operator.add, "-": operator.sub, "*": operator.mul, "/": operator.truediv,
       "**": operator.pow}

def op(name):
    """Decorator factory: register the function under `name` and return it unchanged."""
    def register(func):
        OPS[name] = func
        return func
    return register

@op("hypot")
def hypot(a, b):
    return math.hypot(a, b)

def calc(line):
    left, name, right = line.split()
    try:
        func = OPS[name]
    except KeyError:
        raise ValueError(f"unknown operator {name!r}; choose from {sorted(OPS)}") from None
    return func(float(left), float(right))

for line in ["3 + 4", "2 ** 10", "3 hypot 4", "1 ? 2"]:
    try:
        print(f"{line!r:>12} = {calc(line)}")
    except ValueError as err:
        print(f"{line!r:>12} -> ValueError: {err}")
print(hypot(6, 8))
```
```text
     '3 + 4' = 7.0
   '2 ** 10' = 1024.0
 '3 hypot 4' = 5.0
     '1 ? 2' -> ValueError: unknown operator '?'; choose from ['*', '**', '+', '-', '/', 'hypot']
10.0
```

`hypot(6, 8)` still works as a plain function because `register` returned it, and `from None` hides the chained `KeyError`. The same shape (name to function) is `set_defaults(func=...)` in Example 4.6.

### Example 4.2: A counter closure, with `nonlocal`, and a look inside

`step` rebinds `count`, a variable of the enclosing call, so it declares `nonlocal count`.

```python
def make_counter(start=0):
    count = start
    def step(inc=1):
        nonlocal count
        count += inc
        return count
    return step

a, b = make_counter(), make_counter(100)
print(a(), a(), a(10), b(), a())
cell = a.__closure__[0]
print(a.__code__.co_freevars, cell.cell_contents, b.__closure__[0].cell_contents, cell is b.__closure__[0])
```
```text
1 2 12 101 13
('count',) 13 101 False
```

`a` and `b` have independent cells, one per call of `make_counter`; `1 2 12 101 13` is `a` counting by 1, 1, 10, then `b` from 100, then `a` again. Without `nonlocal`, `count += 1` makes `count` a local of `step` that is read before it has a value, and the first call raises `UnboundLocalError`.

### Example 4.3: `@timer` on the PY-1 root finders

Wrapping does not edit the solvers: `timer(bisect)` is a new function that records each duration in `.timings`. Run each solver 200 times at tolerance $10^{-10}$ and report the median.

```python
import math
import statistics

from numtools import bisect, newton, secant, timer

CASES = [("cos x - x", lambda x: math.cos(x) - x, lambda x: -math.sin(x) - 1, 0.0, 1.0),
         ("x^2 - 2", lambda x: x * x - 2, lambda x: 2 * x, 1.0, 2.0)]
timed = {"bisect": timer(bisect), "newton": timer(newton), "secant": timer(secant)}
call = {"bisect": lambda s, f, df, a, b: s(f, a, b, tol=1e-10),
        "newton": lambda s, f, df, a, b: s(f, df, b, tol=1e-10),
        "secant": lambda s, f, df, a, b: s(f, a, b, tol=1e-10)}
print(f"{'function':>13} {'solver':>7} {'iterations':>10} {'median us':>10}")
for cname, f, df, a, b in CASES:
    for sname, solver in timed.items():
        solver.timings.clear()
        for _ in range(200):
            root, iterations = call[sname](solver, f, df, a, b)
        print(f"{cname:>13} {sname:>7} {iterations:>10} {statistics.median(solver.timings) * 1e6:>10.1f}")
print("bisection bound:", math.ceil(math.log2(1e10)))
```
```text
     function  solver iterations  median us
    cos x - x  bisect         34        9.4
    cos x - x  newton          5        1.7
    cos x - x  secant          6        1.8
      x^2 - 2  bisect         34        8.9
      x^2 - 2  newton          5        1.4
      x^2 - 2  secant          7        2.0
bisection bound: 34
```

Bisection needs 34 iterations, the bound $\lceil \log_2((b-a)/\mathrm{tol}) \rceil$ of PY-1 eq. (3.9); Newton and secant need a handful, and time follows iterations. Microseconds depend on your machine, the ordering does not; `project.py` Part 3 asserts the counts and the ordering.

### Example 4.4: `lru_cache` on Fibonacci, with `cache_info()`

```python
import functools

@functools.lru_cache(maxsize=None)
def fib(n):
    return n if n < 2 else fib(n - 1) + fib(n - 2)

print(fib(30), fib.cache_info())
print(fib(100), fib.cache_info())
fib.cache_clear()
try:
    fib(2000)
except RecursionError:
    print("cold fib(2000): RecursionError")
fib.cache_clear()
for k in range(2001):               # warm the cache bottom-up: recursion depth stays at 2
    fib(k)
print(len(str(fib(2000))), "digits;", fib.cache_info())
```
```text
832040 CacheInfo(hits=28, misses=31, maxsize=None, currsize=31)
354224848179261915075 CacheInfo(hits=99, misses=101, maxsize=None, currsize=101)
cold fib(2000): RecursionError
418 digits; CacheInfo(hits=3999, misses=2001, maxsize=None, currsize=2001)
```

For `fib(30)`: the 31 misses are the distinct arguments $0, \dots, 30$, and 28 hits are the repeated subcalls, $59 = 2 \cdot 30 - 1$ calls against 2,692,537 uncached (eq. (3.1)). A cold `fib(2000)` recurses 2000 deep and fails: a cache bounds the *work*, not the *depth*. Calling `fib(k)` for increasing `k` first means each call finds its two predecessors cached and the stack never grows.

### Example 4.5: A regular expression parses `2024-03-01T12:00:00 sensor_7 23.5`

Name the three fields with groups, compile once, and match the whole line.

```python
import re
from collections import defaultdict
from datetime import datetime
from statistics import mean

LINE = re.compile(r"""
    (?P<ts>\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})   # ISO timestamp
    \s+
    sensor_(?P<id>\d+)                            # sensor number
    \s+
    (?P<value>-?\d+(?:\.\d+)?)                    # signed decimal
    """, re.VERBOSE)

line = "2024-03-01T12:00:00 sensor_7 23.5"
m = LINE.fullmatch(line)
print(m.groupdict())
print((datetime.fromisoformat(m["ts"]), int(m["id"]), float(m["value"])))

log = """2024-03-01T12:00:00 sensor_7 23.5
2024-03-01T12:00:10 sensor_7 23.9
2024-03-01T12:00:10 sensor_12 -4.25
garbage line
2024-03-01T12:00:20 sensor_7 ERR
2024-03-01T12:00:30 sensor_12 -4.0""".splitlines()
values, rejected = defaultdict(list), []
for number, text in enumerate(log, start=1):
    match = LINE.fullmatch(text)
    if match is None:
        rejected.append(number)
    else:
        values[int(match["id"])].append(float(match["value"]))
print("rejected lines:", rejected)
print({sensor: mean(v) for sensor, v in values.items()})
```
```text
{'ts': '2024-03-01T12:00:00', 'id': '7', 'value': '23.5'}
(datetime.datetime(2024, 3, 1, 12, 0), 7, 23.5)
rejected lines: [4, 5]
{7: 23.7, 12: -4.125}
```

`(?P<name>...)` is a named group, `\d{4}` is exactly four digits, `(?:...)` groups without capturing, `-?` makes the sign optional, and `re.VERBOSE` allows layout and comments. `fullmatch` rejects `... ERR`. For CSV or JSON use `csv` and `json` (PY-2 §3.10) instead.

### Example 4.6: The command-line tool, with `--help`

From the module folder. The expression is compiled once into `f(x)`; `bisect` needs a bracket, `newton` a starting point (its derivative is a central difference).

```console
$ python3 -m numtools --help
usage: python -m numtools [-h] [--version] {bisect,newton} ...

Find a root of EXPR (an expression in x) with the PY-1 solvers.

positional arguments:
  {bisect,newton}
    bisect         bisection on a bracket [a, b]
    newton         Newton's method from x0 (derivative by central differences)

options:
  -h, --help       show this help message and exit
  --version        show program's version number and exit
$ python3 -m numtools bisect "cos(x)-x" 0 1 -v
root = 0.739085133216
iterations = 40
f(root) = -8.477e-13
$ python3 -m numtools newton "x*x-2" 1.5 --tol 1e-14 -v
root = 1.41421356237
iterations = 5
f(root) = -4.441e-16
$ python3 -m numtools bisect "x*x+1" -1 1
error: f(a) and f(b) must have opposite signs
$ python3 -m numtools bisect "cos(x)-x" 0
usage: python -m numtools bisect [-h] [--tol TOL] [-v] expr a b
python -m numtools bisect: error: the following arguments are required: b
```

argparse's own errors (a missing argument, a non-number) print a usage line and exit with status 2; errors from the solver are caught in `main` and printed as `error: ...` with the same status. Because `main(argv)` takes the argument list, `project.py` Part 1 tests it by calling `main(["bisect", "x*x-2", "1", "2", "--tol", "1e-6"])` in-process and capturing the output, with no subprocess.

Caution: `make_function` calls `eval` on the text you type, with the built-ins removed. That is acceptable for a tool run on your own input; it is **not** a sandbox. Never pass untrusted text to `eval`; use a small parser or SymPy's `sympify` (PY-8).

### Example 4.7 (from LA-6): the cofactor determinant, and factorial growth

LA-6 eq. (3.12) as a recursion: expand along row 1, each minor is a smaller determinant, the base case is $1 \times 1$. A counter records the products $a_{1j} M_{1j}$, which LA-6 §3.5 counts as $T(n) = n!\sum_{k=1}^{n-1} 1/k!$.

```python
import math
import time

mults = 0

def det_cofactor(A):
    global mults
    n = len(A)
    if n == 1:
        return A[0][0]
    total = 0
    for j in range(n):
        minor = [row[:j] + row[j + 1:] for row in A[1:]]
        mults += 1
        total += (-1) ** j * A[0][j] * det_cofactor(minor)
    return total

print(det_cofactor([[0, 2, 1], [1, 1, 3], [2, -1, 1]]), "(LA-6 Example 4.1 gets 7)")
previous = None
for n in range(2, 8):
    A = [[(i * 7 + j * 3) % 5 + 1 + (i == j) for j in range(n)] for i in range(n)]
    mults = 0
    start = time.perf_counter()
    det_cofactor(A)
    t = time.perf_counter() - start
    assert mults == sum(math.factorial(n) // math.factorial(k) for k in range(1, n))
    note = f" (x{t / previous:.1f})" if previous else ""
    print(f"n = {n}: {mults:>6,} multiplications, {t:.5f} s{note}")
    previous = t
```
```text
7 (LA-6 Example 4.1 gets 7)
n = 2:      2 multiplications, 0.00001 s
n = 3:      9 multiplications, 0.00002 s (x2.6)
n = 4:     40 multiplications, 0.00003 s (x2.2)
n = 5:    205 multiplications, 0.00023 s (x6.7)
n = 6:  1,236 multiplications, 0.00102 s (x4.4)
n = 7:  8,659 multiplications, 0.00734 s (x7.2)
```

The count matches LA-6's closed form exactly (the `assert`), and the time ratio climbs toward $n$ (noisily, at small sizes). At $n = 20$ the count is about $4 \times 10^{18}$, while elimination (what `numpy.linalg.det` does) needs about $\tfrac23 n^3$ operations. The minors are *not* repeated subproblems, so a cache would not help: **a recursion that branches $n$ ways at every level is a factorial algorithm**, and only a different algorithm fixes it. `project.py` Part 4 compares with NumPy and plots the growth.

## 5. Why it matters and where it connects

- **PY-4.** `@property`, `@classmethod` and `@dataclass` are decorators (§3.5); a class with state and one method is `make_counter` in other clothes.
- **PY-5.** `@pytest.mark.parametrize` is a decorator with arguments; `retry` is the simplest error-recovery wrapper; `timer`'s `try/finally` previews context managers.
- **PY-6 and PY-8.** `quad`, `root_scalar`, `minimize` and `solve_ivp` take a function first, with parameters fixed by `partial` or `args=`. Example 4.7's cost is why `numpy.linalg.det` factorizes (LA-6 §3.5).
- **PY-7, PY-9.** One callback per curve in a loop is the late-binding trap of §3.4; `df.apply(func)`, `groupby(...).agg(func)` and `str.extract` (Example 4.5's pattern) are higher-order.
- **PY-10.** `lru_cache` against a hand-built cache is a profiling question; packaging `numtools` is a `pyproject.toml` and `pip install -e .`; Example 4.6 is the template for every analysis script you ship.

## 6. Exercises

See [exercises.md](exercises.md). A1 to A6 are predict-the-output drills; do all of them. B1 to B7 are the module; B1 (`counted`, `trace`), B3 (a decorator with arguments), B4 (recursion limit), B5 (a package) and B6 (an `argparse` tool) are essential. C1 to C3 are a decorator for functions and methods, a measured cache-hit-rate study, and an iterative rewrite with a performance comparison.

## 7. Coding project

`project.py` imports `numtools` from this folder (inserting the folder into `sys.path` first, §3.8) and checks everything above: the import system and CLI through fresh interpreters; the timing table of Example 4.3 (counts and ordering asserted); naive against memoized Fibonacci on a log axis, with eq. (3.1) and the $n + 1$ body executions asserted exactly and the time ratio between $n = 25$ and $n = 20$ within a factor of 2 of $\varphi^5$; the cofactor determinant against `numpy.linalg.det` for $n = 2 \dots 8$ with a growth plot; merge sort and both binary searches against `sorted` and `bisect_left`; `re` parsing of 1000 generated lines; and every number of §4.

**Build it yourself:** before opening `numtools/`, write `timer`, `memoize` and `retry(n)` and test them against the package's versions; then write `cli.py`'s `main` from the four steps of §3.9. Then "Your turn": a call counter, an `lru_cache` hit rate, a recursion rewritten with an explicit stack, a new CLI subcommand.

## 8. Where you'll struggle

- **The late-binding closure.** You will write a loop that builds callbacks (a handler per button, a curve per iteration of a plotting loop) and every callback will use the last value, within a month. Fix: default argument, `partial` or a factory (§3.4); diagnosis: "is this name looked up now, or when the function runs?".
- **Decorators are baffling until you have written one by hand.** Then they are four lines. Do not skip build-it-yourself (1). Typical mistakes: no `return` of the wrapped result (the function now returns `None`), no `*args, **kwargs`, and `@retry` where `@retry(3)` is needed (the function is passed as `n`, and the error points elsewhere). **`functools.wraps` is not optional:** without it every decorated function is `wrapper` in tracebacks, `help()` and test reports.
- **Recursion depth of about 1000** means Python is the wrong tool for deep recursion; a cache does not change the limit. Learn the explicit-stack rewrite. The opposite mistake is recursion that is correct and exponential (naive Fibonacci) or factorial (cofactors): count the calls on a small input first.
- **Binary search off-by-ones** (`len(a)` or `len(a) - 1`, `<` or `<=`, `mid + 1` or `mid`): all eight combinations compile and most are wrong. State the invariant in a comment (§3.7).
- **"Module not found" when the file is right there** is about where you ran Python from. Print `sys.path[0]`; run `python -m package.module` from the project root, or install the package. Related: your own file named `random.py`, `statistics.py` or `math.py` shadows the standard module, because the script's directory is first on `sys.path`; the symptom is an `AttributeError` on a standard-library function. Rename it and delete `__pycache__`. And running a file inside a package directly fails with "attempted relative import with no known parent package" (§3.8); use `-m`.
- **`from module import *`** hides where names come from and shadows built-ins (§3.8). **`random.seed` at import time** resets the global generator for the whole process (§3.9); seed in `main`.

## 9. Check yourself

1. What does `def` do at run time, and what does `del f` delete?
2. What do the lambdas in `[lambda: i for i in range(3)]` return when called, and why? Name two fixes.
3. What are `__closure__` and a cell? When must you write `nonlocal`?
4. What does `@timer` above `def f` expand to? Why must the wrapper take `*args, **kwargs` and return the result?
5. What does `functools.wraps` do? Why does `@retry(3)` need three levels of function?
6. Which functions can be memoized, and what must their arguments be?
7. State the recurrence for the calls of naive `fib(n)`, its solution, and the growth factor per unit of $n$.
8. What is the default recursion limit, and does a cache change it?
9. State the invariant of the binary search and why the loop terminates.
10. What does `import m` do, in order? How often does a module's top-level code run?
11. What is `sys.path[0]` for `python script.py` and for `python -m pkg.mod`?
12. Why is `from m import *` banned here, and why must a library never call `random.seed` at import time?

## 10. Pair with

- **Ramalho, *Fluent Python* 2e, ch. 7 and 9.** The deepest account of §3.1 to §3.6, including `singledispatch` and class-based decorators (PY-4). Read ch. 9 after writing `timer` by hand.
- **Downey, *Think Python* 3e, ch. 5 to 6.** The gentlest recursion and stack-diagram practice.
- **Python Tutorial §4.8 and §6.** Short and authoritative on argument kinds and the module search path.
- **Slatkin, *Effective Python* 2e, items 19 to 26 and 86 to 90.** Opinionated rules for the argument and closure traps and for module layout.
- **The `functools` and `itertools` documentation.** Where `lru_cache`, `partial`, `wraps` and `reduce` live.
- **Python Packaging User Guide, "Packaging Python projects" (first page).** Where `numtools` goes next.
- **Berkeley CS61A lectures on higher-order functions and recursion.** Environment diagrams, drawn better than ASCII allows.
