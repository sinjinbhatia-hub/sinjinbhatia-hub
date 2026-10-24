# PY-1: Objects, Names, Numbers, and Control Flow

**Time:** ~8 hours · **Prerequisites:** none (start here, week 1, alongside F1) · **Unlocks:** PY-2, and every `project.py` in the curriculum

> **Conventions in this module.** Python 3.11+; standard library only (NumPy appears only in the project); PEP 8 names. Blocks starting with `>>>` are REPL transcripts; a script block is followed by its output, and its first line `# file: name.py` is the name to save it under (later scripts import earlier ones). Every output was produced by running the code on CPython 3.13; where 3.11 or 3.12 differ, the text says so. Paths and timings differ on your machine. **Floating point:** a "double" is the IEEE 754 64-bit format Python calls `float`. $\varepsilon = 2^{-52} \approx 2.22\times10^{-16}$ is the gap between $1.0$ and the next double (`sys.float_info.epsilon`); $u = 2^{-53}$ is the largest relative rounding error (LA-13 §3.1 calls this $\varepsilon_{\text{mach}}$). The relative error of $\tilde x \approx x$ is $|\tilde x - x|/|x|$.

## 0. Where this sits

You can already make Python do things; this module is about predicting what it will do. It needs nothing earlier. It makes possible PY-2 (collections, on the aliasing model of §3.2), PY-3 (closures, on the scope rules of §3.11), PY-5 (the traceback of §4.1 becomes exceptions), PY-6 (views and copies are §3.2 for arrays) and LA-13 §3.2 (cancellation for matrices; §3.5 here is the scalar case).

## 1. The problem that forced this idea

The problem is the gap between code that runs and code you can predict. Here are three programs that run, and that most people who "know Python" mispredict:

```pycon
>>> a = [1, 2]
>>> b = a
>>> b.append(3)
>>> a
[1, 2, 3]
>>> 0.1 + 0.2 == 0.3
False
>>> def add_item(x, bucket=[]):
...     bucket.append(x)
...     return bucket
...
>>> add_item(1)
[1]
>>> add_item(2)
[1, 2]
```

Changing `b` changed `a`. One tenth plus two tenths is not three tenths. A function with an "empty" default remembered an earlier call. None of this is a bug: each follows from a rule the interpreter applies without exception, and without the rule you try things until the output looks right. A NumPy slice that aliases its parent (PY-6), a class attribute shared by every instance (PY-4) and a long sum that depends on how it was accumulated (§3.6) are the same puzzles in other clothes.

This module closes the gap with four mental models: names tied to objects (§2, §3.2); calls as stacks of frames, whose printout, the traceback, is the most-used debugging skill in the course (§3.1, §4.1); floats as points on a grid (§3.4–3.6); and control flow and scope as rules you can state (§3.9–3.11). You build each before you see the library tool that packages it.

## 2. Intuition

### Names are labels, not boxes

In C or Fortran a variable is a box: `x = 5` puts 5 in the box `x`, and `y = x` copies it into a second box. Python has no boxes. An **object** lives in memory and has an **identity** (`id(obj)`, fixed for its lifetime), a **type** (fixed) and a **value** (fixed for an `int`, changeable for a `list`). A **name** is an entry in a table (a *namespace*) that points at an object. After `a = [1, 2]; b = a; c = list(a)`:

```text
   names (a namespace)              objects
   -------------------              ------------------------------------------
        a ────────┐
                  ├────────────▶  @1  list   [ ──▶ int 1 , ──▶ int 2 ]
        b ────────┘
        c ─────────────────────▶  @2  list   [ ──▶ int 1 , ──▶ int 2 ]
```

(`@1`, `@2` stand for ids.) `b = a` copied an arrow, not a list. `list(a)` built a new list whose slots point at the *same* two integers. Now `a.append(3)` changes object `@1` and every name tied to it sees that; `c` does not. But `a = [9]` leaves `@1` alone: it moves the arrow `a`, and `b` still points at the old list. The whole model is four rules:

- **R1.** `name = expression` evaluates the expression to an object and ties the name to it. It never copies.
- **R2.** *Mutation* changes an object (`a.append(3)`, `a[0] = 7`); every name tied to it sees the change. *Rebinding* moves one arrow (`a = [9]`); no other name moves.
- **R3.** Calling a function ties the parameter names to the argument objects, exactly as R1 does.
- **R4.** `def` is a statement that runs once and creates a function object; default values are evaluated *then*, once.

R4 is the mutable-default puzzle. §3.2 has you apply R1–R4 to eight snippets.

### The call stack

Each call gets a **frame**: a fresh namespace for its parameters and local names, plus where to resume in the caller. Frames pile up. If `main` calls `summarize`, which calls `parse`:

```text
      ┌────────────────────────────────────┐
      │ parse        token='.'  values=[]  │  ← newest; the error happens here
      ├────────────────────────────────────┤
      │ summarize    lines='1.0,2.0,3.0'   │
      ├────────────────────────────────────┤
      │ main         data='1.0,2.0,3.0'    │
      ├────────────────────────────────────┤
      │ <module>     (the file's top level)│  ← oldest
      └────────────────────────────────────┘
```

A traceback prints this stack oldest first, so the newest frame, where the exception was raised, is at the bottom. That is why you read tracebacks from the bottom up.

### A float is a point on a grid

A double stores a sign, an 11-bit exponent and a 52-bit fraction. Between consecutive powers of two there are always $2^{52}$ doubles, spread over an interval that doubles in width each time, so the grid gets coarser away from zero:

```text
   [1, 2):  |·|·|·|·|·|·|·|·|·|·|·|·|·|·|    spacing 2^-52
   [2, 4):  |  |  |  |  |  |  |  |  |  |    spacing 2^-51
   [4, 8):  |     |     |     |     |        spacing 2^-50
```

Every real number is replaced by the nearest grid point, so the *absolute* error grows with the size of the number while the *relative* error stays near $10^{-16}$. You can already guess three results that §3.4–3.6 prove: `0.1` is not on the grid, so `0.1 + 0.2` is not `0.3`; adding a tiny number to a huge one can do nothing; and subtracting two nearly equal numbers wipes out the digits they shared, since the answer is small but the input errors were sized by the large inputs.

## 3. Formal development

Each subsection states a rule and has you predict or build something before you meet the library version. **Do the predictions on paper**; pasting code without predicting turns this module into reading, and reading does not close the gap.

### 3.1 Running Python three ways

```text
$ python3                  # the REPL: type a statement, see the value of an expression
$ python3 report.py        # run a file top to bottom, then exit
$ python3 -i report.py     # run the file, then stay in a REPL with the file's names defined
```

The REPL prints the `repr` of each expression whose value is not `None`; a script prints only what you `print`. Use the REPL for one-line questions, files for anything you will run twice, and `-i` when a script misbehaves: you can inspect its variables afterwards, even after a crash. When a script crashes Python prints a **traceback**. Read it in three moves, practised in §4.1:

1. **Last line first**: exception type and message, *what* went wrong.
2. **The line above it**: the statement that raised, *where*.
3. If that statement looks fine, walk **up** the frames, each showing who called the one below. The bug is in the frame where your intention and your code disagree, often not the bottom one.

### 3.2 Everything is an object

```pycon
>>> type(3), type(3.0), type("3"), type([3]), type(len)
(<class 'int'>, <class 'float'>, <class 'str'>, <class 'list'>, <class 'builtin_function_or_method'>)
>>> a = [1]
>>> b = [1]
>>> a == b, a is b, id(a) == id(b)
(True, False, False)
>>> c = a
>>> c is a
True
>>> int("256") is int("256"), int("257") is int("257")
(True, False)
```

`==` asks whether two objects have equal *values*; `is` asks whether two names are tied to the *same object*. The last line shows why you never use `is` on numbers or strings: CPython shares one object per small integer (−5 to 256), an implementation detail. The one idiom that uses `is` is `x is None`, since there is exactly one `None`. An object is **immutable** if its value cannot change after creation, **mutable** if it can:

| Immutable | Mutable |
|---|---|
| `int`, `float`, `complex`, `bool`, `None`, `str`, `bytes`, `tuple`, `frozenset` | `list`, `dict`, `set`, `bytearray`, instances of classes you write (PY-4) |

Immutable does not mean everything reachable is frozen: a tuple holds fixed *references*, and a list inside it can still be mutated.

**Build-it-yourself sequence 1: the names-and-values model.**

▶ **Try it first.** Below are eight snippets in one file. For each, write down what its `print` shows, using only R1–R4. Do not run anything yet.

```python
# file: aliasing_quiz.py
# 1
a = [1, 2]
b = a
b.append(3)
print("1:", a)

# 2
a = [1, 2]
b = a
b = b + [3]
print("2:", a)

# 3
a = [1, 2]
b = a
b += [3]
print("3:", a)

# 4
x = 5
y = x
y += 1
print("4:", x, y)

# 5
grid = [[0] * 2] * 2
grid[0][0] = 9
print("5:", grid)


# 6
def grow(items):
    items.append(1)


def regrow(items):
    items = items + [1]


lst = []
grow(lst)
regrow(lst)
print("6:", lst)

# 7
t = (1, [2, 3])
t[1].append(4)
print("7:", t)


# 8
def collect(x, bucket=[]):
    bucket.append(x)
    return bucket


collect(1)
print("8:", collect(2))
```

```text
1: [1, 2, 3]
2: [1, 2]
3: [1, 2, 3]
4: 5 6
5: [[9, 0], [9, 0]]
6: [1]
7: (1, [2, 3, 4])
8: [1, 2]
```

■ **What just happened.** Score yourself, then read the reasons.

1. `b = a` ties a second name to the list (R1); `append` mutates it (R2).
2. `b + [3]` builds a *new* list and `b = ...` moves only the arrow `b`; `a`'s list is untouched.
3. Looks like 2, behaves like 1: `b += [3]` extends the list in place (`list.__iadd__` mutates and returns the same object). On a list `+=` is a mutation; `b = b + [3]` is a rebinding.
4. Integers are immutable, so `y += 1` computes `6` and rebinds `y`; the object `5` cannot change.
5. `[0] * 2` makes one inner list; the outer `* 2` copies the *reference* twice: two arrows, one row. Build rows with `[[0] * 2 for _ in range(2)]`.
6. R3: in `grow`, `items` is tied to the caller's list, so `append` mutates it. In `regrow`, `items = items + [1]` rebinds the *local* name and the caller never sees it. ("Pass by value" and "pass by reference" both mislead: the parameter is assigned the argument object.)
7. The tuple still holds the same two references; the list one of them points at changed. `t[0] = 9` would be a `TypeError`.
8. R4. `def` ran once and stored the list in the function object (`collect.__defaults__` is `([1, 2],)`); every call that omits `bucket` ties it to that list. The first call looks fine, which is why this bug survives testing. Example 4.6 fixes it.

Scoring 4/8 is normal and is why this sequence exists; redo the misses in two days. When you do want a copy: `list(a)`, `a[:]` and `a.copy()` make a new outer list whose slots point at the same elements (a **shallow** copy); `copy.deepcopy(a)` copies recursively. Example 4.2 shows both.

### 3.3 Integers

Python's `int` has arbitrary precision: it grows to hold the value, so there is no overflow.

```pycon
>>> 2**100
1267650600228229401496703205376
>>> 10**20 + 1, float(10**20) + 1
(100000000000000000001, 1e+20)
>>> 7 / 2, 7 // 2, 7 % 2
(3.5, 3, 1)
```

`/` always returns a float (`4 / 2` is `2.0`). `//` is floor division and `%` the matching remainder.

▶ **Try it first.** Predict `-7 // 2`, `-7 % 2`, `7 // -2`, `7 % -2`. Then find a rule that makes $a = b\,(a \mathbin{//} b) + a \mathbin{\%} b$ true for every sign.

**Derivation (the floor rule).** For integers $a, b \ne 0$ let $q = \lfloor a/b \rfloor$ (`a // b`) and $r = a - bq$ (`a % b`), so

$$a = bq + r. \tag{3.1}$$

1. By the definition of the floor, $q \le a/b < q + 1$.
2. If $b > 0$, multiply by $b$: $bq \le a < bq + b$. Subtract $bq$: $0 \le r < b$.
3. If $b < 0$, multiplying flips both inequalities: $bq \ge a > bq + b$. Subtract $bq$: $0 \ge r > b$.
4. So $r$ is zero or has the sign of the *divisor*, and $|r| < |b|$. ∎

```pycon
>>> import math
>>> -7 // 2, -7 % 2, 7 // -2, 7 % -2
(-4, 1, -4, -1)
>>> divmod(-7, 2), math.fmod(-7, 2)
((-4, 1), -1.0)
>>> all(a == 3 * (a // 3) + a % 3 for a in range(-20, 20))
True
```

■ **What just happened.** $\lfloor -3.5 \rfloor = -4$, so `-7 // 2` is $-4$ and $r = -7 - 2(-4) = 1$. The payoff: for positive $n$, `i % n` is always in `range(n)`, so `items[i % len(items)]` wraps correctly for negative `i`. C truncates toward zero and gives remainders with the sign of the dividend; `math.fmod` does that. `divmod(a, b)` returns `(a // b, a % b)`.

The bit operators `&`, `|`, `^` (and, or, exclusive or), `~x` (which is `-x - 1`) and the shifts `<<`, `>>` act on the binary expansion:

```pycon
>>> bin(10), 12 & 10, 12 | 10, 12 ^ 10, 1 << 10, -7 >> 1
('0b1010', 8, 14, 6, 1024, -4)
>>> 40 & (40 - 1), 32 & (32 - 1)
(32, 0)
```

`x >> k` equals `x // 2**k`, floor included. The trick on the second line clears the lowest set bit, so `n & (n - 1)` is `0` exactly for powers of two.

### 3.4 Floats: IEEE 754 doubles in plain words

A double is 64 bits: a sign bit $s$, an 11-bit exponent field $E$, and a 52-bit fraction field $F$ (an integer from $0$ to $2^{52}-1$). For $1 \le E \le 2046$ (a *normal* number) the value is

$$x = (-1)^s\left(1 + \frac{F}{2^{52}}\right) 2^{E-1023}. \tag{3.2}$$

The "1 +" is an implicit leading bit (53 significant bits, 52 stored). $E = 0$ encodes zero and the *subnormals* (no implicit 1; they fill the gap down to $5\times10^{-324}$); $E = 2047$ encodes `inf` ($F = 0$) and `nan` ($F \ne 0$).

For example `1.0` has $E = 1023$ and $F = 0$; halving a number subtracts one from $E$. `0.1` has $E = 1019$ (that is, $2^{-4}$) and a fraction pattern `1001 1001 ...` that repeats and is cut off after 52 bits, which matters below. Exercise A5 has you print these fields with `struct.pack`; the project decodes them back through (3.2) and asserts the values quoted here.

**Build-it-yourself sequence 2: machine epsilon.**

▶ **Try it first.** Predict each boolean.

```pycon
>>> 1 + 2**-52 == 1
False
>>> 1 + 2**-53 == 1
True
>>> (1 + 2**-53) + 2**-53 == 1
True
>>> 1 + (2**-53 + 2**-53) == 1
False
```

The last two differ only in parentheses: floating-point addition is not associative.

1. **Spacing.** Take $x \in [2^e, 2^{e+1})$. By (3.2) neighbouring doubles differ by one unit of $F$, which is $2^e\cdot 2^{-52}$:

   $$\operatorname{ulp}(x) = 2^{e-52}. \tag{3.3}$$

   (`math.ulp(x)` computes it.) On $[1, 2)$ it is $\varepsilon = 2^{-52}$.
2. **Rounding.** An exact result between two doubles becomes the nearer one; an exact tie goes to the one whose last fraction bit is $0$ (*round half to even*).
3. **Halve until it vanishes.**

```python
# file: epsilon.py
import sys

eps = 1.0
while 1.0 + eps / 2 != 1.0:
    eps /= 2
print(eps, eps == sys.float_info.epsilon, eps == 2**-52)
```

```text
2.220446049250313e-16 True True
```

4. **Why it stops there.** `1 + 2**-52` is the next double above $1$; `1 + 2**-53` is exactly halfway to it, a tie, and round-half-even picks $1$. So `eps` ends at $2^{-52}$, the *gap* above $1$. This answers the quiz: `(1 + 2**-53) + 2**-53` rounds back to $1$ twice, while `2**-53 + 2**-53` equals $2^{-52}$ before it is added, and survives.
5. **The rounding model.** If $x \in [2^e, 2^{e+1})$ is replaced by its nearest double, the error is at most half a spacing, $2^{e-53}$ by (3.3). Since $|x| \ge 2^e$ the *relative* error is at most $u = 2^{-53}$:

   $$\mathrm{fl}(x) = x(1+\delta), \qquad |\delta| \le u. \tag{3.4}$$

   The standard requires $+, -, \times, \div, \sqrt{\ }$ to return the correctly rounded exact result, so $\mathrm{fl}(x \circ y) = (x\circ y)(1+\delta)$ with the same bound.

■ **What just happened.** `sys.float_info.epsilon` is the grid spacing at $1$; the worst relative rounding error is half of it. Spacing grows with magnitude and eventually exceeds $1$:

```pycon
>>> import math
>>> math.ulp(1.0), math.ulp(1000.0), math.ulp(1e16)
(2.220446049250313e-16, 1.1368683772161603e-13, 2.0)
>>> 2.0**53 + 1 == 2.0**53, float(2**53 + 1) == float(2**53)
(True, True)
>>> math.nextafter(1.0, 2.0), math.nextafter(1.0, 0.0)
(1.0000000000000002, 0.9999999999999999)
```

Above $2^{53} \approx 9.007\times10^{15}$ not every integer is a double, so counts and indices must never be floats.

**Why `0.1` is not 0.1.** A fraction in lowest terms has a finite binary expansion exactly when its denominator is a power of 2. $0.1 = 1/10$ and $10 = 2\cdot5$, so its expansion is the repeating $0.0001\,1001\,1001\ldots_2$. A double keeps 53 significant bits, so Python stores the nearest grid point, slightly *above* $0.1$. Example 4.3 prints it exactly and traces `0.1 + 0.2` to a single tie.

**Comparing floats.** Computed results carry rounding error, so `==` asks for more than they deliver. `math.isclose(a, b, rel_tol=1e-9, abs_tol=0.0)` tests

$$|a - b| \le \max\bigl(\texttt{rel\_tol}\cdot\max(|a|,|b|),\ \texttt{abs\_tol}\bigr). \tag{3.5}$$

```pycon
>>> import math
>>> 0.1 + 0.2 == 0.3, math.isclose(0.1 + 0.2, 0.3)
(False, True)
>>> math.isclose(1e-20, 0.0), math.isclose(1e-20, 0.0, abs_tol=1e-12)
(False, True)
>>> math.isclose(1e10 + 1, 1e10), round(2.675, 2)
(True, 2.67)
```

The relative tolerance scales with the numbers: a difference of $1$ is "close" at $10^{10}$. Against zero it can never pass (the right side of (3.5) is $0$), so give `abs_tol`, taken from the scale of your problem (say `1e-12 * typical_size`), not from habit. Do not reach for `round()`: it acts on the stored value, so `round(2.675, 2)` is `2.67` because the stored number is slightly below $2.675$. Rounding moves the boundary without removing it.

**Overflow, `inf`, `nan`.**

```pycon
>>> import math
>>> inf, nan = float("inf"), float("nan")
>>> 1e308 * 10, 5e-324 / 2, inf - inf
(inf, 0.0, nan)
>>> nan == nan, nan != nan, nan < 1, math.isnan(nan)
(False, True, False, True)
>>> nan in [nan], float("nan") in [float("nan")]
(True, False)
>>> 2.0 ** 10000
Traceback (most recent call last):
  ...
OverflowError: (34, 'Numerical result out of range')
```

Past the largest double ($\approx 1.8\times10^{308}$) you get `inf`; operations with no meaningful answer give `nan`. IEEE 754 declares `nan` unordered: every comparison with it is false except `!=`, since two "no number came out" results are not known to be equal. So `x != x` detects `nan` (use `math.isnan`). Containers test identity before equality, hence `nan in [nan]`. One `nan` silently breaks `==`, `sorted`, `max` and `min`; check for it where data enters. Python is inconsistent about overflow: bare `*` and `+` return `inf` quietly, while `**` and `math.exp` raise.

### 3.5 Catastrophic cancellation: the quadratic formula

The roots of $ax^2 + bx + c = 0$ are $\dfrac{-b \pm \sqrt{b^2 - 4ac}}{2a}$. Take $a = c = 1$, $b = 10^8$: the product of the roots is $c/a = 1$ and their sum is $-b$, so they are near $-10^8$ and $-10^{-8}$.

**Build-it-yourself sequence 3.**

▶ **Try it first.** Predict how many correct digits the textbook formula gives for the small root. Hint: its numerator is $-b + \sqrt{b^2-4ac}$, the sum of $-10^8$ and a number slightly below $10^8$.

1. **Measure.** Write the naive formula and a reference using 60-digit decimal arithmetic (`decimal` can carry far more digits than a double). Save this as `quadratic.py`:

```python
# file: quadratic.py
import math
from decimal import Decimal, localcontext


def roots_naive(a, b, c):
    """Roots of a*x**2 + b*x + c by the textbook formula."""
    d = math.sqrt(b * b - 4 * a * c)
    return (-b - d) / (2 * a), (-b + d) / (2 * a)


def roots_exact(a, b, c):
    """Reference roots from 60-digit decimal arithmetic."""
    with localcontext() as ctx:
        ctx.prec = 60
        a, b, c = Decimal(a), Decimal(b), Decimal(c)
        d = (b * b - 4 * a * c).sqrt()
        return (-b - d) / (2 * a), (-b + d) / (2 * a)


def rel_err(x, ref):
    """Relative error of the float x against the Decimal reference ref."""
    with localcontext() as ctx:
        ctx.prec = 60
        return float(abs((Decimal(x) - ref) / ref))


if __name__ == "__main__":
    big_ref, small_ref = roots_exact(1.0, 1e8, 1.0)
    big, small = roots_naive(1.0, 1e8, 1.0)
    errs = rel_err(big, big_ref), rel_err(small, small_ref)
    print(f"naive roots: {big!r}  {small!r}")
    print(f"relative errors: {errs[0]:.2e}  {errs[1]:.2e}")
    print(f"sqrt(b*b - 4ac) = {math.sqrt(1e16 - 4)!r}")
    print(f"grid spacing there = {math.ulp(1e8):.3e}")
```

```text
naive roots: -100000000.0  -7.450580596923828e-09
relative errors: 1.00e-16  2.55e-01
sqrt(b*b - 4ac) = 99999999.99999999
grid spacing there = 1.490e-08
```

   The guard `if __name__ == "__main__":` runs the demo only when the file is executed, not when another file imports from it (step 5 does; PY-3 explains).
2. **Diagnose.** The big root is right; the small one is wrong in its *first* digit ($-7.45\times10^{-9}$, not $-10^{-8}$). The computed square root is a double near $10^8$, where the grid spacing ($1.49\times10^{-8}$) exceeds the answer. The subtraction `-b + d` is itself exact; it exposes the rounding already committed in `d`.
3. **Quantify (derivation).** Suppose we hold inexact copies $\tilde x = x(1+\delta_1)$ and $\tilde y = y(1+\delta_2)$ with $|\delta_i| \le u$ by (3.4). Then $(\tilde x - \tilde y) - (x - y) = x\delta_1 - y\delta_2$, of size at most $(|x|+|y|)\,u$. Divide by $|x - y|$:

   $$\frac{|(\tilde x - \tilde y) - (x-y)|}{|x-y|} \le \frac{|x| + |y|}{|x - y|}\, u. \tag{3.6}$$

   The factor explodes when $x \approx y$. Here $|x| + |y| \approx 2\times10^8$ and $|x-y| \approx 2\times10^{-8}$, a factor near $10^{16} \approx 1/u$: the bound allows total loss of all digits. This is **catastrophic cancellation**: subtraction creates no new error, it promotes the old error from the last digit to the first. (LA-13 §3.2 derives the same bound as its equation (3.4) and extends it to matrices.)
4. **Rationalize (derivation).** For $b > 0$ the dangerous root is $x_+ = (-b + d)/(2a)$ with $d = \sqrt{b^2-4ac}$. Multiply top and bottom by $(-b - d)$. The numerator becomes $b^2 - d^2 = 4ac$, so

   $$x_+ = \frac{4ac}{2a(-b-d)} = \frac{2c}{-b-d}.$$

   The denominator adds two negatives: no cancellation. The other root already does. For $b < 0$ the roles swap; choosing the sign that adds covers both cases:

   $$q = -\tfrac12\bigl(b + \operatorname{sign}(b)\,d\bigr), \qquad x_1 = \frac{q}{a}, \qquad x_2 = \frac{c}{q}. \tag{3.7}$$

   (Vieta's $x_1x_2 = c/a$ holds by construction.)
5. **Verify.**

```python
# file: quadratic_stable.py
import math
from quadratic import roots_exact, rel_err


def roots_stable(a, b, c):
    """Roots of a*x**2 + b*x + c without cancellation, eq. (3.7)."""
    d = math.sqrt(b * b - 4 * a * c)
    q = -(b + math.copysign(d, b)) / 2
    return q / a, c / q


big_ref, small_ref = roots_exact(1.0, 1e8, 1.0)
big, small = roots_stable(1.0, 1e8, 1.0)
print(f"stable roots: {big!r}  {small!r}")
print(f"relative errors: {rel_err(big, big_ref):.2e}  {rel_err(small, small_ref):.2e}")
```

```text
stable roots: -100000000.0  -1e-08
relative errors: 1.00e-16  7.91e-17
```

■ **What just happened.** One rewrite took the small root's relative error from $0.25$ to $8\times10^{-17}$, below $u$: equal formulas in exact arithmetic, very different on the grid. The rule: before subtracting two quantities, ask whether they can be nearly equal, and if so look for a form without the subtraction (conjugates, trig identities, `math.expm1` and `math.log1p` for $e^x - 1$ and $\ln(1+x)$ at small $x$). Limits: `b * b` overflows for $|b| > 10^{154}$, and a nearly double root ($b^2 \approx 4ac$) is intrinsically ill-conditioned, which no rewrite cures.

### 3.6 Summation error and Kahan's trick

Adding many floats is the same problem in slow motion: each `total + x` rounds to the grid at the size of `total`, which is coarse once `total` is large.

**Build-it-yourself sequence 6.**

▶ **Try it first.** Sum $0.1$ ten million times with a plain loop. The correct total for the stored double `0.1` is $1000000.0000000000555\ldots$, which is `1000000.0` on the grid. Predict the size of the error: about $10^{-10}$, $10^{-6}$ or $10^{-4}$?

1. **Why it errs.** Near $10^6$ the spacing is `math.ulp(1e6)` $= 1.16\times10^{-10}$, so one addition can be off by $5.8\times10^{-11}$. The same $0.1$ rounds the same way every time, so errors do not cancel: $10^7$ additions can reach $5.8\times10^{-4}$.
2. **The lost part is recoverable (derivation).** Let the running sum be $s$, the next term $y$ with $|s| \ge |y|$, and $t = \mathrm{fl}(s + y)$. Then $t - s$ is computed exactly (this needs $|s| \ge |y|$) and is the part of $y$ that made it into $t$; and $y - (t - s)$ is also computed exactly and is precisely what was rounded away:

   $$s + y = t + e, \qquad e = y - (t - s) \text{ computed exactly.} \tag{3.8}$$

   This is *Fast2Sum* (Dekker, 1971, which has the proof; the project checks (3.8) against exact rational arithmetic). The rounding error of an addition can be computed, with three extra operations.
3. **Carry it forward.** Keep the previous addition's error in $c = -e = (t - s) - y$ and correct the next term before adding: $y' = x - c = x + e_{\text{previous}}$. That is Kahan's algorithm:

```python
# file: kahan.py
import math


def naive_sum(xs):
    total = 0.0
    for x in xs:
        total += x
    return total


def kahan_sum(xs):
    """Compensated summation: carry the rounding error of each addition forward."""
    total = 0.0
    c = 0.0                   # minus the rounding error of the last addition
    for x in xs:
        y = x - c             # put the lost part back into the next term
        t = total + y         # big + small: low bits of y are lost here
        c = (t - total) - y   # (added) - (meant to add) = -(lost part)
        total = t
    return total


xs = [0.1] * 10**7
ref = math.fsum(xs)
for name, func in (("naive", naive_sum), ("kahan", kahan_sum)):
    s = func(xs)
    print(f"{name}: {s!r:>18}  abs error {abs(s - ref):.3e}")
print(f"fsum : {ref!r:>18}")
```

```text
naive:  999999.9998389754  abs error 1.610e-04
kahan:          1000000.0  abs error 0.000e+00
fsum :          1000000.0
```

■ **What just happened.** The naive loop is off by $1.6\times10^{-4}$, a relative $10^{-10}$, about a million times worse than one rounding. Kahan's loop matches `math.fsum`, which returns the correctly rounded exact sum of the doubles. Kahan's bound is about $2u\sum|x_i|$ plus a term of order $n\,u^2\sum|x_i|$; the naive bound grows like $n\,u\sum|x_i|$ (Goldberg 1991, Theorem 8 in the appendix; stated, not proved here). The method fails if a term can exceed the running total; Neumaier's variant repairs that. **In practice** use `math.fsum`. Since Python 3.12 the built-in `sum` of floats uses Neumaier's variant, so `sum(xs)` is accurate here on 3.12+ and equals `naive_sum` on 3.11. NumPy's `np.sum` uses pairwise summation, with error growing like $\log n$ (PY-6; the project measures it).

### 3.7 Complex numbers, booleans, and short-circuit logic

Complex numbers are built in: `3 + 4j` is a pair of doubles (mathematics: F1 §3.8).

```pycon
>>> import cmath, math
>>> z = 3 + 4j
>>> z.real, z.imag, abs(z), z.conjugate(), z * 1j
(3.0, 4.0, 5.0, (3-4j), (-4+3j))
>>> cmath.exp(1j * math.pi), (-8) ** (1 / 3)
((-1+1.2246467991473532e-16j), (1.0000000000000002+1.7320508075688772j))
>>> math.sqrt(-1)
Traceback (most recent call last):
  ...
ValueError: math domain error
```

`math` takes real numbers only; `cmath` has the complex versions. `exp(iπ)` is $-1 + 1.22\times10^{-16}i$ because the double nearest $\pi$ misses $\pi$ by $1.22\times10^{-16}$. And `(-8) ** (1 / 3)` is the principal cube root, a complex number, not $-2$.

`bool` is a subclass of `int` (`True + True` is `2`). **Falsy:** `None`, `False`, zero of any numeric type, and empty containers (`""`, `[]`, `()`, `{}`, `set()`). **Everything else is truthy**, including `nan`, `"0"` and `[0]`. `and` and `or` return one of their *operands*, not a bool, and short-circuit.

```pycon
>>> [bool(v) for v in (0, 0.0, "", [], {}, None)]
[False, False, False, False, False, False]
>>> [bool(v) for v in (0.1, " ", [0], "False")]
[True, True, True, True]
>>> 0 or 10, 5 or 10, 0 and 10, 5 and 10
(10, 5, 0, 10)
>>> n = 0
>>> n != 0 and 10 / n > 1
False
```

`a or b` is `a` if `a` is truthy, else `b`; `a and b` is `a` if `a` is falsy, else `b`. The right side is skipped when the left decides, so `n != 0 and 10 / n > 1` is a safe guard. The idiom `x or default` is a trap when `0` or `""` is legitimate: `timeout = timeout or 30` turns an explicit `timeout=0` into `30`. Write `if timeout is None:`.

### 3.8 Strings

A `str` is an immutable sequence of Unicode code points. Indexing starts at 0, negative indices count from the end, and a slice `s[i:j]` is half-open: indices $i, \ldots, j-1$.

```pycon
>>> s = "Hello, Python"
>>> len(s), s[0], s[-1], s[7:], s[:5], s[::-1]
(13, 'H', 'n', 'Python', 'Hello', 'nohtyP ,olleH')
>>> s.split(", "), "-".join(["a", "b", "c"]), "  pad ".strip()
(['Hello', 'Python'], 'a-b-c', 'pad')
>>> s.replace("l", "L"), s.startswith("Hell"), s.upper()
('HeLLo, Python', True, 'HELLO, PYTHON')
>>> s
'Hello, Python'
>>> s[0] = "J"
Traceback (most recent call last):
  ...
TypeError: 'str' object does not support item assignment
```

Every method returns a *new* string and leaves `s` alone, because strings are immutable. The classic bug is writing `s.strip()` and expecting `s` to change. Use `"-".join(parts)` to glue many strings.

**f-string formatting.** `f"{value:spec}"` uses a mini-language, `[[fill]align][sign][0][width][,][.precision][type]`: alignment `<` `>` `^`, `+` forces a sign, `0` pads with zeros, `,` groups thousands. Types: `f` fixed, `e` scientific, `g` general (that many *significant* digits, switching to exponent form when needed), `%` percent, `d` integer, `x` hex, `b` binary. `{x=}` prints the expression and its value.

```pycon
>>> x = 1234.56789
>>> f"{x:.2f}|{x:10.2f}|{x:<10.2f}|{x:^10.2f}|{x:010.2f}|{x:+.1f}"
'1234.57|   1234.57|1234.57   | 1234.57  |0001234.57|+1234.6'
>>> f"{x:,.1f}|{x:.3e}|{x:.3g}|{0.256:.1%}|{255:#x}|{5:08b}"
'1,234.6|1.235e+03|1.23e+03|25.6%|0xff|00000101'
>>> f"{x=:.1f}|{'ab':*^8}"
'x=1234.6|***ab***'
```

Example 4.4 builds a table with `f`, `e` and `g` side by side.

### 3.9 Control flow

`if/elif/else` runs the first branch whose condition is truthy. Indentation is the syntax. `while` repeats while its condition is truthy; `for x in iterable` takes items one at a time. `break` leaves a loop, `continue` skips to the next pass, and an `else:` clause on a loop runs only if the loop finished *without* `break`.

`range(start, stop, step)` is **half-open**: it includes `start` and excludes `stop`. Three reasons, all one fact: `len(range(a, b))` is `b - a`; `range(a, b)` followed by `range(b, c)` is `range(a, c)` with no overlap or gap; `range(n)` is exactly the $n$ valid indices $0, \ldots, n-1$. Slices follow suit: `s[:k] + s[k:] == s` for every `k`.

```pycon
>>> list(range(3, 8)), len(range(3, 8)), list(range(8, 3, -2))
([3, 4, 5, 6, 7], 5, [8, 6, 4])
>>> list(range(0, 3)) + list(range(3, 7)) == list(range(0, 7))
True
>>> for i in range(3):
...     pass
...
>>> i
2
```

The loop variable outlives the loop. The **accumulator pattern** is the workhorse of numerical code: initialize with the identity of the operation (`0` for a sum, `1` for a product, `[]` to collect, `-math.inf` for a maximum), update inside the loop, use afterwards.

```python
# file: control_flow.py
# accumulator and continue: sum of squares of the odd numbers in [1, 10)
total = 0
for k in range(1, 10):
    if k % 2 == 0:
        continue
    total += k * k
print("sum of odd squares:", total)

# while with break and else: steps of the 3n+1 process from 27, with a cap
n, steps = 27, 0
while n != 1:
    n = 3 * n + 1 if n % 2 else n // 2
    steps += 1
    if steps > 1000:
        print("gave up")
        break
else:
    print("reached 1 after", steps, "steps")   # runs only if no break happened

# nested loops with for/else: primes below 30
primes = []
for n in range(2, 30):
    for d in range(2, n):
        if n % d == 0:
            break          # n is composite
    else:
        primes.append(n)   # the inner loop never broke: n is prime
print(primes)
```

```text
sum of odd squares: 165
reached 1 after 111 steps
[2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
```

### 3.10 Functions

```python
# file: functions_demo.py
def scale(values, factor=1.0, *, offset=0.0):
    """Return a new list: each value times factor, plus offset."""
    return [factor * v + offset for v in values]


def normalize_in_place(values):
    """Divide every entry by the largest one. Mutates values and returns None."""
    top = max(values)
    for i in range(len(values)):
        values[i] /= top


data = [1.0, 2.0, 4.0]
print(scale(data, 2.0, offset=1.0), scale(data, offset=-1.0, factor=3.0), data)
print(normalize_in_place(data), data)
```

```text
[3.0, 5.0, 9.0] [2.0, 5.0, 11.0] [1.0, 2.0, 4.0]
None [0.25, 0.5, 1.0]
```

- A **parameter** is a name in the `def` line; an **argument** is the object passed in the call (R3). Arguments match by position or by name (`factor=3.0`). Parameters after a bare `*` are **keyword-only**: `scale(data, 2.0, 1.0)` is an error, which keeps call sites readable.
- A **default** (`factor=1.0`) is evaluated once, when `def` runs (R4): harmless for numbers and strings, a trap for lists and dicts (Example 4.6).
- `return` ends the function. Reaching the end, or a bare `return`, returns `None`; `return a, b` returns a tuple. A first-statement string is the **docstring**, shown by `help()` and `__doc__`.
- A function is **pure** if its result depends only on its arguments and it changes nothing else (no mutated arguments, printing or globals). `scale` is pure; `normalize_in_place` is not. Prefer pure functions: you can test them without setup. When mutation is the point, follow the standard library and return `None`; this prevents `x = my_list.sort()`, which sets `x` to `None`.

Functions are objects (R1): `f = scale` ties a second name to the function, and a function can be an argument. `lambda x: x * x` is a nameless one-expression function. §3.12 relies on this.

### 3.11 Scope: where names live

Python resolves a name inside a function by searching **L**ocal (this frame), **E**nclosing (outer functions), **G**lobal (module top level), **B**uilt-in: **LEGB**. The surprising part concerns *assignment*: if a name is assigned anywhere in a function body, it is local to the **whole** body, decided before the function runs.

▶ **Try it first.** Predict the three results.

```pycon
>>> x = 10
>>> def show():
...     return x
...
>>> show()
10
>>> def bump():
...     x = x + 1
...     return x
...
>>> bump()
Traceback (most recent call last):
  ...
UnboundLocalError: cannot access local variable 'x' where it is not associated with a value
>>> def shadow():
...     x = 99
...     return x
...
>>> shadow(), x
(99, 10)
```

■ **What just happened.** `show` never assigns `x`, so lookup falls through to the global. `bump` assigns `x`, so `x` is local throughout `bump`, and `x + 1` reads a local with no value yet: `UnboundLocalError` (older versions word the message "local variable 'x' referenced before assignment"). `shadow` creates a new local `x` (R1); the global is untouched. To rebind a module-level name you must say so:

```python
# file: scope_demo.py
counter = 0
log = []


def bump_global():
    global counter        # without this line: UnboundLocalError
    counter += 1
    log.append("bumped")  # mutation needs no declaration (R2)


def make_counter():
    count = 0

    def step():
        nonlocal count    # rebind the name in the enclosing function
        count += 1
        return count

    return step


bump_global()
bump_global()
print(counter, log)
step = make_counter()
print(step(), step(), step())
```

```text
2 ['bumped', 'bumped']
1 2 3
```

`global` rebinds a module-level name; `nonlocal` rebinds a name in the nearest enclosing function. Only *rebinding* needs them; mutating `log` does not. Avoid `global`: the function then depends on state that is not in its signature, so you cannot test it alone. Pass values in, return results out. `step` is a **closure**: it remembers `count` after `make_counter` returned (PY-3 builds decorators on this).

### 3.12 Root finding and polynomial evaluation: your first numerical algorithms

**Build-it-yourself sequence 4: bisection, then Newton.**

▶ **Try it first.** $f(x) = \cos x - x$ has $f(0) = 1 > 0$ and $f(1) = \cos 1 - 1 < 0$, so by the intermediate value theorem there is a root in $(0, 1)$. Sketch how to shrink the interval to width $10^{-6}$ using only evaluations of $f$, and predict how many halvings it takes.

1. **The idea.** If $f$ is continuous on $[a, b]$ with $f(a)f(b) < 0$, a root lies in $(a, b)$. Evaluate at the midpoint $m$; keep the half that still has a sign change. After $k$ halvings the width is $(b-a)/2^k$.
2. **The count.** Return the midpoint, which is within half the current width of a root. After $k$ iterations that half-width is $(b-a)/2^{k}$, so stopping at `tol` requires $2^k \ge (b-a)/\text{tol}$:

   $$k = \left\lceil \log_2 \frac{b-a}{\text{tol}} \right\rceil. \tag{3.9}$$

   Each iteration buys one binary digit: $10^{-6}$ on $[0,1]$ takes $\lceil 19.93 \rceil = 20$.
3. **The code.** A tolerance alone is not enough: if `tol` is below the grid spacing near the root, the width stops shrinking and the loop never ends. So there is an iteration cap, and `f` is a parameter. (`raise` stops a function with an exception and `try`/`except` catches it; PY-5 covers both.)

```python
# file: roots.py
import math


def bisect(f, a, b, tol=1e-12, max_iter=200):
    """Find a root of f in [a, b] by bisection; f(a) and f(b) must differ in sign.

    Returns (midpoint, iterations). The midpoint is within tol of a root.
    """
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        raise ValueError("f(a) and f(b) must have opposite signs")
    for k in range(1, max_iter + 1):
        m = a + (b - a) / 2          # not (a + b) / 2: cannot overflow
        fm = f(m)
        if fm == 0 or (b - a) / 2 <= tol:
            return m, k
        if fa * fm < 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    raise RuntimeError("bisection did not reach the tolerance")


def newton(f, df, x0, tol=1e-12, max_iter=50):
    """Find a root of f near x0 by Newton's method; df is the derivative of f.

    Returns (x, iterations). Stops when the step is at most tol.
    """
    x = x0
    for k in range(1, max_iter + 1):
        step = f(x) / df(x)
        x -= step
        if abs(step) <= tol:
            return x, k
    raise RuntimeError("Newton's method did not converge")


def f(x):
    return math.cos(x) - x


def df(x):
    return -math.sin(x) - 1


if __name__ == "__main__":
    root, k = bisect(f, 0.0, 1.0, tol=1e-6)
    expected = math.ceil(math.log2(1 / 1e-6))
    print(f"bisect: {root!r} in {k} iterations; (3.9) says {expected}")
    try:
        bisect(lambda x: x * x - 2, 1.0, 2.0, tol=1e-20)
    except RuntimeError as err:
        print("tol below the grid spacing:", err)
    x, r = 1.0, 0.7390851332151607
    for n in range(5):
        print(f"Newton n={n}  x={x:.16f}  error={abs(x - r):.2e}")
        x -= f(x) / df(x)
```

```text
bisect: 0.7390851974487305 in 20 iterations; (3.9) says 20
tol below the grid spacing: bisection did not reach the tolerance
Newton n=0  x=1.0000000000000000  error=2.61e-01
Newton n=1  x=0.7503638678402439  error=1.13e-02
Newton n=2  x=0.7391128909113617  error=2.78e-05
Newton n=3  x=0.7390851333852840  error=1.70e-10
Newton n=4  x=0.7390851332151607  error=0.00e+00
```

4. **Newton's method (derivation).** Replace $f$ near the guess $x_n$ by its tangent line, $f(x) \approx f(x_n) + f'(x_n)(x - x_n)$, and take the tangent's root as the next guess:

   $$x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}. \tag{3.10}$$

   Both `f` and its derivative `df` are *arguments*: functions are values, and the root finder does not care what they compute.
5. **Why it is fast (derivation).** Let $r$ be a root with $f'(r) \ne 0$ and $e_n = x_n - r$. Taylor's theorem gives, for some $\xi$ between $x_n$ and $r$, $0 = f(r) = f(x_n) - f'(x_n)e_n + \tfrac12 f''(\xi)e_n^2$. Dividing by $f'(x_n)$: $f(x_n)/f'(x_n) = e_n - f''(\xi)e_n^2/(2f'(x_n))$. Since $e_{n+1} = e_n - f(x_n)/f'(x_n)$,

   $$e_{n+1} = \frac{f''(\xi)}{2f'(x_n)}\,e_n^2 \quad\Longrightarrow\quad |e_{n+1}| \approx \left|\frac{f''(r)}{2f'(r)}\right| |e_n|^2. \tag{3.11}$$

■ **What just happened.** Bisection gains one bit per evaluation and cannot fail once the bracket is valid. Newton squares its error each step, so correct digits roughly *double*: $1.1\times10^{-2}$, $2.8\times10^{-5}$, $1.7\times10^{-10}$, then $0$ (the double nearest the root); bisection needs 20 steps for six digits. The price is a derivative, a decent start and $f'(x_n) \ne 0$. Production solvers (`scipy.optimize.brentq`, PY-8) combine both. The `tol=1e-20` demo shows the cap earning its keep: near $\sqrt2$ the grid spacing is $2.2\times10^{-16}$, so the bracket cannot shrink below it and the cap, not the tolerance, ends the loop.

**Build-it-yourself sequence 5: Horner's rule.**

▶ **Try it first.** Evaluate $p(x) = 2x^3 - 6x^2 + 2x - 1$ at $x = 3$. Then factor $x$ out of the first terms, then out of the first two groups, until no powers remain. Count your multiplications.

1. **The nested form.** For $p(x) = a_nx^n + \cdots + a_1x + a_0$, repeated factoring gives

   $$p(x) = \bigl(\cdots\bigl((a_n x + a_{n-1})x + a_{n-2}\bigr)x + \cdots + a_1\bigr)x + a_0. \tag{3.12}$$

   For the example, $((2x - 6)x + 2)x - 1$.
2. **The count.** Term by term, with $x^k$ built by $k-1$ multiplications and one more for the coefficient, term $k$ costs $k$ multiplications: $\sum_{k=1}^n k = n(n+1)/2$ in all, plus $n$ additions. The nested form uses exactly $n$ multiplications and $n$ additions: $5050$ against $100$ at $n = 100$. The project counts both and asserts these numbers.
3. **The code.**

```python
# file: horner.py
def horner(coeffs, x):
    """Evaluate coeffs = [a_n, ..., a_0] by the nested form (3.12)."""
    result = coeffs[0]
    for a in coeffs[1:]:
        result = result * x + a
    return result


coeffs = [2, -6, 2, -1]
print(horner(coeffs, 3), sum(a * 3 ** (3 - i) for i, a in enumerate(coeffs)))
```

```text
5 5
```

■ **What just happened.** The cost went from quadratic to linear in the degree. Libraries evaluate polynomials this way (`numpy.polyval`), and differentiating (3.12) as you go gives $p'(x)$ almost free, so a polynomial can feed Newton's method.

### 3.13 PEP 8 in ten rules

PEP 8 is the community style guide; it feels like bureaucracy until you read someone else's code.

1. Indent with 4 spaces; never mix tabs and spaces.
2. Keep lines to 79 characters (or one agreed limit, enforced by a formatter).
3. `snake_case` functions and variables, `CapWords` classes, `UPPER_CASE` constants; short formula names (`x`, `dx`, `eps`) are fine.
4. Never reuse a built-in name (`list`, `sum`, `max`, `id`, `input`) for a variable.
5. Two blank lines between top-level definitions, one between methods; spaces around operators and after commas, none around `=` in keyword arguments (`f(x=1)`).
6. Imports at the top, one per line: standard library, third party, yours. Never `from module import *`.
7. `is None`, not `== None`; `if items:`, not `len(items) > 0`.
8. Catch specific exceptions; never a bare `except:` (PY-5).
9. Comments say *why*, not what.
10. Public functions get docstrings.

Formatters and linters (PY-10) enforce most of this.

## 4. Worked examples

### Example 4.1: Reading a four-frame traceback

This program should print the mean of each comma-separated line. It crashes.

```python
# file: report.py
def parse(line):
    values = []
    for token in line.split(","):
        values.append(float(token))
    return values


def mean(values):
    return sum(values) / len(values)


def summarize(lines):
    means = []
    for line in lines:
        means.append(mean(parse(line)))
    return means


def main():
    data = "1.0,2.0,3.0"
    print(summarize(data))


if __name__ == "__main__":
    main()
```

```text
Traceback (most recent call last):
  File "/tmp/py1/report.py", line 25, in <module>
    main()
    ~~~~^^
  File "/tmp/py1/report.py", line 21, in main
    print(summarize(data))
          ~~~~~~~~~^^^^^^
  File "/tmp/py1/report.py", line 15, in summarize
    means.append(mean(parse(line)))
                      ~~~~~^^^^^^
  File "/tmp/py1/report.py", line 4, in parse
    values.append(float(token))
                  ~~~~~^^^^^^^
ValueError: could not convert string to float: '.'
```

Four frames: `<module>`, `main`, `summarize`, `parse`. (Python 3.11+ draws marks under the failing expression; their exact shape varies by version, and your paths differ.) Read from the bottom:

1. **What:** `ValueError: could not convert string to float: '.'`.
2. **Where:** `values.append(float(token))` in `parse`. That line is fine for a token like `"2.0"`; the surprise is a lone `.`, which is no field of `"1.0,2.0,3.0"`.
3. **Walk up.** `summarize` called `parse` with a one-character `line`, so `for line in lines` is iterating over *characters*: `'1'` parsed, then `'.'` failed. `lines` is a string.
4. **Walk up again.** `main` passes `data`, a `str`. The bug is in the second frame from the bottom: the crash is only where the bad value finally hurt. To confirm, `print(repr(line))` in `summarize`.

The fix is in `main`: pass a list of lines. We check it from a second file (the `__main__` guard lets us import):

```python
# file: report_fixed.py
from report import summarize

print(summarize(["1.0,2.0,3.0", "4.5,5.5"]))
```

```text
[2.0, 5.0]
```

### Example 4.2: The aliasing diagram, run

```python
# file: aliasing_diagram.py
import copy

a = [1, 2]
b = a
c = list(a)
print("a is b:", a is b, "| a is c:", a is c, "| a == c:", a == c)
print("elements shared:", a[0] is c[0])
a.append(3)
print("after a.append(3): a =", a, " b =", b, " c =", c)
a = [9]
print("after a = [9]:     a =", a, " b =", b)

m = [[1, 2], [3]]
shallow, deep = list(m), copy.deepcopy(m)
m[0].append(99)
print("m =", m, " shallow =", shallow, " deep =", deep)
```

```text
a is b: True | a is c: False | a == c: True
elements shared: True
after a.append(3): a = [1, 2, 3]  b = [1, 2, 3]  c = [1, 2]
after a = [9]:     a = [9]  b = [1, 2, 3]
m = [[1, 2, 99], [3]]  shallow = [[1, 2, 99], [3]]  deep = [[1, 2], [3]]
```

`a is c` is false though `a == c`: equal values, different objects. Appending through `a` shows through `b`, not `c`; rebinding `a` moves one arrow. The shallow copy shares the inner list, so a change to `m[0]` shows in `shallow` and not in `deep`.

### Example 4.3: Why `0.1 + 0.2 != 0.3`, exactly

```python
# file: point_one.py
import math
from decimal import Decimal
from fractions import Fraction

x, y, z = 0.1, 0.2, 0.3
for name, v in (("0.1", x), ("0.2", y), ("0.3", z), ("0.1+0.2", x + y)):
    print(f"{name:>8} = {Decimal(v)}   {v.hex()}")
print("0.1 as a ratio:", x.as_integer_ratio())

exact = Fraction(x) + Fraction(y)       # exact sum of the two stored values
lo, hi = Fraction(z), Fraction(x + y)   # the stored 0.3 and the computed sum
print("computed sum is one ulp above 0.3:", hi - lo == Fraction(math.ulp(z)))
print("exact sum is exactly halfway between them:", exact == (lo + hi) / 2)
```

```text
     0.1 = 0.1000000000000000055511151231257827021181583404541015625   0x1.999999999999ap-4
     0.2 = 0.200000000000000011102230246251565404236316680908203125   0x1.999999999999ap-3
     0.3 = 0.299999999999999988897769753748434595763683319091796875   0x1.3333333333333p-2
 0.1+0.2 = 0.3000000000000000444089209850062616169452667236328125   0x1.3333333333334p-2
0.1 as a ratio: (3602879701896397, 36028797018963968)
computed sum is one ulp above 0.3: True
exact sum is exactly halfway between them: True
```

The stored `0.1` is slightly above $0.1$, as §3.4 predicted, with denominator $2^{55}$. `0.2` is exactly `2 * 0.1` in binary (same fraction in the hex, exponent one higher), so the exact sum is three times the stored `0.1`: an odd numerator over $2^{55}$. On $[0.25, 0.5)$ the doubles are spaced $2^{-54}$, so that sum lies *exactly halfway* between two doubles, as the program confirms, and round-half-even picks the one ending in `...4`. The stored `0.3` ends in `...3`.

### Example 4.4: A table with f-string specs

```python
# file: table.py
import math
import sys

rows = [
    ("pi", math.pi),
    ("speed of light, m/s", 299_792_458),
    ("Avogadro, 1/mol", 6.02214076e23),
    ("Planck, J s", 6.62607015e-34),
    ("machine epsilon", sys.float_info.epsilon),
]
print(f"{'quantity':<20}{'fixed ,.4f':>38}{'sci .4e':>13}{'general .6g':>14}")
print("-" * 85)
for name, value in rows:
    print(f"{name:<20}{value:>38,.4f}{value:>13.4e}{value:>14.6g}")
```

```text
quantity                                        fixed ,.4f      sci .4e   general .6g
-------------------------------------------------------------------------------------
pi                                                  3.1416   3.1416e+00       3.14159
speed of light, m/s                       299,792,458.0000   2.9979e+08   2.99792e+08
Avogadro, 1/mol       602,214,075,999,999,987,023,872.0000   6.0221e+23   6.02214e+23
Planck, J s                                         0.0000   6.6261e-34   6.62607e-34
machine epsilon                                     0.0000   2.2204e-16   2.22045e-16
```

Fixed notation suits $\pi$ and fails on the tiny constants: Planck's constant prints as `0.0000`. For Avogadro it prints 24 integer digits, of which about 9 mean anything; the rest are the exact decimal expansion of the nearest double. Scientific keeps four digits at every scale; `g` keeps six *significant* digits and picks the notation per value. `<20` left-aligns; the numeric columns right-align.

### Example 4.5: Bisection on $\cos x = x$, with the iteration table

```python
# file: bisect_table.py
import math

from roots import bisect

calls = []


def f(x):
    calls.append(x)             # record every point where bisect evaluates f
    return math.cos(x) - x


root, k = bisect(f, 0.0, 1.0, tol=1e-3)
print(f"{'k':>2} {'midpoint':>10} {'f(midpoint)':>12} {'half-width':>11}")
for i, m in enumerate(calls[2:], start=1):   # calls[0], calls[1] are f(a), f(b)
    print(f"{i:>2} {m:>10.6f} {math.cos(m) - m:>12.6f} {0.5**i:>11.2e}")
expected = math.ceil(math.log2(1 / 1e-3))
print(f"returned {root:.6f} after {k} iterations; formula (3.9): {expected}")
```

```text
 k   midpoint  f(midpoint)  half-width
 1   0.500000     0.377583    5.00e-01
 2   0.750000    -0.018311    2.50e-01
 3   0.625000     0.185963    1.25e-01
 4   0.687500     0.085335    6.25e-02
 5   0.718750     0.033879    3.12e-02
 6   0.734375     0.007875    1.56e-02
 7   0.742188    -0.005196    7.81e-03
 8   0.738281     0.001345    3.91e-03
 9   0.740234    -0.001924    1.95e-03
10   0.739258    -0.000289    9.77e-04
returned 0.739258 after 10 iterations; formula (3.9): 10
```

The sign of `f(midpoint)` tells `bisect` which half to keep, and the half-width halves exactly. The count 10 is $\lceil\log_2 10^3\rceil$ from (3.9). Because `bisect` receives `f` as a value, we could wrap `f` to record its calls.

### Example 4.6: The mutable default, and the `None` sentinel

```python
# file: mutable_default.py
def record(value, log=[]):
    log.append(value)
    return log


def record_fixed(value, log=None):
    if log is None:
        log = []
    log.append(value)
    return log


run1, run2 = record(1.5), record(2.5)
print("buggy:", run1, run2, run1 is run2, record.__defaults__)
run3, run4 = record_fixed(1.5), record_fixed(2.5)
print("fixed:", run3, run4, run3 is run4, record_fixed.__defaults__)
mine = []
print("caller's own list is used:", record_fixed(7.0, mine) is mine)
```

```text
buggy: [1.5, 2.5] [1.5, 2.5] True ([1.5, 2.5],)
fixed: [1.5] [2.5] False (None,)
caller's own list is used: True
```

By R4, `def record` ran once and stored one list in `record.__defaults__`; every call that omits `log` appends to it, so `run1` and `run2` are the same object. The fix evaluates `[]` *inside* the body, once per call. Test `is None`, not `if not log`: the latter would replace a caller's own empty list and never fill it. (The same mechanism used on purpose is a persistent cache; PY-3.)

## 5. Why it matters and where it connects

- **PY-2:** list and dict methods mutate (R2); dictionary keys must be immutable (§3.2 table). **PY-4:** a class-level list shared by all instances is §3.2 snippet 8 in a class.
- **PY-3:** closures capture *names* (the `count` of §3.11); decorators take and return functions, as `bisect` takes `f`.
- **PY-5:** the traceback of §4.1 is an exception propagating up the frames of §2; PY-5 adds `try`/`except` and `pdb`.
- **PY-6:** a NumPy slice is a view, aliasing by design: `b = a[:2]; b[0] = 9` changes `a`, unlike a list slice. `np.finfo(float).eps` is the `epsilon` of §3.4; `np.sum` is the pairwise sum of §3.6.
- **PY-8:** `scipy.optimize.brentq` and `newton` are the production root finders; `math.fsum` is the exact `kahan_sum`.
- **LA-13 §3.2:** the factor (3.6) becomes the condition number of a linear system.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 predict outputs and fix small bugs. B1–B7 are the module; B1 (aliasing audit), B3 (secant method), B5 (series summation with a stopping rule) and B7 (`project.py` extension) are essential. C1–C3: a grid audit, a design question about `x or default`, and a measured study of summation order.

## 7. Coding project

`project.py` is a floating-point lab (standard library, plus NumPy for three comparisons):

1. **Epsilon and the grid:** the halving loop; an `ulp` table; (3.2), (3.3), (3.4) checked against exact `Fraction` arithmetic.
2. **`0.1 + 0.2`:** the `Decimal`, `hex` and tie-break facts of Example 4.3; `isclose`, `inf`, `nan`.
3. **Quadratic formula:** naive against stable (3.7) for $b = 10^2, \ldots, 10^{10}$, errors against a `Decimal` reference (stable $< 10^{-15}$, naive $> 10^{-8}$ at $b = 10^8$); the bound (3.6).
4. **Summation:** Fast2Sum (3.8) on random pairs; naive, Kahan, `fsum` and NumPy on $10^7$ terms of $0.1$.
5. **Root finders:** bisection and Newton on three functions; counts against (3.9); the Newton error ratio against (3.11).
6. **Horner:** operation counts ($n(n+1)/2$ against $n$) and agreement of values.
7. **Report generator:** the table of Example 4.4, plus assertions for the traceback, aliasing, scope and default-argument examples.

**Build it yourself first:** write the epsilon loop, `kahan_sum`, `bisect` and `horner`, then compare. **Your turn** (end of the script): the secant method; a series sum with a stopping rule, and when it stops early; the largest $n$ for which `float(10**n)` is exact.

## 8. Where you'll struggle

- **Knowing `b = a` is not a copy, and writing it anyway.** Most people who have read Batchelder's essay still do this within the month. The fix is procedural: before mutating a list, dict or array you received, ask "who else holds this?", and copy at the boundary if the answer is the caller.
- **Floating point.** The usual response is to sprinkle `round()` or `abs(a - b) < 1e-9` until the test passes, which hides the problem. Ask instead: what is the grid spacing here, what is the scale of what I compare against, does my algorithm subtract nearly equal numbers? The fix is a tolerance taken from the problem, or a rewritten formula.
- **The mutable default works the first time.** `log=[]` passes a test that calls the function once; the bug appears on the second call, perhaps in another test. Make the `None` sentinel a reflex.
- **Reading tracebacks top-down.** The top frame is the oldest and least informative. Read the last line, then the line above, then walk up only until you find the frame whose *arguments* are wrong.
- **PEP 8 feels like bureaucracy** until you read someone else's code. Adopt it from the first file.
- **Loop `else`.** Legal; it means "no `break` happened", which almost no reader remembers. Know it, rarely write it: a helper with an early `return` is clearer.
- **`UnboundLocalError`** surprises everyone once, because the rule is about assignment anywhere in the body, not execution order. A4 drills it.

## 9. Check yourself

1. After `a = [1, 2]; b = a; b = b + [3]`, what is `a`? What if the third statement is `b += [3]`?
2. State R1–R4 in one sentence each. Which explains the mutable default?
3. Why is `-7 // 2` equal to `-4`, and what identity does `%` preserve?
4. How many bits are the sign, exponent and fraction of a double? Why is the relative spacing of doubles between $2^{-53}$ and $2^{-52}$?
5. Why does `1 + 2**-53 == 1`, and what does the loop of sequence 2 find? Why is `0.1` not stored exactly (one sentence about denominators)?
6. When does `math.isclose(x, 0.0)` pass with default arguments? What do you pass instead?
7. What is catastrophic cancellation? Derive the amplification factor (3.6) in three lines, and say how rationalizing removes it from the quadratic formula.
8. What does `c` hold in Kahan summation, and why is it subtracted from the next term?
9. How many iterations does bisection need for tolerance $10^{-9}$ on $[0, 4]$? Why is Newton faster, and when does it fail?
10. Why does `x = x + 1` inside a function raise `UnboundLocalError` when `x` is a module-level name, while `items.append(1)` on a module-level list needs no `global`?

## 10. Pair with

- **Ned Batchelder, "Facts and myths about Python names and values"** (PyCon 2015 talk and essay). The best 30 minutes on §2 and §3.2; it has the pictures that ASCII cannot carry.
- **Downey, *Think Python* 3e, ch. 1-8** (free online). Gentler and slower on variables, conditionals, functions and loops; extra practice for §3.9-3.10.
- **The official Python Tutorial, §3-4** (docs.python.org/3/tutorial). The short authoritative statement of syntax for numbers, strings, `if`, `for`, `range`, `def`.
- **Goldberg, "What every computer scientist should know about floating-point arithmetic," §1-2.** Reference: the rounding material behind §3.4, and the Kahan theorem of §3.6.
- **MIT 6.100A, lectures 1-4** (OCW). A live-coding version of the first half of this module.
