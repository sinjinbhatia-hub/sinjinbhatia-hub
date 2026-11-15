# PY-2 Solutions

Every code block below was run; the ones with `assert`s passed them, and outputs shown are the real ones. Timing-based solutions assert slopes with wide margins, so they pass on any ordinary machine but print numbers that differ from yours.

## A. Warm-up

**A1.** Walk the indices. `s[7:2:-2]` visits 7, 5, 3 (stops before 2). `s[-1:-4:-1]` visits 9, 8, 7. `s[::4]` visits 0, 4, 8. `s[:-7]` stops before index 3. `s[3:3]` is empty. `s[8:100]` clips the stop to 10.

```python
s = "0123456789"
print(s[7:2:-2], s[-1:-4:-1], s[::4], s[:-7], repr(s[3:3]), s[8:100])
lst = list(range(6))
lst[1:4] = []            # delete positions 1, 2, 3
print(lst)
del lst[::2]             # delete positions 0 and 2 of [0, 4, 5]
print(lst)
```
```text
753 987 048 012 '' 89
[0, 4, 5]
[4]
```

**A2.** `b = a[:]` is a shallow copy, so `b` has its own row of arrows but shares the inner list with `a`; `c = a` is the same list as `a`.

```python
a = [1, [2, 3]]
b = a[:]
c = a
b[1].append(4)           # inner list shared by a and b
b[0] = 99                # replaces an arrow in b's own row only
c.append(5)              # c is a: appends to a's row
print(a, b, c)
print(a is c, a[1] is b[1])
```
```text
[1, [2, 3, 4], 5] [99, [2, 3, 4]] [1, [2, 3, 4], 5]
True True
```

**A3.** (a) `list.sort()` returns `None`; assigning it destroys the list and then `None[0]` fails. (b) Removing from a list while iterating shifts later items left, so the loop skips the neighbor of each removed item; build a new list. (c) A generator is consumed by the first `sum`; make a list, or build a fresh generator for each use.

```python
scores = [3, 1, 2]
scores = scores.sort()
print(scores)            # None
try:
    scores[0]
except TypeError as err:
    print("TypeError:", err)

nums = [2, 4, 6, 7, 8]
for x in nums:
    if x % 2 == 0:
        nums.remove(x)
print(nums)              # buggy: skipped items, 4 survived

squares = (x * x for x in range(4))
print(sum(squares), sum(squares))   # 14 0

# Repairs
scores = [3, 1, 2]
scores.sort()
assert scores[0] == 1 and min([3, 1, 2]) == 1
assert [x for x in [2, 4, 6, 7, 8] if x % 2] == [7]
squares = [x * x for x in range(4)]
assert (sum(squares), sum(squares)) == (14, 14)
```
```text
None
TypeError: 'NoneType' object is not subscriptable
[4, 7]
14 0
```

**A4.** The first comprehension reads `for a ... for b ...` as nested loops with `a` outermost. A set display of lengths dedupes. `next(g)` consumes 0, so `list(g)` is `[1, 2]` and the second `list(g)` is empty. `zip` stops at the shorter input.

```python
print([(a, b) for a in "xy" for b in range(2)])
print({len(w) for w in ["aa", "b", "cc"]})
g = (n for n in range(3))
print(next(g), list(g), list(g))
print(list(zip("abc", [1, 2])), list(enumerate("ab", start=5)))
```
```text
[('x', 0), ('x', 1), ('y', 0), ('y', 1)]
{1, 2}
0 [1, 2] []
[('a', 1), ('b', 2)] [(5, 'a'), (6, 'b')]
```

**A5.**

```python
from collections import defaultdict


def invert(d):
    out = defaultdict(list)
    for key, value in d.items():
        out[value].append(key)
    return dict(out)


assert invert({"a": 1, "b": 2, "c": 1}) == {1: ["a", "c"], 2: ["b"]}
assert invert({}) == {}
try:
    invert({"a": [1, 2]})
except TypeError as err:
    print("TypeError:", err)
```
```text
TypeError: unhashable type: 'list'
```

A list value would have to become a dict key, and lists are unhashable (their contents, hence their hash, can change). Convert with `tuple(value)` first if you need such an inversion.

**A6.** `1.0 == 1`, so assigning `d[1]` after `d[1.0]` updates the existing key (which stays `1.0`, the first one inserted) with the new value. `0.1 + 0.2` is not `0.3` in floating point, so the first set has two elements. `1`, `True` and `1.0` are one key. Set equality ignores order. The last line raises.

```python
d = {(1, 2): "tuple"}
d[1.0] = "float"
d[1] = "int"
print(d, len(d))
print(len({0.1 + 0.2, 0.3}), len({1, True, 1.0}), {1, 2} == {2, 1})
d[[1, 2]] = 0
```
```text
{(1, 2): 'tuple', 1.0: 'int'} 2
2 1 True
Traceback (most recent call last):
  File "demo.py", line 6, in <module>
    d[[1, 2]] = 0
    ~^^^^^^^^
TypeError: unhashable type: 'list'
```

## B. Core

**B1.** `yield from` delegates to another iterable. `my_islice` stops right after yielding its last item, without asking for one more, which matters for infinite or shared sources. `my_accumulate` starts from the first item, like `itertools.accumulate`, so it also works for strings and floats. In strict `my_zip`, a short iterator after some items were pulled is a mismatch; if the first one ran out, every other must be exhausted too, so we probe them once.

```python
from itertools import accumulate, chain, count, islice

_END = object()


def my_chain(*iterables):
    for it in iterables:
        yield from it


def my_islice(iterable, stop):
    if stop <= 0:
        return
    for i, item in enumerate(iterable):
        yield item
        if i + 1 >= stop:
            return


def my_accumulate(iterable):
    it = iter(iterable)
    try:
        total = next(it)
    except StopIteration:
        return
    yield total
    for x in it:
        total += x
        yield total


def my_zip(*iterables, strict=False):
    iterators = [iter(it) for it in iterables]
    while iterators:
        row = []
        for it in iterators:
            item = next(it, _END)
            if item is _END:
                break
            row.append(item)
        if len(row) < len(iterators):
            if strict and (row or any(next(it, _END) is not _END for it in iterators[1:])):
                raise ValueError("zip() arguments have different lengths")
            return
        yield tuple(row)


assert list(my_chain([1, 2], "ab", (3.5,))) == list(chain([1, 2], "ab", (3.5,)))
assert list(my_islice("abcdef", 3)) == list(islice("abcdef", 3)) == ["a", "b", "c"]
assert list(my_islice([1, 2], 10)) == [1, 2] and list(my_islice([1, 2], 0)) == []
assert list(my_accumulate([3, 1, 4, 1, 5])) == list(accumulate([3, 1, 4, 1, 5]))
assert list(my_accumulate(["a", "b", "c"])) == ["a", "ab", "abc"] and list(my_accumulate([])) == []
assert list(my_islice(count(), 5)) == [0, 1, 2, 3, 4]
assert list(my_islice(my_accumulate(count(1)), 4)) == [1, 3, 6, 10]
assert list(my_zip("abc", [1, 2])) == list(zip("abc", [1, 2]))
assert list(my_zip("ab", [1, 2], strict=True)) == [("a", 1), ("b", 2)]
for a, b in (("abc", [1, 2]), ("a", [1, 2]), ("", [1])):
    for z in (my_zip, zip):
        try:
            list(z(a, b, strict=True))
        except ValueError:
            pass
        else:
            raise AssertionError((z, a, b))
```

**B2.** Generators are pulled from the right: `islice` asks `running_mean` for an item, which asks `numbers`, which asks `lines`, which reads one line from the file. Nothing asks for a third number, so the file is never read past `2.5`. The `traced` wrapper records exactly what flowed through.

```python
import tempfile
from itertools import islice
from pathlib import Path


def lines(path):
    with Path(path).open(encoding="utf-8") as f:
        for line in f:
            yield line.strip()


def numbers(line_iter):
    for line in line_iter:
        if not line or line.startswith("#"):
            continue
        try:
            yield float(line)
        except ValueError:
            continue


def running_mean(values):
    total = 0.0
    for n, v in enumerate(values, start=1):
        total += v
        yield total / n


def traced(it, log):
    for x in it:
        log.append(x)
        yield x


with tempfile.TemporaryDirectory() as tmp:
    p = Path(tmp) / "x.txt"
    p.write_text("# header\n1.5\n\n2.5\nbad\n3.5\n", encoding="utf-8")
    assert list(numbers(lines(p))) == [1.5, 2.5, 3.5]
    assert list(running_mean(numbers(lines(p)))) == [1.5, 2.0, 2.5]
    seen = []
    first_two = list(islice(running_mean(traced(numbers(lines(p)), seen)), 2))
    assert first_two == [1.5, 2.0] and seen == [1.5, 2.5]      # 3.5 was never parsed
```

**B3.** `DictReader` takes any iterable of lines, so `LOG.splitlines()` works. The `float(...)` conversion fails on `oops`; we count and skip. The `defaultdict` pass needs no sorting; `groupby` needs the sort.

```python
import csv
from collections import defaultdict
from itertools import groupby
from operator import itemgetter

LOG = """t,sensor,value
0.0,a,1.0
0.1,b,10.0
0.2,a,3.0
0.3,b,oops
0.4,a,5.0
0.5,b,20.0
"""

rows, bad = [], 0
for r in csv.DictReader(LOG.splitlines()):
    try:
        rows.append((r["sensor"], float(r["value"])))
    except ValueError:
        bad += 1

count, total, lo, hi = defaultdict(int), defaultdict(float), {}, {}
for name, v in rows:
    count[name] += 1
    total[name] += v
    lo[name] = min(lo.get(name, v), v)
    hi[name] = max(hi.get(name, v), v)
summary = {s: (count[s], total[s] / count[s], lo[s], hi[s]) for s in sorted(count)}

check = {}
for name, grp in groupby(sorted(rows, key=itemgetter(0)), key=itemgetter(0)):
    vals = [v for _, v in grp]
    check[name] = (len(vals), sum(vals) / len(vals), min(vals), max(vals))

assert bad == 1
assert summary == check == {"a": (3, 3.0, 1.0, 5.0), "b": (2, 15.0, 10.0, 20.0)}
```

**B4.** `poly_compose` builds $q^k$ incrementally: for each power $k$ of $p$ in increasing order it multiplies the running power of $q$ by $q$ once per step of the gap, then adds $c_k q^k$. `Fraction` keeps the integral exact so the round trip is an exact equality, not an approximation. (Here $p(q(x)) = 1 + 2(x-1) + 3(x-1)^2 = 2 - 4x + 3x^2$.)

```python
from collections import defaultdict
from fractions import Fraction


def poly_clean(p):
    return {k: c for k, c in p.items() if c != 0}


def poly_add(p, q):
    r = dict(p)
    for k, c in q.items():
        r[k] = r.get(k, 0) + c
    return poly_clean(r)


def poly_mul(p, q):
    r = defaultdict(int)
    for i, a in p.items():
        for j, b in q.items():
            r[i + j] += a * b
    return poly_clean(r)


def poly_eval(p, x):
    return sum(c * x ** k for k, c in p.items())


def poly_deriv(p):
    return {k - 1: k * c for k, c in p.items() if k > 0}


def poly_integ(p):
    return {k + 1: Fraction(c) / (k + 1) for k, c in p.items()}


def poly_compose(p, q):
    result, power, prev = {}, {0: 1}, 0            # power holds q**prev
    for k in sorted(p):
        for _ in range(k - prev):
            power = poly_mul(power, q)
        prev = k
        result = poly_add(result, {e: p[k] * c for e, c in power.items()})
    return result


p, q = {0: 1, 1: 2, 2: 3}, {0: -1, 1: 1}
assert poly_deriv(poly_integ(p)) == p
assert poly_integ({2: 3}) == {3: 1}
assert poly_compose(p, q) == {0: 2, 1: -4, 2: 3}
assert poly_compose(p, {1: 1}) == p and poly_compose({}, q) == {}
for x in range(-3, 4):
    assert poly_eval(poly_compose(p, q), x) == poly_eval(p, poly_eval(q, x))
```

**B5.** Each operation hashes once to pick a bucket and then compares within the chain. With growth disabled and keys that are all multiples of 64, every key lands in bucket 0: one chain of 64 items, and a lookup of the last key does 64 comparisons, which is a list scan. Real dicts avoid this by mixing the higher bits of the hash into the probe sequence, so keys that agree in their low bits do not follow the same path; a table that only looked at the low bits would degrade on such patterns, and Python's integer keys are exactly the kind that have them (multiples of powers of two).

```python
import random


class ToyDict:
    def __init__(self, nbuckets=8, grow=True):
        self.buckets = [[] for _ in range(nbuckets)]
        self.n, self.grow = 0, grow

    def _chain(self, key):
        return self.buckets[hash(key) % len(self.buckets)]

    def put(self, key, value):
        chain = self._chain(key)
        for i, (k, _) in enumerate(chain):
            if k == key:
                chain[i] = (key, value)
                return
        chain.append((key, value))
        self.n += 1
        if self.grow and self.n > len(self.buckets):
            old = [kv for ch in self.buckets for kv in ch]
            self.buckets = [[] for _ in range(2 * len(self.buckets))]
            for k, v in old:
                self._chain(k).append((k, v))

    def get(self, key, default=None):
        for k, v in self._chain(key):
            if k == key:
                return v
        return default

    def delete(self, key):
        chain = self._chain(key)
        for i, (k, _) in enumerate(chain):
            if k == key:
                del chain[i]
                self.n -= 1
                return
        raise KeyError(key)

    def __len__(self):
        return self.n


rng = random.Random(0)
toy, ref = ToyDict(), {}
for _ in range(5000):
    op, key = rng.choice(["put", "put", "get", "delete"]), rng.randrange(300)
    if op == "put":
        toy.put(key, key * 2)
        ref[key] = key * 2
    elif op == "get":
        assert toy.get(key, "none") == ref.get(key, "none")
    elif key in ref:
        toy.delete(key)
        del ref[key]
    else:
        try:
            toy.delete(key)
        except KeyError:
            pass
        else:
            raise AssertionError("delete of a missing key must raise")
    assert len(toy) == len(ref)
assert all(toy.get(k) == v for k, v in ref.items())

fixed = ToyDict(64, grow=False)
for key in range(0, 64 * 64, 64):
    fixed.put(key, key)
assert max(len(ch) for ch in fixed.buckets) == 64 and len(fixed.buckets[0]) == 64
```

**B6.** Three stable passes, least significant key first. The department pass uses `reverse=True`, which still preserves the order among equal departments.

```python
from operator import itemgetter

staff = [("Ana", "eng", 95), ("Raj", "ops", 70), ("Lee", "eng", 95), ("Kim", "ops", 82),
         ("Bo", "eng", 80), ("Zed", "ops", 70)]
out = sorted(staff, key=itemgetter(0))                          # name ascending
out = sorted(out, key=itemgetter(2))                            # salary ascending
out = sorted(out, key=itemgetter(1), reverse=True)              # department descending
for row in out:
    print(row)

for (n1, d1, s1), (n2, d2, s2) in zip(out, out[1:]):
    assert d1 > d2 or (d1 == d2 and (s1 < s2 or (s1 == s2 and n1 <= n2))), (n1, n2)
```
```text
('Raj', 'ops', 70)
('Zed', 'ops', 70)
('Kim', 'ops', 82)
('Bo', 'eng', 80)
('Ana', 'eng', 95)
('Lee', 'eng', 95)
```

**B7.** Add the two entries to `EXPERIMENTS`, extend the asserts, and give the figure a third panel. The three edits to `project.py` Part 2 are:

1. After the `"linked: prepend"` entry add
   `"list queue: append+pop(0)": (lambda n: list(range(n)), lambda c: (c.append(0), c.pop(0)), lambda n: 100),`
   and
   `"deque queue: append+popleft": (lambda n: deque(range(n)), lambda c: (c.append(0), c.popleft()), lambda n: 100),`
2. After the existing asserts add
   `assert 0.7 < slopes["list queue: append+pop(0)"] < 1.5 and abs(slopes["deque queue: append+popleft"]) < 0.3`
3. Change `plt.subplots(1, 2, figsize=(10, 4))` to `plt.subplots(1, 3, figsize=(14, 4))`, make the second panel's keys `list(EXPERIMENTS)[2:5]`, and add a third tuple `(axes[2], "FIFO queue", list(EXPERIMENTS)[5:])` to the panel loop.

The block below (run it from this module's folder) applies those three edits to a copy of `project.py` and runs the copy, so you can see that they work and what the new slopes are.

```python
import re
import subprocess
import sys
import tempfile
from pathlib import Path

src = Path("project.py").read_text(encoding="utf-8")
edits = [
    ('    "linked: prepend": (linked, lambda h: (0, h), lambda n: 100),\n}',
     '    "linked: prepend": (linked, lambda h: (0, h), lambda n: 100),\n'
     '    "list queue: append+pop(0)": (lambda n: list(range(n)), lambda c: (c.append(0), c.pop(0)), lambda n: 100),\n'
     '    "deque queue: append+popleft": (lambda n: deque(range(n)), lambda c: (c.append(0), c.popleft()), lambda n: 100),\n}'),
    ('assert times["list: -1 in c"][-1] > 100 * times["set: -1 in c"][-1]\n',
     'assert times["list: -1 in c"][-1] > 100 * times["set: -1 in c"][-1]\n'
     'assert 0.7 < slopes["list queue: append+pop(0)"] < 1.5 and abs(slopes["deque queue: append+popleft"]) < 0.3\n'),
    ('plt.subplots(1, 2, figsize=(10, 4))', 'plt.subplots(1, 3, figsize=(14, 4))'),
    ('"Insert at the front, then remove it", list(EXPERIMENTS)[2:])):',
     '"Insert at the front, then remove it", list(EXPERIMENTS)[2:5]),\n'
     '                        (axes[2], "FIFO queue", list(EXPERIMENTS)[5:])):'),
]
for old, new in edits:
    assert src.count(old) == 1, old
    src = src.replace(old, new)
with tempfile.TemporaryDirectory() as tmp:
    copy = Path(tmp) / "project.py"
    copy.write_text(src, encoding="utf-8")
    run = subprocess.run([sys.executable, str(copy)], capture_output=True, text=True, cwd=tmp,
                         env={"MPLBACKEND": "Agg", "PATH": "/usr/bin:/bin"})
    assert run.returncode == 0, run.stderr
    print("\n".join(line for line in run.stdout.splitlines() if "queue" in line or "All checks" in line))
```
```text
  list queue: append+pop(0)    slope  1.11   time at n=10^6:       438355 ns
  deque queue: append+popleft  slope  0.02   time at n=10^6:           99 ns
All checks passed.
```

The list queue's slope is near 1 because `pop(0)` shifts every remaining arrow; the deque's is near 0. (Fitted over $n = 10^3 \dots 10^6$ the list exponent often lands between 1.1 and 1.4, for the cache reason given in §3.1, which is why the bound is 1.5 rather than 1.3.) This is the reason `collections.deque` exists.

**B8.** Work for `freq_lists` is about (number of words) × (average position of a hit) ≈ $n V/2$, so time should grow like $V$ (slope 1); `freq_dict` does one hash lookup per word, so slope 0. Measured slopes land near 1 for the list version (a little under it when $V$ is small, where the fixed per-word interpreter overhead is comparable to the search) and near 0 for the dict; the dict's slope comes out slightly negative or slightly positive from run to run ($-0.05$ below, $+0.07$ on a rerun), which is timing noise on a cost that does not depend on $V$.

```python
import math
import random
import statistics
import time


def freq_lists(words):
    seen, counts = [], []
    for w in words:
        for i, s in enumerate(seen):
            if s == w:
                counts[i] += 1
                break
        else:
            seen.append(w)
            counts.append(1)
    return dict(zip(seen, counts))


def freq_dict(words):
    counts = {}
    for w in words:
        counts[w] = counts.get(w, 0) + 1
    return counts


def best_time(fn, words, repeat=3):
    best = math.inf
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn(words)
        best = min(best, time.perf_counter() - t0)
    return best


vs = [50, 100, 200, 400, 800]
rng = random.Random(0)
t_list, t_dict = [], []
for v in vs:
    vocab = [f"w{i}" for i in range(v)]
    words = rng.choices(vocab, k=20_000)
    assert freq_lists(words) == freq_dict(words)
    t_list.append(best_time(freq_lists, words))
    t_dict.append(best_time(freq_dict, words))
logv = [math.log(v) for v in vs]
slope_list = statistics.linear_regression(logv, [math.log(t) for t in t_list]).slope
slope_dict = statistics.linear_regression(logv, [math.log(t) for t in t_dict]).slope
print(f"slope vs vocabulary size: freq_lists {slope_list:.2f}, freq_dict {slope_dict:.2f}")
assert 0.7 < slope_list < 1.2 and slope_dict < 0.3
```
```text
slope vs vocabulary size: freq_lists 1.01, freq_dict -0.05
```


## C. Challenge

**C1.** `delete` that empties the slot cuts the probe path: a later key that collided with the deleted one is behind an empty slot, and `get` stops at the first empty slot. The tombstone keeps the path intact; `get` steps over it, and `put` may reuse the first tombstone it passes (after confirming the key is not further along the path). Because tombstones count as occupied, `used` (live plus tombstones) drives the resize, and a resize drops them. Measured probes follow the classical linear-probing estimate $\tfrac12\bigl(1 + \tfrac{1}{1-\alpha}\bigr)$ for a successful search, and the curve bends sharply upward as $\alpha \to 1$, which is why tables resize at about two-thirds full.

```python
import random

_EMPTY, _TOMB = object(), object()


class OpenTable:
    """Linear probing in one flat list of slots; deletions leave tombstones."""

    def __init__(self, nslots=8, tombstones=True, max_load=2 / 3):
        self.slots = [_EMPTY] * nslots
        self.live = self.used = self.probes = 0          # used = live + tombstones
        self.tombstones, self.max_load = tombstones, max_load

    def _path(self, key):
        n = len(self.slots)
        start = hash(key) % n
        return ((start + step) % n for step in range(n))

    def get(self, key, default=None):
        for i in self._path(key):
            slot = self.slots[i]
            self.probes += 1
            if slot is _EMPTY:
                return default
            if slot is not _TOMB and slot[0] == key:
                return slot[1]
        return default

    def put(self, key, value):
        free = None                                      # first tombstone on the path
        for i in self._path(key):                        # an EMPTY slot always exists (load < 1)
            slot = self.slots[i]
            if slot is _EMPTY:
                break
            if slot is _TOMB:
                free = i if free is None else free
            elif slot[0] == key:
                self.slots[i] = (key, value)
                return
        target = i if free is None else free
        if self.slots[target] is _EMPTY:
            self.used += 1
        self.slots[target] = (key, value)
        self.live += 1
        if self.used > self.max_load * len(self.slots):
            self._resize()

    def delete(self, key):
        for i in self._path(key):
            slot = self.slots[i]
            if slot is _EMPTY:
                break
            if slot is not _TOMB and slot[0] == key:
                self.slots[i] = _TOMB if self.tombstones else _EMPTY
                self.live -= 1
                if not self.tombstones:
                    self.used -= 1
                return
        raise KeyError(key)

    def _resize(self):
        items = [s for s in self.slots if s is not _EMPTY and s is not _TOMB]
        n = len(self.slots) * (2 if self.live > len(self.slots) // 3 else 1)
        self.slots, self.live, self.used = [_EMPTY] * n, 0, 0
        for k, v in items:
            self.put(k, v)

    def __len__(self):
        return self.live


for flag in (False, True):                               # 3 and 11 both start at slot 3 of 8
    t = OpenTable(tombstones=flag)
    t.put(3, "a")
    t.put(11, "b")
    t.delete(3)
    print(f"tombstones={flag}: get(11) ->", t.get(11, "LOST"))

rng = random.Random(0)
table, ref = OpenTable(), {}
for _ in range(10_000):
    op, key = rng.choice(["put", "put", "get", "delete"]), rng.randrange(400)
    if op == "put":
        table.put(key, -key)
        ref[key] = -key
    elif op == "get":
        assert table.get(key, "none") == ref.get(key, "none")
    elif key in ref:
        table.delete(key)
        del ref[key]
    else:
        try:
            table.delete(key)
        except KeyError:
            pass
        else:
            raise AssertionError("delete of a missing key must raise")
    assert len(table) == len(ref)
assert all(table.get(k) == v for k, v in ref.items())

for alpha in (0.2, 0.4, 0.6, 0.8):
    size = 4096
    tab = OpenTable(size, max_load=0.95)                 # no resizing during the experiment
    keys = random.Random(1).sample(range(10**9), int(alpha * size))
    for k in keys:
        tab.put(k, 1)
    tab.probes = 0
    for k in keys:
        tab.get(k)
    avg, theory = tab.probes / len(keys), 0.5 * (1 + 1 / (1 - alpha))
    print(f"load {alpha}: {avg:.3f} probes per get (classical estimate {theory:.3f})")
    assert abs(avg - theory) < 0.2 * theory
```
```text
tombstones=False: get(11) -> LOST
tombstones=True: get(11) -> b
load 0.2: 1.132 probes per get (classical estimate 1.125)
load 0.4: 1.296 probes per get (classical estimate 1.333)
load 0.6: 1.740 probes per get (classical estimate 1.750)
load 0.8: 3.197 probes per get (classical estimate 3.000)
```

**C2.** All three caches are checked against each other first. The measured picture: the list version's cost per access grows steeply with $k$ because `remove` and `pop(0)` shift or scan $O(k)$ items (§3.1 table). The plain-dict version is cheap per hit (a `pop` and a re-insert are $O(1)$), but its eviction `next(iter(d))` must skip the holes left by earlier evictions at the front of the entries array, so its cost per access also creeps up with $k$ until the dict rebuilds itself. `OrderedDict` keeps its entries in a linked structure, so `move_to_end` and `popitem(last=False)` are $O(1)$ whatever $k$. The second block isolates the effect: delete the oldest key and add a new one, repeatedly.

```python
import math
import random
import statistics
import time
from collections import OrderedDict


class LRUDict:
    def __init__(self, capacity):
        self.capacity, self.data = capacity, {}

    def get(self, key, default=None):
        if key not in self.data:
            return default
        value = self.data.pop(key)                       # re-insert: now the newest key
        self.data[key] = value
        return value

    def put(self, key, value):
        if key in self.data:
            del self.data[key]
        elif len(self.data) >= self.capacity:
            del self.data[next(iter(self.data))]         # the first key is the least recently used
        self.data[key] = value


class LRUList:
    def __init__(self, capacity):
        self.capacity, self.order, self.values = capacity, [], {}

    def get(self, key, default=None):
        if key not in self.values:
            return default
        self.order.remove(key)                           # O(k) scan
        self.order.append(key)
        return self.values[key]

    def put(self, key, value):
        if key in self.values:
            self.order.remove(key)
        elif len(self.values) >= self.capacity:
            del self.values[self.order.pop(0)]           # O(k) shift
        self.order.append(key)
        self.values[key] = value


class LRUOrdered:
    def __init__(self, capacity):
        self.capacity, self.data = capacity, OrderedDict()

    def get(self, key, default=None):
        if key not in self.data:
            return default
        self.data.move_to_end(key)
        return self.data[key]

    def put(self, key, value):
        if key in self.data:
            self.data.move_to_end(key)
        elif len(self.data) >= self.capacity:
            self.data.popitem(last=False)
        self.data[key] = value


def run(cache, ops):
    hits = 0
    for key in ops:
        if cache.get(key) is None:
            cache.put(key, key)
        else:
            hits += 1
    return hits


def slope(ks, ts):
    return statistics.linear_regression([math.log(k) for k in ks], [math.log(t) for t in ts]).slope


r = random.Random(0)
small = [r.randrange(30) for _ in range(5000)]
hits = [run(c(10), small) for c in (LRUDict, LRUList, LRUOrdered)]
assert len(set(hits)) == 1                                # identical behavior

ks, kinds = [100, 1000, 10_000], (LRUDict, LRUList, LRUOrdered)
per_access = {c.__name__: [] for c in kinds}
for k in ks:
    r = random.Random(k)
    ops = [r.randrange(3 * k) for _ in range(60_000)]       # about 2/3 of accesses miss
    for c in kinds:
        elapsed = []
        for _ in range(2):
            cache, t0 = c(k), time.perf_counter()
            run(cache, ops)
            elapsed.append(time.perf_counter() - t0)
        per_access[c.__name__].append(min(elapsed) / len(ops))
for name, ts in per_access.items():
    print(f"{name:11s} " + "  ".join(f"{t * 1e6:6.2f} us" for t in ts) + f"   slope {slope(ks, ts):.2f}")
assert slope(ks, per_access["LRUOrdered"]) < 0.3
assert slope(ks, per_access["LRUList"]) > 0.5
assert slope(ks, per_access["LRUDict"]) > slope(ks, per_access["LRUOrdered"]) + 0.1

# Isolate the eviction idiom: repeatedly drop the oldest key and add a new one.
ks2, plain, ordered = [1_000, 10_000, 100_000], [], []
for k in ks2:
    d, od, new = dict.fromkeys(range(k)), OrderedDict.fromkeys(range(k)), k
    t0 = time.perf_counter()
    for _ in range(50_000):
        del d[next(iter(d))]
        d[new] = None
        new += 1
    plain.append((time.perf_counter() - t0) / 50_000)
    new, t0 = k, time.perf_counter()
    for _ in range(50_000):
        od.popitem(last=False)
        od[new] = None
        new += 1
    ordered.append((time.perf_counter() - t0) / 50_000)
print("evict-oldest, plain dict:", "  ".join(f"{t * 1e6:6.2f} us" for t in plain), f"  slope {slope(ks2, plain):.2f}")
print("evict-oldest, OrderedDict:", "  ".join(f"{t * 1e6:6.2f} us" for t in ordered), f"  slope {slope(ks2, ordered):.2f}")
assert slope(ks2, plain) > 0.3 and slope(ks2, ordered) < 0.35
```
```text
LRUDict       0.40 us    0.61 us    2.07 us   slope 0.36
LRUList       0.65 us    2.56 us   18.96 us   slope 0.73
LRUOrdered    0.39 us    0.42 us    0.54 us   slope 0.07
evict-oldest, plain dict:   0.87 us    4.31 us   18.68 us   slope 0.67
evict-oldest, OrderedDict:   0.33 us    0.36 us    0.34 us   slope 0.01
```

So "the dict remembers insertion order" gives you an LRU quickly, but not an $O(1)$ one in CPython; when the eviction path matters, use `OrderedDict`. Check this kind of claim with a measurement, as here.

**C3.** (a) The double loop does $n(n-1)/2$ pair checks, slope about 2; the dict version does one lookup and one insertion per element, slope about 1. (b) A single `-1 in xs` is one scan; building `set(xs)` costs about as much as several scans, so converting only pays off after a few lookups. The break-even $q$ is a small constant, which is the honest reading of big-$O$: it says what wins eventually and how the cost grows, not where the crossover is. For one lookup the list wins; for $q$ lookups it is $q \cdot O(n)$ against $O(n) + q \cdot O(1)$.

```python
import math
import random
import statistics
import time


def two_sum_loop(xs, target):
    for i in range(len(xs)):
        for j in range(i + 1, len(xs)):
            if xs[i] + xs[j] == target:
                return i, j
    return None


def two_sum_dict(xs, target):
    seen = {}
    for j, x in enumerate(xs):
        if target - x in seen:
            return seen[target - x], j
        seen[x] = j
    return None


rng = random.Random(0)
for _ in range(200):
    xs = rng.sample(range(-50, 50), rng.randrange(2, 40))
    target = rng.randrange(-60, 60)
    a, b = two_sum_loop(xs, target), two_sum_dict(xs, target)
    assert (a is None) == (b is None)
    for pair in (a, b):
        if pair is not None:
            i, j = pair
            assert i != j and xs[i] + xs[j] == target


def best(fn, *args, repeat=3):
    out = math.inf
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn(*args)
        out = min(out, time.perf_counter() - t0)
    return out


ns = [500, 1000, 2000, 4000]
t_loop = [best(two_sum_loop, list(range(n)), -1) for n in ns]
t_dict = [best(two_sum_dict, list(range(n)), -1) for n in ns]
logn = [math.log(n) for n in ns]
s_loop = statistics.linear_regression(logn, [math.log(t) for t in t_loop]).slope
s_dict = statistics.linear_regression(logn, [math.log(t) for t in t_dict]).slope
print(f"slopes: double loop {s_loop:.2f}, dict {s_dict:.2f}")
assert 1.6 < s_loop < 2.2 and 0.8 < s_dict < 1.4

xs = list(range(1000))
scan = best(lambda: [-1 in xs for _ in range(200)]) / 200      # one list scan
build = best(lambda: [set(xs) for _ in range(200)]) / 200      # one set construction
q = math.ceil(build / scan)
print(f"scan {scan * 1e6:.1f} us, build {build * 1e6:.1f} us, break-even at q = {q} lookups")
assert 1 <= q < 50
```
```text
slopes: double loop 1.93, dict 1.11
scan 7.9 us, build 34.2 us, break-even at q = 5 lookups
```
