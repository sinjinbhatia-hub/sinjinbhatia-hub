# PY-2: Collections, Iteration, and the Data Model

**Time:** ~9 hours · **Prerequisites:** PY-1 · **Unlocks:** PY-3, PY-6, PY-9

> **Conventions in this module.** Python 3.11+; every output below was produced by a real run under Python 3.13. `>>>` blocks are interactive transcripts; other blocks are scripts followed by their real output. $n$ is the number of elements in a container. "$O(1)$" and "$O(n)$" describe the cost of one operation as $n$ grows; for dicts and sets "$O(1)$" means *average case* (§3.4 says when that fails), and for `list.append` it means *amortized*. Timings come from one machine in one run: trust the growth and the ratios, not the nanoseconds. Containers are named by their literal syntax: list `[]`, tuple `()`, dict `{k: v}`, set `{x}`. Standard library only; NumPy appears only in the last cell of `project.py`.

## 0. Where this sits

PY-1 gave you names, objects, mutability and floats. This module is about the objects that hold other objects, and about the one protocol, `for`, that walks them. Everything in PY-6 (arrays) and PY-9 (tables) assumes it: NumPy slicing is the slicing of §3.1 with more axes, pandas `groupby` is §4.4 with a nicer face, and every `for` loop in the rest of the curriculum runs the machinery of §3.7. PY-3 turns the `key=` function of §3.1 into a general tool, and PY-4 shows how a class opts into the protocols of §3.7 (the "data model").

## 1. The problem that forced this idea

A numerical program is mostly bookkeeping: collections of numbers, lookups, grouping, filtering, reading a file of measurements. The container you pick decides whether the program is three clear lines or a mess of index arithmetic that does quadratic work.

The forcing problem: ten thousand readings, each tagged with one of 500 sensor names, and you need the mean per sensor. With lists you keep parallel lists of names and sums and search for the name each time. With a dict you never search.

```python
import random

rng = random.Random(0)
names = [f"sensor_{i:03d}" for i in range(500)]
readings = [(rng.choice(names), rng.random()) for _ in range(10_000)]

comparisons, seen, sums, counts = 0, [], [], []       # version 1: parallel lists, linear search
for name, value in readings:
    for i in range(len(seen)):
        comparisons += 1
        if seen[i] == name:
            sums[i] += value
            counts[i] += 1
            break
    else:                                              # not found: first time we see this sensor
        seen.append(name)
        sums.append(value)
        counts.append(1)
means_lists = {n: s / c for n, s, c in zip(seen, sums, counts)}

totals = {}                                            # version 2: dict from name to (sum, count)
for name, value in readings:
    s, c = totals.get(name, (0.0, 0))
    totals[name] = (s + value, c + 1)
means_dict = {n: s / c for n, (s, c) in totals.items()}

print(f"list version: {comparisons:,} string comparisons for {len(readings):,} readings")
print(f"same answer from the dict version: {means_lists == means_dict}")
```
```text
list version: 2,448,233 string comparisons for 10,000 readings
same answer from the dict version: True
```

The list version does about $nk/2$ comparisons for $n$ readings and $k$ sensors. The dict version does one hash lookup per reading, whatever $k$ is. That is the case for hash tables: **replace a search with arithmetic on the key.**

The second forcing problem is size. A log with 10⁸ lines does not fit in a list, but a loop over it needs only one line at a time. Iterators and generators are the language feature for "one at a time, on demand"; §3.7–3.8 build them from scratch.

## 2. Intuition

**A list is a row of arrows.** The list object stores references (arrows) to its items, in consecutive memory, plus a length and some spare room:

```
lst = [10, 'a', 3.5]

 lst ──► [ length 3 | spare 1 ]        objects live elsewhere
         [ ref 0 ] ───────────────────► 10
         [ ref 1 ] ───────────────────► 'a'
         [ ref 2 ] ───────────────────► 3.5
         [ spare ]
```

Item $i$ is at address start + 8·$i$ bytes, so `lst[i]` is $O(1)$. Appending writes into the spare slot (occasionally the row moves to a bigger block). Inserting at the front shifts every arrow one place right: $O(n)$. `x in lst` follows each arrow and compares: $O(n)$.

**A dict is a table addressed by a fingerprint of the key.**

```
key ──hash()──► integer ──reduce to table size──► slot ──► (key, value)
```

Compute where the key must be and look there. An empty slot means absent; an equal key means found; a different key (a collision) means try the next slot in a fixed, key-determined sequence. The work does not grow with the number of items as long as the table stays sparse. The price: a key's fingerprint must never change while it sits in the table, so keys must be immutable.

**A set is a dict with keys and no values:** the same idea, with the same average costs.

**An iterator is a bookmark.** It supports one operation: "give me the next item, or say there are none." A list is *iterable*: it can hand you a fresh bookmark (`iter(lst)`). A **generator** is a function that can pause; it computes the next item when asked, remembers where it stopped, and stores nothing else. That is how a generator can "contain" a billion items.

Guess before reading on: which of `lst[i]`, `lst.append(x)`, `lst.insert(0, x)`, `x in lst`, `x in d`, `x in s` are $O(1)$? §3.1 and §3.4 check you.

## 3. Formal development

### 3.1 Lists, slicing, and sorting

**Indexing and slicing.** For a sequence of length $n$, index $i$ means position $i$ if $0 \le i < n$ and $n + i$ if $-n \le i < 0$; anything else raises `IndexError`. The slice `s[i:j:k]` selects positions $i, i+k, i+2k, \dots$ for as long as they stay strictly before $j$ (if $k > 0$) or strictly after $j$ (if $k < 0$). Omitted endpoints take defaults that **depend on the sign of $k$**: for $k > 0$, $i = 0$ and $j = n$; for $k < 0$, $i = n - 1$ and $j$ = "immediately before position 0". Out-of-range endpoints are clipped, so slices never raise `IndexError`. The method `slice(i, j, k).indices(n)` returns the resolved triple, which is the quickest way to see what a slice means.

▶ **Try it first.** Predict the output of each line for `s = "abcdefghij"` (indices 0 to 9). The outputs below appear in order, one per expression line (tracebacks included); assignments print nothing.

```python
s = "abcdefghij"
s[2:7]
s[-3:]
s[1:8:3]
s[5:2:-1]
s[::-1]
s[2:5:-1]
slice(None, None, -1).indices(10)
```

■ **What just happened.**

```text
'cdefg'
'hij'
'beh'
'fed'
'jihgfedcba'
''
(9, -1, -1)
```

`s[5:2:-1]` walks from position 5 down to, but not including, position 2: `fed`. `s[2:5:-1]` is empty because walking *down* from 2 never reaches "after 5". The last line shows the rule: with a negative stride the default start is 9 and the default stop is $-1$, which here means "before position 0", not "the last element". Example 4.1 draws the picture.

Slices also appear on the left of `=`. **Slice assignment** replaces a stretch of the list with the items of any iterable, and the lengths need not match:

```pycon
>>> lst = list(range(10))
>>> lst[2:5] = ["x"]
>>> lst
[0, 1, 'x', 5, 6, 7, 8, 9]
>>> lst[1:1] = [7, 7]
>>> lst
[0, 7, 7, 1, 'x', 5, 6, 7, 8, 9]
>>> del lst[::2]
>>> lst
[7, 1, 5, 7, 9]
```

**Methods and what they cost.** The costs follow from the picture in §2.

| Operation | Cost | Why |
|---|---|---|
| `lst[i]`, `lst[i] = x`, `len(lst)` | $O(1)$ | address arithmetic; length is stored |
| `lst.append(x)`, `lst.pop()` | $O(1)$ amortized | write into spare room; occasionally move the row |
| `lst.insert(i, x)`, `lst.pop(i)`, `del lst[i]` | $O(n - i)$ | shift every later arrow; `insert(0, x)` is $O(n)$ |
| `x in lst`, `lst.index(x)`, `lst.count(x)`, `lst.remove(x)` | $O(n)$ | scan, comparing with `==` |
| `lst[i:j]` | $O(j - i)$ | copies the arrows (not the objects) |
| `lst.extend(other)`, `lst + other` | $O(\text{len(other)})$, $O(n + m)$ | `+` builds a new list |
| `lst.sort()`, `sorted(lst)` | $O(n \log n)$ | Timsort; close to $O(n)$ if already nearly sorted |
| `lst.reverse()` | $O(n)$ | |

"Amortized $O(1)$" for `append` rests on over-allocation: when the spare room runs out, CPython allocates a block somewhat larger (about 12% more, plus a few slots) and copies the arrows once. Because each reallocation buys room for a *proportional* number of future appends, the average cost per append stays constant. The growth factor is an implementation detail; any geometric factor gives the same conclusion.

```python
import sys

lst = []
last, resizes = sys.getsizeof(lst), []
for n in range(100):
    lst.append(n)
    if sys.getsizeof(lst) != last:
        last = sys.getsizeof(lst)
        resizes.append(len(lst))
print("the list moved to a bigger block when its length reached", resizes)
```
```text
the list moved to a bigger block when its length reached [1, 5, 9, 17, 25, 33, 41, 53, 65, 77, 93]
```

Eleven reallocations in 100 appends, and the gaps between them widen: that is geometric growth.

**Sorting.** `lst.sort()` sorts in place and returns `None`; `sorted(iterable)` returns a new list from any iterable. Both take `key=` (a function applied once per element; elements are compared by their keys) and `reverse=True`, and both are **stable**: equal keys keep their input order, even with `reverse=True`. Stability makes multi-key sorting work: sort by the secondary key first and the primary second, or sort once with a tuple key `(primary, secondary)`, since tuples compare element by element. Example 4.3 does both.

**`collections.deque`** is the container for "both ends": `appendleft`, `append`, `popleft`, `pop` are all $O(1)$, indexing the middle is $O(n)$, and `deque(maxlen=k)` discards from the opposite end automatically (a sliding window).

```pycon
>>> from collections import deque
>>> window = deque(range(6), maxlen=3)
>>> window
deque([3, 4, 5], maxlen=3)
>>> window.appendleft(9)
>>> window
deque([9, 3, 4], maxlen=3)
```

▶ **Build it yourself (1): a linked list, and why the front is cheap.** A linked list stores each item together with a reference to the *rest* of the list. Put an item in front by making a new cell that points at the old head. Nothing else moves. Predict first: if you prepend 10⁶ items to a linked list and to a Python list, which one slows down as it grows?

1. Represent the empty list as `None` and a cell as a pair `(value, rest)`. Build `[1, 2, 3]` as `(1, (2, (3, None)))`.
2. Prepend is `head = (x, head)`: one allocation, independent of length. Because tuples are immutable, two lists can *share* a tail safely.
3. Walking is the price: to reach item $i$ you follow $i$ links, so indexing is $O(i)$.
4. Make the cells mutable (dicts) and you can splice into the middle in $O(1)$ once you hold the cell: `node["next"] = {"value": x, "next": node["next"]}`.

```python
def to_list(head):
    out = []
    while head is not None:
        value, head = head
        out.append(value)
    return out


def nth(head, i):
    for _ in range(i):
        head = head[1]
    return head[0]


tail = (2, (3, None))
a, b = (1, tail), (9, tail)
print(to_list(a), to_list(b), "shared tail:", a[1] is b[1])
print("item 2 of a:", nth(a, 2))

# Mutable cells: splice 'x' in after the first cell, in O(1).
node = {"value": 1, "next": {"value": 3, "next": None}}
node["next"] = {"value": "x", "next": node["next"]}
cell, vals = node, []
while cell is not None:
    vals.append(cell["value"])
    cell = cell["next"]
print(vals)
```
```text
[1, 2, 3] [9, 2, 3] shared tail: True
item 2 of a: 3
[1, 'x', 3]
```

Now the measurement. For each size $n$ we time one *insert at the front followed by a removal from the front*, so the container stays at size $n$.

```python
from collections import deque
from timeit import repeat

print(f"{'n':>9} {'list':>9} {'deque':>9} {'linked':>9}   microseconds per front insert + remove")
for n in (10**3, 10**4, 10**5, 10**6):
    head = None
    for i in range(n):
        head = (i, head)
    t = {}
    for name, stmt, c in (("list", "c.insert(0, 0); c.pop(0)", list(range(n))),
                          ("deque", "c.appendleft(0); c.popleft()", deque(range(n))),
                          ("linked", "(0, c)", head)):
        t[name] = min(repeat(stmt, globals={"c": c}, number=100, repeat=5)) / 100 * 1e6
    print(f"{n:>9,} {t['list']:>9.3f} {t['deque']:>9.3f} {t['linked']:>9.3f}")
```

■ **What just happened.**

```text
        n      list     deque    linked   microseconds per front insert + remove
    1,000     0.359     0.046     0.031
   10,000     3.362     0.046     0.032
  100,000    42.372     0.047     0.031
1,000,000   711.136     0.047     0.031
```

The list column grows in proportion to $n$ (the arrows really are being shifted); the deque and linked columns stay flat. A deque is a linked list of fixed-size blocks, which is why both ends are cheap. `project.py` Part 2 fits the slopes on a log-log plot. In real code reach for `deque`; build linked lists to understand why.

### 3.2 Tuples

A tuple is an immutable sequence: its length and its references cannot change. It is the right type for a fixed record, a dict key, and "return several values". **Packing** builds a tuple (`point = 3, 4`) and **unpacking** takes one apart (`x, y = point`). The right-hand side is evaluated completely before any name on the left is bound, which is why the swap idiom `a, b = b, a` works. A starred name collects the rest:

```pycon
>>> first, *middle, last = [1, 2, 3, 4, 5]
>>> first, middle, last
(1, [2, 3, 4], 5)
>>> (x, y), z = (1, 2), 3
>>> x + y + z
6
>>> a, b = 1, 2
>>> a, b = b, a
>>> a, b
(2, 1)
>>> from collections import namedtuple
>>> Point = namedtuple("Point", ["x", "y"])
>>> p = Point(3, 4)
>>> p.x, p[1], p == (3, 4), hash(p) == hash((3, 4))
(3, 4, True, True)
```

A `namedtuple` is a tuple with field names. Immutability is **shallow**: it freezes which objects the tuple refers to, not their contents.

▶ **Try it first.** Predict what happens to `t` and whether the last line raises an error.

```python
t = ([1, 2], 3)
t[0].append(99)
t
t[0] += [4]
t
```

■ **What just happened.**

```text
([1, 2, 99], 3)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: 'tuple' object does not support item assignment
([1, 2, 99, 4], 3)
```

`t[0] += [4]` first mutates the list in place (this succeeds), then tries to assign the result back into `t[0]` (this fails). The error and the mutation both happen. Mutable objects inside tuples are legal and usually a design smell.

### 3.3 Copying: shallow and deep

PY-1's rule still holds: assignment never copies. To copy a list you ask explicitly. All of `list(a)`, `a[:]`, `a.copy()`, and `copy.copy(a)` make a **shallow copy**: a new row of arrows pointing at the *same* objects. For flat lists of numbers or strings (immutable items) that is a full copy. For nested containers it is not.

▶ **Try it first.** Predict the printed `a` and `b`.

```python
a = [[1, 2], [3, 4]]
b = list(a)
b[0].append(99)
b.append([5])
print("a =", a)
print("b =", b)
print("rows shared:", a[0] is b[0], " lists shared:", a is b)
```

■ **What just happened.**

```text
a = [[1, 2, 99], [3, 4]]
b = [[1, 2, 99], [3, 4], [5]]
rows shared: True  lists shared: False
```

`b.append([5])` changed only `b`'s row of arrows. `b[0].append(99)` followed an arrow that `a` also holds, so it changed what both see. `copy.deepcopy(a)` copies recursively with a memo of what it has already copied, so shared references stay shared and cycles terminate; it costs time proportional to the whole structure. Example 4.2 shows the most common version of this bug, `[[0] * n] * m`.

### 3.4 Dicts and the hashing story

A dict maps keys to values. Keys must be **hashable**; values can be anything.

```pycon
>>> d = {"a": 1, "b": 2}
>>> d["c"] = 3
>>> d.get("z", 0), d.get("z"), "a" in d, len(d)
(0, None, True, 3)
>>> d.setdefault("z", []).append(1)
>>> d
{'a': 1, 'b': 2, 'c': 3, 'z': [1]}
>>> for k, v in d.items():
...     print(k, v)
...
a 1
b 2
c 3
z [1]
>>> list(d), list(d.values())
(['a', 'b', 'c', 'z'], [1, 2, 3, [1]])
>>> d2 = {"a": 100, "y": 5}
>>> d | d2
{'a': 100, 'b': 2, 'c': 3, 'z': [1], 'y': 5}
>>> {**d, **d2} == d | d2
True
>>> {v: k for k, v in {"x": 1, "y": 2}.items()}
{1: 'x', 2: 'y'}
>>> dict.fromkeys(["p", "q"], 0)
{'p': 0, 'q': 0}
>>> d["nope"]
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
KeyError: 'nope'
```

Four facts to hold on to. **(i)** `d[k]` raises `KeyError` for a missing key; `d.get(k, default)` never does. **(ii)** Iterating a dict yields its *keys* (`for k in d`, not `for k in d.keys()`); `.values()` and `.items()` are live views. **(iii)** Since Python 3.7 a dict iterates in insertion order, and that is a language guarantee; re-assigning an existing key keeps its position, deleting and re-inserting moves it to the end. **(iv)** In `d | d2` the right operand wins on duplicate keys, and changing a dict's size while iterating over it raises `RuntimeError`.

Two specialized dicts remove the most common boilerplate. `defaultdict(factory)` calls `factory()` to create a missing value on first access, so `groups[key].append(x)` needs no `if key not in groups`. `Counter(iterable)` counts hashable items and adds `most_common(k)`.

▶ **Try it first.** Predict the dict, its length, and the five hashes.

```python
d = {1: "int", 1.0: "float", True: "bool"}
d, len(d)
hash(1), hash(1.0), hash(True), hash(-1), hash(-2)
```

■ **What just happened.**

```text
({1: 'bool'}, 1)
(1, 1, 1, -2, -2)
```

`1 == 1.0 == True`, so they are the *same key*: the first inserted key (`1`) stays, the last value (`'bool'`) wins. And `hash(-1) == hash(-2)`, an example of a collision (CPython reserves -1 as an error code internally).

**The hashing story, in numbered steps.**

1. `hash(key)` computes an integer fingerprint from the key's *value*. Small ints hash to themselves; strings hash by their characters, with a per-process random salt (so the numbers differ between runs, deliberately).
2. **The contract:** if `a == b` then `hash(a) == hash(b)`. The converse is not required. Two unequal keys with the same hash (or the same slot) are a **collision**.
3. The table has $m$ slots, a power of two. The first slot comes from the low bits of the hash. Empty: the key is absent. Occupied: same hash and an equal key means found; otherwise probe another slot, chosen by a fixed rule that depends on the hash, and repeat. (CPython calls this *open addressing*.)
4. Insertion uses the first empty slot on that same probe path. When the table is about two-thirds full CPython builds a bigger one and re-inserts everything, which keeps the expected number of probes bounded independent of $n$: **average $O(1)$**.
5. The guarantee is only *average*. Many keys with the same hash make every lookup a scan, $O(n)$; per-process random string hashes defend against adversaries who would try.
6. **Why keys must be immutable.** The slot was chosen from the hash at insertion time. If the key's value, and so its hash, changes later, the entry sits where the new hash never looks.

Step 6 deserves a demonstration. This class hashes by value, which is legal but a trap when the value is mutable (classes are PY-4; read `self` as "this object").

```python
class Box:
    def __init__(self, v):
        self.v = v
    def __eq__(self, other):
        return self.v == other.v
    def __hash__(self):
        return hash(self.v)

b = Box(1)
d = {b: "found"}
b.v = 2
len(d), Box(1) in d, Box(2) in d, b in d
[k.v for k in d]
```

■ **What just happened.**

```text
(1, False, False, False)
[2]
```

The entry is still listed, and `len` is 1, but no lookup reaches it: `Box(1)` finds the slot, but the stored key no longer equals it; `Box(2)` and `b` hash to a different slot. This is why `list`, `dict` and `set` are unhashable (`hash([1])` raises `TypeError`), why a tuple is hashable only if everything in it is, and why `frozenset` exists.

Hash, slot, equality check, collision: you can watch the first two steps directly. In an 8-slot table, the keys 3, 11 and 19 all land in slot 3:

```pycon
>>> [hash(k) % 8 for k in (3, 11, 19, 4)]
[3, 3, 3, 4]
```

A lookup of 19 would find slot 3 occupied by 3, fail the equality check, and probe on. Exercise B5 builds a working toy table with buckets, and C1 the open-addressing version.

▶ **Build it yourself (2): word frequency three ways.** Count the words in a text with (a) nested loops over two lists, (b) a dict with `get`, (c) `Counter`. Predict which is slowest and how its time scales.

```python
import inspect
import random
import time
from collections import Counter


def freq_lists(words):                 # deliberately naive: a linear search per word
    seen, counts = [], []
    for w in words:
        for i in range(len(seen)):
            if seen[i] == w:
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


def freq_counter(words):
    return Counter(words)


text = "the cat saw the dog and the dog saw the cat"
print(Counter(text.split()).most_common(3))

rng = random.Random(2)
vocab = [f"w{i}" for i in range(400)]
words = rng.choices(vocab, [1 / (i + 1) for i in range(400)], k=40_000)   # Zipf-like frequencies
for fn in (freq_lists, freq_dict, freq_counter):
    t0 = time.perf_counter()
    result = fn(words)
    dt = time.perf_counter() - t0
    n_lines = len(inspect.getsource(fn).strip().splitlines())
    print(f"{fn.__name__:13s} {n_lines:2d} lines  {dt:8.4f} s   total {sum(result.values())}")
```
```text
[('the', 4), ('cat', 2), ('saw', 2)]
freq_lists    11 lines    0.1340 s   total 40000
freq_dict      5 lines    0.0034 s   total 40000
freq_counter   2 lines    0.0023 s   total 40000
```

■ **What just happened.** The three agree on the answer. The nested-loop version is far slower (work $\approx n \times$ distinct words, so it degrades with both) and takes 11 lines against 5 and 2. `Counter` is the same loop as `freq_dict` written in C, so keep `Counter` for counting and the `get` loop for anything else you accumulate per key.

▶ **Build it yourself (3): a sparse polynomial as a dict.** A polynomial $p(x) = \sum_k c_k x^k$ is its list of nonzero $(k, c_k)$ pairs, so store `{power: coefficient}`. The polynomial $x^{10000} + 1$ costs two entries instead of 10,001. Before reading, write `poly_add` on paper: for every power in either input, add the coefficients. Which dict method returns 0 for a missing power?

1. **Add:** copy `p`, then for each `(k, c)` in `q` set `r[k] = r.get(k, 0) + c`. Drop zero coefficients so that equal polynomials are equal dicts.
2. **Multiply:** $c_i x^i \cdot d_j x^j = c_i d_j\, x^{i+j}$, so loop over all pairs and accumulate into `r[i + j]` (a `defaultdict(int)` removes the `get`).
3. **Evaluate:** Horner's rule (PY-1) with gaps. Walk the powers from highest to lowest and multiply the running value by $x^{\text{gap}}$ between terms.
4. **Differentiate:** $c x^k \mapsto k c\, x^{k-1}$, and constants disappear.

```python
from collections import defaultdict


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
    if not p:
        return 0
    result, prev = 0, max(p)
    for k in sorted(p, reverse=True):
        result = result * x ** (prev - k) + p[k]     # x**gap between consecutive powers
        prev = k
    return result * x ** prev


def poly_deriv(p):
    return {k - 1: k * c for k, c in p.items() if k > 0}


p, q = {0: 1, 1: 2, 2: 3}, {0: -1, 1: 1}             # 1 + 2x + 3x^2 and x - 1
print("p*q      =", poly_mul(p, q))
print("p+q      =", poly_add(p, q))
print("p(2)     =", poly_eval(p, 2), "   q(5) =", poly_eval(q, 5))
print("p'       =", poly_deriv(p))
big = {10_000: 1, 0: 1}
print("big*big  =", poly_mul(big, big))
lhs = poly_deriv(poly_mul(p, q))                       # product rule: (pq)' = p'q + pq'
rhs = poly_add(poly_mul(poly_deriv(p), q), poly_mul(p, poly_deriv(q)))
print("(pq)'    =", lhs, " product rule holds:", lhs == rhs)
```
```text
p*q      = {0: -1, 1: -1, 2: -1, 3: 3}
p+q      = {1: 3, 2: 3}
p(2)     = 17    q(5) = 4
p'       = {0: 2, 1: 6}
big*big  = {20000: 1, 10000: 2, 0: 1}
(pq)'    = {0: -1, 1: -2, 2: 9}  product rule holds: True
```

■ **What just happened.** $(1 + 2x + 3x^2)(x - 1) = -1 - x - x^2 + 3x^3$, matching `{0: -1, 1: -1, 2: -1, 3: 3}`; $p(2) = 1 + 4 + 12 = 17$; $p' = 2 + 6x$. The `big*big` line is the point of the representation: squaring $x^{10000} + 1$ touches three entries, not 20,001. The final line verifies the calculus identity $(pq)' = p'q + pq'$ exactly, since the coefficients are integers. `project.py` Part 9 checks the module against `numpy.polynomial.Polynomial` on random inputs.

### 3.5 Sets, and the list-versus-set experiment

A set is an unordered collection of distinct hashable items. Build it with `{1, 2, 3}`, `set(iterable)`, or a comprehension; the empty set is `set()` because `{}` is an empty dict. Algebra: `a | b` union, `a & b` intersection, `a - b` difference, `a ^ b` symmetric difference, `a <= b` subset, `a.isdisjoint(b)`. `s.add(x)`; `s.discard(x)` is silent if `x` is absent while `s.remove(x)` raises `KeyError`. A `frozenset` is an immutable (hence hashable) set, usable as a dict key or a set element.

```pycon
>>> a, b = {1, 2, 3, 4}, {3, 4, 5}
>>> a | b, a & b, a - b, a ^ b
({1, 2, 3, 4, 5}, {3, 4}, {1, 2}, {1, 2, 5})
>>> {x % 3 for x in range(10)}
{0, 1, 2}
>>> xs = [3, 1, 3, 2, 1]
>>> list(set(xs)), list(dict.fromkeys(xs))
([1, 2, 3], [3, 1, 2])
>>> {frozenset({1, 2}): "pair"}[frozenset({2, 1})]
'pair'
```

`list(set(xs))` removes duplicates with no ordering promise; `list(dict.fromkeys(xs))` removes them and keeps first-seen order. For strings, set order depends on the per-process hash salt:

```python
import subprocess
import sys

code = "print(list({'delta', 'alpha', 'gamma', 'beta'}))"
for seed in ("1", "2", "3"):
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                         env={"PYTHONHASHSEED": seed}).stdout.strip()
    print(f"PYTHONHASHSEED={seed}: {out}")
```
```text
PYTHONHASHSEED=1: ['beta', 'delta', 'gamma', 'alpha']
PYTHONHASHSEED=2: ['delta', 'gamma', 'alpha', 'beta']
PYTHONHASHSEED=3: ['gamma', 'beta', 'delta', 'alpha']
```

Small ints often *look* sorted when printed from a set; that is an accident of `hash(n) == n` and the table size, not a promise. Sort when you need order: `sorted(s)`.

▶ **Build it yourself (4): the list-versus-set experiment.** For $n = 10^3, 10^4, 10^5, 10^6$, time `-1 in container` for a list and a set of `range(n)` (the item is absent, which is the list's worst case). Predict the ratio of the $n = 10^6$ time to the $n = 10^3$ time for each, then fit the exponent $s$ in $t \approx C n^{s}$ as the slope of $\log t$ against $\log n$.

```python
import math
import statistics
from timeit import repeat

ns = [10**3, 10**4, 10**5, 10**6]
t_list, t_set = [], []
for n in ns:
    reps = max(2, 2_000_000 // n)
    t_list.append(min(repeat("-1 in c", globals={"c": list(range(n))}, number=reps, repeat=5)) / reps)
    t_set.append(min(repeat("-1 in c", globals={"c": set(range(n))}, number=100_000, repeat=5)) / 100_000)
    print(f"n={n:>9,}  list {t_list[-1] * 1e6:>10.3f} us   set {t_set[-1] * 1e6:>7.3f} us   ratio {t_list[-1] / t_set[-1]:>9,.0f}")
logn = [math.log(n) for n in ns]
for label, ts in (("list", t_list), ("set", t_set)):
    print(f"slope, {label}:", round(statistics.linear_regression(logn, [math.log(t) for t in ts]).slope, 2))
```

■ **What just happened.**

```text
n=    1,000  list      8.023 us   set   0.025 us   ratio       317
n=   10,000  list     79.602 us   set   0.027 us   ratio     2,937
n=  100,000  list    849.062 us   set   0.026 us   ratio    32,691
n=1,000,000  list   9173.432 us   set   0.026 us   ratio   348,975
slope, list: 1.02
slope, set: 0.0
```

A thousand times more data costs the list a thousand times more time (slope near 1: $O(n)$) and costs the set the same time (slope near 0: $O(1)$). The consequence for your code: **a membership test inside a loop wants a set**. `if x in big_list` inside a loop over $m$ items is $O(mn)$; convert to a set once, outside the loop, and it is $O(m + n)$. Sets need hashable items, lose order and duplicates, and use more memory per item, so use them where membership dominates.

### 3.6 Comprehensions

A comprehension is a loop that builds a container: `[expr for x in xs if cond]` (list), `{k: v for ...}` (dict), `{x for ...}` (set), and with parentheses a lazy **generator expression**. Each `for` is a loop and each `if` a filter, in the same order as the equivalent statements; the front expression is the loop body. The loop variable stays inside the comprehension.

▶ **Try it first.** Predict all six results, then run.

```python
rows = [[1, 2, 3], [4, 5, 6]]
[x for row in rows for x in row]
[[x for x in row] for row in rows]
[(i, j) for i in range(3) for j in range(i)]
[row[j] for row in rows for j in range(3) if j != 1]
[[row[j] for row in rows] for j in range(3)]
[x if x % 2 else -x for x in range(5)]
```

■ **What just happened.**

```text
[1, 2, 3, 4, 5, 6]
[[1, 2, 3], [4, 5, 6]]
[(1, 0), (2, 0), (2, 1)]
[1, 3, 4, 6]
[[1, 4], [2, 5], [3, 6]]
[0, 1, -2, 3, -4]
```

The third line settles the order question: `for i in range(3)` is the outer loop and the inner `range(i)` uses `i`, as in nested statements. The fifth is a transpose: the inner comprehension builds one output row. A trailing `if` filters; `x if c else y` at the front chooses a value. Write the loop when you need more than two `for` clauses, an `elif` chain, a `break`, or a side effect.

### 3.7 The iteration protocol and the data model

An **iterable** is anything `iter(x)` accepts. An **iterator** is an object that returns the next item from `next(it)` and raises `StopIteration` when none remain; `iter(it)` returns `it` itself. A list is iterable but not an iterator; `iter(lst)` makes a new bookmark each time, which is why you can loop over a list twice.

▶ **Build it yourself (5): desugar `for`.** Predict each line first.

```python
it = iter([10, 20])
next(it), next(it)
next(it)
list(it)
next([10, 20])
iter(it) is it
```

■ **What just happened.**

```text
(10, 20)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
StopIteration
[]
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: 'list' object is not an iterator
True
```

Once an iterator is exhausted it stays exhausted, and a list is not an iterator. Now write the loop `for x in data: total += x` without using `for`:

1. Ask for the bookmark once: `it = iter(data)`.
2. Loop forever, asking for the next item.
3. Catch `StopIteration` and leave the loop; that exception is the "no more" signal, not an error.

```python
def for_sum(data):
    total, it = 0, iter(data)
    while True:
        try:
            x = next(it)
        except StopIteration:
            break
        total += x
    return total


print(for_sum([1, 2, 3]), for_sum((4, 5)), for_sum({6: "a", 7: "b"}), for_sum(range(8)),
      for_sum(x * x for x in range(4)), for_sum({1, 2, 3}))
```
```text
6 9 13 28 14 6
```

■ **What just happened.** One function sums a list, tuple, dict (its *keys*), range, generator and set, because `for` only calls `iter` and `next`. Files, `zip` and `enumerate` objects and generators are iterators too, and their items need not exist until asked for. That contract makes laziness (§3.8) possible.

**The data model.** `len(x)`, `x[i]`, `i in x`, `for v in x` and `hash(x)` are not magic for lists and dicts. Each is a protocol: Python calls a specially named method on the object, so any type that defines it takes part.

| You write | Python calls | So it works on |
|---|---|---|
| `len(x)` | `x.__len__()` | list, str, dict, set, ... |
| `x[i]`, `x[i:j]` | `x.__getitem__(i)`; for a slice, `i` is a `slice` object | sequences, dicts, NumPy arrays |
| `i in x` | `x.__contains__(i)`; if absent, Python iterates | list, set, dict, range |
| `for v in x` | `iter(x)` calls `x.__iter__()`, then `__next__` repeatedly | everything iterable |
| `hash(x)`, `d[x]` | `x.__hash__()` and `x.__eq__(y)` | dict keys, set items |
| `sorted(xs)` | `a.__lt__(b)` | numbers, strings, tuples |

```pycon
>>> lst = [10, 20, 30, 40]
>>> lst.__len__(), lst.__getitem__(1), lst.__getitem__(slice(1, 3))
(4, 20, [20, 30])
>>> type(lst).__contains__(lst, 30)
True
```

PY-4 writes these methods for your own classes. NumPy arrays and pandas Series work because their authors did exactly that, which is why this module's indexing syntax carries into PY-6 and PY-9.

### 3.8 Generators and itertools

A function containing `yield` is a **generator function**. Calling it runs nothing; it returns a generator, an iterator whose `next` runs the body until the next `yield`, freezes there, and resumes on the following `next`. When the body ends (or `return`s), the generator raises `StopIteration`.

▶ **Try it first.** At which line of output does `start` appear?

```python
def countdown(n):
    print("start")
    while n > 0:
        yield n
        n -= 1
    print("done")


g = countdown(2)
print("created")
print(next(g))
print(next(g))
print(list(g))
```

■ **What just happened.**

```text
created
start
2
1
done
[]
```

`start` appears at the first `next`, not at creation. `list(g)` ran the rest of the body, printing `done`, and found no more items. A **generator expression** is the same idea in one line, and a generator can be infinite because it holds only its current state:

```pycon
>>> import sys
>>> from itertools import islice
>>> sys.getsizeof([x * x for x in range(10**6)]), sys.getsizeof((x * x for x in range(10**6)))
(8448728, 200)
>>> sum(x * x for x in range(10**6))
333332833333500000
>>> def fib():
...     a, b = 0, 1
...     while True:
...         yield a
...         a, b = b, a + b
...
>>> list(islice(fib(), 10))
[0, 1, 1, 2, 3, 5, 8, 13, 21, 34]
>>> g = (c for c in "ab")
>>> list(g), list(g)
(['a', 'b'], [])
```

The generator object is about 200 bytes whatever the range, and the second `list(g)` is empty: **a generator is consumed once.**

▶ **Build it yourself (6): `my_range`, `my_enumerate`, `my_zip`.** Write each as a generator, then compare with the built-in on several inputs. Hints: `range` has two modes (one argument means `stop`); `enumerate` counts as it goes; `zip` pulls one item from each iterator per round and stops at the first exhausted one.

```python
def my_range(start, stop=None, step=1):
    if stop is None:
        start, stop = 0, start
    i = start
    while (i < stop) if step > 0 else (i > stop):
        yield i
        i += step


def my_enumerate(iterable, start=0):
    i = start
    for item in iterable:
        yield i, item
        i += 1


def my_zip(*iterables):
    iterators = [iter(it) for it in iterables]
    while iterators:
        row = []
        for it in iterators:
            try:
                row.append(next(it))
            except StopIteration:
                return                    # shortest input ends the whole zip
        yield tuple(row)


print(list(my_range(5)) == list(range(5)), list(my_range(2, 11, 3)) == list(range(2, 11, 3)),
      list(my_range(5, 0, -2)) == list(range(5, 0, -2)))
print(list(my_enumerate("abc", 1)) == list(enumerate("abc", 1)))
print(list(my_zip("abc", [1, 2])) == list(zip("abc", [1, 2])), list(my_zip("abc", [1, 2])))
```
```text
True True True
True
True [('a', 1), ('b', 2)]
```

■ **What just happened.** Three built-ins that look like language features are short generators. (The real `range` is a lazy *sequence* with `len` and indexing, so it can be re-iterated; the real `zip` also has `strict=True`, below.) `project.py` Part 6 tests all three against the built-ins on 200 random cases.

**`itertools` and the lazy built-ins.** These compose iterators without building lists:

```pycon
>>> from itertools import count, islice, chain, product, combinations, accumulate, groupby
>>> list(islice(count(10, 5), 4))
[10, 15, 20, 25]
>>> list(chain([1, 2], "ab"))
[1, 2, 'a', 'b']
>>> list(product("ab", [0, 1]))
[('a', 0), ('a', 1), ('b', 0), ('b', 1)]
>>> list(combinations("abc", 2))
[('a', 'b'), ('a', 'c'), ('b', 'c')]
>>> list(accumulate([3, 1, 4, 1, 5])), list(accumulate([3, 1, 4, 1, 5], max))
([3, 4, 8, 9, 14], [3, 3, 4, 4, 5])
>>> [(k, len(list(g))) for k, g in groupby("aabccc")]
[('a', 2), ('b', 1), ('c', 3)]
>>> words = ["pear", "fig", "banana"]
>>> list(enumerate(words, start=1))
[(1, 'pear'), (2, 'fig'), (3, 'banana')]
>>> list(zip(words, range(10)))
[('pear', 0), ('fig', 1), ('banana', 2)]
>>> list(zip(words, range(10), strict=True))
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
ValueError: zip() argument 2 is longer than argument 1
>>> list(reversed(words)), max(words, key=len), max([], default=None)
(['banana', 'fig', 'pear'], 'banana', None)
>>> any(w.startswith("b") for w in words), all(len(w) > 3 for w in words)
(True, False)
>>> any(n > 3 for n in count())
True
```

`groupby` groups only *consecutive* equal keys, so it needs sorted input (Example 4.4). `any` and `all` stop at the first decisive item, so the last line terminates on an infinite stream. `zip` stops silently at the shortest input; pass `strict=True` when unequal lengths would be a bug.

### 3.9 Strings, revisited

A `str` is an immutable sequence of Unicode code points; slicing, `in` and iteration work as for lists. A `bytes` object is an immutable sequence of integers 0 to 255; `.encode(enc)` converts text to bytes and `.decode(enc)` back, and the two types never mix silently. A text-mode file decodes for you with a platform-dependent default encoding, so always write `encoding="utf-8"`.

```pycon
>>> "a,b,,c".split(","), "  a  b ".split(), "-".join(["2025", "01", "07"])
(['a', 'b', '', 'c'], ['a', 'b'], '2025-01-07')
>>> "Hello, world!".translate(str.maketrans("", "", "!,"))
'Hello world'
>>> "lo" in "hello", "hello".find("l")
(True, 2)
>>> len("é"), "é".encode("utf-8"), len("é".encode("utf-8")), b"\xc3\xa9".decode("utf-8")
(1, b'\xc3\xa9', 2, 'é')
```

`join` builds a long string in $O(n)$ total, whereas repeated `s += piece` can be quadratic (CPython sometimes optimizes it; do not rely on that). `str.translate` with `str.maketrans(from, to, delete)` replaces or deletes characters in one pass.

### 3.10 Files, paths, CSV and JSON

`open(path, mode, encoding=...)` returns a file object. Mode `"r"` (default) reads text, `"w"` truncates and writes, `"a"` appends, and `"b"` adds bytes with no decoding. Always use `with open(...) as f:`, which closes the file even if an exception escapes. A text file is an iterator over its lines, so `for line in f` reads lazily in constant memory, while `f.read()` and `f.readlines()` load everything. `pathlib.Path` handles paths: `/` joins, and `exists`, `glob`, `read_text`, `write_text`, `parent` / `parents`, `name`, `stem`, `suffix`, `mkdir(parents=True, exist_ok=True)` replace most of `os.path`.

```python
import csv
import tempfile
from pathlib import Path

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    path = root / "data" / "run1.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("name,reading\nbaro,1013.2\nlux,402.5\nbaro,1011.9\n", encoding="utf-8")
    print(path.exists(), path.name, path.stem, path.suffix, path.parent.name,
          [p.name for p in root.glob("**/*.csv")])
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))                  # each row is a dict; every value is a str
    print(rows[0], type(rows[0]["reading"]).__name__)
    with path.open(encoding="utf-8") as f:              # lazy: the header, then one line at a time
        header = next(f).strip()
        total = sum(float(line.split(",")[1]) for line in f)
    print(header, round(total, 1))
    with (root / "out.csv").open("w", newline="", encoding="utf-8") as f:
        csv.writer(f).writerows([["name", "mean"], ["baro", 1012.55]])
    print((root / "out.csv").read_text(encoding="utf-8"), end="")
```
```text
True run1.csv run1 .csv data ['run1.csv']
{'name': 'baro', 'reading': '1013.2'} str
name,reading 2427.6
name,mean
baro,1012.55
```

Open CSV files with `newline=""` (the `csv` module handles line endings). `csv.reader` yields lists, `DictReader` dicts keyed by the header, and **everything comes back as `str`**: converting types is your job. **JSON** (`dump`/`load` for files, `dumps`/`loads` for strings) maps object ↔ `dict`, array ↔ `list`, string ↔ `str`, number ↔ `int` or `float`, `true`/`false`/`null` ↔ `True`/`False`/`None`. Not preserved: tuples become lists, non-string keys become strings, and `nan` is written as the non-standard token `NaN` (Example 4.6).

### 3.11 Choosing a container

| You need | Use | Because |
|---|---|---|
| An ordered, growable sequence indexed by position | `list` | $O(1)$ index and append |
| A fixed record, a dict key, several return values | `tuple` / `namedtuple` | immutable, hashable |
| Lookup by key | `dict` | $O(1)$ average |
| Uniqueness or fast membership, order irrelevant | `set` | $O(1)$ average `in` |
| Counting / grouping by key | `Counter` / `defaultdict(list)` | removes the "first time" branch |
| Cheap insert and remove at both ends (queue, window) | `deque` | $O(1)$ at both ends |
| Remove duplicates but keep order | `dict.fromkeys(xs)` | dicts are ordered |
| Data too big to hold, or a stream | generator / iterator | $O(1)$ memory |
| Millions of numbers | NumPy array (PY-6) | contiguous raw values |

Ask in order: do I look things up by something other than position (dict, set)? Do I need duplicates or order (list)? Do I work at the front (deque)? Must I hold all of it at once (maybe not: generator)?

## 4. Worked examples

### Example 4.1: Slicing with a negative stride, and reversing a string

```
index:      0   1   2   3   4   5   6   7   8   9
s:          a   b   c   d   e   f   g   h   i   j
negative: -10  -9  -8  -7  -6  -5  -4  -3  -2  -1
```

`s[5:2:-1]`: start at position 5 (`f`), step $-1$, stop *before* position 2: positions 5, 4, 3, so `fed`. `s[-2:-6:-2]`: positions 8, 6, and 4 is excluded, so `ig`. `s[::-1]` starts at the default 9 and runs down to position 0.

```python
s = "abcdefghij"
print(s[5:2:-1], s[-2:-6:-2], s[8:2:-3], s[::-1])
print("level" == "level"[::-1], "plane" == "plane"[::-1])
lst = [1, 2, 3]
print(lst[::-1], list(reversed(lst)), lst.reverse(), lst)
```
```text
fed ig if jihgfedcba
True False
[3, 2, 1] [3, 2, 1] None [3, 2, 1]
```

`lst[::-1]` makes a reversed copy; `reversed(lst)` is a lazy iterator over the original; `lst.reverse()` reverses in place and returns `None`. In PY-6 `a[::-1]` on a NumPy array is a *view*, not a copy.

### Example 4.2: The `[[0] * 3] * 3` trap

`[0] * 3` is fine: it copies three references to the immutable `0`. The outer `* 3` copies three references to *the same inner list*.

```
grid = [[0]*3]*3                   grid = [[0]*3 for _ in range(3)]

grid ─► [ ● ● ● ]                  grid ─► [ ●   ●   ● ]
          └─┴─┴──► [1, 0, 0]                 │   │   │
                                              ▼   ▼   ▼
                                          [1,0,0] [0,0,0] [0,0,0]
```

```python
bad = [[0] * 3] * 3
bad[0][0] = 1
good = [[0] * 3 for _ in range(3)]
good[0][0] = 1
print(bad, bad[0] is bad[2])
print(good, good[0] is good[2])
```
```text
[[1, 0, 0], [1, 0, 0], [1, 0, 0]] True
[[1, 0, 0], [0, 0, 0], [0, 0, 0]] False
```

The comprehension evaluates `[0] * 3` once per row, so each row is a fresh list.

### Example 4.3: Sort records by score descending, then name

```python
from operator import itemgetter

records = [("Ana", 82), ("Raj", 91), ("Lee", 82), ("Kim", 91), ("Bo", 75), ("Zed", 82)]
print(sorted(records, key=itemgetter(1), reverse=True))
by_tuple = sorted(records, key=lambda r: (-r[1], r[0]))
two_pass = sorted(sorted(records, key=itemgetter(0)), key=itemgetter(1), reverse=True)
print(by_tuple)
print(by_tuple == two_pass)
```
```text
[('Raj', 91), ('Kim', 91), ('Ana', 82), ('Lee', 82), ('Zed', 82), ('Bo', 75)]
[('Kim', 91), ('Raj', 91), ('Ana', 82), ('Lee', 82), ('Zed', 82), ('Bo', 75)]
True
```

The first line sorts by score only: ties keep their input order (Raj before Kim, Ana-Lee-Zed) because the sort is stable, even with `reverse=True`. The tuple key `(-score, name)` sorts descending on score by negation and ascending on name. Where negation is impossible (a string key you want descending), use two passes: secondary key first, then primary with `reverse=True`; stability keeps the secondary order inside ties.

### Example 4.4: `groupby` on sorted records versus `defaultdict(list)`

```python
from collections import defaultdict
from itertools import groupby
from operator import itemgetter

readings = [("b", 3.0), ("a", 1.0), ("b", 5.0), ("a", 2.0), ("c", 7.0), ("a", 4.0)]

groups = defaultdict(list)
for key, value in readings:
    groups[key].append(value)
print(dict(groups))

by_key = itemgetter(0)
print({k: [v for _, v in g] for k, g in groupby(sorted(readings, key=by_key), key=by_key)})
print([(k, [v for _, v in g]) for k, g in groupby(readings, key=by_key)])   # forgot to sort
spent = list(groupby(readings, key=by_key))
print([list(g) for _, g in spent])                                          # groups are spent
```
```text
{'b': [3.0, 5.0], 'a': [1.0, 2.0, 4.0], 'c': [7.0]}
{'a': [1.0, 2.0, 4.0], 'b': [3.0, 5.0], 'c': [7.0]}
[('b', [3.0]), ('a', [1.0]), ('b', [5.0]), ('a', [2.0]), ('c', [7.0]), ('a', [4.0])]
[[], [], [], [], [], []]
```

`defaultdict` makes one $O(n)$ pass, needs no sort, and keeps keys in first-seen order, but it holds everything. `groupby` needs sorted input ($O(n \log n)$), and in exchange it streams. Unsorted input gives six groups with `a` and `b` repeated; and each inner group is a view into the same underlying stream, so materializing `list(groupby(...))` first leaves every group empty. Consume a group before asking for the next.

### Example 4.5: A pipeline over a million lines, without holding the file

Create a file with $10^6$ lines `i,v` where $v = (7919 i) \bmod 1000$, keep the lines with $v \ge 900$, and add them up. Since $\gcd(7919, 1000) = 1$, each residue 0 to 999 occurs exactly $1000$ times, so the answer is known in closed form: $100\,000$ kept lines and $1000 \times (900 + \dots + 999) = 94\,950\,000$.

```python
import tempfile
import tracemalloc
from pathlib import Path


def read_lines(path):
    with path.open(encoding="utf-8") as f:
        for line in f:
            yield line


def streaming(path):
    n = total = 0
    for line in read_lines(path):                     # holds nothing
        v = int(line.split(",")[1])
        if v >= 900:
            n, total = n + 1, total + v
    return n, total


def eager(path):
    lines = path.read_text(encoding="utf-8").splitlines()     # the whole file, then a list of lines
    values = [int(line.split(",")[1]) for line in lines]
    kept = [v for v in values if v >= 900]
    return len(kept), sum(kept)


def peak(fn, *args):
    tracemalloc.start()
    result = fn(*args)
    _, top = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return result, top / 1024


with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp) / "big.txt"
    path.write_text("".join(f"{i},{i * 7919 % 1000}\n" for i in range(10**6)), encoding="utf-8")
    for fn in (eager, streaming):
        result, kib = peak(fn, path)
        print(f"{fn.__name__:15s} {result}  peak {kib:>10,.0f} KiB")
```
```text
eager           (100000, 94950000)  peak     87,189 KiB
streaming       (100000, 94950000)  peak         22 KiB
```

Both give the closed-form answer. The streaming loop's peak memory is a few KiB whatever the file size; `eager` holds the whole file as text, a list of lines, and two more lists. Generators let the working set be one line instead of the file.

### Example 4.6: A JSON round trip that is not the identity

```python
import json

data = {"run": 7, "t": (0.0, 0.1), "x": 0.1 + 0.2, "tiny": 1e-7, "big": 1e22, 3: "three"}
text = json.dumps(data)
back = json.loads(text)
print(text)
print(back)
print("floats exact:", back["x"] == data["x"] and back["tiny"] == data["tiny"] and back["big"] == data["big"])
print("equal overall:", back == data, " tuple became", type(back["t"]).__name__, " key 3 became", repr(list(back)[-1]))

nan_text = json.dumps({"v": float("nan")})
v = json.loads(nan_text)["v"]
print(nan_text, v != v, {"v": float("nan")} == json.loads(nan_text))
```
```text
{"run": 7, "t": [0.0, 0.1], "x": 0.30000000000000004, "tiny": 1e-07, "big": 1e+22, "3": "three"}
{'run': 7, 't': [0.0, 0.1], 'x': 0.30000000000000004, 'tiny': 1e-07, 'big': 1e+22, '3': 'three'}
floats exact: True
equal overall: False  tuple became list  key 3 became '3'
{"v": NaN} True False
```

The floats print differently than you typed them (`1e-07`, `1e+22`, `0.30000000000000004`) but are bit-identical after loading, because Python writes the shortest decimal string that reads back to the same float. What changed is the *structure*: the tuple is now a list and `3` is now `"3"`. A `nan` is written as `NaN` (not valid JSON for strict parsers), is not equal to itself, and so the "same" data compares unequal after a round trip. `json.dumps(..., allow_nan=False)` turns that into an error at write time.

### Example 4.7: A parameter grid, lazily (from your ODE course)

A damped oscillator $m\ddot x + c\dot x + kx = 0$ has damping ratio $\zeta = c / (2\sqrt{km})$. It oscillates with decaying amplitude if $\zeta < 1$, is critically damped at $\zeta = 1$, and decays without oscillating if $\zeta > 1$. Classify every combination of three masses, four dampers and three springs, without building the 36-point grid as a list. (The same code sweeps $10^6$ points: `product` holds one tuple at a time. CM-6 returns to oscillators in many degrees of freedom.)

```python
import math
from collections import Counter
from itertools import islice, product


def regime(m, c, k):
    zeta = c / (2 * math.sqrt(k * m))
    return "critical" if math.isclose(zeta, 1.0) else "under" if zeta < 1 else "over"


masses, dampers, springs = (0.5, 1.0, 2.0), (0.5, 1.0, 2.0, 4.0), (1.0, 2.0, 8.0)
print(list(islice(product(masses, dampers, springs), 3)))
print(Counter(regime(*p) for p in product(masses, dampers, springs)))
```
```text
[(0.5, 0.5, 1.0), (0.5, 0.5, 2.0), (0.5, 0.5, 8.0)]
Counter({'under': 26, 'over': 6, 'critical': 4})
```

`math.isclose` (PY-1) guards the critical case against rounding. The `*p` unpacks each tuple into three arguments.

## 5. Why it matters and where it connects

- **PY-3** generalizes `key=`: functions are values, so `sorted(xs, key=f)`, dispatch dicts and `functools` follow, and generators gain decorators. Its `numtools` package reuses the sparse polynomial of §3.4.
- **PY-4** implements the data-model table of §3.7 for your own classes (`__getitem__`, `__iter__`, `__eq__`/`__hash__`).
- **PY-6**: NumPy indexing is §3.1 with an index per axis; `a[::-1]` and `a[1:3]` are *views* sharing memory, the array analogue of the shared inner lists in Example 4.2. Boolean masks replace the filtering comprehension of §3.6.
- **PY-9**: pandas `groupby` is Example 4.4, and PY-9's hand-written groupby is `defaultdict(list)`. A `Series` is a dict-like label-to-value map with an array underneath.
- **PY-10** fixes slow programs by changing the data structure (§3.5's set, `deque`) before touching constants.
- PY-9's `read_csv` is the industrial version of §3.10's `csv.DictReader`; the cleaning problems (strings that should be numbers, missing values) are the same.

## 6. Exercises

See [exercises.md](exercises.md). A1–A6 are predict-the-output and short-function drills; do all of them. B1–B8 are the module; B1 (iterator tools as generators), B3 (log parsing with `defaultdict`), B5 (a toy hash table) and B7 (extend `project.py`'s timing experiment) are essential. C1–C3 (an open-addressing table, an LRU cache and what measuring it reveals about dict order, and an $O(n^2)$-to-$O(n)$ rewrite with measurements) go further.

## 7. Coding project

`project.py` generates a synthetic sensor log of $10^5$ lines (seeded `random`, 0.5% corrupt) and builds:

1. **A generator pipeline** (`read_lines` → `parse` → `summarize`) giving per-sensor count, mean, min and max with `defaultdict`/`Counter`, asserted against a second computation by sorting and `groupby`.
2. **A complexity lab:** list vs set membership for $n = 10^3 \dots 10^6$, and front insertion for list, deque and linked list, with log-log slopes fitted and asserted (about 1 for the list cases, about 0 for the others) and a figure saved.
3. **The sparse polynomial** (§3.4), tested against `numpy.polynomial.Polynomial` on random inputs.
4. **CSV and JSON writers** for the summaries with exact round-trip asserts; shallow vs deep copy asserted.
5. **`my_range`, `my_enumerate`, `my_zip`** and a `for` loop desugared by hand, against the built-ins on random cases.
6. **Every number in §4** asserted, including the 10⁶-line pipeline's closed-form answer and memory ratio, and the three word-frequency functions timed.

**Build it yourself:** before opening the file, write `parse` (yield `(t, sensor, value)`, collect rejected lines) and `summarize` using only `defaultdict`, `Counter` and `dict.get`. Then do "Your turn": a `Counter`-based anagram grouper, a lazy `itertools.product` grid, and a loop-to-comprehension rewrite (and back). Exercise B7 extends the timing lab.

## 8. Where you'll struggle

- **Slice endpoints with negative strides.** Everyone writes `s[5:2:-1]` expecting `s[2:5]` reversed once. A negative stride reverses the direction of *both* endpoints and flips the defaults. Write out the visited indices, or print `slice(i, j, k).indices(len(s))`.
- **`a = a.sort()`.** Methods that mutate in place (`sort`, `append`, `reverse`, `extend`, `insert`) return `None`, so the assignment destroys your list. Use `sorted(a)` for a new list.
- **Shallow copies of nested data.** The bug that costs a day, because it surfaces only when a second function reads the structure you "copied". Usual suspects: `[[0]*n]*m`, `list(a)` of a list of lists, and `dict.fromkeys(keys, [])` (one list shared by every key). Rule: if the elements are mutable, copy them too.
- **Order.** Dict iteration is insertion order *by guarantee since 3.7*; set order is not defined at all. Small-integer sets often print sorted and teach the wrong lesson. If order matters, use a list, a dict, or `sorted()`.
- **Generators are consumed once.** A second pass silently yields nothing, with no error. The same holds for `map`, `filter`, `zip`, `enumerate`, `reversed` and file objects, and none of them support `len()` or indexing. For two passes, make a list or call the generator function again.
- **`zip` truncates silently.** Lists of different lengths give a result of the shorter length and no complaint. Use `strict=True` unless truncation is intended.
- **Comprehension nesting order.** Clauses read left to right like the loops they replace (`for row in rows for x in row`), which nobody believes until they test it. When in doubt, write the loops first.
- **Changing a container while iterating it.** Deleting from a list inside `for x in lst` skips elements; growing a dict raises `RuntimeError`. Build a new container or iterate over `list(lst)`.
- **`groupby` without sorting,** and `in` on a list inside a loop: quietly wrong, and quietly quadratic. Example 4.4 and the §3.5 experiment are the antidotes.

## 9. Check yourself

1. Why is `lst.insert(0, x)` $O(n)$ but `lst.append(x)` amortized $O(1)$? Draw the list.
2. What are the default start and stop of `s[::-1]`, and why is `s[2:5:-1]` empty?
3. Why does `sorted(records, key=itemgetter(1), reverse=True)` keep tied records in input order?
4. How do a shallow copy and `copy.deepcopy` differ? Give a nested-list example where the results differ.
5. State the contract between `__eq__` and `__hash__`. Why does `{1: "a", 1.0: "b"}` have one key?
6. Why are lists unhashable? What goes wrong if a key's hash changes after insertion?
7. When is a dict lookup not $O(1)$?
8. What does `for x in obj:` call, and what ends the loop?
9. Why does a generator use almost no memory, and what happens if you iterate it twice?
10. Which of `set`, `list`, `deque`, `dict` do you choose for membership tests in a loop; a FIFO queue; counting words; ordered unique items?
11. Name two things that change in a JSON round trip, and what CSV does to every value.

## 10. Pair with

- **Downey, *Think Python* 3e, ch. 9–14.** The gentlest walk through lists, dictionaries, tuples and files; use it for extra practice at a slower pace.
- **Ramalho, *Fluent Python* 2e, ch. 2–3 and 17.** The deepest treatment of sequences, dict and set internals, and iterators and generators. Read ch. 3 after §3.4 and ch. 17 after §3.8.
- **Python Tutorial §5 (data structures) and §7 (input and output)**, docs.python.org. Short and authoritative; the place to check a method name.
- **Python docs, "Sorting HOW TO" and the `itertools` recipes.** The recipes are a catalog of what lazy pipelines can do.
- **Beazley, *Python Distilled*, ch. 2–3.** Compact reference prose for the same material.
- **Berkeley CS61A lectures on sequences and iterators.** Environment diagrams for §3.7.
