# PY-2 Exercises

Conventions as in the notes: Python 3.11+, standard library only (B7 edits `project.py`, which also uses matplotlib). Write your prediction down before you run anything. "Write a function" means: write it, then test it with `assert`s you chose, including an edge case.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** (Slicing.) With `s = "0123456789"`, predict `s[7:2:-2]`, `s[-1:-4:-1]`, `s[::4]`, `s[:-7]`, `s[3:3]`, and `s[8:100]`. Then with `lst = list(range(6))`, predict `lst` after `lst[1:4] = []`, and after a further `del lst[::2]`.

**A2.** (Aliasing and copies.) Predict the three printed lists, then draw the arrows.

```python
a = [1, [2, 3]]
b = a[:]
c = a
b[1].append(4)
b[0] = 99
c.append(5)
print(a, b, c)
```

**A3.** (Fix the bug.) Each snippet misbehaves. Say what it does and repair it.

```python
# (a) goal: print the smallest score
scores = [3, 1, 2]
scores = scores.sort()
print(scores[0])

# (b) goal: keep only the odd numbers
nums = [2, 4, 6, 7, 8]
for x in nums:
    if x % 2 == 0:
        nums.remove(x)
print(nums)

# (c) goal: print "14 14"
squares = (x * x for x in range(4))
print(sum(squares), sum(squares))
```

**A4.** (Predict.)

```python
print([(a, b) for a in "xy" for b in range(2)])
print({len(w) for w in ["aa", "b", "cc"]})
g = (n for n in range(3))
print(next(g), list(g), list(g))
print(list(zip("abc", [1, 2])), list(enumerate("ab", start=5)))
```

**A5.** (Write a function.) `invert(d)` maps each value of `d` to the list of keys that had it, built with `defaultdict`: `invert({"a": 1, "b": 2, "c": 1})` is `{1: ["a", "c"], 2: ["b"]}`. What happens, and why, if a value of `d` is a list?

**A6.** (Hashing.) Predict the output of the three lines, including the error message.

```python
d = {(1, 2): "tuple"}
d[1.0] = "float"
d[1] = "int"
print(d, len(d))
print(len({0.1 + 0.2, 0.3}), len({1, True, 1.0}), {1, 2} == {2, 1})
d[[1, 2]] = 0
```

## B. Core (do all; these are the module)

**B1.** (Iterator tools.) Write generators `my_chain(*iterables)`, `my_islice(iterable, stop)`, `my_accumulate(iterable)` (running sums), and add a `strict=False` keyword to the `my_zip` of §3.8 that raises `ValueError` on unequal lengths when true. Test each against `itertools` or the built-in on finite inputs, and test `my_islice` and `my_accumulate` on the infinite `itertools.count()`.

**B2.** (A lazy pipeline.) Write three generators: `lines(path)` yields stripped lines of a file; `numbers(lines)` skips blank lines, lines starting with `#`, and lines that are not floats, yielding floats; `running_mean(values)` yields the mean of everything seen so far. Test on a file containing `# header`, `1.5`, a blank line, `2.5`, `bad`, `3.5`. Then use `itertools.islice` to take the first two running means. Why does the pipeline not read the third number until you ask for it?

**B3.** (Per-key summaries.) This text has a header, one malformed row, and a repeated sensor:

```python
LOG = """t,sensor,value
0.0,a,1.0
0.1,b,10.0
0.2,a,3.0
0.3,b,oops
0.4,a,5.0
0.5,b,20.0
"""
```

Parse it with `csv.DictReader` over `LOG.splitlines()`, and compute per-sensor count, mean, min and max with a `defaultdict`, skipping and counting the malformed row. Verify against a second computation using `sorted` and `groupby`.

**B4.** (Sparse polynomials.) Extend §3.4: write `poly_integ(p)` (antiderivative with zero constant term) and `poly_compose(p, q)` returning $p(q(x))$. Check that `poly_deriv(poly_integ(p)) == p`, that `poly_eval(poly_compose(p, q), x) == poly_eval(p, poly_eval(q, x))` at several integer $x$, and that composing with `{1: 1}` (the polynomial $x$) changes nothing.

**B5.** (A toy hash table.) Implement `ToyDict` with separate chaining: a list of buckets, each a list of `(key, value)` pairs; methods `put`, `get(key, default=None)`, `delete`, and `__len__`; double the bucket count when the number of items exceeds the number of buckets. Test against a real `dict` on 5,000 random operations. Then insert the keys `0, 64, 128, ...` (64 of them) into a fresh 64-bucket table with growth disabled and report the longest chain. What does that do to lookup cost, and what would a table have to do differently to avoid it?

**B6.** (Multi-key sorting.) Sort these records by department **descending**, then salary ascending, then name ascending. A tuple key cannot negate a string, so use three stable passes. Then verify the result by checking every adjacent pair directly.

```python
staff = [("Ana", "eng", 95), ("Raj", "ops", 70), ("Lee", "eng", 95), ("Kim", "ops", 82),
         ("Bo", "eng", 80), ("Zed", "ops", 70)]
```

**B7.** (Extend `project.py`.) A FIFO queue appends at the back and removes from the front. Add two entries to `EXPERIMENTS` in Part 2, `"list queue: append+pop(0)"` and `"deque queue: append+popleft"`, each timing one append followed by one removal from the front at constant size $n$. Fit the slopes, assert them (about 1 for the list, about 0 for the deque), and plot them in a third panel.

**B8.** (Measuring scaling.) With the word-frequency functions of §3.4, fix 20,000 words drawn uniformly from a vocabulary of size $V$ and time `freq_lists` and `freq_dict` for $V = 50, 100, 200, 400, 800$. Fit the log-log slope of time against $V$ for each. Predict both slopes first, using the work estimate in §3.4, and explain any gap between prediction and measurement.

## C. Challenge (do at least one)

**C1.** (Open addressing.) Implement a hash table with one flat list of slots and linear probing: on a collision try slots $i+1, i+2, \dots$ (wrapping around). Support `put`, `get` and `delete`. Show with a failing example why `delete` cannot reset the slot to empty, and fix it with a *tombstone* marker that lookups skip and inserts may reuse. Resize when occupied slots (live plus tombstones) exceed two-thirds of the table. Test against `dict` on 10,000 random operations. Then count probes: at load factors 0.2, 0.4, 0.6 and 0.8 (no resizing), report the average number of probes per successful `get` and describe the shape of the curve.

**C2.** (An LRU cache.) A least-recently-used cache of capacity $k$ evicts the entry used longest ago. (a) Write `LRUDict` with a plain `dict` and its insertion order: `get` removes and re-inserts the key to mark it newest; `put` evicts `next(iter(self.data))` when full. (b) Write `LRUList`: a list of keys, oldest first, plus a dict of values. (c) Write `LRUOrdered` using `collections.OrderedDict`, which has `move_to_end(key)` and `popitem(last=False)`. Check all three agree on a random access sequence. Time the cost per access for $k = 100, 1000, 10000$ with a 2-in-3 miss rate and fit the slopes. Hint for the explanation: CPython stores a dict's entries in one array in insertion order, and deleting an entry leaves a hole until the dict is next rebuilt.

**C3.** (From $O(n^2)$ to $O(n)$, and when it pays off.) Given a list of $n$ distinct integers and a target, find two positions whose values sum to the target. (a) Write the double-loop version and a one-pass dict version, assert they agree on 200 random instances, and time both (target absent, so every pair is checked) for $n = 500, 1000, 2000, 4000$. Fit the slopes. (b) For $n = 1000$, time one membership test `-1 in xs` on the list against building `set(xs)`, and compute the smallest number of lookups $q$ for which building the set first is cheaper. State in a sentence what this says about big-$O$ and small $q$.
