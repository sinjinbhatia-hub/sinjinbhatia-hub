# %% [markdown]
# # PY-2 Project: A sensor log, five containers, and a lazy pipeline
#
# You generate a synthetic sensor log (10^5 lines), parse it with a generator
# pipeline, summarise it per sensor with a `defaultdict`, and check the answer
# against an independent computation (sorted + groupby). Then you measure the
# costs the notes claim (list `in` is O(n), set `in` is O(1), front insertion
# into a list is O(n) but O(1) into a deque or linked list) and fit the slopes
# on a log-log plot. You also build a sparse polynomial as a dict, round-trip
# data through CSV and JSON, demonstrate shallow vs deep copy, rebuild
# `range`/`enumerate`/`zip` as generators, and assert every number quoted in
# the notes' worked examples. NumPy appears only in the last cross-check.
#
# Run headless: `MPLBACKEND=Agg python3 project.py`. Figures and data go to ./figures/.

# %%
import copy
import csv
import itertools
import json
import math
import random
import statistics
import tempfile
import time
import tracemalloc
from collections import Counter, defaultdict, deque
from operator import itemgetter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
DATA = FIG / "data"
DATA.mkdir(exist_ok=True)

# %% [markdown]
# ## Part 1: A synthetic log and a generator pipeline (notes §3.8, §4.5)
#
# Each line is `timestamp,sensor,reading`. About 0.5% of lines are corrupt
# (reading `ERR` or a missing field), because real logs are. The pipeline is
# three generators chained: read lines -> parse -> summarise. No stage ever
# holds the whole file.

# %%
PROFILE = {  # sensor -> (mean, standard deviation)
    "baro": (1013.0, 3.0), "humid": (55.0, 8.0), "temp_a": (21.0, 2.0),
    "temp_b": (19.0, 2.5), "wind": (4.5, 1.5), "soil": (30.0, 5.0),
    "lux": (400.0, 90.0), "co2": (415.0, 20.0),
}
SENSORS = list(PROFILE)


def make_log(path: Path, n_lines: int, seed: int = 0) -> int:
    """Write n_lines of `t,sensor,reading` to path; return how many are corrupt."""
    rng = random.Random(seed)
    n_bad, t = 0, 0.0
    with path.open("w", encoding="utf-8") as f:
        for _ in range(n_lines):
            t += rng.expovariate(10.0)  # about 0.1 s between lines
            name = rng.choice(SENSORS)
            mu, sd = PROFILE[name]
            u = rng.random()
            if u < 0.0025:
                f.write(f"{t:.3f},{name},ERR\n")
                n_bad += 1
            elif u < 0.005:
                f.write(f"{t:.3f},{name}\n")  # truncated line
                n_bad += 1
            else:
                f.write(f"{t:.3f},{name},{rng.gauss(mu, sd):.3f}\n")
    return n_bad


def read_lines(path: Path):
    """Stage 1: yield stripped lines lazily (constant memory)."""
    with path.open(encoding="utf-8") as f:
        for line in f:
            yield line.rstrip("\n")


def parse(lines, rejects: list):
    """Stage 2: yield (t, sensor, value) tuples; append malformed lines to rejects."""
    for line in lines:
        parts = line.split(",")
        if len(parts) != 3:
            rejects.append(line)
            continue
        t, name, value = parts
        try:
            yield float(t), name, float(value)
        except ValueError:
            rejects.append(line)


def summarize(records) -> dict:
    """Stage 3: one pass, per-sensor count/mean/min/max using dict accumulators."""
    count, total = Counter(), defaultdict(float)
    lo, hi = {}, {}
    for _, name, v in records:
        count[name] += 1
        total[name] += v
        lo[name] = min(lo.get(name, v), v)
        hi[name] = max(hi.get(name, v), v)
    return {s: {"count": count[s], "mean": total[s] / count[s], "min": lo[s], "max": hi[s]}
            for s in sorted(count)}


def summarize_groupby(records: list) -> dict:
    """Independent check: sort by sensor, then itertools.groupby."""
    out = {}
    for name, grp in itertools.groupby(sorted(records, key=itemgetter(1)), key=itemgetter(1)):
        vals = [v for _, _, v in grp]
        out[name] = {"count": len(vals), "mean": math.fsum(vals) / len(vals),
                     "min": min(vals), "max": max(vals)}
    return out


LOG = DATA / "sensors.log"
n_bad_written = make_log(LOG, 100_000)
rejects: list = []
records = list(parse(read_lines(LOG), rejects))
assert len(records) + len(rejects) == 100_000
assert len(rejects) == n_bad_written
assert all(r.endswith("ERR") or r.count(",") == 1 for r in rejects)

summary = summarize(iter(records))
check = summarize_groupby(records)
assert summary.keys() == check.keys() == set(SENSORS)
for s in SENSORS:
    assert summary[s]["count"] == check[s]["count"]
    assert summary[s]["min"] == check[s]["min"] and summary[s]["max"] == check[s]["max"]
    assert math.isclose(summary[s]["mean"], check[s]["mean"], rel_tol=1e-12)
    mu, sd = PROFILE[s]  # sample mean of ~12k draws is within 5 standard errors of mu
    assert abs(summary[s]["mean"] - mu) < 5 * sd / math.sqrt(summary[s]["count"])
print(f"Part 1: {len(records)} good lines, {len(rejects)} rejected; "
      f"per-sensor summary agrees with sorted+groupby.")
for s, row in summary.items():
    print(f"  {s:7s} n={row['count']:5d} mean={row['mean']:8.3f} "
          f"min={row['min']:8.3f} max={row['max']:8.3f}")

# %% [markdown]
# ## Part 2: Measuring big-O (notes §3.1, §3.5, build sequence 2)
#
# Cost of `-1 in container` (absent key = list's worst case) for n = 10^3..10^6,
# and cost of one insert-at-front + pop-at-front pair at constant size n.
# The slope of log(time) against log(n) is the exponent: 1 means O(n), 0 means O(1).

# %%
def per_op(container, op, reps: int, repeat: int = 5) -> float:
    """Seconds per call of op(container): the minimum over `repeat` timed runs of `reps` calls."""
    best = math.inf
    for _ in range(repeat):
        t0 = time.perf_counter()
        for _ in range(reps):
            op(container)
        best = min(best, time.perf_counter() - t0)
    return best / reps


def loglog_slope(ns, ts) -> float:
    """Least-squares slope of log t against log n (the empirical exponent)."""
    return statistics.linear_regression([math.log(n) for n in ns], [math.log(t) for t in ts]).slope


def linked(n: int):
    """Linked list as nested tuples (value, rest), built by n prepends."""
    head = None
    for i in range(n):
        head = (i, head)
    return head


NS = [10**3, 10**4, 10**5, 10**6]
EXPERIMENTS = {  # label -> (build container of size n, one operation, calls per timing)
    "list: -1 in c": (lambda n: list(range(n)), lambda c: -1 in c, lambda n: max(2, 2_000_000 // n)),
    "set: -1 in c": (lambda n: set(range(n)), lambda c: -1 in c, lambda n: 100_000),
    "list: insert(0)+pop(0)": (lambda n: list(range(n)), lambda c: (c.insert(0, 0), c.pop(0)), lambda n: 100),
    "deque: appendleft+popleft": (lambda n: deque(range(n)), lambda c: (c.appendleft(0), c.popleft()), lambda n: 100),
    "linked: prepend": (linked, lambda h: (0, h), lambda n: 100),
}
times = {k: [per_op(make(n), op, reps(n)) for n in NS] for k, (make, op, reps) in EXPERIMENTS.items()}
slopes = {k: loglog_slope(NS, ts) for k, ts in times.items()}
print("Part 2: fitted exponents (log-log slope of time against n):")
for k, sl in slopes.items():
    print(f"  {k:28s} slope {sl:5.2f}   time at n=10^6: {times[k][-1] * 1e9:12.0f} ns")
assert 0.8 < slopes["list: -1 in c"] < 1.25                      # O(n)
assert -0.15 < slopes["set: -1 in c"] < 0.35                     # O(1); cache misses add a little at 10^6
assert 0.7 < slopes["list: insert(0)+pop(0)"] < 1.3              # O(n)
assert abs(slopes["deque: appendleft+popleft"]) < 0.3 and abs(slopes["linked: prepend"]) < 0.3
assert times["list: -1 in c"][-1] > 100 * times["set: -1 in c"][-1]

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for ax, title, keys in ((axes[0], "Membership test", list(EXPERIMENTS)[:2]),
                        (axes[1], "Insert at the front, then remove it", list(EXPERIMENTS)[2:])):
    for k in keys:
        ax.loglog(NS, times[k], "o-", label=f"{k}, slope {slopes[k]:.2f}")
    ax.set(xlabel="n (elements)", ylabel="seconds per operation", title=title)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIG / "complexity.png", dpi=120)
plt.close(fig)

# %% [markdown]
# ## Part 3: A sparse polynomial as a dict (notes §3.4, build sequence 5)
#
# A polynomial `sum_k c_k x^k` is stored as `{k: c_k}` with no zero entries,
# so x^10000 + 1 costs two entries, not 10001.

# %%
def poly_clean(p: dict) -> dict:
    """Drop zero coefficients so equal polynomials are equal dicts."""
    return {k: c for k, c in p.items() if c != 0}


def poly_add(p: dict, q: dict) -> dict:
    """(p + q): coefficients add power by power."""
    r = dict(p)
    for k, c in q.items():
        r[k] = r.get(k, 0) + c
    return poly_clean(r)


def poly_mul(p: dict, q: dict) -> dict:
    """(p * q): c_i x^i times d_j x^j contributes c_i d_j to power i + j."""
    r = defaultdict(int)
    for i, a in p.items():
        for j, b in q.items():
            r[i + j] += a * b
    return poly_clean(r)


def poly_eval(p: dict, x: float) -> float:
    """Sparse Horner: walk powers downward, multiplying by x**gap between terms."""
    if not p:
        return 0
    result, prev = 0, max(p)
    for k in sorted(p, reverse=True):
        result = result * x ** (prev - k) + p[k]
        prev = k
    return result * x ** prev


def poly_deriv(p: dict) -> dict:
    """d/dx: c x^k -> k c x^(k-1); constants vanish."""
    return {k - 1: k * c for k, c in p.items() if k > 0}


def poly_integ(p: dict) -> dict:
    """Antiderivative with zero constant: c x^k -> c/(k+1) x^(k+1)."""
    return {k + 1: c / (k + 1) for k, c in p.items()}


def test_sparse_poly() -> None:
    p, q = {0: 1, 1: 2, 2: 3}, {0: -1, 1: 1}  # 1 + 2x + 3x^2 and x - 1
    assert poly_mul(p, q) == {0: -1, 1: -1, 2: -1, 3: 3}
    assert poly_add(p, {2: -3}) == {0: 1, 1: 2}          # the x^2 term cancels exactly
    assert poly_eval(p, 2) == 17 and poly_eval(q, 5) == 4
    assert poly_deriv(p) == {0: 2, 1: 6}
    lhs = poly_deriv(poly_mul(p, q))                      # product rule: (pq)' = p'q + pq'
    rhs = poly_add(poly_mul(poly_deriv(p), q), poly_mul(p, poly_deriv(q)))
    assert lhs == rhs == {0: -1, 1: -2, 2: 9}
    assert poly_deriv(poly_integ({0: 3.0, 4: 5.0})) == {0: 3.0, 4: 5.0}
    big = {10_000: 1, 0: 1}
    assert poly_mul(big, big) == {20_000: 1, 10_000: 2, 0: 1} and len(poly_mul(big, big)) == 3


test_sparse_poly()
print("Part 3: sparse polynomial tests passed.")

# %% [markdown]
# ## Part 4: CSV and JSON round trips (notes §3.10, §4.6)

# %%
csv_path, json_path = DATA / "summary.csv", DATA / "summary.json"
with csv_path.open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["sensor", "count", "mean", "min", "max"])
    w.writeheader()
    for s, row in summary.items():
        w.writerow({"sensor": s, **row})
json_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

with csv_path.open(newline="", encoding="utf-8") as f:
    from_csv = {r["sensor"]: {"count": int(r["count"]), "mean": float(r["mean"]),
                              "min": float(r["min"]), "max": float(r["max"])}
                for r in csv.DictReader(f)}
from_json = json.loads(json_path.read_text(encoding="utf-8"))
assert from_csv == summary and from_json == summary  # floats round-trip exactly via repr
print("Part 4: CSV and JSON round trips are exact.")

# %% [markdown]
# ## Part 5: Shallow vs deep copy (notes §3.3, §4.2)

# %%
orig = {"name": "run1", "grid": [[1, 2], [3, 4]]}
shallow, deep = copy.copy(orig), copy.deepcopy(orig)
orig["grid"][0][0] = 99
assert shallow["grid"] is orig["grid"] and shallow["grid"][0][0] == 99
assert deep["grid"] is not orig["grid"] and deep["grid"][0][0] == 1
trap = [[0] * 3] * 3
trap[0][0] = 1
assert trap == [[1, 0, 0]] * 3 and trap[0] is trap[2]
fixed = [[0] * 3 for _ in range(3)]
fixed[0][0] = 1
assert fixed == [[1, 0, 0], [0, 0, 0], [0, 0, 0]]
print("Part 5: copy demonstrations hold.")

# %% [markdown]
# ## Part 6: The iteration protocol (notes §3.7-3.8, build sequences 3 and 6)

# %%
def my_range(*args):
    """range(stop) / range(start, stop[, step]) as a generator."""
    start, stop, step = (0, *args, 1)[-3:] if len(args) == 1 else (*args, 1)[:3]
    if step == 0:
        raise ValueError("step must not be zero")
    i = start
    while (i < stop) if step > 0 else (i > stop):
        yield i
        i += step


def my_enumerate(iterable, start=0):
    """Yield (index, item) pairs."""
    i = start
    for item in iterable:
        yield i, item
        i += 1


_END = object()  # a sentinel no iterator can yield


def my_zip(*iterables, strict=False):
    """Stop at the shortest input; with strict=True raise ValueError on unequal lengths."""
    iterators = [iter(it) for it in iterables]
    while iterators:
        row = []
        for it in iterators:
            item = next(it, _END)
            if item is _END:
                break
            row.append(item)
        if len(row) < len(iterators):  # someone ran out
            if strict and (row or any(next(it, _END) is not _END for it in iterators[1:])):
                raise ValueError("zip() arguments have different lengths")
            return
        yield tuple(row)


def for_by_hand(iterable):
    """`for x in iterable: out.append(x)` desugared into iter/next/StopIteration."""
    out, it = [], iter(iterable)
    while True:
        try:
            x = next(it)
        except StopIteration:
            break
        out.append(x)
    return out


rnd = random.Random(1)
for _ in range(200):
    a, b, c = rnd.randint(-20, 20), rnd.randint(-20, 20), rnd.choice([-3, -2, -1, 1, 2, 5])
    assert list(my_range(a)) == list(range(a))
    assert list(my_range(a, b)) == list(range(a, b))
    assert list(my_range(a, b, c)) == list(range(a, b, c))
    xs, ys = rnd.sample(range(100), rnd.randint(0, 6)), "abcde"[: rnd.randint(0, 5)]
    assert list(my_enumerate(xs, 7)) == list(enumerate(xs, 7))
    assert list(my_zip(xs, ys)) == list(zip(xs, ys))
    assert for_by_hand(xs) == xs and for_by_hand(iter(ys)) == list(ys)
    if len(xs) == len(ys):
        assert list(my_zip(xs, ys, strict=True)) == list(zip(xs, ys, strict=True))
    else:
        for z in (my_zip, zip):
            try:
                list(z(xs, ys, strict=True))
            except ValueError:
                pass
            else:
                raise AssertionError("strict zip must raise on unequal lengths")
gen = (x * x for x in range(4))
assert list(gen) == [0, 1, 4, 9] and list(gen) == []  # generators are consumed once
print("Part 6: my_range, my_enumerate, my_zip, for_by_hand match the built-ins.")

# %% [markdown]
# ## Part 7: Every number in the notes' worked examples (notes §4)

# %%
# 4.1 slicing
s = "abcdefghij"
assert (s[2:7], s[:3], s[-3:], s[::2], s[1:8:3]) == ("cdefg", "abc", "hij", "acegi", "beh")
assert (s[5:2:-1], s[::-1], s[-2:-6:-2], s[2:5:-1]) == ("fed", "jihgfedcba", "ig", "")
assert slice(None, None, -1).indices(10) == (9, -1, -1) and slice(None, None, 1).indices(10) == (0, 10, 1)
lst10 = list(range(10))
lst10[2:5] = ["x"]  # slice assignment may change the length
assert lst10 == [0, 1, "x", 5, 6, 7, 8, 9]

# 4.3 sorting records: score descending, then name
records_ns = [("Ana", 82), ("Raj", 91), ("Lee", 82), ("Kim", 91), ("Bo", 75), ("Zed", 82)]
by_tuple = sorted(records_ns, key=lambda r: (-r[1], r[0]))
two_pass = sorted(sorted(records_ns, key=itemgetter(0)), key=itemgetter(1), reverse=True)
assert by_tuple == two_pass == [("Kim", 91), ("Raj", 91), ("Ana", 82), ("Lee", 82), ("Zed", 82), ("Bo", 75)]
assert sorted(records_ns, key=itemgetter(1), reverse=True)[:2] == [("Raj", 91), ("Kim", 91)]  # stable

# 4.4 groupby vs defaultdict
readings = [("b", 3.0), ("a", 1.0), ("b", 5.0), ("a", 2.0), ("c", 7.0), ("a", 4.0)]
dd = defaultdict(list)
for k, v in readings:
    dd[k].append(v)
assert dict(dd) == {"b": [3.0, 5.0], "a": [1.0, 2.0, 4.0], "c": [7.0]}
gb = {k: [v for _, v in g] for k, g in itertools.groupby(sorted(readings, key=itemgetter(0)), key=itemgetter(0))}
assert gb == {"a": [1.0, 2.0, 4.0], "b": [3.0, 5.0], "c": [7.0]}
unsorted_groups = [(k, [v for _, v in g]) for k, g in itertools.groupby(readings, key=itemgetter(0))]
assert [k for k, _ in unsorted_groups] == ["b", "a", "b", "a", "c", "a"]
stale = list(itertools.groupby(readings, key=itemgetter(0)))
assert [list(g) for _, g in stale][:-1] == [[]] * 5 and list(stale[-1][1]) == []  # inner iterators are spent

# 4.5 a 10^6-line file, filtered and summed without holding it in memory
N_BIG, THRESH = 1_000_000, 900
big = DATA / "big.txt"
with big.open("w", encoding="utf-8") as f:
    f.writelines(f"{i},{(i * 7919) % 1000}\n" for i in range(N_BIG))


def peak_kib(fn):
    tracemalloc.start()
    result = fn()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return result, peak / 1024


def lazy_total():
    values = (int(line.split(",")[1]) for line in read_lines(big))
    kept_sum, kept_n = 0, 0
    for v in values:
        if v >= THRESH:
            kept_sum += v
            kept_n += 1
    return kept_n, kept_sum


def eager_total():
    lines = big.read_text(encoding="utf-8").splitlines()          # whole file as one str, then a list of str
    values = [int(line.split(",")[1]) for line in lines]          # a second list
    kept = [v for v in values if v >= THRESH]                     # a third
    return len(kept), sum(kept)


(lazy_res, lazy_kib), (eager_res, eager_kib) = peak_kib(lazy_total), peak_kib(eager_total)
expected = (100_000, 1000 * sum(range(900, 1000)))  # each residue 0..999 occurs 1000 times
assert lazy_res == eager_res == expected == (100_000, 94_950_000)
assert eager_kib > 100 * lazy_kib
print(f"Part 7: 10^6-line file -> {lazy_res}; peak memory lazy {lazy_kib:.0f} KiB vs eager {eager_kib:.0f} KiB")

# 4.6 JSON round trip
data = {"run": 7, "t": (0.0, 0.1), "x": 0.1 + 0.2, "tiny": 1e-7, "big": 1e22, 3: "three", "ok": None}
text = json.dumps(data)
assert text == ('{"run": 7, "t": [0.0, 0.1], "x": 0.30000000000000004, "tiny": 1e-07, '
                '"big": 1e+22, "3": "three", "ok": null}')
back = json.loads(text)
assert back["x"] == data["x"] and back["tiny"] == data["tiny"] and back["big"] == data["big"]
assert back["t"] == [0.0, 0.1] and back["t"] != data["t"] and 3 not in back and back["3"] == "three"
nan_text = json.dumps({"v": float("nan")})
assert nan_text == '{"v": NaN}'
nan_back = json.loads(nan_text)["v"]
assert nan_back != nan_back                      # nan != nan, even to itself
assert {"v": float("nan")} != json.loads(nan_text)  # so the "same" data is not equal after a round trip
try:
    json.dumps(float("nan"), allow_nan=False)
except ValueError:
    pass
else:
    raise AssertionError("allow_nan=False must raise")

# 4.8 a parameter grid, lazily
masses, dampers, springs = (0.5, 1.0, 2.0), (0.5, 1.0, 2.0, 4.0), (1.0, 2.0, 8.0)
grid = itertools.product(masses, dampers, springs)


def regime(m: float, c: float, k: float) -> str:
    """Damping ratio zeta = c / (2 sqrt(k m)) decides under-, critically, or over-damped."""
    zeta = c / (2 * math.sqrt(k * m))
    return "critical" if math.isclose(zeta, 1.0) else ("under" if zeta < 1 else "over")


assert list(itertools.islice(itertools.product(masses, dampers, springs), 3)) == [
    (0.5, 0.5, 1.0), (0.5, 0.5, 2.0), (0.5, 0.5, 8.0)]
regimes = Counter(regime(*p) for p in grid)
assert sum(regimes.values()) == 36 == len(masses) * len(dampers) * len(springs)
print(f"Part 7: regimes over the 36-point grid: {dict(regimes)}")
assert dict(regimes) == {"under": 26, "over": 6, "critical": 4}

# %% [markdown]
# ## Part 8: Word frequency three ways (notes §3.4, build sequence 4)

# %%
def freq_lists(words):
    """Nested loops over two parallel lists: a linear search per word, O(n * distinct)."""
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
    """dict + get: one hash lookup per word."""
    counts = {}
    for w in words:
        counts[w] = counts.get(w, 0) + 1
    return counts


def freq_counter(words):
    """Counter: the same loop, written in C."""
    return Counter(words)


rng_w = random.Random(2)
vocab = [f"w{i}" for i in range(400)]
weights = [1 / (i + 1) for i in range(400)]  # Zipf-like: a few common words, a long tail
words = rng_w.choices(vocab, weights, k=40_000)
timings = {}
results = {}
for fn in (freq_lists, freq_dict, freq_counter):
    t0 = time.perf_counter()
    results[fn.__name__] = fn(words)
    timings[fn.__name__] = time.perf_counter() - t0
assert results["freq_lists"] == results["freq_dict"] == dict(results["freq_counter"])
assert sum(results["freq_dict"].values()) == 40_000
assert timings["freq_lists"] > 3 * timings["freq_dict"]
print("Part 8: seconds -> " + ", ".join(f"{k}: {v:.4f}" for k, v in timings.items()))

# %% [markdown]
# ## Part 9: Cross-checks against NumPy (the only place NumPy appears)

# %%
import numpy as np
from numpy.polynomial import Polynomial

names = np.array([r[1] for r in records])
vals = np.array([r[2] for r in records])
for s in SENSORS:
    v = vals[names == s]
    assert len(v) == summary[s]["count"]
    assert np.isclose(v.mean(), summary[s]["mean"], rtol=1e-12)
    assert v.min() == summary[s]["min"] and v.max() == summary[s]["max"]


def dense(p: dict) -> list:
    """Dict {power: coeff} -> dense coefficient list, lowest power first."""
    return [p.get(k, 0) for k in range(max(p) + 1)] if p else [0]


rng_np = np.random.default_rng(0)
for _ in range(50):
    pa = {int(k): int(c) for k, c in zip(rng_np.integers(0, 12, 4), rng_np.integers(-5, 6, 4)) if c}
    pb = {int(k): int(c) for k, c in zip(rng_np.integers(0, 12, 3), rng_np.integers(-5, 6, 3)) if c}
    Pa, Pb = Polynomial(dense(pa)), Polynomial(dense(pb))
    x0 = float(rng_np.uniform(-1.5, 1.5))
    assert np.allclose(dense(poly_mul(pa, pb)), (Pa * Pb).coef, atol=1e-9) or not (pa and pb)
    assert np.allclose(dense(poly_add(pa, pb)), (Pa + Pb).coef, atol=1e-9) or poly_add(pa, pb) == {}
    assert math.isclose(poly_eval(pa, x0), Pa(x0), rel_tol=1e-9, abs_tol=1e-9)
    if len(pa) > 1 or (pa and max(pa) > 0):
        assert np.allclose(dense(poly_deriv(pa)), Pa.deriv().coef, atol=1e-9)
print("Part 9: NumPy agrees with the dict-based statistics and polynomials.")

# %% [markdown]
# ## Your turn
#
# 1. **Anagram grouper.** Write `group_anagrams(words)` returning lists of words that are
#    anagrams of each other, using a `Counter` of letters turned into a hashable key.
# 2. **Lazy grid.** Use `itertools.product` to enumerate a three-parameter grid without
#    building it as a list, and count how many points satisfy a condition.
# 3. **Rewrite.** Turn one nested-loop solution into a comprehension and one comprehension
#    back into a loop. Decide which version is clearer, and write one sentence why.
# 4. Exercise B7: add `list.pop(0)` vs `deque.popleft()` to the Part 2 experiment.
#
# Reference solutions follow. Try first.

# %%
def group_anagrams(words):
    """Group words whose letter multisets are equal; key = frozenset of (letter, count) pairs."""
    groups = defaultdict(list)
    for w in words:
        groups[frozenset(Counter(w).items())].append(w)
    return list(groups.values())


assert group_anagrams(["listen", "silent", "enlist", "google", "gogole", "cat", "act", "tac", "dog"]) == [
    ["listen", "silent", "enlist"], ["google", "gogole"], ["cat", "act", "tac"], ["dog"]]

pairs_loop = []  # nested loops: all (i, j) with i < j <= 6 and i + j prime
for i in range(1, 7):
    for j in range(i + 1, 7):
        if i + j in (2, 3, 5, 7, 11):
            pairs_loop.append((i, j))
pairs_comp = [(i, j) for i in range(1, 7) for j in range(i + 1, 7) if i + j in (2, 3, 5, 7, 11)]
assert pairs_loop == pairs_comp and len(pairs_comp) == 7

# A comprehension with two conditions and a nested ternary reads worse than the loop:
labels_comp = ["small" if x < 3 else "mid" if x < 6 else "big" for x in range(8) if x % 2 == 0 or x > 5]
labels_loop = []
for x in range(8):
    if x % 2 == 0 or x > 5:
        if x < 3:
            labels_loop.append("small")
        elif x < 6:
            labels_loop.append("mid")
        else:
            labels_loop.append("big")
assert labels_comp == labels_loop  # here the loop wins: elif chains do not fit in an expression

# %%
if __name__ == "__main__":
    print("All checks passed.")
