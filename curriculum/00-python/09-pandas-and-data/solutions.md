# PY-9 Solutions

Every block below was run, in order, under Python 3.13, pandas 3.0.5 and NumPy 2.5.3, and its output pasted from the run; its asserts passed. Each solution starts from a fresh session in the module folder (so `import project` works) after this preamble:

```python
import contextlib
import io
import time
import warnings

import numpy as np
import pandas as pd


def best_of(f, repeat: int = 3) -> float:
    """Smallest wall time of `repeat` calls of f(), in seconds."""
    times = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        f()
        times.append(time.perf_counter() - t0)
    return min(times)
```

Timings are from one run; only their ratios transfer.

## A. Warm-up

**A1.**

```python
s = pd.Series([10, 20, 30, 40, 50], index=[5, 4, 3, 2, 1])
print(s.loc[4:2].tolist(), s.iloc[4:2].tolist(), s.iloc[1:3].tolist())
print(s[1], s[1:3].tolist(), s.loc[[1, 5]].tolist(), s.sort_index().loc[2:4].tolist())
```
```text
[20, 30, 40] [] [20, 30]
50 [20, 30] [50, 10] [40, 30, 20]
```

- `s.loc[4:2]` runs from the *position* of label 4 to the position of label 2, both included: three values. The labels decrease, but a label slice follows the row order, not the numeric order of the labels.
- `s.iloc[4:2]` is positions 4 up to (not including) 2: empty, as for a list.
- `s.iloc[1:3]` is positions 1 and 2.
- `s[1]` with an integer index is a **label** lookup: the value labelled 1 is 50.
- `s[1:3]` with integers is a **position** slice: 20, 30. Bare `[]` changes meaning between a scalar and a slice; this is why the notes never use it on rows.
- `s.loc[[1, 5]]` returns the labels in the order asked: 50, 10.
- After `sort_index()` the labels run 1..5, and `loc[2:4]` is labels 2, 3, 4: 40, 30, 20.

**A2.**

```python
df = pd.DataFrame({"a": [1.0, 2.0], "b": [3.0, 4.0]}, index=["x", "y"])
print(df - df.mean())
print(df - df.mean(axis=1))


def center_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Subtract each row's mean from that row: align the row means with the row index (axis=0)."""
    return df.sub(df.mean(axis=1), axis=0)


centered = center_rows(df)
print(centered)
assert np.allclose(centered.sum(axis=1), 0.0)
```
```text
     a    b
x -0.5 -0.5
y  0.5  0.5
    a   b   x   y
x NaN NaN NaN NaN
y NaN NaN NaN NaN
     a    b
x -1.0  1.0
y -1.0  1.0
```

`df.mean()` is indexed by the columns `a, b`; subtracting it aligns its labels with the columns, so each *column* is centred. `df.mean(axis=1)` is indexed by the rows `x, y`, but `DataFrame - Series` still aligns the Series with the *columns*. The union of `{a, b}` and `{x, y}` gives four columns, and no label is on both sides, so everything is `NaN`: a 2×4 frame of nothing, without an error. `sub(..., axis=0)` says "align with the rows", which centres each row.

**A3.**

```python
csv_text = "machine,value\nA,1.0\nA,\nB,3.0\nB,500.0\n"
df = pd.read_csv(io.StringIO(csv_text))
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    df[df["machine"] == "B"]["value"] = df["value"] / 10
    df["value"].fillna(0, inplace=True)
df["big"] = df.apply(lambda row: row["value"] > 2, axis=1)
print([w.category.__name__ for w in caught])
print(df.to_dict("list"))

fixed = (
    pd.read_csv(io.StringIO(csv_text))
    # machine B logs tenths: rescale its rows only
    .assign(value=lambda d: d["value"].where(d["machine"] != "B", d["value"] / 10))
    # vectorized comparison; a missing reading gives a missing answer, not False
    .assign(big=lambda d: d["value"].gt(2).astype("boolean").mask(d["value"].isna()))
)
print(fixed["value"].tolist(), fixed["big"].tolist(), fixed["big"].dtype)
assert fixed["value"].tolist()[2:] == [0.3, 50.0] and fixed["value"].isna().sum() == 1
```
```text
['ChainedAssignmentError', 'ChainedAssignmentError']
{'machine': ['A', 'A', 'B', 'B'], 'value': [1.0, nan, 3.0, 500.0], 'big': [False, False, True, True]}
[1.0, nan, 0.3, 50.0] [False, <NA>, False, True] boolean
```

Line 2 is chained assignment: `df[mask]` is a new frame, and the assignment changes only that, so pandas 3 warns and `df` keeps 3.0 and 500.0. Line 3 is chained too (`df["value"]` then an in-place method on it), so it also does nothing; and even if it worked, filling a missing temperature with 0 invents a reading of zero, which every later mean would include. Line 4 runs a Python function per row, and `NaN > 2` silently answers `False`. The chain rescales with one `where`, compares with one `gt`, and keeps the unknown as `<NA>` in a nullable `boolean` column.

**A4.**

```python
d = pd.DataFrame({"g": ["a", "a", "b", None, "b", "b"], "v": [1.0, 1.0, np.nan, 2.0, 3.0, np.nan]})
print(d.groupby("g").size().to_dict(), d.groupby("g")["v"].count().to_dict(), d.groupby("g")["v"].nunique().to_dict())
print(d.groupby("g", dropna=False).size().to_dict(), d["g"].value_counts(dropna=False).to_dict())
```
```text
{'a': 2, 'b': 3} {'a': 2, 'b': 1} {'a': 1, 'b': 1}
{'a': 2, 'b': 3, nan: 1} {'b': 3, 'a': 2, nan: 1}
```

`size` counts rows (a: 2, b: 3); `count` counts non-missing values (b has two `NaN`s, so 1); `nunique` counts distinct non-missing values (a's two readings are both 1.0). The row with a missing key is dropped by `groupby` unless `dropna=False`, and `value_counts` keeps it only with `dropna=False`. "Readings per sensor" is `size`; "valid readings per sensor" is `count`. Report the difference too: it is the missing count.

**A5.**

```python
r = pd.Series([1.0, 2.0, 3.0, 4.0], index=pd.to_datetime(["2024-03-04 08:05", "2024-03-04 08:55",
                                                         "2024-03-04 09:00", "2024-03-04 11:30"]))
outs = {"resample mean": r.resample("1h").mean(), "resample size": r.resample("1h").size(),
        "floor + groupby": r.groupby(r.index.floor("h")).mean(),
        "closed/label right": r.resample("1h", closed="right", label="right").mean()}
for name, out in outs.items():
    print(f"{name:18s}", {k.strftime("%H:%M"): v for k, v in out.items()})
```
```text
resample mean      {'08:00': 1.5, '09:00': 3.0, '10:00': nan, '11:00': 4.0}
resample size      {'08:00': 2, '09:00': 1, '10:00': 0, '11:00': 1}
floor + groupby    {'08:00': 1.5, '09:00': 3.0, '11:00': 4.0}
closed/label right {'09:00': 2.0, '10:00': nan, '11:00': nan, '12:00': 4.0}
```

Default bins are $[h, h+1)$ labelled by $h$: 08:05 and 08:55 average to 1.5 at 08:00, 09:00 starts its own bin, 10:00 is empty and appears as `NaN` (size 0), 11:30 is in 11:00. `floor` + `groupby` gives the same means but no row for the empty hour. With `closed="right"` the bins are $(h-1, h]$ labelled by the right edge: the reading at exactly 09:00 moves into the bin with 08:05 and 08:55 (mean 2.0, labelled 09:00), and two empty bins appear before 11:30's bin, labelled 12:00.

## B. Core

**B1.**

```python
rng = np.random.default_rng(9)
n = 200_000
log = pd.DataFrame({"machine": rng.choice(["A", "B", "C"], n, p=[0.5, 0.3, 0.2]),
                    "sensor": rng.choice(["temp_1", "temp_2", "temp_3", "temp_4"], n)})
log["value"] = (log["machine"].map({"A": 60.0, "B": 65.0, "C": 58.0})
                + log["sensor"].str[-1].astype(int) - 1 + rng.normal(0, 2.0, n))


def group_stats(k1: pd.Series, k2: pd.Series, values: np.ndarray) -> pd.DataFrame:
    """count, mean, std (ddof 1, two-pass), min, max per (k1, k2) pair, by integer codes and scatter operations."""
    c1, u1 = pd.factorize(k1, sort=True)
    c2, u2 = pd.factorize(k2, sort=True)
    code = c1 * len(u2) + c2                                   # one integer per (k1, k2) pair
    G = len(u1) * len(u2)
    cnt = np.bincount(code, minlength=G)
    mean = np.bincount(code, weights=values, minlength=G) / np.maximum(cnt, 1)
    ss = np.bincount(code, weights=(values - mean[code]) ** 2, minlength=G)    # second pass
    std = np.sqrt(ss / np.maximum(cnt - 1, 1))
    lo, hi = np.full(G, np.inf), np.full(G, -np.inf)
    np.minimum.at(lo, code, values)
    np.maximum.at(hi, code, values)
    out = pd.DataFrame({"count": cnt, "mean": mean, "std": std, "min": lo, "max": hi},
                       index=pd.MultiIndex.from_product([u1, u2], names=[k1.name, k2.name]))
    return out.loc[out["count"] > 0]


mine = group_stats(log["machine"], log["sensor"], log["value"].to_numpy())
ref = log.groupby(["machine", "sensor"])["value"].agg(["count", "mean", "std", "min", "max"])
pd.testing.assert_frame_equal(mine, ref, check_exact=False, rtol=1e-9)
t_mine = best_of(lambda: group_stats(log["machine"], log["sensor"], log["value"].to_numpy()))
t_ref = best_of(lambda: log.groupby(["machine", "sensor"])["value"].agg(["count", "mean", "std", "min", "max"]))
print(f"by hand {t_mine * 1e3:.1f} ms, groupby {t_ref * 1e3:.1f} ms; {len(mine)} groups agree to 1e-9")
print(mine.iloc[:2].round(4))

x = 1e9 + rng.normal(0, 1.0, 10_000)                          # e.g. epoch seconds with 1 s jitter
one_pass = (np.sum(x**2) - np.sum(x) ** 2 / len(x)) / (len(x) - 1)
two_pass = np.sum((x - x.mean()) ** 2) / (len(x) - 1)
print(f"variance of 1e9 + N(0, 1): one-pass {one_pass:.4g}, two-pass {two_pass:.6f}")
```
```text
by hand 19.7 ms, groupby 22.1 ms; 12 groups agree to 1e-9
                count     mean     std      min      max
machine sensor
A       temp_1  25081  60.0049  2.0003  51.3416  70.3789
        temp_2  25104  61.0179  2.0041  52.3700  69.5645
variance of 1e9 + N(0, 1): one-pass 0, two-pass 1.017263
```

The codes `c1 * len(u2) + c2` number the pairs like the cells of a `len(u1) × len(u2)` array, so one `bincount` per statistic does all the groups at once; `np.minimum.at` is the unbuffered scatter for min and max (PY-6 §3.4). I used the two-pass formula: first the group means, then $\sum (x_i - \bar x_g)^2$ through `mean[code]`. The one-pass formula $\big(\sum x^2 - (\sum x)^2/n\big)/(n-1)$ subtracts two numbers of size $n\bar x^2$ to get one of size $n s^2$; it loses about $\log_{10}(\bar x^2/s^2)$ of the 16 digits. For temperatures (mean 60, sd 2) that is 3 digits, harmless; for $10^9 \pm 1$, such as timestamps in seconds, it is all 16 digits and more, and the result above is garbage. The hand version is competitive here because pandas' `agg` with a list runs one pass per statistic too.

**B2.**

```python
def merge_by_hand(left: pd.DataFrame, right: pd.DataFrame, on: str, how: str = "inner") -> list[dict]:
    """Relational join on one key column by nested loops; `_merge` as with indicator=True."""
    L, R = left.to_dict("records"), right.to_dict("records")
    no_left = {c: None for c in left.columns if c != on}
    no_right = {c: None for c in right.columns if c != on}
    out, matched_right = [], set()
    for row in L:
        hits = [j for j, r in enumerate(R) if r[on] == row[on]]
        matched_right.update(hits)
        out += [{**row, **R[j], "_merge": "both"} for j in hits]
        if not hits and how in ("left", "outer"):
            out.append({**row, **no_right, "_merge": "left_only"})
    if how in ("right", "outer"):
        out += [{**no_left, **r, "_merge": "right_only"} for j, r in enumerate(R) if j not in matched_right]
    return out


def check_validate(left: pd.DataFrame, right: pd.DataFrame, on: str, kind: str) -> None:
    """Raise ValueError exactly when pd.merge(..., validate=kind) would."""
    if kind in ("one_to_one", "one_to_many") and not left[on].is_unique:
        raise ValueError(f"left keys not unique: not {kind}")
    if kind in ("one_to_one", "many_to_one") and not right[on].is_unique:
        raise ValueError(f"right keys not unique: not {kind}")


def canonical(rows) -> list[tuple]:
    """Rows as sorted tuples, with every kind of missing value as None, for order-free comparison."""
    df = pd.DataFrame(rows, columns=["k", "x", "y", "_merge"]).astype(object)
    tuples = [tuple(None if pd.isna(v) else v for v in r) for r in df.itertuples(index=False)]
    return sorted(tuples, key=repr)


rng = np.random.default_rng(12)
kinds = ("one_to_one", "one_to_many", "many_to_one", "many_to_many")
n_raised = 0
for trial in range(300):
    nl, nr = rng.integers(0, 7, size=2)
    left = pd.DataFrame({"k": rng.choice(list("abc"), nl), "x": rng.integers(0, 100, nl)})
    right = pd.DataFrame({"k": rng.choice(list("abc"), nr), "y": rng.integers(0, 100, nr)})
    for how in ("inner", "left", "right", "outer"):
        ref = pd.merge(left, right, on="k", how=how, indicator=True).astype({"_merge": str})
        assert canonical(merge_by_hand(left, right, "k", how)) == canonical(ref), (trial, how)
    for kind in kinds:
        try:
            pd.merge(left, right, on="k", validate=kind)
            pandas_raised = False
        except pd.errors.MergeError:
            pandas_raised = True
        try:
            check_validate(left, right, "k", kind)
            mine_raised = False
        except ValueError:
            mine_raised = True
        assert mine_raised == pandas_raised, (trial, kind)
        n_raised += pandas_raised
print(f"300 random pairs x 4 joins and x 4 validate kinds agree; validate raised {n_raised} times")

left = pd.DataFrame({"k": list("aab"), "x": [1, 2, 3]})
right = pd.DataFrame({"k": list("aaac"), "y": [4, 5, 6, 7]})
for how in ("inner", "left", "right", "outer"):
    print(how, len(merge_by_hand(left, right, "k", how)), len(pd.merge(left, right, on="k", how=how)))
```
```text
300 random pairs x 4 joins and x 4 validate kinds agree; validate raised 595 times
inner 6 6
left 7 7
right 7 7
outer 8 8
```

The comparison sorts both results because row order is not part of what a join defines: pandas keeps the left frame's order for `inner` and `left`, the right frame's for `right`, and sorts the keys for `outer`, while my loop always follows the left frame; the *multiset* of rows is the definition. `MergeError` is a subclass of `ValueError`, so the two checks raise the same kind of exception. The row counts follow from the nested loop: with $m_L(k)$ and $m_R(k)$ the multiplicities of key $k$,

$$|\text{inner}|=\sum_k m_L(k)\,m_R(k),\qquad |\text{left}|=|\text{inner}|+\sum_{k:\,m_R(k)=0} m_L(k),$$

and symmetrically for `right`; `outer` adds both unmatched sums. For the last pair, $a$ appears 2 and 3 times, so inner has $2\cdot3=6$ rows; `left` adds the unmatched `b` (7), `right` adds `c` (7), `outer` both (8).

**B3.**

```python
wide = pd.DataFrame({"date": pd.date_range("2024-03-04", periods=3, freq="D"),
                     "A_temp_mean": [53.0, 53.1, 53.4], "A_temp_max": [62.3, 62.8, 61.6],
                     "B_temp_mean": [54.1, 54.6, 54.2], "B_temp_max": [63.2, 63.0, 62.7],
                     "B_vib_mean": [0.94, 0.95, 0.93], "B_vib_max": [1.41, 1.38, 1.40]})
long = (
    wide.melt(id_vars="date", var_name="col", value_name="value")              # headers are values: melt
    .pipe(lambda d: d.join(d["col"].str.split("_", expand=True)                  # three variables in one string
                            .set_axis(["machine", "kind", "stat"], axis=1)))
    .loc[:, ["date", "machine", "kind", "stat", "value"]]
)
tidy = (
    long.pivot(index=["date", "machine", "kind"], columns="stat", values="value")   # stat is two variables
    .reset_index()
    .rename_axis(columns=None)
    .loc[:, ["date", "machine", "kind", "mean", "max"]]
)
print(long.shape, tidy.shape)
print(tidy.head(3))
rebuilt = (
    tidy.melt(id_vars=["date", "machine", "kind"], var_name="stat", value_name="value")
    .assign(col=lambda d: d["machine"] + "_" + d["kind"] + "_" + d["stat"])
    .pivot(index="date", columns="col", values="value")
    .loc[:, wide.columns[1:]]
    .reset_index()
    .rename_axis(columns=None)
)
pd.testing.assert_frame_equal(rebuilt, wide)
print("wide -> long -> tidy -> wide: identical")
```
```text
(18, 5) (9, 5)
        date machine  kind   mean    max
0 2024-03-04       A  temp  53.00  62.30
1 2024-03-04       B  temp  54.10  63.20
2 2024-03-04       B   vib   0.94   1.41
wide -> long -> tidy -> wide: identical
```

Three messes at once: the column headers are values (machine, kind and statistic: mess 1), each header packs three variables into one string (mess 2), and the statistic is a variable stored across columns while `mean` and `max` are really two measured variables of one observation (mess 3). The fixes are `melt`, `str.split(expand=True)`, and a `pivot` of the `stat` column. The tidy table has one row per (date, machine, sensor kind): 3 dates × 3 machine-kind pairs = 9. `pivot` is safe here because every (date, machine, kind, stat) occurs once; with duplicates you would need `pivot_table` and a decision about how to aggregate.

**B4.**

```python
rng = np.random.default_rng(3)
day = pd.Timestamp("2024-03-04")


def jittered_sensor(step: str, gain: float, offset: float, sigma: float) -> pd.Series:
    """offset + gain * L(t) + noise at a regular grid jittered by U(-2, 2) min, rounded to the second."""
    base = pd.date_range(day, day + pd.Timedelta("1D"), freq=step, inclusive="left")
    t = (base + pd.to_timedelta(rng.uniform(-120, 120, len(base)), unit="s")).round("s")
    hours = (t - day) / pd.Timedelta("1h")
    load = 0.5 + 0.3 * np.sin(2 * np.pi * hours / 24)
    return pd.Series(offset + gain * load + rng.normal(0, sigma, len(t)), index=t).sort_index()


temp = jittered_sensor("10min", 25.0, 40.0, 0.5).rename("temp")
pres = jittered_sensor("15min", 3.0, 2.0, 0.05).rename("pres")
exact = pd.concat([temp, pres], axis=1, join="inner")
binned = pd.concat([temp.resample("30min").mean(), pres.resample("30min").mean()], axis=1).dropna()
asof = pd.merge_asof(pres.rename_axis("time").reset_index(), temp.rename_axis("time").reset_index(),
                     on="time", direction="nearest", tolerance=pd.Timedelta("5min")).dropna()
for name, tab in (("exact join", exact), ("30-min resample", binned), ("merge_asof, 5 min", asof)):
    r = tab["temp"].corr(tab["pres"]) if len(tab) > 2 else float("nan")
    print(f"{name:18s} {len(tab):4d} pairs, corr {r:.4f}")
```
```text
exact join            2 pairs, corr nan
30-min resample      48 pairs, corr 0.9975
merge_asof, 5 min    89 pairs, corr 0.9936
```

With jitter, two sensors almost never share a timestamp to the second: the exact join found 2 chance coincidences out of 96 pressure readings, too few for a correlation (and without jitter it would keep only the stamps on the half hour). Resampling to a common grid uses every reading and averages noise away (three temperatures and two pressures per bin), so the correlation rises; use it when the question is about the slow signal and you want a regular table for regression or plotting (`project.py` Part 7). `merge_asof` pairs each pressure reading with the nearest temperature within a tolerance, keeping individual readings and their timing; use it for event-level questions (what was the temperature when this pressure reading, alarm or sample was taken?), and choose the tolerance so that no pair is staler than the process allows. `resample` does not care about row order; `merge_asof` needs both tables sorted on the key and refuses unsorted ones.

**B5.**

```python
with contextlib.redirect_stdout(io.StringIO()):
    import project as P                        # runs the project's checks once, quietly


def add_z_scores(clean: pd.DataFrame) -> pd.DataFrame:
    """Global and per-day z-scores within (machine, sensor), by groupby().transform (notes 3.8)."""
    g = clean.groupby(["machine", "sensor"], observed=True)["value"]
    gd = clean.groupby(["machine", "sensor", pd.Grouper(key="time", freq="D")], observed=True)["value"]
    return clean.assign(z_global=(clean["value"] - g.transform("mean")) / g.transform("std"),
                        z_daily=(clean["value"] - gd.transform("mean")) / gd.transform("std"))


rows = []
for size in (30.0, 15.0, 8.0):
    P.BURST = (*P.BURST[:4], size)
    rng = np.random.default_rng(P.SEED)
    raw, planted = P.plant_problems(P.make_truth(rng), rng)
    clean = add_z_scores(P.flag_rolling_mad(P.clean_log(raw)))
    burst = (clean["machine"] == "B") & (clean["sensor"] == "temp_3") & clean["time"].isin(planted["burst_times"])
    rules = {"global |z|>4": clean["z_global"].abs() > 4, "daily |z|>4": clean["z_daily"].abs() > 4,
             "Hampel k=10": clean["outlier"]}
    for rule, flag in rules.items():
        rows.append({"burst °C": size, "rule": rule, "burst rows flagged": int((flag & burst).sum()),
                     "other rows flagged": int((flag & ~burst).sum())})
    print(f"burst {size:4.0f}: z of burst rows, global {clean.loc[burst, 'z_global'].mean():.2f}, "
          f"daily {clean.loc[burst, 'z_daily'].mean():.2f}, Hampel {clean.loc[burst, 'robust_z'].mean():.2f}")
table = pd.DataFrame(rows).pivot(index="rule", columns="burst °C", values=["burst rows flagged", "other rows flagged"])
print(table)
```
```text
burst   30: z of burst rows, global 6.05, daily 4.14, Hampel 25.69
burst   15: z of burst rows, global 3.74, daily 3.43, Hampel 13.17
burst    8: z of burst rows, global 2.58, daily 2.74, Hampel 7.33
             burst rows flagged           other rows flagged
burst °C                   8.0  15.0 30.0               8.0  15.0 30.0
rule
Hampel k=10                   0    6    6                  0    0    0
daily |z|>4                   0    0    6                  0    0    0
global |z|>4                  0    0    6                  0    0    0
```

To add this to `project.py`, put `add_z_scores` after Part 2 and the loop in a new cell; `plant_problems` reads `BURST` at call time, so rebinding the module constant is enough. The table is the answer. The daily z-score of the 30 °C burst is *smaller* than the global one because the burst sits inside the group whose mean and standard deviation it is measured against: six readings of +30 °C in a day of 144 raise that day's standard deviation far more than they raise the two-week one, so they hide themselves. This is called **masking**, and it is why a mean and a standard deviation are the wrong yardstick for the outliers that contaminate them. The Hampel rule uses a median and a MAD, which six points out of 25 barely move, and it measures each reading against its own four-hour neighbourhood, so the daily cycle does not inflate the yardstick either. Its threshold of 10 was set from the clean data (largest clean robust $z$ 6.33, which Part 2 prints), so it raises no false alarms and catches the 15 °C burst that both z-scores miss. The 8 °C burst reaches a robust $z$ of only 7.3, barely above the clean maximum: no threshold of this rule separates it reliably, and at that size you need a better model of the normal signal (a shorter window, or the daily cycle subtracted first).

**B6.**

```python
rng = np.random.default_rng(6)
n = 100_000
big = pd.DataFrame({"value": rng.normal(300, 20, n), "unit": rng.choice(["degC", "degF", "K"], n)})
SCALE, OFFSET = {"degC": 1.0, "degF": 5 / 9, "K": 1.0}, {"degC": 0.0, "degF": -32.0, "K": -273.15}


def to_c(value: float, unit: str) -> float:
    return (value + OFFSET[unit]) * SCALE[unit]


ways = {
    "iterrows": lambda: [to_c(r["value"], r["unit"]) for _, r in big.iterrows()],
    "apply(axis=1)": lambda: big.apply(lambda r: to_c(r["value"], r["unit"]), axis=1),
    "itertuples": lambda: [to_c(r.value, r.unit) for r in big.itertuples(index=False)],
    "loop over to_numpy()": lambda: [to_c(v, u) for v, u in zip(big["value"].to_numpy(), big["unit"].to_numpy())],
    "np.select": lambda: np.select([big["unit"] == "degF", big["unit"] == "K"],
                                   [(big["value"] - 32) * 5 / 9, big["value"] - 273.15], big["value"]),
    "map scale/offset": lambda: (big["value"] + big["unit"].map(OFFSET)) * big["unit"].map(SCALE),
}
reference = np.asarray(ways["np.select"]())
times = {}
for name, f in ways.items():
    assert np.allclose(np.asarray(f(), dtype=float), reference)
    times[name] = best_of(f, repeat=1 if name in ("iterrows", "apply(axis=1)") else 3)
for name, sec in times.items():
    print(f"{name:22s} {sec * 1e3:9.2f} ms   {sec / times['map scale/offset']:8.0f}x the map version")

as_cat = big.astype({"unit": "category"})
for label, frame in (("unit as str", big), ("unit as category", as_cat)):
    mem = frame["unit"].memory_usage(deep=True) / 1e6
    t = best_of(lambda: frame.groupby("unit", observed=True)["value"].mean())
    print(f"{label:17s} {mem:6.2f} MB, groupby mean {t * 1e3:.2f} ms")
```
```text
iterrows                 2068.24 ms        258x the map version
apply(axis=1)             372.04 ms         46x the map version
itertuples                 51.44 ms          6x the map version
loop over to_numpy()       17.85 ms          2x the map version
np.select                  12.49 ms          2x the map version
map scale/offset            8.00 ms          1x the map version
unit as str         5.20 MB, groupby mean 6.32 ms
unit as category    0.10 MB, groupby mean 1.54 ms
```

Every row-wise method builds or looks up Python objects per row; the two vectorized ones do a few passes in compiled code. `iterrows` is the worst because it builds a full `Series` per row (and would silently upcast if the columns had mixed numeric dtypes); `apply(axis=1)` does the same with less interpreter overhead; `itertuples` and the plain loop over arrays avoid the `Series` and are roughly 40 and 100 times faster than `iterrows`, but still make a Python call per row. At $10^7$ rows every time grows about 100-fold (all are linear), so the *ratios* stay, while the absolute cost changes what is feasible: the `iterrows` version would take minutes, the vectorized ones one or two seconds. The memory saving of `category` grows with $n$ (codes are 1 byte per row against a pointer plus a string object), and the groupby speed-up stays, since the codes are already the integer labels that §3.8's NumPy version had to compute.

## C. Challenge

**C1.**

```python
true_utc = pd.date_range("2024-10-25 22:00", "2024-10-27 23:00", freq="15min", tz="UTC", inclusive="left")
written = true_utc.tz_convert("Europe/Berlin").strftime("%Y-%m-%d %H:%M")      # what the logger writes
naive = pd.DatetimeIndex(pd.to_datetime(written))
print(len(naive), "stamps,", naive.duplicated().sum(), "repeated:", sorted({t.strftime("%H:%M") for t in naive[naive.duplicated()]}))
try:
    naive.tz_localize("Europe/Berlin")
except Exception as err:                                                      # the exact class varies by version
    print(type(err).__name__, "|", err)
repaired = naive.tz_localize("Europe/Berlin", ambiguous="infer")
assert (repaired == true_utc).all()
ones = pd.Series(1, index=repaired)
print("per local day:", ones.resample("D").size().to_dict())
print("per UTC day:  ", ones.tz_convert("UTC").resample("D").size().to_dict())

march = pd.DatetimeIndex(["2024-03-31 01:45", "2024-03-31 02:30", "2024-03-31 03:00"])
for policy in ("shift_forward", "NaT"):
    print(policy, march.tz_localize("Europe/Berlin", nonexistent=policy).strftime("%H:%M %z").tolist())
```
```text
196 stamps, 4 repeated: ['02:00', '02:15', '02:30', '02:45']
ValueError | Cannot infer dst time from 2024-10-27 02:00:00, try using the 'ambiguous' argument
per local day: {Timestamp('2024-10-26 00:00:00+0200', tz='Europe/Berlin'): 96, Timestamp('2024-10-27 00:00:00+0200', tz='Europe/Berlin'): 100}
per UTC day:   {Timestamp('2024-10-25 00:00:00+0000', tz='UTC'): 8, Timestamp('2024-10-26 00:00:00+0000', tz='UTC'): 96, Timestamp('2024-10-27 00:00:00+0000', tz='UTC'): 92}
shift_forward ['01:45 +0100', '03:00 +0200', '03:00 +0200']
NaT ['01:45 +0100', nan, '03:00 +0200']
```

The October hour from 02:00 to 02:59 local happens twice, so a naive logger writes four stamps twice and `tz_localize` refuses to guess. `ambiguous="infer"` uses the order of the stamps: the first 02:15 must be summer time and the second winter time, because the stamps would otherwise run backwards; this works only for data still in logging order, and the repaired index equals the true instants. Local days hold 96 and 100 readings (24 and 25 hours); UTC days hold whatever the 22:00 UTC start leaves, and every full UTC day has 96. In March, local 02:00–02:59 does not exist, so an honest logger never writes it; `nonexistent=` is for timestamps you *construct* (a schedule written as "every day at 02:30"), and you choose whether they move forward to 03:00 or become `NaT`.

Design note. The logger should write UTC with an explicit offset in ISO 8601 (`2024-10-27T01:15:00+00:00`) or integer epoch seconds, never naive local time. UTC is monotonic and every hour has 60 minutes, so sorting, differencing and resampling are correct with no special cases. Local time is a display concern: convert at the edge, with `tz_convert`, for a report's axis labels. If a plant insists on local time, it must at least write the offset (`+02:00`/`+01:00`), which makes every stamp unambiguous. And a repair by `ambiguous="infer"` depends on row order, which the first `sort_values` or `drop_duplicates` destroys.

**C2.**

```python
import tracemalloc
from pathlib import Path

rng = np.random.default_rng(5)
seconds = 250_000                                       # about 2.9 days at 1 s, 8 series -> 2,000,000 rows
pairs = [(m, s) for m in ("A", "B") for s in ("temp_1", "temp_2", "temp_3", "temp_4")]
t = np.tile(np.arange(seconds), len(pairs))
pair = np.repeat(np.arange(len(pairs)), seconds)
df = pd.DataFrame({"time": pd.Timestamp("2024-03-04") + pd.to_timedelta(t, unit="s"),
                   "machine": np.array([m for m, _ in pairs])[pair],
                   "sensor": np.array([s for _, s in pairs])[pair],
                   "value": 60 + 5 * np.sin(2 * np.pi * t / 86_400) + rng.normal(0, 1, len(t))})
dups = df.iloc[rng.choice(len(df), len(df) // 100, replace=False)]
df = pd.concat([df, dups]).sample(frac=1.0, random_state=5)              # duplicates anywhere in the file
path = Path("c2_big.csv")
df.to_csv(path, index=False)
print(f"{len(df):,} rows, {path.stat().st_size / 1e6:.0f} MB on disk")
SENSOR_CODE = {f"{m}|{s}": i for i, (m, s) in enumerate(pairs)}
T0 = pd.Timestamp("2024-03-04")


def chunked_daily_means(path: Path, chunksize: int = 200_000, dedupe_across_chunks: bool = True) -> pd.DataFrame:
    """Exact per-(machine, day) sum, count and mean; duplicates found with a bitmap over the dense key (second, series)."""
    seen = np.zeros(seconds * len(pairs), dtype=bool)
    tot, cnt = None, None
    for chunk in pd.read_csv(path, chunksize=chunksize, parse_dates=["time"]):
        if not dedupe_across_chunks:
            seen[:] = False
        key = (((chunk["time"] - T0) // pd.Timedelta("1s")).to_numpy() * len(pairs)
               + (chunk["machine"] + "|" + chunk["sensor"]).map(SENSOR_CODE).to_numpy())
        first_in_chunk = ~pd.Series(key).duplicated().to_numpy()
        keep = first_in_chunk & ~seen[key]
        seen[key] = True
        part = chunk.loc[keep].groupby(["machine", chunk.loc[keep, "time"].dt.floor("D")])["value"].agg(["sum", "count"])
        tot = part if tot is None else tot.add(part, fill_value=0)
    return tot.assign(mean=tot["sum"] / tot["count"])


def full_daily_means(path: Path) -> pd.Series:
    d = pd.read_csv(path, parse_dates=["time"]).drop_duplicates()
    return d.groupby(["machine", d["time"].dt.floor("D")])["value"].mean().rename("mean")


results = {}
for name, f in (("chunked", lambda: chunked_daily_means(path)["mean"]), ("full read", lambda: full_daily_means(path))):
    tracemalloc.start()
    t0 = time.perf_counter()
    results[name] = f()
    elapsed = time.perf_counter() - t0
    peak = tracemalloc.get_traced_memory()[1] / 1e6
    tracemalloc.stop()
    print(f"{name:10s} {elapsed:5.1f} s, peak traced memory {peak:6.0f} MB")
pd.testing.assert_series_equal(results["chunked"], results["full read"], check_exact=False, rtol=1e-12)
naive = chunked_daily_means(path, dedupe_across_chunks=False)
print(results["full read"].round(4).to_dict())
print(f"within-chunk dedupe only: {int(naive['count'].sum()) - len(df) + len(dups):,} of {len(dups):,} duplicates survive "
      f"({-(-len(df) // 200_000)} chunks); max error in a daily mean {(naive['mean'] - results['full read']).abs().max():.2e}")
path.unlink()
```
```text
2,020,000 rows, 95 MB on disk
chunked      9.3 s, peak traced memory     40 MB
full read    5.2 s, peak traced memory    230 MB
{('A', Timestamp('2024-03-04 00:00:00')): 60.0008, ('A', Timestamp('2024-03-05 00:00:00')): 60.0039, ('A', Timestamp('2024-03-06 00:00:00')): 60.1915, ('B', Timestamp('2024-03-04 00:00:00')): 60.0015, ('B', Timestamp('2024-03-05 00:00:00')): 59.9963, ('B', Timestamp('2024-03-06 00:00:00')): 60.1954}
within-chunk dedupe only: 18,047 of 20,000 duplicates survive (11 chunks); max error in a daily mean 7.97e-04
```

The chunked version holds one chunk (200,000 rows) plus a bitmap of 2,000,000 booleans (2 MB) and the running sums; the full read holds the whole table, so its peak grows with the file while the chunked peak does not. Sums and counts combine across chunks exactly; means would not. Deduplication is the hard part: a duplicate may land in a different chunk from its original, so "seen" must persist across chunks. A Python `set` of 2,000,000 tuples would cost hundreds of MB; because the key (second, series) is dense, a boolean array indexed by it is exact and tiny. Deduplicating within chunks only keeps every cross-chunk copy, which is most of them (a copy lands in its original's chunk with probability about 1/11), so each daily mean is computed with the wrong weights, here off in the fourth decimal, and every count is inflated.

**C3.**

```python
def rolling_time_mean(index: pd.DatetimeIndex, x: np.ndarray, window: pd.Timedelta, closed: str = "right") -> np.ndarray:
    """Trailing time-window mean by binary search and cumulative sums, O(n log n) -> O(n) after the search.

    closed="right": window (t - w, t];  closed="both": [t - w, t].
    """
    t = index.as_unit("ns").asi8
    w = window // pd.Timedelta("1ns")
    first = np.searchsorted(t, t - w, side="right" if closed == "right" else "left")
    cs = np.concatenate([[0.0], np.cumsum(x)])
    i = np.arange(len(t))
    return (cs[i + 1] - cs[first]) / (i + 1 - first)


rng = np.random.default_rng(8)
idx = pd.Timestamp("2024-03-04") + pd.to_timedelta(np.cumsum(rng.exponential(60, 2000)).round(), unit="s")
s = pd.Series(rng.normal(0, 1, 2000), index=pd.DatetimeIndex(idx))
for closed in ("right", "both"):
    assert np.allclose(rolling_time_mean(s.index, s.to_numpy(), pd.Timedelta("1h"), closed),
                       s.rolling("1h", closed=closed).mean().to_numpy(), atol=1e-12)
hit = np.isin(s.index - pd.Timedelta("1h"), s.index).sum()
print(f"both closings agree with pandas; {hit} readings have a neighbour exactly 1 h earlier, where they differ")
default_is_right = np.allclose(rolling_time_mean(s.index, s.to_numpy(), pd.Timedelta("1h")), s.rolling("1h").mean())
print("default rolling('1h') is (t - 1h, t]:", default_is_right)

n = 1_000_000
idx_big = pd.Timestamp("2024-03-04") + pd.to_timedelta(np.cumsum(rng.exponential(1.0, n)), unit="s")
s_big = pd.Series(rng.normal(0, 1, n), index=pd.DatetimeIndex(idx_big))
t_hand = best_of(lambda: rolling_time_mean(s_big.index, s_big.to_numpy(), pd.Timedelta("1h")))
t_pd = best_of(lambda: s_big.rolling("1h").mean())
err = np.abs(rolling_time_mean(s_big.index, s_big.to_numpy(), pd.Timedelta("1h")) - s_big.rolling("1h").mean()).max()
print(f"n = 10^6: searchsorted + cumsum {t_hand * 1e3:.0f} ms, pandas {t_pd * 1e3:.0f} ms, max difference {err:.1e}")
shifted = s_big + 60.0                                         # same data around 60, like a temperature
err60 = np.abs(rolling_time_mean(shifted.index, shifted.to_numpy(), pd.Timedelta("1h")) - shifted.rolling("1h").mean()).max()
print(f"the same data plus 60: max difference {err60:.1e}")
```
```text
both closings agree with pandas; 26 readings have a neighbour exactly 1 h earlier, where they differ
default rolling('1h') is (t - 1h, t]: True
n = 10^6: searchsorted + cumsum 50 ms, pandas 49 ms, max difference 7.3e-15
the same data plus 60: max difference 1.2e-10
```

The window for row $i$ starts at the first stamp inside it; `searchsorted(t, t - w, side="right")` finds the first $t_j > t_i - w$, which is the half-open window $(t_i - w, t_i]$, and `side="left"` finds the first $t_j \ge t_i - w$, the closed one. Agreement with pandas in both cases, on data where some readings sit exactly one hour after another, shows that the default for a time-based window is `closed="right"`, i.e. $(t-1\,\mathrm h, t]$. The cumulative sum makes each window sum one subtraction, so the cost is the $O(n\log n)$ binary search plus $O(n)$ arithmetic. Its weakness is rounding: each window sum is the difference of two running sums, which reach $6\times10^7$ when $10^6$ readings sit near 60, so the rounding error of the running sum, not of the window, sets the accuracy; shifting the same data by 60 multiplies the error by about $10^4$ (last line). pandas' rolling mean adds and subtracts values as the window moves, with compensated (Kahan) summation, so it does not degrade this way.
