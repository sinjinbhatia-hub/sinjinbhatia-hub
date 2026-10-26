# PY-9: pandas — Tables, Alignment, Groupby, Merge, Time Series, and Tidy Data

**Time:** ~8 hours · **Prerequisites:** PY-6, PY-7 (F4 helpful) · **Unlocks:** ST-1 onward (every statistics project handles real tables), ML-1, ST-7 (time series)

> **Conventions in this module.** Every output below came from a real run (Python 3.13, pandas 3.0.5, NumPy 2.5.3; the first block prints them). Blocks run top to bottom in one session, so names persist; the `text` block after a block is its real output, and timings are from one run (only their ratios transfer). `import pandas as pd`. Select with `df["col"]`, `df.loc[rows, cols]` (labels) and `df.iloc[rows, cols]` (positions) only; **never assign through chained indexing** `df[a][b] = ...` (§3.3). Multi-step transforms are one method chain inside parentheses, one step per line. Time is `pd.Timestamp`/`DatetimeIndex`; frequency strings are pandas 3's lowercase `"10min"`, `"h"`, `"D"` (`"H"` is an error now). Random numbers come from `rng = np.random.default_rng(seed)`. **pandas 3 differs from most tutorials online in three ways that matter here:** Copy-on-Write is always on, so the `SettingWithCopyWarning` of pandas 1–2 no longer exists (§3.3); text columns get the dtype `str` instead of `object`; parsed timestamps are `datetime64[us]`, not `[ns]`. A *tidy* table has one variable per column and one observation per row (§3.14).

## 0. Where this sits

PY-2 gave you the dict and `defaultdict(list)` grouping (its Example 4.4), PY-6 the array and the view/copy distinction (§3.4), PY-7 the figure. pandas glues an index of labels onto PY-6's arrays and gives you the table operations that statistics runs on. After this module ST-1 through ST-9 read, clean and summarize real tables, ML-1 builds its train/test splits from them, ST-7 starts from a `DatetimeIndex`, and ST-3's regressions start from §3.14's clean long table.

## 1. The problem that forced this idea

Real measurements arrive as tables with labels: a timestamp, a machine, a sensor, a reading, a unit. Here are ten rows of a plant's logger export, exactly as they come:

```python
import io
import sys
import time
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd

pd.set_option("display.width", 100)
print(sys.version.split()[0], pd.__version__, np.__version__)
raw_csv = """time,machine,sensor,value,unit
2024-03-04 08:00,A, Temp_1,71.2,degC
2024-03-04 08:00,A,vib_1,0.31,mm/s
2024-03-04 08:10,A,temp_1,71.9,degC
2024-03-04 08:10,A, Temp_1,71.9,degC
2024-03-04 08:20,A,temp_1,999,degC
2024-03-04 08:00,B,TEMP-1,68.0,degC
2024-03-04 08:10,B,temp_1,,degC
2024-03-04 08:00,C,temp_1,160.3,degF
2024-03-04 08:10,C,temp_1,161.6,degF
2024-03-04 08:10,C,vib_1,0.29,mm/s
"""
raw = pd.read_csv(io.StringIO(raw_csv))
print(raw)
print(raw["sensor"].unique().tolist(), raw.dtypes.astype(str).to_dict())
```
```text
3.13.16 3.0.5 2.5.3
               time machine   sensor   value  unit
0  2024-03-04 08:00       A   Temp_1   71.20  degC
1  2024-03-04 08:00       A    vib_1    0.31  mm/s
2  2024-03-04 08:10       A   temp_1   71.90  degC
3  2024-03-04 08:10       A   Temp_1   71.90  degC
4  2024-03-04 08:20       A   temp_1  999.00  degC
5  2024-03-04 08:00       B   TEMP-1   68.00  degC
6  2024-03-04 08:10       B   temp_1     NaN  degC
7  2024-03-04 08:00       C   temp_1  160.30  degF
8  2024-03-04 08:10       C   temp_1  161.60  degF
9  2024-03-04 08:10       C    vib_1    0.29  mm/s
[' Temp_1', 'vib_1', 'temp_1', 'TEMP-1'] {'time': 'str', 'machine': 'str', 'sensor': 'str', 'value': 'float64', 'unit': 'str'}
```

Count the problems. One sensor is spelled three ways (`" Temp_1"`, `"temp_1"`, `"TEMP-1"`). Row 3 is a retransmitted copy of row 2, spelled like row 0. `999` is the logger's "no reading" sentinel. Row 6 has an empty cell, which `read_csv` turned into `NaN`. Machine C logs Fahrenheit. And `time` came in as text, not as time. Multiply this by 50,000 rows from 12 sensors on 3 machines over two weeks, with gaps, and ask three questions: what is each machine's daily mean temperature; how correlated are a temperature sensor sampled every 10 minutes and a pressure sensor sampled every 15; and give me a clean table to regress one on the other (ST-3).

With NumPy you would keep five parallel arrays in step by hand through every filter, sort and deduplication, and match two sampling rates with `searchsorted`. A dict of lists has labels but no vectorized arithmetic, alignment or grouping. pandas has all three: columns that are arrays, rows that carry labels, and operations that line labels up before they compute. The price is that you must learn to think in aligned, labelled columns, and that pandas has more API than any library in this course. This module is the subset you will use every week.

## 2. Intuition

**A Series is an array with labels glued on.** A **DataFrame** is a dict of Series that share one index:

```text
Series s (name "temp")        DataFrame df: one shared index, a dict of columns
 index | value                index | time              machine  value
 ------+------                ------+-------------------------------------
  "a"  | 71.2                     0 | 2024-03-04 08:00  "A"      71.2
  "b"  | 71.9                     1 | 2024-03-04 08:10  "A"      71.9
  "c"  | 68.0                     2 | 2024-03-04 08:00  "B"      68.0
 values: one float64 array (PY-6)  each column: one array with one dtype
 index:  labels, used to find and to align rows
```

Two ways to name a row: by **label** (`"b"`, or `2024-03-04 08:10`) and by **position** (the second row). pandas keeps them strictly apart: `.loc` speaks labels, `.iloc` speaks positions. The default index `0, 1, 2, ...` makes them look the same, until you sort or filter and they stop being the same.

**Every operation aligns on labels first.** Add two Series and pandas matches labels, not positions:

```text
  s1:   a 1     b 2     c 3
  s2:           b 10    c 20    d 30
  s1 + s2:  a ?     b 12    c 23    d ?
```

▶ **Guess:** what goes where the `?` are? A label present on one side only has nothing to add to. pandas's answer is "missing" (`NaN`), and the index of the result is the union of the two indexes. NumPy would have added position by position, or refused because the shapes differ.

**`groupby` is split, apply, combine.** "Mean temperature per machine" means: split the rows by the value of `machine`, apply `mean` to each piece, and combine the answers into a new Series indexed by machine:

```text
 machine value          split           apply mean      combine
  A      71.2   ->  A: [71.2, 71.9]  ->  A: 71.55   ->  machine
  A      71.9       B: [68.0]           B: 68.0         A  71.55
  B      68.0                                           B  68.00
```

That is PY-2's `defaultdict(list)` with a nicer face, and §3.8 builds it three ways.

**Tidy tables.** The logger export above has one *observation* (one reading) per row and one *variable* (time, machine, sensor, value, unit) per column. That shape is called **tidy**, and it is the shape every pandas operation wants: `groupby` splits on a column, `merge` matches on a column, plotting picks columns. A table with one column per sensor (`temp_1`, `vib_1`, ...) stores a *variable* (which sensor) in the column headers; that is fine to look at and painful to compute with. Most "pandas is confusing" is "my table is not tidy" (§3.14).

## 3. Formal development

### 3.1 Series and DataFrame

A `Series` is a 1-D array of values plus an `Index` of the same length; a `DataFrame` is an ordered set of columns, each a Series, sharing one row index, plus a column index. You build them from a list or array (default index `0..n-1`), a dict (keys become labels or column names), a list of records, or a 2-D array:

```python
s = pd.Series([71.2, 71.9, 68.0], index=["a", "b", "c"], name="temp")
df = pd.DataFrame({"machine": ["A", "A", "B"], "temp": [71.2, 71.9, 68.0], "ok": [True, True, False]})
df_rec = pd.DataFrame.from_records([("A", 71.2, True), ("A", 71.9, True), ("B", 68.0, False)],
                                   columns=["machine", "temp", "ok"])
df_arr = pd.DataFrame(np.arange(6.0).reshape(3, 2), columns=["x", "y"])
assert df.equals(df_rec)
print(s.index.tolist(), s.to_numpy(), s.dtype, s.name)
print(df.shape, df.columns.tolist(), df.index)
print(df.dtypes.astype(str).to_dict())
print(df_arr.describe().loc[["count", "mean", "std", "max"]])
df.info()
```
```text
['a', 'b', 'c'] [71.2 71.9 68. ] float64 temp
(3, 3) ['machine', 'temp', 'ok'] RangeIndex(start=0, stop=3, step=1)
{'machine': 'str', 'temp': 'float64', 'ok': 'bool'}
         x    y
count  3.0  3.0
mean   2.0  3.0
std    2.0  2.0
max    4.0  5.0
<class 'pandas.DataFrame'>
RangeIndex: 3 entries, 0 to 2
Data columns (total 3 columns):
 #   Column   Non-Null Count  Dtype
---  ------   --------------  -----
 0   machine  3 non-null      str
 1   temp     3 non-null      float64
 2   ok       3 non-null      bool
dtypes: bool(1), float64(1), str(1)
memory usage: 183.0 bytes
```

`describe` gives count, mean, standard deviation (with $n-1$), minimum, quartiles and maximum per numeric column; `info` gives the non-null counts and the dtypes, and it is the first thing to run on any table you did not make (on `raw` it would show the missing value in `value` and the text `time`). `head(n)` and `tail(n)` show the ends. Each column is one array with one dtype: `int64`, `float64`, `bool`, `str` (pandas 3's text dtype), `datetime64[us]`, `category` (§3.4), or the nullable `Int64`/`boolean` (§3.6).

### 3.2 Selecting rows and columns

`df["col"]` is a column (a Series), `df[["a", "b"]]` is a DataFrame of those columns, and `df[mask]` with a boolean Series keeps the rows where the mask is `True`. Everything else goes through two indexers: **`.loc[rows, cols]` takes labels, label slices and boolean masks; `.iloc[rows, cols]` takes integer positions and position slices.** `.at[label, col]` and `.iat[i, j]` are the fast scalar versions.

▶ **Try it first.** Predict the five lines, paying attention to where each slice stops.

```python
t = pd.DataFrame({"temp": [71.2, 71.9, 68.0, 69.4, 70.1], "vib": [0.31, 0.35, 0.28, 0.30, 0.33]},
                 index=["a", "b", "c", "d", "e"])
print(t.loc["b":"d", "temp"].tolist())
print(t.iloc[1:3, 0].tolist())
print(t.loc[t["temp"] > 70, ["vib"]])
print(t.at["c", "vib"], t.iat[2, 1])
print(type(t["temp"]).__name__, type(t[["temp"]]).__name__)
```
```text
[71.9, 68.0, 69.4]
[71.9, 68.0]
    vib
a  0.31
b  0.35
e  0.33
0.28 0.28
Series DataFrame
```

■ **What just happened.** **A label slice includes both ends**, unlike every other slice in Python: `"b":"d"` is three rows, `1:3` is two. It has to be so: labels need not be integers or consecutive, so "stop just before `"d"`" would require knowing which label comes before `"d"`. A label slice runs from the *position* of the first label to the position of the second, so on an unsorted index it follows the current row order (Example 4.1). The `t.loc[mask, ["vib"]]` form, with a list of columns, returns a DataFrame; a single column name returns a Series.

### 3.3 Views, copies, and chained assignment

PY-6 §3.4 showed that a NumPy slice is a view and a mask selection is a copy. Through pandas 2, an indexing result was sometimes a view and sometimes a copy, depending on dtypes and the kind of indexer, and assigning into it might or might not change the parent. pandas warned you with `SettingWithCopyWarning` because it could not tell which. pandas 3 removed the ambiguity with **Copy-on-Write**: every result of indexing *behaves as a copy*. Internally it still shares memory, and the copy is made only when one side is written to. So:

1. `sub = df[df["temp"] > 70]; sub["x"] = 0` changes `sub` and never `df`. Silent and correct.
2. `df["temp"][mask] = 0` (chained assignment: a `__getitem__`, then a `__setitem__` on its result) **never** changes `df`; pandas 3 emits `ChainedAssignmentError` to tell you.
3. `df.loc[mask, "temp"] = 0` is one `__setitem__` on `df` itself, and is the only way to modify part of a frame.

Example 4.2 triggers the error and fixes it. Two related rules: `inplace=True` saves nothing under Copy-on-Write and breaks method chains, so do not use it; and `.copy()` is now only needed when you want to cut the memory link to a big parent explicitly.

### 3.4 Creating and modifying columns

Column arithmetic is vectorized and aligned: `df["F"] = df["temp"] * 9 / 5 + 32`. In a chain, `assign(name=func)` adds or replaces columns and returns a new frame, calling `func` on the frame *at that point of the chain*. Conditional columns come from `Series.where(cond, other)` (keep where `cond`, else `other`), its mirror `mask(cond, other)`, `np.where(cond, a, b)` for two choices and `np.select(conds, choices, default)` for several:

```python
t2 = t.assign(temp_F=lambda d: d["temp"] * 9 / 5 + 32,
              level=lambda d: np.select([d["temp"] > 71, d["temp"] < 69], ["high", "low"], "normal"))
print(t2)
print(t["vib"].where(t["vib"] < 0.33).tolist(), t["vib"].mask(t["vib"] < 0.33, 0.0).tolist())
print(t["temp"].astype("int64").tolist(), t2["level"].astype("category").cat.codes.tolist())
```
```text
   temp   vib  temp_F   level
a  71.2  0.31  160.16    high
b  71.9  0.35  161.42    high
c  68.0  0.28  154.40     low
d  69.4  0.30  156.92  normal
e  70.1  0.33  158.18  normal
[0.31, nan, 0.28, 0.3, nan] [0.0, 0.35, 0.0, 0.0, 0.33]
[71, 71, 68, 69, 70] [0, 0, 1, 2, 2]
```

`astype` converts dtypes (float to int truncates; it raises on `NaN`). **Categorical dtype** stores each distinct string once, plus a small integer code per row. For a column with few distinct values it saves most of the memory and speeds up grouping:

```python
labels = pd.Series(np.random.default_rng(0).choice(["temp_1", "vib_1", "pres_1"], 100_000))
as_cat = labels.astype("category")
print(labels.dtype, labels.memory_usage(deep=True), "bytes;", as_cat.dtype, as_cat.memory_usage(deep=True), "bytes")
print(as_cat.cat.categories.tolist(), as_cat.cat.codes.dtype)
```
```text
str 5466642 bytes; category 100296 bytes
['pres_1', 'temp_1', 'vib_1'] int8
```

A 55-fold saving. Do not make a column categorical if most values are distinct (IDs, free text): you pay for the categories *and* the codes. And a categorical with a value outside its categories cannot be assigned until you add the category.

### 3.5 Alignment (Build 2)

**Semantics.** For a binary operation on two Series, pandas forms the union of the two indexes (sorted when it can be), places each operand's values at their labels, fills absent labels with `NaN`, and operates elementwise. For a DataFrame and a Series, the Series' index aligns with the DataFrame's *columns*, which is why row-wise operations need `df.sub(series, axis=0)`.

▶ **Try it first.** Write `add_aligned(x, y)` for two dicts: the union of keys, the sum where both have the key, `nan` otherwise. Then predict `s1 + s2`, `s1.add(s2, fill_value=0)` and the plain-array sum.

```python
s1 = pd.Series([1.0, 2.0, 3.0], index=["a", "b", "c"])
s2 = pd.Series([10.0, 20.0, 30.0], index=["b", "c", "d"])


def add_aligned(x: dict, y: dict) -> dict:
    """x + y on the union of labels; nan where a label is on one side only."""
    return {k: x[k] + y[k] if k in x and k in y else float("nan") for k in sorted(x.keys() | y.keys())}


by_hand = add_aligned(s1.to_dict(), s2.to_dict())
print(by_hand)
print((s1 + s2).to_dict())
assert pd.Series(by_hand).equals(s1 + s2)
print(s1.add(s2, fill_value=0).to_dict(), (s1 + s2).fillna(0).to_dict())
left, right = s1.align(s2)
print(left.tolist(), right.tolist())
print(s1.reindex(["a", "b", "z"]).tolist(), (s1.to_numpy() + s2.to_numpy()).tolist())
```
```text
{'a': nan, 'b': 12.0, 'c': 23.0, 'd': nan}
{'a': nan, 'b': 12.0, 'c': 23.0, 'd': nan}
{'a': 1.0, 'b': 12.0, 'c': 23.0, 'd': 30.0} {'a': 0.0, 'b': 12.0, 'c': 23.0, 'd': 0.0}
[1.0, 2.0, 3.0, nan] [nan, 10.0, 20.0, 30.0]
[1.0, 2.0, nan] [11.0, 22.0, 33.0]
```

■ **What just happened.** Your dict comprehension is the semantics. `fill_value=0` treats a missing side as 0 *before* adding (both missing stays `NaN`); `fillna(0)` *after* the sum replaces the `NaN`s by 0, a different answer. `align` returns both operands reindexed to the union, which is what `+` does internally. `reindex` conforms to a given index, inserting `NaN` for new labels. The plain arrays ignored the labels and added position by position, giving a number for `a + b` that means nothing. **Alignment is the reason to use pandas, and the reason a result can be full of `NaN` without an error**: after any binary operation on data from two sources, count the `NaN`s.

### 3.6 Missing data

pandas marks a missing float as `NaN` (IEEE not-a-number). Because `NaN != NaN`, you can never find missing values with `==`; use `isna()`/`notna()`. A missing value forces an integer column to float (third line below); the nullable dtypes `Int64`, `Float64`, `boolean` keep the type and use the separate marker `pd.NA`. Reductions skip missing values by default (`skipna=True`).

```python
m = pd.Series([1.0, np.nan, 3.0, np.nan, 5.0])
print(np.nan == np.nan, (m == np.nan).sum(), m.isna().sum(), m.mean())
print(m.fillna(0).tolist(), m.ffill().tolist(), m.interpolate().tolist())
mi = pd.Series([1, None, 3], dtype="Int64")
print(mi.tolist(), mi.sum(), pd.Series([1, None, 3]).dtype)
d = pd.DataFrame({"x": [1.0, np.nan, np.nan], "y": [1.0, 2.0, np.nan]})
print(len(d.dropna()), len(d.dropna(how="all")), len(d.dropna(subset=["y"])))
```
```text
False 0 2 3.0
[1.0, 0.0, 3.0, 0.0, 5.0] [1.0, 1.0, 3.0, 3.0, 5.0] [1.0, 2.0, 3.0, 4.0, 5.0]
[1, <NA>, 3] 4 float64
1 2 2
```

`dropna()` drops a row with any missing cell, `how="all"` only rows missing everywhere, `subset=` looks at the named columns only. `fillna(value)` fills a constant, `ffill()`/`bfill()` carry the last/next valid value, `interpolate()` draws straight lines (by position by default; by time with `method="time"`, Example 4.6). **Every fill is a modelling choice.** Missing-data mechanisms: values are *missing completely at random* (MCAR) if missingness is independent of everything, *missing at random* (MAR) if it depends only on observed columns (the logger drops readings when the machine is off, and "off" is recorded), and *missing not at random* (MNAR) if it depends on the missing value itself (a sensor saturates and reports nothing exactly when the reading is extreme). Dropping rows is harmless only under MCAR; under MNAR every mean you compute is biased, and no pandas method can tell you which case you are in. That is a selection argument, of the kind ST-9 makes with DAGs.

### 3.7 Sorting, strings, and the cost of `apply`

`sort_values(["a", "b"], ascending=[False, True])` sorts by several keys (stable, so ties keep their order); `sort_index()` sorts by labels; `nlargest(k, "col")` is a partial sort. String methods live behind the `.str` accessor and run on the whole column: `strip`, `lower`, `replace`, `contains`, `startswith`, `split(..., expand=True)` (one column per piece), and `extract(regex)` with named groups becoming columns (PY-3 Example 4.5 wrote the regex by hand):

```python
names = pd.Series([" Temp_1", "TEMP-1", "vib_2 ", "Pres_10"])
norm = names.str.strip().str.lower().str.replace("-", "_")
print(norm.tolist(), norm.str.contains("temp").tolist())
print(norm.str.extract(r"(?P<kind>[a-z]+)_(?P<num>\d+)").assign(num=lambda d: d["num"].astype(int)))
print(norm.str.split("_", expand=True).to_numpy().tolist())
print(t.sort_values(["vib", "temp"], ascending=[False, True]).index.tolist(), t.nlargest(2, "temp").index.tolist())
```
```text
['temp_1', 'temp_1', 'vib_2', 'pres_10'] [True, True, False, False]
   kind  num
0  temp    1
1  temp    1
2   vib    2
3  pres   10
[['temp', '1'], ['temp', '1'], ['vib', '2'], ['pres', '10']]
['b', 'e', 'a', 'd', 'c'] ['b', 'a']
```

**Three ways to transform a column, and their prices.** `Series.map(dict_or_func)` maps each value; `df.apply(func, axis=1)` calls a Python function on every row, built as a Series; `iterrows()` yields `(label, row Series)` pairs to a `for` loop. The vectorized expression does the same work in compiled loops. Here is "convert Fahrenheit rows to Celsius" on 50,000 rows:

▶ **Try it first.** Guess the ratio between the slowest and the fastest.

```python
rng = np.random.default_rng(1)
n = 50_000
big = pd.DataFrame({"value": rng.normal(100, 10, n), "unit": rng.choice(["degC", "degF"], n)})


def to_celsius_row(row: pd.Series) -> float:
    return (row["value"] - 32) * 5 / 9 if row["unit"] == "degF" else row["value"]


def best_of(f, repeat: int = 3) -> tuple[float, object]:
    """Smallest wall time of `repeat` calls, and the last result."""
    times = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        out = f()
        times.append(time.perf_counter() - t0)
    return min(times), out


ways = {
    "iterrows loop": lambda: pd.Series([to_celsius_row(r) for _, r in big.iterrows()]),
    "apply(axis=1)": lambda: big.apply(to_celsius_row, axis=1),
    "map + arithmetic": lambda: (big["value"] + big["unit"].map({"degC": 0.0, "degF": -32.0}))
                                * big["unit"].map({"degC": 1.0, "degF": 5 / 9}),
    "where (vectorized)": lambda: big["value"].where(big["unit"] == "degC", (big["value"] - 32) * 5 / 9),
}
results = {name: best_of(f, repeat=1 if "iter" in name else 3) for name, f in ways.items()}
for name, (sec, out) in results.items():
    assert np.allclose(out.to_numpy(), results["where (vectorized)"][1].to_numpy())
    print(f"{name:20s} {sec * 1e3:8.1f} ms")
```
```text
iterrows loop           940.9 ms
apply(axis=1)           172.8 ms
map + arithmetic          6.2 ms
where (vectorized)        4.8 ms
```

`apply(axis=1)` builds a Series per row; `iterrows` does the same and also runs your loop in the interpreter. Both are the slow path that *feels* natural because it reads like the row-by-row description of the task. The rule: if you are writing `lambda row:`, look for the column expression (`where`, `np.select`, `map`, `.str`, `.dt`, a merge). `apply` on a *group* (§3.8) or on a column with a function that is itself vectorized is fine.

### 3.8 Groupby (Build 1)

**Semantics.** `df.groupby(keys)` splits rows by the distinct values of one or more key columns (sorted, missing keys dropped unless `dropna=False`). Then:

- **aggregate** (`.mean()`, `.agg(...)`): one row per group;
- **transform** (`.transform(...)`): one value per *original row*, the group result broadcast back, aligned on the original index;
- **filter** (`.filter(func)`): keep or drop whole groups.

▶ **Try it first.** Write the per-machine mean three ways before reading: with a `defaultdict(list)` (PY-2 Example 4.4), with NumPy (labels to integer codes by `np.unique(..., return_inverse=True)`, then a scatter-add), then with pandas. Why does NumPy need `np.add.at` rather than `sums[codes] += values` (PY-6 §3.4)?

```python
rng = np.random.default_rng(9)
n = 200_000
log = pd.DataFrame({"machine": rng.choice(["A", "B", "C"], n, p=[0.5, 0.3, 0.2]),
                    "sensor": rng.choice(["temp_1", "temp_2", "temp_3", "temp_4"], n)})
log["value"] = (log["machine"].map({"A": 60.0, "B": 65.0, "C": 58.0})
                + log["sensor"].str[-1].astype(int) - 1 + rng.normal(0, 2.0, n))


def mean_by_key_dict(keys, values) -> dict:
    """Split into a dict of lists, apply mean, combine into a dict."""
    groups = defaultdict(list)
    for k, v in zip(keys, values):
        groups[k].append(v)
    return {k: sum(vs) / len(vs) for k, vs in sorted(groups.items())}


def mean_by_key_numpy(keys: np.ndarray, values: np.ndarray) -> dict:
    """Labels -> codes 0..G-1 with np.unique; unbuffered scatter-add with np.add.at; counts by bincount."""
    uniq, codes = np.unique(keys, return_inverse=True)
    sums = np.zeros(len(uniq))
    np.add.at(sums, codes, values)
    return dict(zip(uniq.tolist(), (sums / np.bincount(codes)).tolist()))


keys, vals = log["machine"].to_numpy(), log["value"].to_numpy()
t_dict, r_dict = best_of(lambda: mean_by_key_dict(keys.tolist(), vals.tolist()))
t_np, r_np = best_of(lambda: mean_by_key_numpy(keys, vals))
t_pd, r_pd = best_of(lambda: log.groupby("machine")["value"].mean())
codes = pd.Categorical(keys).codes                       # labels -> int8 codes, once
t_codes, r_codes = best_of(lambda: np.bincount(codes, weights=vals) / np.bincount(codes))
assert r_dict.keys() == r_np.keys() == r_pd.to_dict().keys()
for r in (list(r_dict.values()), list(r_np.values()), r_codes):
    assert np.allclose(r, r_pd.to_numpy())
print({k: round(v, 4) for k, v in r_pd.items()})
print(f"dict of lists {t_dict * 1e3:.1f} ms | np.unique + np.add.at {t_np * 1e3:.1f} ms | groupby {t_pd * 1e3:.1f} ms"
      f" | bincount on ready-made codes {t_codes * 1e3:.2f} ms")
bad = np.zeros(2)
bad[[0, 0, 1]] += 1
print("sums[codes] += 1 with codes [0, 0, 1]:", bad)
```
```text
{'A': 61.5004, 'B': 66.4979, 'C': 59.5014}
dict of lists 21.6 ms | np.unique + np.add.at 56.6 ms | groupby 10.3 ms | bincount on ready-made codes 0.95 ms
sums[codes] += 1 with codes [0, 0, 1]: [1. 1.]
```

■ **What just happened.** All three are split–apply–combine; they differ in where the loop runs. The surprise is that NumPy came last: `np.unique` on an array of 200,000 Python strings is a sort of objects and takes nearly all of its time, while the scatter-add on integer codes (last timing) costs under a millisecond. Turning labels into integers is the expensive step, which is why `groupby` hashes them instead of sorting and why a `category` column (§3.4), which stores the codes, groups faster. `np.add.at` accumulates repeated indices; `sums[codes] += v` does not, because the fancy-index form reads all the old values first, so a repeated code adds once (last line). The means agree to rounding; only the order of floating-point summation differs.

**Aggregations.** `.agg` takes a function name, a list of names (columns become a MultiIndex), or a dict `{column: names}`; **named aggregation** `agg(new_name=(column, func))` gives flat, named output columns (Example 4.3). `size()` counts rows per group; `count()` counts non-missing values per column. Grouping by two keys gives a **MultiIndex** (rows labelled by pairs); `reset_index()` turns the levels back into ordinary columns, and `unstack()` moves the inner level to the columns.

```python
g = log.groupby(["machine", "sensor"])["value"]
per_pair = g.agg(["mean", "std", "size"]).round(3)
print(per_pair.head(3))
print(per_pair.index.nlevels, per_pair.index[0], "->", per_pair.reset_index().columns.tolist())
z = (log["value"] - g.transform("mean")) / g.transform("std")
zs = z.groupby([log["machine"], log["sensor"]]).agg(["mean", "std"])
print(len(z) == len(log), float(zs["mean"].abs().max()) < 1e-12, np.allclose(zs["std"], 1.0))
keep = log.groupby("sensor").filter(lambda grp: grp["value"].mean() > 62.0)
print(sorted(keep["sensor"].unique()), len(keep))
tiny = pd.DataFrame({"g": ["a", "a", None, "b"], "v": [1.0, np.nan, 3.0, 4.0]})
print(tiny.groupby("g")["v"].size().to_dict(), tiny.groupby("g")["v"].count().to_dict(),
      tiny.groupby("g", dropna=False)["v"].size().to_dict())
```
```text
                  mean    std   size
machine sensor
A       temp_1  60.005  2.000  25081
        temp_2  61.018  2.004  25104
        temp_3  62.006  1.993  25148
2 ('A', 'temp_1') -> ['machine', 'sensor', 'mean', 'std', 'size']
True True True
['temp_2', 'temp_3', 'temp_4'] 149916
{'a': 2, 'b': 1} {'a': 1, 'b': 1} {'a': 2, 'b': 1, nan: 1}
```

`transform` is how you compute "relative to my group" quantities without a merge: the group-wise z-score above has mean 0 and standard deviation 1 within every group, and one value per original row, aligned on the original index. Iterating `for key, grp in df.groupby(...)` works and is occasionally the clearest code (one plot per group), but it is a Python loop over groups.

### 3.9 Reshaping (Build 4)

**Long and wide.** A long table has one row per (time, sensor) reading; a wide table has one row per time and one column per sensor. `melt` goes wide to long (the tidy move); `pivot(index, columns, values)` goes long to wide and **raises if an (index, column) pair occurs twice**; `pivot_table` goes long to wide *and aggregates* duplicates with `aggfunc` (default `"mean"`). `stack` moves columns into the innermost index level and `unstack` the reverse; `crosstab` counts combinations.

▶ **Try it first.** For a wide table of 3 times × 3 sensors, predict the shape of the melted table, and write the pivot that undoes it.

```python
wide = pd.DataFrame({"time": pd.date_range("2024-03-04 08:00", periods=3, freq="10min"),
                     "temp_1": [71.2, 71.9, 72.3], "vib_1": [0.31, 0.35, 0.30], "pres_1": [2.1, 2.2, 2.2]})
long = wide.melt(id_vars="time", var_name="sensor", value_name="value")
print(long.shape)
print(long.iloc[[0, 3, 8]])
back = (long.pivot(index="time", columns="sensor", values="value")
        .loc[:, ["temp_1", "vib_1", "pres_1"]]      # pivot sorts columns; restore the order
        .reset_index()
        .rename_axis(columns=None))
assert back.equals(wide)
via_unstack = long.set_index(["time", "sensor"])["value"].unstack()
assert via_unstack.equals(long.pivot(index="time", columns="sensor", values="value"))
dup = pd.concat([long, long.iloc[[0]].assign(value=71.4)])
try:
    dup.pivot(index="time", columns="sensor", values="value")
except ValueError as err:
    print("pivot:", err)
print(dup.pivot_table(index="time", columns="sensor", values="value", aggfunc="mean").iloc[0].round(2).to_dict())
print(pd.crosstab(log["machine"], log["sensor"]).iloc[:, :2])
```
```text
(9, 3)
                 time  sensor  value
0 2024-03-04 08:00:00  temp_1  71.20
3 2024-03-04 08:00:00   vib_1   0.31
8 2024-03-04 08:20:00  pres_1   2.20
pivot: Index contains duplicate entries, cannot reshape
{'pres_1': 2.1, 'temp_1': 71.3, 'vib_1': 0.31}
sensor   temp_1  temp_2
machine
A         25081   25104
B         15030   14957
C          9973    9939
```

■ **What just happened.** Melting 3 rows × 3 sensor columns gave 9 rows of (time, sensor, value): the column headers became values of a `sensor` column. `pivot` undid it exactly once the column order and the leftover column-axis name were restored; `set_index` + `unstack` is the same operation. With a duplicate (time, sensor) pair, `pivot` refused and `pivot_table` averaged 71.2 and 71.4.

### 3.10 Combining tables (Build 3)

`pd.concat(list, axis=0)` stacks frames vertically, aligning columns by name (`keys=[...]` adds an outer index level saying which frame each row came from); `axis=1` places them side by side, aligning rows by index. **`merge`** is the relational join: rows of `left` and `right` are paired wherever their key columns are equal.

▶ **Try it first.** Write an inner join with two nested loops. Then predict the number of rows when the key `"A"` appears 2 times on the left and 3 times on the right.

```python
def inner_join(left: list[dict], right: list[dict], on: str) -> list[dict]:
    """Every pair (l, r) with l[on] == r[on], as one merged record."""
    return [{**l, **r} for l in left for r in right if l[on] == r[on]]


machines = pd.DataFrame({"machine": ["A", "B", "C"], "line": [1, 1, 2]})
vendors = pd.DataFrame({"machine": ["A", "B", "D"], "vendor": ["Acme", "Bolt", "Dyna"]})
hand = pd.DataFrame(inner_join(machines.to_dict("records"), vendors.to_dict("records"), "machine"))
assert hand.equals(machines.merge(vendors, on="machine", how="inner"))
for how in ("inner", "left", "right", "outer"):
    out = machines.merge(vendors, on="machine", how=how)
    print(f"{how:5s} {len(out)} rows, machines {out['machine'].tolist()}, vendor {out['vendor'].tolist()}")
shifts = pd.DataFrame({"machine": ["A", "A"], "shift": ["day", "night"]})
alarms = pd.DataFrame({"machine": ["A", "A", "A"], "alarm": ["E1", "E7", "E9"]})
print(len(shifts.merge(alarms, on="machine")), len(inner_join(shifts.to_dict("records"), alarms.to_dict("records"), "machine")))
try:
    shifts.merge(alarms, on="machine", validate="one_to_one")
except pd.errors.MergeError as err:
    print("MergeError:", str(err).splitlines()[0])
```
```text
inner 2 rows, machines ['A', 'B'], vendor ['Acme', 'Bolt']
left  3 rows, machines ['A', 'B', 'C'], vendor ['Acme', 'Bolt', nan]
right 3 rows, machines ['A', 'B', 'D'], vendor ['Acme', 'Bolt', 'Dyna']
outer 4 rows, machines ['A', 'B', 'C', 'D'], vendor ['Acme', 'Bolt', nan, 'Dyna']
6 6
MergeError: Merge keys are not unique in either left or right dataset; not a one-to-one merge.
```

```text
 machines        vendors           inner          left               outer
 A  1            A  Acme           A 1 Acme       A 1 Acme           A 1   Acme
 B  1            B  Bolt           B 1 Bolt       B 1 Bolt           B 1   Bolt
 C  2            D  Dyna                          C 2 NaN            C 2   NaN
                                                                     D NaN Dyna
 shifts (A twice) x alarms (A three times) on "machine"  ->  2 x 3 = 6 rows
```

■ **What just happened.** Your nested loop *is* the definition of a join: every matching pair is a row. With duplicate keys on both sides that is a Cartesian product within the key, the **many-to-many explosion**, and pandas does it silently. `left` keeps every left row (no match: `NaN` on the right's columns), `right` the mirror, `outer` both. **`validate=` turns your assumption into a check**: `"one_to_one"`, `"one_to_many"` (left keys unique), `"many_to_one"` (right keys unique, the usual case for attaching metadata to readings). `indicator=True` adds a `_merge` column saying `both`, `left_only` or `right_only` (Example 4.5). `df.join(other)` is `merge` on the index. Use `on=[...]` for several key columns and `left_on=`/`right_on=` when names differ.

### 3.11 Time series (Builds 5 and 6)

`pd.to_datetime` parses text (it infers one format from the first element; mixed formats need `format="mixed"` or `format="ISO8601"`, and `errors="coerce"` turns failures into `NaT`, the missing time). A `DatetimeIndex` makes a Series a time series: `loc` accepts **partial strings** (`ts.loc["2024-03"]` is all of March; `ts.loc["2024-03-04 08:00":"2024-03-04 09:00"]` includes both ends, §3.2), and `.dt` gives fields (`.dt.hour`, `.dt.floor("h")`, `.dt.date`). `pd.date_range(start, periods=, freq=)` makes regular indexes.

**`resample(freq)`** is a `groupby` on time bins. For `"min"`, `"h"` and `"D"` bins are closed on the left and labelled by their left edge: the 08:00 bin of `"h"` holds $08{:}00 \le t < 09{:}00$. Unlike `groupby`, `resample` emits **every** bin in the range, empty ones included (as `NaN` for `mean`, 0 for `size`).

▶ **Try it first (Build 5).** Implement the hourly mean with `index.floor("h")` and `groupby`, then compare with `resample("1h").mean()`.

```python
rng = np.random.default_rng(4)
idx = pd.date_range("2024-03-04 00:00", periods=144, freq="10min")
temp = pd.Series(65 + 6 * np.sin(2 * np.pi * np.arange(144) / 144) + rng.normal(0, 0.5, 144), index=idx, name="temp")
hourly_hand = temp.groupby(temp.index.floor("h")).mean()
hourly = temp.resample("1h").mean()
assert hourly_hand.index.equals(hourly.index) and np.allclose(hourly_hand, hourly)
print(len(temp), "->", len(hourly), "bins; first bin:", hourly.index[0], round(hourly.iloc[0], 4),
      "=", round(temp.iloc[:6].mean(), 4))
print(pd.to_datetime(pd.Series(["2024-03-04 08:00", "2024-03-04 08:10", "n/a"]), errors="coerce").tolist())
print(temp.resample("6h").ohlc().round(2).head(2))
print(temp.loc["2024-03-04 08:00":"2024-03-04 08:30"].size, temp.loc["2024-03"].size)
```
```text
144 -> 24 bins; first bin: 2024-03-04 00:00:00 65.6389 = 65.6389
[Timestamp('2024-03-04 08:00:00'), Timestamp('2024-03-04 08:10:00'), NaT]
                      open   high    low  close
2024-03-04 00:00:00  64.67  71.50  64.67  71.50
2024-03-04 06:00:00  71.48  71.48  65.15  65.19
4 144
```

`ohlc` gives the first, maximum, minimum and last value per bin. Other aggregations work as in `groupby` (`.agg`, `.size()`, `.max()`).

**Rolling windows.** `rolling(w)` slides a window of `w` rows; `rolling("1h")` a window of one hour *of time*, which is what you want for irregular sampling. The result at row $i$ uses rows $i-w+1, \dots, i$ (trailing) unless `center=True`, and is `NaN` until `min_periods` (default `w`) values are present. `ewm(halflife=...)` is the exponentially weighted mean, $y_i = (1-\alpha) y_{i-1} + \alpha x_i$ with `adjust=False`.

▶ **Try it first (Build 6).** PY-6 §3.11's `sliding_window_view(x, w)` gives an `(n - w + 1, w)` view; its row means are a rolling mean. How long is that array, and where do its values sit in `rolling(w).mean()`?

```python
from numpy.lib.stride_tricks import sliding_window_view

w = 6
by_view = sliding_window_view(temp.to_numpy(), w).mean(axis=1)
by_pandas = temp.rolling(w).mean()
assert len(by_view) == len(temp) - w + 1 and np.allclose(by_view, by_pandas.to_numpy()[w - 1:])
print(by_pandas.isna().sum(), temp.rolling(w, center=True).mean().isna().sum(),
      temp.rolling(w, min_periods=1).mean().isna().sum())
assert np.allclose(temp.rolling("1h").mean().iloc[w - 1:], by_pandas.iloc[w - 1:])
alpha = 0.3
ew, prev = [], temp.iloc[0]
for x in temp:
    prev = (1 - alpha) * prev + alpha * x
    ew.append(prev)
assert np.allclose(ew, temp.ewm(alpha=alpha, adjust=False).mean())
print(temp.shift(1).iloc[:3].round(2).tolist(), temp.diff().iloc[:3].round(2).tolist(),
      temp.pct_change().iloc[1:3].round(5).tolist())
```
```text
5 5 0
[nan, 64.67, 65.17] [nan, 0.5, 1.18] [0.00773, 0.01811]
```

■ **What just happened.** The view has $n-w+1 = 139$ windows; pandas returns $n$ values and puts window $j$ at row $j+w-1$ (its right end), so the first $w-1 = 5$ are `NaN`. Centering moves the label to the middle and leaves `NaN` at both ends; `min_periods=1` fills the start with shorter windows. On a regular grid `rolling("1h")` equals `rolling(6)` once a full hour is available (before that the time window averages what it has). `shift(1)` moves values one row later (lag), `diff()` is $x_i - x_{i-1}$ and `pct_change()` is $x_i/x_{i-1}-1$.

**Gaps and `asfreq`.** `asfreq("10min")` conforms a series to an exact grid, inserting `NaN` rows at missing stamps without aggregating (Example 4.6).

**Time zones, in one paragraph.** A naive timestamp has no zone; `tz_localize("Europe/Berlin")` declares one, `tz_convert("UTC")` converts. Store and compute in UTC; convert for display. Local clocks have days of 23 and 25 hours, so a daily resample of local time carries a different number of hours per day:

```python
for start, end in (("2024-03-30", "2024-04-02"), ("2024-10-26", "2024-10-29")):
    local = pd.date_range(start, end, freq="h", tz="Europe/Berlin", inclusive="left")
    print(start, "to", end, "local time:", len(local), "hours; per day", pd.Series(1, index=local).resample("D").size().tolist())
```
```text
2024-03-30 to 2024-04-02 local time: 71 hours; per day [24, 23, 24]
2024-10-26 to 2024-10-29 local time: 73 hours; per day [24, 25, 24]
```

A naive logger that writes local time without a zone duplicates an hour every October and skips one every March, and nothing in pandas will complain until your daily sums are off by 4%.

### 3.12 Reading and writing

`read_csv` is the workhorse, and its arguments are the first cleaning step: `usecols=` (read only some columns), `dtype={...}` (no guessing; `"category"` saves memory at read time), `parse_dates=[...]`, `na_values=[...]` (extra strings or numbers that mean missing, such as the sentinel), `skipinitialspace=True`, and `chunksize=` (iterate over pieces of a file too big for memory).

```python
r2 = pd.read_csv(io.StringIO(raw_csv), usecols=["time", "machine", "sensor", "value"],
                 dtype={"machine": "category"}, parse_dates=["time"], na_values=[999])
print(r2.dtypes.astype(str).to_dict(), int(r2["value"].isna().sum()))
tot, cnt = defaultdict(float), defaultdict(int)
for chunk in pd.read_csv(io.StringIO(raw_csv), chunksize=4, na_values=[999]):
    for unit, sm, ct in chunk.groupby("unit")["value"].agg(["sum", "count"]).itertuples():
        tot[unit] += sm                               # a loop over the groups of a chunk, not its rows
        cnt[unit] += ct
by_chunks = {u: round(tot[u] / cnt[u], 3) for u in sorted(tot)}
print(by_chunks, by_chunks == pd.read_csv(io.StringIO(raw_csv), na_values=[999]).groupby("unit")["value"].mean().round(3).to_dict())
buf = io.StringIO()
r2.to_csv(buf, index=False)
again = pd.read_csv(io.StringIO(buf.getvalue()), dtype={"machine": "category"}, parse_dates=["time"])
pd.testing.assert_frame_equal(again, r2)
print("CSV round trip OK")
```
```text
{'time': 'datetime64[us]', 'machine': 'category', 'sensor': 'str', 'value': 'float64'} 2
{'degC': 70.75, 'degF': 160.95, 'mm/s': 0.3} True
CSV round trip OK
```

`na_values=[999]` turned the sentinel into `NaN` at read time (two missing values now). The chunked loop reads four rows at a time and carries per-unit sums and counts, which add across chunks; per-chunk *means* do not. The CSV round trip needs the dtypes restated, because CSV stores text: categories, dtypes and time zones are lost. **Parquet** (`to_parquet`/`read_parquet`) is a binary, columnar format that keeps dtypes and compresses well; it needs the `pyarrow` package, which this run does not have, so it is described and not run here. `read_excel` and `to_excel` need `openpyxl`; `read_json(path, orient="records")` reads a list of objects, and `lines=True` reads one JSON object per line.

### 3.13 Plotting and performance

`df.plot(kind=..., ax=ax)` is PY-7's Matplotlib with defaults read from the frame: the index becomes the x axis, each column a line, column names the legend. Always pass `ax=` so the plot lands in a figure you made with the OO API, and label the axes yourself:

```python
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(9, 3))
temp.resample("1h").mean().plot(ax=axes[0], label="hourly mean")
temp.rolling(6, center=True).mean().plot(ax=axes[0], label="rolling 1 h")
axes[0].set(xlabel="time", ylabel="temperature (°C)")
axes[0].legend()
log.groupby("machine")["value"].plot(kind="hist", bins=40, alpha=0.5, ax=axes[1], legend=True)
axes[1].set(xlabel="temperature (°C)", ylabel="count")
print(len(axes[0].lines), len(axes[1].patches), [txt.get_text() for txt in axes[1].get_legend().get_texts()])
plt.close(fig)
```
```text
2 120 ['A', 'B', 'C']
```

`groupby(...).plot` draws one artist set per group on the same axes. Drop to Matplotlib (PY-7) when you need anything the frame does not know: error bands with `fill_between`, twin axes, annotations, a figure built for a paper. One trap (pandas 3.0 with Matplotlib 3.11): on a time axis that `df.plot` drew, a later `ax.plot` warns about incompatible date converters, and under `sharex=True` a frame whose index has no regular frequency can be drawn in different x units from one that has. Keep each time axis to `df.plot` plus `fill_between`, with `asfreq` on every index, or to Matplotlib alone (`project.py` Part 6).

**Performance rules,** in order of payoff. (1) Vectorize: column expressions, `.str`, `.dt`, `where`, `merge`, `groupby` (§3.7's table). (2) Never `iterrows` over a large frame; if a loop is truly needed, loop over `.to_numpy()` arrays or `itertuples()`, and in PY-10 compile it. (3) Use `category` for repeated strings and the narrowest numeric dtype that is safe. (4) Read only the columns you need. `df.query("value > 70 and machine == 'A'")` and `df.eval("F = value * 9 / 5 + 32")` take expression strings; they read well in chains and can be faster on very large frames, but the column names inside a string are invisible to your editor and linter, so prefer plain boolean masks in library code.

### 3.14 Tidy data, the cleaning pipeline (Build 7), and pandas versus NumPy

**The three rules** (Wickham 2014): every variable is a column; every observation is a row; every type of observational unit is a table. **Five common messes and their fixes:**

1. *Column headers are values, not variable names* (`temp_1`, `vib_1` as columns): `melt`.
2. *Several variables stored in one column* (`"A_temp_1"`): `.str.split(expand=True)` or `.str.extract`.
3. *Variables stored in both rows and columns* (a `kind` column with entries `min` and `max`, and days as columns): `melt`, then `pivot` the `kind` column out.
4. *Several types of observational unit in one table* (machine location repeated on every reading): split into a readings table and a machines table; `merge` (`validate="many_to_one"`) when needed.
5. *One type of unit spread over several tables* (one file per machine): `concat` with `keys=`.

▶ **Try it first (Build 7).** Clean §1's `raw` into a tidy table with one spelling per sensor, no duplicates, no sentinel, one unit, typed time. Decide the order of the steps: which step can *create* new duplicates?

```python
clean = (
    pd.read_csv(io.StringIO(raw_csv), parse_dates=["time"])
    # 1. one spelling per sensor: trim, lower-case, hyphen -> underscore
    .assign(sensor=lambda d: d["sensor"].str.strip().str.lower().str.replace("-", "_"))
    # 2. drop retransmitted rows; after step 1, so differently spelled copies also match
    .drop_duplicates()
    # 3. the logger's sentinel 999 means "no reading"
    .assign(value=lambda d: d["value"].mask(d["value"] >= 999))
    # 4. one unit: Fahrenheit rows to Celsius, then relabel the unit
    .assign(value=lambda d: d["value"].where(d["unit"] != "degF", (d["value"] - 32) * 5 / 9),
            unit=lambda d: d["unit"].replace({"degF": "degC"}))
    # 5. rows without a value cannot be used for means
    .dropna(subset=["value"])
    # 6. a stable order and a fresh 0..n-1 index
    .sort_values(["machine", "sensor", "time"])
    .reset_index(drop=True)
)
print(clean.round({"value": 2}))
print(len(raw), "rows;", len(raw.drop_duplicates()), "after drop_duplicates on raw;", len(clean), "after the chain")
print(clean.query("sensor == 'temp_1'").groupby("machine")["value"].mean().round(2).to_dict())
```
```text
                 time machine  sensor  value  unit
0 2024-03-04 08:00:00       A  temp_1  71.20  degC
1 2024-03-04 08:10:00       A  temp_1  71.90  degC
2 2024-03-04 08:00:00       A   vib_1   0.31  mm/s
3 2024-03-04 08:00:00       B  temp_1  68.00  degC
4 2024-03-04 08:00:00       C  temp_1  71.28  degC
5 2024-03-04 08:10:00       C  temp_1  72.00  degC
6 2024-03-04 08:10:00       C   vib_1   0.29  mm/s
10 rows; 10 after drop_duplicates on raw; 7 after the chain
{'A': 71.55, 'B': 68.0, 'C': 71.64}
```

■ **What just happened.** Ten rows became seven. Name normalization must come before `drop_duplicates`: on the raw names the retransmitted row is not a duplicate (second line). The chain reads top to bottom like the description of the cleaning, every intermediate is a new frame (cheap under Copy-on-Write), and you can comment out any line to debug. The last line answers §1's first question on this toy table; `project.py` does it on 50,000 rows.

**pandas or NumPy?** Labels, mixed dtypes, missing values, grouping, joining, time: pandas. Homogeneous numerical arrays, linear algebra, simulations, anything in an inner loop: NumPy. Cross with `.to_numpy()`: for a Series or a single-dtype frame it is a *read-only* view (Copy-on-Write protecting the frame; `.copy()` it before writing into it), for mixed dtypes a slow `object` copy, so select the numeric columns first. Come back with `pd.DataFrame(array, index=..., columns=...)`. ST-3's least squares is `np.linalg.lstsq(X, y)` on `df[cols].to_numpy()`, after pandas has made the table clean.

## 4. Worked examples

### Example 4.1: `loc` versus `iloc` on a non-default integer index

After sorting or filtering, an integer index stops matching positions. This is where `loc` and `iloc` bite.

```python
e = pd.DataFrame({"temp": [71.2, 69.9, 68.0, 72.4, 70.1]}, index=[10, 20, 30, 40, 50]).sort_values("temp")
print(e.index.tolist())
print("e.loc[20:40]  ->", e.loc[20:40, "temp"].tolist())
print("e.iloc[1:3]   ->", e.iloc[1:3]["temp"].tolist())
print("e.iloc[0]     ->", e.iloc[0].to_dict(), "| e.loc[30] ->", e.loc[30].to_dict())
try:
    e.loc[0]
except KeyError as err:
    print("e.loc[0]      -> KeyError", err)
s_int = pd.Series([1.0, 2.0, 3.0], index=[2, 0, 1])
print("s_int[0] =", s_int[0], "| s_int[0:2] =", s_int[0:2].tolist(), "| s_int.iloc[0] =", s_int.iloc[0])
```
```text
[30, 20, 50, 10, 40]
e.loc[20:40]  -> [69.9, 70.1, 71.2, 72.4]
e.iloc[1:3]   -> [69.9, 70.1]
e.iloc[0]     -> {'temp': 68.0} | e.loc[30] -> {'temp': 68.0}
e.loc[0]      -> KeyError 0
s_int[0] = 2.0 | s_int[0:2] = [1.0, 2.0] | s_int.iloc[0] = 1.0
```

`e.loc[20:40]` runs from the row *labelled* 20 to the row labelled 40 in current order, inclusive: four rows, including 50 and 10, which are not "between 20 and 40" numerically. `e.iloc[1:3]` is positions 1 and 2. `e.loc[0]` fails because no row is labelled 0. Worst of all, plain `[]` on an integer-indexed Series treats a scalar as a label (`s_int[0]` is 2.0, the value labelled 0) but a slice as positions (`s_int[0:2]` is the first two rows). Never use bare `[]` on rows: say `loc` or `iloc`.

### Example 4.2: chained assignment, triggered and fixed

Task: cap temperatures above 71 at 71.

```python
cap = pd.DataFrame({"temp": [71.2, 69.9, 72.4], "vib": [0.31, 0.28, 0.35]})
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    cap["temp"][cap["temp"] > 71] = 71.0                # chained: a getitem, then a setitem on its result
print(caught[0].category.__name__, "|", str(caught[0].message).splitlines()[0])
print("after chained assignment:", cap["temp"].tolist())
cap.loc[cap["temp"] > 71, "temp"] = 71.0                 # one setitem on cap itself
print("after .loc:              ", cap["temp"].tolist())
hot = cap[cap["vib"] > 0.3]                               # a new frame (behaves as a copy)
hot["vib"] = 0.0
print("cap untouched by hot:    ", cap["vib"].tolist())
```
```text
ChainedAssignmentError | A value is being set on a copy of a DataFrame or Series through chained assignment.
after chained assignment: [71.2, 69.9, 72.4]
after .loc:               [71.0, 69.9, 71.0]
cap untouched by hot:     [0.31, 0.28, 0.35]
```

`cap["temp"]` returned a Series; the `[...] = 71.0` then wrote into that Series, which under Copy-on-Write is detached from `cap` as soon as it is written. pandas 3 catches the pattern and says so; pandas 1–2 printed `SettingWithCopyWarning: A value is trying to be set on a copy of a slice from a DataFrame`, and the outcome depended on the pattern: `df["col"][mask] = v` usually did change `df` (the column was a view), `df[mask]["col"] = v` never did (the mask had made a copy). The fix is always the same: one `.loc[rows, col] = value`. The last three lines are the pattern that used to trigger the warning without being a bug: `hot` is a new table and changing it leaves `cap` alone.

### Example 4.3: a summary table by named aggregation

Per machine on the 200,000-row `log` of §3.8: number of readings, mean, standard deviation, 95th percentile, and the fraction above 64 °C.

```python
summary = (
    log.assign(hot=lambda d: d["value"] > 64.0)
    .groupby("machine")
    .agg(n=("value", "size"), mean=("value", "mean"), sd=("value", "std"),
         p95=("value", lambda v: v.quantile(0.95)), frac_hot=("hot", "mean"))
    .round(3)
)
print(summary)
print("check A:", round(log.loc[log["machine"] == "A", "value"].mean(), 3), int((log["machine"] == "A").sum()))
```
```text
              n    mean     sd     p95  frac_hot
machine
A        100372  61.500  2.285  65.245     0.138
B         59743  66.498  2.291  70.279     0.859
C         39885  59.501  2.290  63.254     0.024
check A: 61.5 100372
```

Each keyword becomes a column: `name=(column, function)`. The means sit 1.5 °C above the machine baselines 60, 65, 58 because the four equally likely sensors add offsets 0, 1, 2, 3 (mean 1.5, variance 1.25), and the `sd` of 2.29 is $\sqrt{2.0^2 + 1.25} = 2.291$: noise variance plus offset variance. `frac_hot` is the mean of a boolean column, which is a proportion. The `lambda` for `p95` runs once per group (three times), so it costs nothing; a `lambda` per row would.

### Example 4.4: `pivot_table` with margins

Mean reading by machine × sensor, with the overall means in the margins:

```python
pt = log.pivot_table(index="machine", columns="sensor", values="value", aggfunc="mean", margins=True)
print(pt.round(3))
shares = log["machine"].value_counts(normalize=True).sort_index()
print(f"grand mean {log['value'].mean():.3f} | share-weighted row means {(shares * pt['All'].iloc[:3]).sum():.3f}"
      f" | plain mean of row means {pt['All'].iloc[:3].mean():.3f} | shares {shares.round(3).tolist()}")
```
```text
sensor   temp_1  temp_2  temp_3  temp_4     All
machine
A        60.005  61.018  62.006  62.974  61.500
B        65.005  66.017  66.980  68.007  66.498
C        58.027  58.998  59.978  60.985  59.501
All      61.111  62.112  63.089  64.069  62.595
grand mean 62.595 | share-weighted row means 62.595 | plain mean of row means 62.500 | shares [0.502, 0.299, 0.199]
```

The `All` row and column are computed from the raw rows, not by averaging the cells. So the bottom-right 62.595 is the grand mean, which weights each machine by its share of readings (half are A's); the plain average of the three machine means is 62.500. When groups have unequal sizes the two differ, and you must say which one you report.

### Example 4.5: a merge with `indicator=True`

A per-machine table of readings meets the machines' metadata, where the maintenance system knows a machine D that never logged, and the logger knows a machine C that maintenance forgot:

```python
per_machine = log.groupby("machine", as_index=False).agg(n=("value", "size"))
meta = pd.DataFrame({"machine": ["A", "B", "D"], "site": ["Hall 1", "Hall 1", "Hall 2"], "year": [2015, 2019, 2023]})
joined = per_machine.merge(meta, on="machine", how="outer", validate="one_to_one", indicator=True)
print(joined)
print(joined["_merge"].value_counts().to_dict())
print("unmatched readings:", joined.loc[joined["_merge"] == "left_only", "n"].sum())
```
```text
  machine         n    site    year      _merge
0       A  100372.0  Hall 1  2015.0        both
1       B   59743.0  Hall 1  2019.0        both
2       C   39885.0     NaN     NaN   left_only
3       D       NaN  Hall 2  2023.0  right_only
{'both': 2, 'left_only': 1, 'right_only': 1}
unmatched readings: 39885.0
```

`_merge` says which side each row came from: C's 39,885 readings have no metadata (`left_only`), D has metadata and no readings (`right_only`). Run this check every time you attach a table; an inner join would have dropped C's readings without a word. Notice that `n` and `year` became `float64`, because the `NaN` the outer join inserted cannot live in an `int64` column (`.astype("Int64")` brings integers back, with `<NA>`).

### Example 4.6: `resample` across a gap, and the `interpolate` fix

A 10-minute temperature log, rising 0.5 °C per reading, loses its connection from 09:00 to 10:29:

```python
stamps = pd.date_range("2024-03-04 08:00", "2024-03-04 11:00", freq="10min")
lost = (stamps >= "2024-03-04 09:00") & (stamps < "2024-03-04 10:30")
gappy = pd.Series(70 + 0.5 * np.arange(len(stamps)), index=stamps).loc[~lost]
print(len(stamps), "->", len(gappy), "readings")
half = gappy.resample("30min").mean()
print({k.strftime("%H:%M"): v for k, v in half.items()})
print(gappy.asfreq("10min").isna().sum(), "NaN rows inserted by asfreq")
print(half.interpolate(method="time").tolist())
print(half.interpolate(method="time", limit=1).isna().sum(), "bins left NaN with limit=1")
```
```text
19 -> 10 readings
{'08:00': 70.5, '08:30': 72.0, '09:00': nan, '09:30': nan, '10:00': nan, '10:30': 78.0, '11:00': 79.0}
9 NaN rows inserted by asfreq
[70.5, 72.0, 73.5, 75.0, 76.5, 78.0, 79.0]
2 bins left NaN with limit=1
```

The 09:00, 09:30 and 10:00 bins have no readings, so their means are `NaN`; `resample` keeps them so the gap is visible, where `groupby` on `floor("30min")` would have produced a series with three bins silently missing. `asfreq` shows the same gap at the original resolution: nine missing stamps. `interpolate(method="time")` draws a straight line in time between the 08:30 and 10:30 bins; the underlying signal *is* a straight line and both neighbouring bins are complete, so the fill is exact (73.5 at 09:00 is the mean of the lost readings 73.0, 73.5, 74.0). `limit=1` fills at most one bin per gap, so a long outage stays visible. Filling a gap is a claim about what happened while you were not looking (§3.6): fine for a 30-minute dropout in a slow temperature, wrong for a vibration burst.

## 5. Why it matters and where it connects

- **ST-1, ST-2:** most estimates in their projects are a `groupby` or a column reduction on a tidy table; ST-2's bootstrap resamples rows with `df.sample(frac=1, replace=True)`. **ST-3:** the design matrix $X$ of its (3.3) is `df[cols].to_numpy()` plus a column of ones, after §3.14's cleaning; categorical predictors become dummy columns (`pd.get_dummies`); the regression in `project.py` previews it, including the attenuation by measurement error that ST-3 §3.11's assumption table warns about. **ST-7, ST-8:** `resample`, `shift`, `diff` and `rolling` are the first steps of every time-series analysis (the lag-$k$ sample autocorrelation is `s.autocorr(k)`). **ST-9:** §3.6's missing-data mechanisms are selection mechanisms.
- **ML-1:** train/test splits, `groupby`-aware splits (no machine in both sets), and feature tables are pandas; then `.to_numpy()` for the model. **PY-10:** profiling §3.7's slow paths; Numba for the loop that genuinely cannot be vectorized. **PY-7:** `df.plot` is Matplotlib, and every PY-7 rule about honest figures applies.

## 6. Exercises

See [exercises.md](exercises.md). A1–A5 are predict-and-fix drills (slices, alignment with a DataFrame, chained assignment, `size` versus `count`, `resample` bins). B1–B6 are the module: B1 (groupby with several statistics by hand), B2 (a `merge` by hand with property tests), B3 (tidying a messy wide table), B5 (the transform-based z-score in `project.py`) are essential. C1–C3: daylight-saving repair, an out-of-core daily summary with `chunksize`, and an irregular-time rolling mean with `searchsorted`.

## 7. Coding project

`project.py` generates a reproducible two-week sensor dataset (12 sensors × 3 machines, 10/15/20-minute sampling, about 50,000 rows) with planted problems: gaps, duplicate rows, an outlier burst, one machine logging Fahrenheit, inconsistent sensor names. It writes `figures/data.csv`, then: (1) cleans it in one method chain, with asserts on every count; (2) flags the burst by a rolling median/MAD (Hampel) rule built on `sliding_window_view`; (3) computes per-machine daily means and checks two machines against a pure-Python computation from the CSV; (4) builds the machine × sensor `pivot_table`; (5) merges the machine metadata with `validate=` and checks the unmatched rows; (6) aligns a 10-minute and a 15-minute sensor by `resample` and checks their correlation against `np.corrcoef`; (7) exports and reloads the tidy table; (8) plots daily means with a rolling band and small multiples per machine; (9) regresses temperature on pressure by `np.linalg.lstsq` and recovers the planted slope; (10) re-checks the numbers of §4. It runs in about 4 s.

**Build it yourself:** before opening the file, write the generator's planted problems and the cleaning chain from the description above. Then "Your turn": a `transform`-based z-score against the Hampel rule (Exercise B5 shows why a mean and a standard deviation hide the outliers they should find); a heatmap of missing data per sensor per day; the daily summary in `query` form; `iterrows` against vectorized on $10^5$ rows.

## 8. Where you'll struggle

- **`loc` slices include the end.** Every other slice in Python excludes it; `iloc` excludes it too. And plain `[]` on an integer index is a label for scalars and a position for slices (Example 4.1). The fix: never index rows with bare `[]`; write `loc` or `iloc`, and after a sort or filter that you index by position, `reset_index(drop=True)`.
- **Chained assignment.** `df[df.x > 0]["y"] = 1` reads like English and does nothing. In pandas 1–2 it raised `SettingWithCopyWarning`, and the common response was to silence the warning; the bug stayed. In pandas 3 it never works. One `.loc[mask, "y"] = 1`, always (Example 4.2). Code copied from a 2019 answer will contain this.
- **The MultiIndex after `groupby`.** Two keys or a list of aggregations give hierarchical labels you will fight with `.loc[("A", "temp_1")]` and `xs`. Use named aggregation (flat column names) and `reset_index()` (keys back as columns) until you need the hierarchy.
- **`merge` multiplies rows on duplicate keys, silently.** A metadata table with one duplicated machine row doubles that machine's readings and every sum computed from them. `validate=` is not optional; neither is comparing `len` before and after.
- **Time zones and daylight saving** will corrupt a resample once in your life: a naive local-time logger produces a duplicated hour in October, a missing hour in March, and 23- and 25-hour days (§3.11). Store UTC; localize at read time with `ambiguous=` and `nonexistent=` decided on purpose (Exercise C1).
- **`apply(lambda row: ...)` is the slow path and feels natural** (§3.7: hundreds of times slower). When you catch yourself writing it, stop and look for the column expression.
- **Untidy tables.** If an operation needs `iterrows`, a list of column names built by string formatting, or a loop over columns, the table is probably wide where it should be long. `melt` first.
- **The API moves.** pandas changes faster than NumPy: `append`, `ix`, `applymap`, `fillna(method=...)` and uppercase frequency aliases are gone or going, and Copy-on-Write changed the semantics of assignment. Read the `FutureWarning`s and `Pandas4Warning`s instead of silencing them; they are telling you what will break next year.

## 9. Check yourself

1. What is `pd.Series([1, 2], index=["a", "b"]) + pd.Series([10], index=["b"])`, and why?
2. How many rows does `df.loc["b":"d"]` return on the index `a..e`, and `df.iloc[1:3]`?
3. Why does `df["temp"][mask] = 0` never change `df` in pandas 3, and what is the one-line fix?
4. What is the difference between `groupby(...).agg("mean")` and `groupby(...).transform("mean")` in the shape of the result?
5. `size()` versus `count()`: which counts missing values?
6. When does `pivot` raise, and what does `pivot_table` do instead?
7. Left has key A twice, right has A three times. How many rows does an inner merge produce? Which `validate` value catches it?
8. Which bin does a reading at 09:00 fall in for `resample("1h")`, and how is that bin labelled?
9. Why does `resample` show `NaN` for an empty bin where `groupby(index.floor(...))` shows nothing?
10. Name the three tidy-data rules and the fix for "column headers are values".
11. Why must sensor names be normalized before `drop_duplicates`?
12. When would you leave pandas for NumPy, and how do you cross?

## 10. Pair with

- **McKinney, *Python for Data Analysis* 3e (free at wesmckinney.com/book), ch. 5–8 and 10–11.** By pandas's author; more examples of every method here, and chapter 11 on time series is the best single reference for `resample` and offsets. Written for pandas 1.x: read it with §3.3's Copy-on-Write rules in mind.
- **pandas docs: "10 minutes to pandas", and the User Guide chapters "Group by", "Merge, join, concatenate and compare", "Reshaping and pivot tables", "Time series / date functionality" and "Copy-on-Write" (pandas.pydata.org).** The authoritative semantics, current for pandas 3.
- **VanderPlas, *Python Data Science Handbook* 2e, ch. 3.** A faster tour with good pictures of alignment and `groupby`.
- **Wickham, "Tidy Data", *Journal of Statistical Software* 59(10), 2014.** The principles of §3.14; short and worth reading in full.
- **pandas "Comparison with SQL".** If you know SQL: `merge` is `JOIN`, `groupby` is `GROUP BY`, `query` is `WHERE`.
