# PY-9 Exercises

Conventions as in the notes: pandas 3, `import pandas as pd`, `rng = np.random.default_rng(seed)`, selection by `loc`/`iloc` only, multi-step transforms as one method chain. Put your answers in a file next to `project.py`; importing `project` re-runs its checks (a few seconds), so copy the helpers you need. **Write each prediction down before you run anything**, and turn every "check" into an `assert`.

## A. Warm-up (do all; 5–10 minutes each)

**A1.** (Labels, positions, and bare `[]`.) For `s = pd.Series([10, 20, 30, 40, 50], index=[5, 4, 3, 2, 1])`, predict `s.loc[4:2]`, `s.iloc[4:2]`, `s.iloc[1:3]`, `s[1]`, `s[1:3]`, `s.loc[[1, 5]]` and `s.sort_index().loc[2:4]`. Run them, and for each one that surprised you write the rule in one line.

**A2.** (Alignment against columns.) For `df = pd.DataFrame({"a": [1.0, 2.0], "b": [3.0, 4.0]}, index=["x", "y"])`, predict the shape and contents of `df - df.mean()`, `df - df.mean(axis=1)` and `df.sub(df.mean(axis=1), axis=0)`. Which one centres each row? Explain the middle one using the alignment rule of notes §3.5, then write `center_rows(df)` and assert that every row of its result sums to zero.

**A3.** (Fix the bug.) A colleague's cleaning script:

```python
df = pd.read_csv(io.StringIO("machine,value\nA,1.0\nA,\nB,3.0\nB,500.0\n"))
df[df["machine"] == "B"]["value"] = df["value"] / 10      # B logs in tenths of a unit
df["value"].fillna(0, inplace=True)                        # treat missing as zero
df["big"] = df.apply(lambda row: row["value"] > 2, axis=1)
```

Predict `df` after each line under pandas 3. Name the three problems (one line does nothing, one does nothing *and* encodes a modelling mistake, one is the slow path), and rewrite the script as one method chain that does what was intended, with missing values left missing.

**A4.** (`size`, `count`, `nunique`, missing keys.) For `d = pd.DataFrame({"g": ["a", "a", "b", None, "b", "b"], "v": [1.0, 1.0, np.nan, 2.0, 3.0, np.nan]})`, predict `d.groupby("g").size()`, `d.groupby("g")["v"].count()`, `d.groupby("g")["v"].nunique()`, `d.groupby("g", dropna=False).size()` and `d["g"].value_counts(dropna=False)`. Which of them would you use to report "readings per sensor", and which "valid readings per sensor"?

**A5.** (Resample bins.) `r` holds the values 1, 2, 3, 4 at 08:05, 08:55, 09:00 and 11:30 on one day. Predict the index and values of `r.resample("1h").mean()`, `r.resample("1h").size()`, `r.groupby(r.index.floor("h")).mean()` and `r.resample("1h", closed="right", label="right").mean()`. Which reading changes bin between the default and `closed="right"`?

## B. Core (do all; these are the module)

**B1.** (Groupby with several statistics, by hand.) On the 200,000-row `log` of notes §3.8, compute per (machine, sensor) the count, mean, sample standard deviation, minimum and maximum *without* `groupby`: turn the two keys into one integer code per row (`pd.factorize` on each key, then combine; or `np.unique` on a structured array), then use `np.bincount` (with `weights=` for sums and sums of squares) and `np.minimum.at`/`np.maximum.at`. Return a DataFrame indexed like `log.groupby(["machine", "sensor"])["value"].agg(["count", "mean", "std", "min", "max"])` and assert it equal to that to $10^{-9}$. Time both. Which formula for the standard deviation did you use, and when would the one-pass "sum of squares" formula lose digits?

**B2.** (A merge, by hand, against pandas.) Write `merge_by_hand(left, right, on, how)` for `how` in `inner`, `left`, `right`, `outer`: it takes two DataFrames, works only on their `.to_dict("records")` lists with loops, and returns a list of dicts with an `_merge` field like `indicator=True` and `None` for the missing side. Then write `check_validate(left, right, on, kind)` that raises `ValueError` exactly when pandas' `validate=kind` would raise, for the four kinds. Property-test both on 300 random pairs of small tables (0–6 rows, keys from `"abc"`, so duplicates are common) by comparing with `pd.merge(..., indicator=True)` after sorting both results. For one pair with duplicates on both sides, state the row count as a formula in the key multiplicities.

**B3.** (Tidy a messy table.) A maintenance spreadsheet arrives wide:

```python
wide = pd.DataFrame({"date": pd.date_range("2024-03-04", periods=3, freq="D"),
                     "A_temp_mean": [53.0, 53.1, 53.4], "A_temp_max": [62.3, 62.8, 61.6],
                     "B_temp_mean": [54.1, 54.6, 54.2], "B_temp_max": [63.2, 63.0, 62.7],
                     "B_vib_mean": [0.94, 0.95, 0.93], "B_vib_max": [1.41, 1.38, 1.40]})
```

Identify which of the five messes of notes §3.14 it contains. Produce (a) the fully long table with columns `date, machine, kind, stat, value` (18 rows), and (b) the tidy table with columns `date, machine, kind, mean, max` (9 rows), each as one method chain. Then (c) rebuild `wide` exactly from (b) and assert equality, column order included.

**B4.** (Two sampling rates, three alignments.) Simulate one day of a latent load $L(t)=0.5+0.3\sin(2\pi t/24\,\mathrm h)$ and two sensors: temperature $40+25L+\varepsilon_T$ every 10 minutes and pressure $2+3L+\varepsilon_P$ every 15 minutes, both with timestamps jittered by up to ±2 minutes (uniform), $\sigma_T=0.5$, $\sigma_P=0.05$, seed 3. Align them (a) by an exact-timestamp join, (b) by resampling both to 30-minute means, (c) with `pd.merge_asof(..., direction="nearest", tolerance=pd.Timedelta("5min"))` (read its docstring). Report the number of aligned pairs and the correlation for each, and explain in three sentences when you would use (b) and when (c).

**B5.** (Extend `project.py`: three outlier rules.) Add to `project.py` a global per-sensor z-score $z=(x-\bar x_g)/s_g$ computed with `groupby(["machine", "sensor"])["value"].transform`, and a per-sensor-per-day version (add `pd.Grouper(key="time", freq="D")` to the keys). For burst sizes 30, 15 and 8 °C (change `BURST`), tabulate for each rule (global $|z|>4$, daily $|z|>4$, Part 2's Hampel rule with $k=10$) how many of the 6 burst rows it flags and how many other rows. Explain why the daily z-score of the 30 °C burst is *smaller* than the global one, and what that says about using a mean and a standard deviation to detect the outliers that contaminate them.

**B6.** (Performance, measured.) On $10^5$ rows with columns `value` (float) and `unit` (`"degC"`, `"degF"` or `"K"`), convert every value to Celsius six ways: an `iterrows` loop, an `itertuples` loop, `apply(axis=1)`, `np.select` on the columns, `map` of per-unit scale and offset, and a loop over `.to_numpy()` arrays in plain Python. Assert they agree, report the times and ratios, and then measure the memory and the `groupby("unit")["value"].mean()` time with `unit` as `str` and as `category`. Which result would change if the frame had $10^7$ rows, and which would not?

## C. Challenge (do at least one)

**C1.** (Daylight saving, repaired.) A logger in Berlin writes naive local time every 15 minutes from 2024-10-26 00:00 to 2024-10-28 00:00 local; during the October change it writes the 02:00–02:45 stamps twice (once in summer time, once in winter time). Build that list of strings honestly (generate the true instants in UTC, convert to local, drop the zone, format). Show that `tz_localize("Europe/Berlin")` raises on it, then repair it with `ambiguous="infer"` and verify against the true instants. Report the number of readings per local calendar day and per UTC day. Then build the March case (2024-03-31) and explain what `nonexistent=` is for. Write a five-line design note: what should the logger write instead, and why?

**C2.** (Out of core.) Write a CSV with 2,000,000 rows of `time, machine, sensor, value` (2 machines, 4 sensors, 1-second sampling over about 3 days, seed 5), with 1% of rows duplicated and placed *anywhere* in the file. Compute exact per-machine daily means by reading it with `chunksize=200_000`, carrying sums and counts, and handling duplicates that land in different chunks. Compare with a full `read_csv` and report both run times and the peak memory of each approach (`tracemalloc`). What would break if you deduplicated within each chunk only?

**C3.** (A time-based rolling mean by hand.) For an irregular `DatetimeIndex` (2,000 readings at exponentially distributed intervals with mean 1 minute, seed 8), implement the trailing one-hour mean of `rolling("1h")` with `np.searchsorted` on the int64 timestamps and a cumulative sum, in $O(n \log n)$ (a two-pointer sweep over the sorted stamps would make it $O(n)$; not required). Determine from your implementation's agreement with pandas whether the window is $(t-1\,\mathrm h,\,t]$ or $[t-1\,\mathrm h,\,t]$, then check `closed="both"`. Time both for $n=10^6$.
