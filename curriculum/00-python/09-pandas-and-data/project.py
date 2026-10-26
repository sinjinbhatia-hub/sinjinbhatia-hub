# %% [markdown]
# # PY-9 Project: Two Weeks of Sensor Data, Cleaned, Summarized, Aligned and Regressed
#
# This script generates a reproducible two-week log from 12 sensors on each of 3 machines
# with planted problems (gaps, duplicates, an outlier burst, one machine in Fahrenheit,
# inconsistent sensor names), writes it to figures/data.csv, and then does what the notes
# teach: a cleaning method chain (notes 3.14), outlier flags by a rolling median/MAD rule
# (3.11), daily summaries by groupby checked against pure Python (3.8), a pivot_table (3.9),
# a validated merge (3.10), two sampling rates aligned by resample (3.11), a CSV round trip
# (3.12), figures (3.13), a least-squares fit on the clean table (ST-3 preview), and a
# re-check of the numbers quoted in the notes' section 4.
#
# Run headless: `MPLBACKEND=Agg python3 project.py` (about 4 s). Files go to ./figures/.

# %%
import csv
import importlib.util
import io
import time
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numpy.lib.stride_tricks import sliding_window_view

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
T_START = time.perf_counter()
SEED = 2024

MACHINES = ["A", "B", "C"]
KINDS = {"temp": ("degC", "10min"), "pres": ("bar", "15min"), "vib": ("mm/s", "20min")}
SENSORS = [f"{kind}_{k}" for kind in KINDS for k in range(1, 5)]           # 12 sensors per machine
START, DAYS = pd.Timestamp("2024-03-04"), 14
GAP = ("A", pd.Timestamp("2024-03-09 02:00"), pd.Timestamp("2024-03-09 08:00"))   # A offline 6 h
DROP_FRAC = 0.01                                                            # random lost rows
N_DUP = 400                                                                 # retransmitted rows
BURST = ("B", "temp_3", pd.Timestamp("2024-03-12 14:00"), 6, 30.0)          # 6 readings, +30 degC
F_MACHINE = "C"                                                             # logs temperature in degF
TEMP0, TEMP_GAIN, PRES0, PRES_GAIN = 40.0, 25.0, 2.0, 3.0                    # temp = 40 + 25 load, pres = 2 + 3 load
TRUE_SLOPE = TEMP_GAIN / PRES_GAIN                                           # planted temp-on-pres slope, 8.333


# %% [markdown]
# ## Part 1: the generator
#
# Each machine has a latent load on a 5-minute master grid: a daily cycle plus a slow AR(1)
# wander. Every sensor samples that load at its own rate (10, 15 or 20 minutes) with its own
# noise. Then the problems are planted in the order a real logger creates them: lost rows,
# the burst, retransmissions, Fahrenheit on machine C, and name spellings that differ by
# machine (and, on C, from row to row, so a retransmitted copy may be spelled differently).

# %%
def machine_load(n: int, rng: np.random.Generator) -> np.ndarray:
    """Load on the 5-min grid: 0.5 + 0.3 sin(2 pi (h - 6) / 24) + AR(1) wander (phi 0.995, sd 0.1)."""
    hours = np.arange(n) * 5 / 60
    eps = rng.normal(0, 0.1 * np.sqrt(1 - 0.995**2), n)
    ar = np.zeros(n)
    for i in range(1, n):                       # a recursion, so a loop (PY-10 would compile it)
        ar[i] = 0.995 * ar[i - 1] + eps[i]
    return 0.5 + 0.3 * np.sin(2 * np.pi * (hours - 6) / 24) + ar


def make_truth(rng: np.random.Generator) -> pd.DataFrame:
    """The tidy table the logger should have written: time, machine, sensor, value, unit."""
    grid = pd.date_range(START, periods=DAYS * 24 * 12, freq="5min")
    parts = []
    for m in MACHINES:
        load = pd.Series(machine_load(len(grid), rng), index=grid)
        for sensor in SENSORS:
            kind, k = sensor.split("_")
            stamps = pd.date_range(START, START + pd.Timedelta(days=DAYS), freq=KINDS[kind][1], inclusive="left")
            x = load.loc[stamps].to_numpy()
            noise = {"temp": 0.5, "pres": 0.05, "vib": 0.03}[kind] * rng.standard_normal(len(x))
            value = {"temp": TEMP0 + TEMP_GAIN * x + 0.5 * (int(k) - 1),
                     "pres": PRES0 + PRES_GAIN * x, "vib": 0.3 + 1.2 * x}[kind] + noise
            parts.append(pd.DataFrame({"time": stamps, "machine": m, "sensor": sensor,
                                       "value": value, "unit": KINDS[kind][0]}))
    return pd.concat(parts, ignore_index=True)


def plant_problems(truth: pd.DataFrame, rng: np.random.Generator) -> tuple[pd.DataFrame, dict]:
    """Return the raw log (messy, shuffled) and a record of everything planted."""
    m, t0, t1 = GAP
    in_gap = (truth["machine"] == m) & truth["time"].between(t0, t1, inclusive="left")
    lost = (rng.random(len(truth)) < DROP_FRAC) & ~in_gap
    log = truth.loc[~in_gap & ~lost].reset_index(drop=True)
    bm, bs, bt, blen, bsize = BURST
    burst_idx = log.index[(log["machine"] == bm) & (log["sensor"] == bs) & (log["time"] >= bt)][:blen]
    log.loc[burst_idx, "value"] += bsize
    burst_times = set(log.loc[burst_idx, "time"])
    log = pd.concat([log, log.iloc[rng.choice(len(log), N_DUP, replace=False)]], ignore_index=True)
    is_f = (log["machine"] == F_MACHINE) & log["sensor"].str.startswith("temp")
    log.loc[is_f, "value"] = log.loc[is_f, "value"] * 9 / 5 + 32
    log.loc[is_f, "unit"] = "degF"
    names = log["sensor"].where(log["machine"] != "B", log["sensor"].str.capitalize())       # Temp_1
    shout = (log["machine"] == "C") & (rng.random(len(log)) < 0.5)
    log["sensor"] = names.where(~shout, names.str.upper().str.replace("_", "-") + " ")      # "TEMP-1 "
    log = log.sample(frac=1.0, random_state=SEED).reset_index(drop=True)                    # loggers interleave
    return log, {"gap": int(in_gap.sum()), "lost": int(lost.sum()), "burst_times": burst_times}


rng = np.random.default_rng(SEED)
truth = make_truth(rng)
raw_mem, planted = plant_problems(truth, rng)
raw_mem.to_csv(FIG / "data.csv", index=False)
raw = pd.read_csv(FIG / "data.csv")                     # from here on, work from the file like a user
print(f"truth {len(truth):,} rows; raw log {len(raw):,} rows ({planted['gap']} in the gap, "
      f"{planted['lost']} lost, {N_DUP} duplicated); {raw['sensor'].nunique()} spellings of 12 sensor names")
assert len(truth) == 3 * 4 * (2016 + 1344 + 1008)
assert len(raw) == len(truth) - planted["gap"] - planted["lost"] + N_DUP
assert raw["sensor"].nunique() > 12 and len(planted["burst_times"]) == BURST[3]


# %% [markdown]
# ## Part 2: the cleaning chain (notes 3.14, Build 7) and the outlier rule (3.11)
#
# One method chain, one comment per step. Names are normalized *before* drop_duplicates, so
# retransmitted copies spelled differently are caught. Then outliers are flagged per sensor by
# a Hampel rule: in the centered window W_i of 25 readings around x_i, robust z_i =
# |x_i - med(W_i)| / (1.4826 median|W_i - med(W_i)|); 1.4826 makes the MAD estimate the standard
# deviation for Gaussian noise. The windows come from sliding_window_view (notes Build 6). Do not
# shortcut it as "rolling median of |x - rolling median|": on a steady ramp each point *is* its
# window's median, those deviations are 0, and that MAD collapses to 0 (try it).

# %%
def clean_log(raw: pd.DataFrame) -> pd.DataFrame:
    """Raw logger table -> tidy, typed, deduplicated table in one unit per kind (degC, bar, mm/s)."""
    return (
        raw
        # 1. one spelling per sensor: trim, lower-case, hyphen -> underscore
        .assign(sensor=lambda d: d["sensor"].str.strip().str.lower().str.replace("-", "_"))
        # 2. drop retransmissions (after step 1, so differently spelled copies match)
        .drop_duplicates()
        # 3. one unit per kind: Fahrenheit rows to Celsius
        .assign(value=lambda d: d["value"].where(d["unit"] != "degF", (d["value"] - 32) * 5 / 9),
                unit=lambda d: d["unit"].replace({"degF": "degC"}))
        # 4. types: real timestamps, the sensor kind as its own variable, categories for repeated strings
        .assign(time=lambda d: pd.to_datetime(d["time"]), kind=lambda d: d["sensor"].str.split("_").str[0])
        .astype({"machine": "category", "sensor": "category", "unit": "category", "kind": "category"})
        # 5. a stable order (rolling windows below need time order within each sensor)
        .sort_values(["machine", "sensor", "time"])
        .reset_index(drop=True)
    )


def hampel_z(v: np.ndarray, window: int) -> np.ndarray:
    """|x_i - med(W_i)| / (1.4826 MAD(W_i)), W_i the centered window of `window` readings (shorter at the ends)."""
    h = window // 2
    win = sliding_window_view(np.concatenate([np.full(h, np.nan), v, np.full(h, np.nan)]), window)   # (n, window)
    med = np.nanmedian(win, axis=1)
    mad = np.nanmedian(np.abs(win - med[:, None]), axis=1)
    return np.abs(v - med) / (1.4826 * mad)


def flag_rolling_mad(df: pd.DataFrame, window: int = 25, k: float = 10.0) -> pd.DataFrame:
    """Add robust_z (hampel_z per machine and sensor, rows in time order) and outlier = robust_z > k."""
    z = df.groupby(["machine", "sensor"], observed=True)["value"].transform(lambda v: hampel_z(v.to_numpy(), window))
    return df.assign(robust_z=z, outlier=z > k)


clean = flag_rolling_mad(clean_log(raw))
naive_dedup = len(raw) - len(raw.drop_duplicates())
print(f"clean: {len(clean):,} rows; drop_duplicates on the raw names would catch only {naive_dedup} of {N_DUP}")
assert len(raw) - len(clean) == N_DUP and naive_dedup < N_DUP
assert sorted(clean["sensor"].cat.categories) == sorted(SENSORS) and set(clean["unit"]) == {"degC", "bar", "mm/s"}
check = clean.merge(truth, on=["time", "machine", "sensor"], how="outer", validate="one_to_one",
                    indicator=True, suffixes=("", "_true"))
assert (check["_merge"] == "left_only").sum() == 0
assert (check["_merge"] == "right_only").sum() == planted["gap"] + planted["lost"]
err = (check["value"] - check["value_true"]).abs()
assert np.isclose(err.nlargest(BURST[3]), BURST[4]).all() and err.nlargest(BURST[3] + 1).iloc[-1] < 1e-9
flagged = clean.loc[clean["outlier"]]
print(f"rolling-MAD rule flags {len(flagged)} rows: {flagged[['machine', 'sensor']].drop_duplicates().to_numpy().tolist()}"
      f", robust z from {flagged['robust_z'].min():.1f} to {flagged['robust_z'].max():.1f}; "
      f"largest unflagged z {clean.loc[~clean['outlier'], 'robust_z'].max():.2f}")
assert set(flagged["time"]) == planted["burst_times"] and set(flagged["machine"]) == {"B"}


# %% [markdown]
# ## Part 3: daily summaries by groupby, checked against pure Python (notes 3.8, 3.11)
#
# Named aggregation per machine and calendar day for the temperature sensors (outliers
# excluded). The check re-reads data.csv with the csv module and redoes the cleaning for
# machines A and C in plain Python, with no pandas at all.

# %%
temp = clean.loc[(clean["kind"] == "temp") & ~clean["outlier"]]
daily = (temp.groupby(["machine", pd.Grouper(key="time", freq="D")], observed=True)["value"]
         .agg(mean="mean", min="min", max="max", n="size")
         .reset_index())
print(daily.head(3).round({"mean": 2, "min": 2, "max": 2}).to_string(index=False))


def daily_temp_means_pure_python(path: Path, machine: str) -> dict[str, float]:
    """Mean temperature (degC) per calendar day for one machine, from the raw CSV, by dicts and sets."""
    seen, acc = set(), defaultdict(lambda: [0.0, 0])
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            name = row["sensor"].strip().lower().replace("-", "_")
            if row["machine"] != machine or not name.startswith("temp") or (row["time"], name) in seen:
                continue
            seen.add((row["time"], name))
            v = float(row["value"])
            acc[row["time"][:10]][0] += (v - 32) * 5 / 9 if row["unit"] == "degF" else v
            acc[row["time"][:10]][1] += 1
    return {day: s / n for day, (s, n) in sorted(acc.items())}


for m in ("A", "C"):
    by_hand = daily_temp_means_pure_python(FIG / "data.csv", m)
    sel = daily.loc[daily["machine"] == m]
    by_pandas = pd.Series(sel["mean"].to_numpy(), index=sel["time"].dt.strftime("%Y-%m-%d"))
    assert list(by_hand) == by_pandas.index.tolist() and np.allclose(list(by_hand.values()), by_pandas)
assert len(daily) == 3 * DAYS and daily.loc[(daily["machine"] == "A") & (daily["time"] == "2024-03-09"), "n"].item() == \
    len(temp.loc[(temp["machine"] == "A") & (temp["time"].dt.day == 9)])
print("daily means of machines A and C agree with the pure-Python computation")


# %% [markdown]
# ## Part 4: pivot_table (3.9) and a validated merge (3.10)

# %%
pt = clean.loc[~clean["outlier"]].pivot_table(index="machine", columns="sensor", values="value",
                                              aggfunc="mean", observed=True)
print(pt[["temp_1", "temp_2", "temp_3", "temp_4", "pres_1", "vib_1"]].round(3))
assert pt.shape == (3, 12)
assert np.allclose(pt, clean.loc[~clean["outlier"]].groupby(["machine", "sensor"], observed=True)["value"].mean().unstack())
offsets = pt[[f"temp_{k}" for k in range(1, 5)]].sub(pt["temp_1"], axis=0)
assert np.allclose(offsets, [[0.0, 0.5, 1.0, 1.5]] * 3, atol=0.05)          # planted offsets 0.5 (k - 1)

meta = pd.DataFrame({"machine": ["A", "B", "C", "D"], "site": ["Hall 1", "Hall 1", "Hall 2", "Hall 2"],
                     "installed": [2015, 2019, 2021, 2024]})
joined = meta.merge(daily.astype({"machine": str}), on="machine", how="left", validate="one_to_many", indicator=True)
counts = joined["_merge"].value_counts()
print(f"merge: {counts.to_dict()}; unmatched metadata: {joined.loc[joined['_merge'] == 'left_only', 'machine'].tolist()}")
assert counts["both"] == len(daily) and counts["left_only"] == 1 and counts["right_only"] == 0
try:
    pd.concat([meta, meta.iloc[[1]]]).merge(daily.astype({"machine": str}), on="machine", validate="one_to_many")
    raise AssertionError("a duplicated metadata row must fail validation")
except pd.errors.MergeError:
    print("validate='one_to_many' rejects a duplicated metadata row")


# %% [markdown]
# ## Part 5: two sampling rates aligned by resample (3.11), and the CSV round trip (3.12)
#
# temp_1 is sampled every 10 minutes and pres_1 every 15. Joining on exact timestamps keeps
# only the stamps they share (every 30 min, minus lost rows); resampling both to 30-minute
# means uses every reading.

# %%
def series_of(df: pd.DataFrame, machine: str, sensor: str) -> pd.Series:
    """One sensor's clean readings as a time-indexed Series, outliers removed."""
    sel = df.loc[(df["machine"] == machine) & (df["sensor"] == sensor) & ~df["outlier"]]
    return sel.set_index("time")["value"]


t1, p1 = series_of(clean, "A", "temp_1"), series_of(clean, "A", "pres_1")
exact = pd.concat({"temp": t1, "pres": p1}, axis=1, join="inner")
aligned = pd.concat({"temp": t1.resample("30min").mean(), "pres": p1.resample("30min").mean()}, axis=1).dropna()
r_pd = aligned["temp"].corr(aligned["pres"])
r_np = np.corrcoef(aligned["temp"].to_numpy(), aligned["pres"].to_numpy())[0, 1]
print(f"machine A: {len(t1)} temp and {len(p1)} pres readings; {len(exact)} shared stamps; {len(aligned)} aligned "
      f"30-min bins; corr pandas {r_pd:.6f}, np.corrcoef {r_np:.6f}")
assert abs(r_pd - r_np) < 1e-12 and r_pd > 0.95

out = FIG / "clean_long.csv"
clean.to_csv(out, index=False)
back = pd.read_csv(out, parse_dates=["time"], dtype={c: "category" for c in ("machine", "sensor", "unit", "kind")})
pd.testing.assert_frame_equal(back, clean)
print(f"CSV round trip of {len(clean):,} rows: identical")
if importlib.util.find_spec("pyarrow") is not None:
    clean.to_parquet(FIG / "clean_long.parquet")
    pd.testing.assert_frame_equal(pd.read_parquet(FIG / "clean_long.parquet"), clean)
    print("Parquet round trip: identical, dtypes kept without being restated")


# %% [markdown]
# ## Part 6: figures (3.13)
#
# (a) One panel per machine: hourly mean temperature, its 24-hour rolling mean (df.plot), a
# +-2 sd band of the hourly values in that window (Matplotlib's fill_between), and the daily
# means as dashed steps. (b) Small multiples: machine x sensor kind, hourly means per sensor,
# flagged outliers in red. Keep each time axis to pandas' plotting (plus fill_between) and
# give every plotted index a regular frequency (asfreq): with pandas 3 and Matplotlib 3.11,
# mixing ax.plot with df.plot, or regular with irregular indexes, on one time axis warns about
# date converters or, with sharex, puts the curves in different units.

# %%
hourly = (temp.groupby(["machine", pd.Grouper(key="time", freq="h")], observed=True)["value"].mean()
          .unstack("machine").asfreq("h"))
roll = hourly.rolling(24, center=True, min_periods=12)
mean24, sd24 = roll.mean(), roll.std()
daily_wide = daily.pivot(index="time", columns="machine", values="mean").asfreq("D")
fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
for ax, m, color in zip(axes, MACHINES, ["C0", "C1", "C2"]):
    hourly[m].plot(ax=ax, lw=0.5, alpha=0.6, color=color, label="hourly mean")
    mean24[m].plot(ax=ax, lw=2, color=color, label="24 h rolling mean")
    ax.fill_between(hourly.index, mean24[m] - 2 * sd24[m], mean24[m] + 2 * sd24[m], color=color, alpha=0.15,
                    label="± 2 sd within the 24 h window")
    daily_wide[m].plot(ax=ax, drawstyle="steps-post", ls="--", lw=1, color="k", label="daily mean")
    ax.set(ylabel=f"{m}: temperature (°C)", xlabel="time")
axes[0].legend(ncol=4, fontsize=8, loc="upper left")
axes[0].set_title("Machine temperature (temp_1..temp_4 averaged); the 6 h gap on A shows as a break")
fig.tight_layout()
fig.savefig(FIG / "daily_rolling_band.png", dpi=120)
plt.close(fig)

fig, axes = plt.subplots(3, 3, figsize=(12, 7), sharex=True)
for i, m in enumerate(MACHINES):
    for j, kind in enumerate(KINDS):
        sub = clean.loc[(clean["machine"] == m) & (clean["kind"] == kind)]
        wide = sub.pivot_table(index=pd.Grouper(key="time", freq="h"), columns="sensor", values="value", observed=True).asfreq("h")
        wide.plot(ax=axes[i, j], lw=0.6, legend=(i == 0))
        bad = sub.loc[sub["outlier"]].set_index("time")["value"]
        if len(bad):
            bad.resample("h").max().reindex(wide.index).plot(ax=axes[i, j], style="r.", ms=7, label="outlier")
        axes[i, j].set(ylabel=f"{m}: {kind} ({KINDS[kind][0].replace('degC', '°C')})", xlabel="time")
fig.tight_layout()
fig.savefig(FIG / "small_multiples.png", dpi=110)
plt.close(fig)
assert (FIG / "daily_rolling_band.png").stat().st_size > 10_000 and (FIG / "small_multiples.png").stat().st_size > 10_000


# %% [markdown]
# ## Part 7: regression of temperature on pressure (ST-3 preview)
#
# The generator plants temp = 40 + 25 load and pres = 2 + 3 load, so temp = 23.33 + 8.333 pres
# without noise. Least squares on the aligned 30-minute table (all three machines) via
# np.linalg.lstsq on .to_numpy() arrays: pandas made the table, NumPy does the algebra.
# Noise in the predictor shrinks the slope by var(signal) / (var(signal) + var(noise)).

# %%
pairs = pd.concat([pd.concat({"temp": series_of(clean, m, "temp_1").resample("30min").mean(),
                              "pres": series_of(clean, m, "pres_1").resample("30min").mean()}, axis=1).dropna()
                   for m in MACHINES], keys=MACHINES, names=["machine", "time"])
X = np.column_stack([np.ones(len(pairs)), pairs["pres"].to_numpy()])
(b0, b1), *_ = np.linalg.lstsq(X, pairs["temp"].to_numpy(), rcond=None)
noise_var = 0.05**2 / 2                                                       # two 15-min readings per bin
shrink = (pairs["pres"].var() - noise_var) / pairs["pres"].var()
print(f"temp = {b0:.3f} + {b1:.4f} pres on {len(pairs)} bins; planted slope {TRUE_SLOPE:.4f}, "
      f"intercept {TEMP0 - TRUE_SLOPE * PRES0:.3f}; predicted attenuation factor {shrink:.4f}")
assert abs(b1 - TRUE_SLOPE) < 0.1 and abs(b1 - TRUE_SLOPE * shrink) < 0.05 and abs(b0 - (TEMP0 - TRUE_SLOPE * PRES0)) < 0.5


# %% [markdown]
# ## Part 8: the numbers quoted in the notes (sections 3.4, 3.11, 3.14 and 4)

# %%
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
toy = (pd.read_csv(io.StringIO(raw_csv), parse_dates=["time"])
       .assign(sensor=lambda d: d["sensor"].str.strip().str.lower().str.replace("-", "_"))
       .drop_duplicates()
       .assign(value=lambda d: d["value"].mask(d["value"] >= 999))
       .assign(value=lambda d: d["value"].where(d["unit"] != "degF", (d["value"] - 32) * 5 / 9))
       .dropna(subset=["value"]))
assert len(toy) == 7 and toy.query("sensor == 'temp_1'").groupby("machine")["value"].mean().round(2).to_dict() == \
    {"A": 71.55, "B": 68.0, "C": 71.64}

labels = pd.Series(np.random.default_rng(0).choice(["temp_1", "vib_1", "pres_1"], 100_000))   # 3.4
assert labels.memory_usage(deep=True) / labels.astype("category").memory_usage(deep=True) > 50
for start, end, per_day in (("2024-03-30", "2024-04-02", [24, 23, 24]), ("2024-10-26", "2024-10-29", [24, 25, 24])):
    local = pd.date_range(start, end, freq="h", tz="Europe/Berlin", inclusive="left")          # 3.11
    assert pd.Series(1, index=local).resample("D").size().tolist() == per_day

e = pd.DataFrame({"temp": [71.2, 69.9, 68.0, 72.4, 70.1]}, index=[10, 20, 30, 40, 50]).sort_values("temp")   # 4.1
s_int = pd.Series([1.0, 2.0, 3.0], index=[2, 0, 1])
assert e.index.tolist() == [30, 20, 50, 10, 40] and e.loc[20:40, "temp"].tolist() == [69.9, 70.1, 71.2, 72.4]
assert e.iloc[1:3]["temp"].tolist() == [69.9, 70.1] and s_int[0] == 2.0 and s_int[0:2].tolist() == [1.0, 2.0]

cap = pd.DataFrame({"temp": [71.2, 69.9, 72.4]})                                                   # 4.2
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    cap["temp"][cap["temp"] > 71] = 71.0          # the anti-pattern, kept only to check its effect
assert caught[0].category is pd.errors.ChainedAssignmentError and cap["temp"].tolist() == [71.2, 69.9, 72.4]
cap.loc[cap["temp"] > 71, "temp"] = 71.0
assert cap["temp"].tolist() == [71.0, 69.9, 71.0]

rng9 = np.random.default_rng(9)                                                                    # 3.8, 4.3-4.5
n9 = 200_000
log = pd.DataFrame({"machine": rng9.choice(["A", "B", "C"], n9, p=[0.5, 0.3, 0.2]),
                    "sensor": rng9.choice(["temp_1", "temp_2", "temp_3", "temp_4"], n9)})
log["value"] = (log["machine"].map({"A": 60.0, "B": 65.0, "C": 58.0})
                + log["sensor"].str[-1].astype(int) - 1 + rng9.normal(0, 2.0, n9))
summary = (log.assign(hot=lambda d: d["value"] > 64.0).groupby("machine")
           .agg(n=("value", "size"), mean=("value", "mean"), sd=("value", "std"),
                p95=("value", lambda v: v.quantile(0.95)), frac_hot=("hot", "mean")))
expected = {"n": [100372, 59743, 39885], "mean": [61.500, 66.498, 59.501], "sd": [2.285, 2.291, 2.290],
            "p95": [65.245, 70.279, 63.254], "frac_hot": [0.138, 0.859, 0.024]}
for col, vals in expected.items():
    assert np.allclose(summary[col], vals, atol=5e-4), col
pt9 = log.pivot_table(index="machine", columns="sensor", values="value", aggfunc="mean", margins=True)
assert round(pt9.loc["All", "All"], 3) == 62.595 and round(pt9["All"].iloc[:3].mean(), 3) == 62.500
assert np.isclose(np.sqrt(2.0**2 + 1.25), 2.291, atol=5e-4)

stamps = pd.date_range("2024-03-04 08:00", "2024-03-04 11:00", freq="10min")                       # 4.6
lost = (stamps >= "2024-03-04 09:00") & (stamps < "2024-03-04 10:30")
gappy = pd.Series(70 + 0.5 * np.arange(len(stamps)), index=stamps).loc[~lost]
half = gappy.resample("30min").mean()
assert half.isna().tolist() == [False, False, True, True, True, False, False] and gappy.asfreq("10min").isna().sum() == 9
assert half.interpolate(method="time").tolist() == [70.5, 72.0, 73.5, 75.0, 76.5, 78.0, 79.0]
assert half.interpolate(method="time", limit=1).isna().sum() == 2
print("notes checks: Build 7 toy table, category memory, DST day lengths, Examples 4.1-4.6 all reproduce")
print(f"Total run time: {time.perf_counter() - T_START:.1f} s")


# %% [markdown]
# ## Your turn
#
# 1. **A transform-based z-score.** Add a column z = (value - group mean) / group sd with
#    groupby(["machine", "sensor"])["value"].transform, and list the rows with |z| > 4. Does it find
#    the burst? Now group by machine, sensor *and calendar day* (pd.Grouper(key="time", freq="D"))
#    and repeat. Why are the per-day z-scores of the burst *smaller*? Compare both with Part 2's
#    Hampel flags (misses, false alarms) for burst sizes 30, 15 and 8 degC (Exercise B5).
# 2. **Missing data as a heatmap.** For every (machine, sensor, day), compute the fraction of the
#    expected readings that are missing (expected = 144, 96 or 72 per day by kind). Show it as an
#    imshow heatmap (rows: machine/sensor, columns: day) with a colorbar. The 6-hour gap on A
#    and the 1% random losses should both be visible; check the gap cells numerically.
# 3. **The daily summary in query form.** Rewrite Part 3's filter with clean.query(...) and compare
#    the result with `daily` using pd.testing.assert_frame_equal. Which version would you rather
#    debug, and why?
# 4. **iterrows against vectorized.** Time the Fahrenheit conversion of Part 2 on 10^5 rows of raw
#    three ways: a loop over iterrows(), apply(axis=1), and the where() expression. Report the
#    ratios, then double the rows and check that every method scales linearly.

# %%
if __name__ == "__main__":
    print("All checks passed.")
