# PY-1 Exercises

Conventions as in the notes: Python 3.11+, standard library only. **Write your prediction down before you run anything.** Every solution in `solutions.md` was run and passes its own asserts. In the Predict exercises the interesting part is the reason, so write one sentence per answer.

## A. Warm-up (do all; 5-10 minutes each)

**A1. Aliasing.** Predict what each `print` shows, citing R1-R4 of the notes.

```python
x = [1, 2, 3]
y = x
y[0] = 99
z = x[:]
z[1] = 77
print(x, y, z)

def add_one(lst):
    lst.append(1)
    return lst

a = [0]
b = add_one(a)
print(a is b, a)

s = "abc"
t = s
t += "d"
print(s, t)

d = {"k": [1]}
e = dict(d)
e["k"].append(2)
e["j"] = 0
print(d, e)
```

**A2. Floor rule.** Predict `-17 // 5`, `-17 % 5`, `17 // -5`, `17 % -5`, `divmod(-17, 5)`, `2 ** -2`, `7 // 2.0`, `-7.5 % 2`, `-3 ** 2`, `(-3) ** 2`. Then write `hour_after(h, delta)` in one line, returning the clock hour in `range(24)` that is `delta` hours (any integer, possibly negative) after hour `h`. It must satisfy `hour_after(22, 5) == 3`, `hour_after(2, -5) == 21` and `hour_after(0, -48) == 0`. Which property of `%` makes the one-liner correct?

**A3. The grid.** Predict each result, and for the `False` ones say which spacing or tie caused it.

```python
print(0.1 + 0.1 + 0.1 == 0.3)
print(1e16 + 1 == 1e16)
print((1e16 + 1) + 1 == 1e16, 1e16 + (1 + 1) == 1e16)
print(0.5 + 0.25 == 0.75)
print(float("inf") == float("inf") + 1)
print(2**53 + 1 == float(2**53 + 1))
```

**A4. Scope.** Predict what this does, then fix it in two ways, one of which uses `global` and one of which does not. Say which you would ship.

```python
total = 0

def add(values):
    for v in values:
        total += v
    return total

print(add([1, 2, 3]))
```

Also predict the output of `names = ["a", "b"]`, then `def rename(): names.append("c")`, then `rename(); print(names)`, and explain why it needs no `global`.

**A5. The bits of a double.** Write `fields(x)` returning `(sign, E, F)` (notes eq. 3.2). Hint: `int.from_bytes(struct.pack(">d", x), "big")` is the 64-bit pattern as an integer; shift and mask it. Predict, then print, `E` for `1.0`, `0.5`, `8.0` and `0.1`, and the fields of `-0.0`. What does `E - 1023` mean for `8.0`?

**A6. Truthiness and formatting.** Predict the output.

```python
print(bool([]), bool([0]), bool(0.0), bool("0"), bool(None))
print(0 or "default", "" or 0 or [], [1] and "yes", None and 1 / 0)
print(f"{3.14159:8.3f}|{42:<6d}|{42:>6d}|{0.00012345:.2e}|{1234567:,}|{0.5:.0%}")
```

Then find the bug in `def connect(timeout=None): timeout = timeout or 30`, say what call exposes it, and fix it.

## B. Core (do all)

**B1. Aliasing audit.** Write `describe_sharing(a, b)` for two lists, returning a dict with `same_object` (bool), `equal` (bool) and `shared` (the indices `i` where `a[i] is b[i]`, for equal-length lists). Use it to audit `rows = [[0] * 3] * 3` against `[[0] * 3 for _ in range(3)]`, and `list(m)`, `m[:]` and `copy.deepcopy(m)` of `m = [[1, 2], [3, 4]]`. Report which construction shares what, and give the rule for when each copy is safe.

**B2. `isclose` from scratch.** Implement `my_isclose(a, b, rel_tol=1e-9, abs_tol=0.0)` from equation (3.5), with correct handling of equal values, infinities and `nan` (what should `my_isclose(inf, inf)` and `my_isclose(inf, 1e308)` return?). Test your function against `math.isclose` on 100,000 seeded random pairs, half of them constructed so the relative difference is within a factor of three of `rel_tol`, and on the special cases. Then show a pair where `my_isclose(a, b)` is true and `my_isclose(a, b, rel_tol=0)` is false, and explain what `abs_tol` should be if the values are voltages in volts around `1e-6`.

**B3. The secant method.** Newton needs `df`. Replace `f'(x_n)` by the slope through the last two points to get $x_{n+1} = x_n - f(x_n)\,\dfrac{x_n - x_{n-1}}{f(x_n) - f(x_{n-1})}$. Implement `secant(f, x0, x1, tol, max_iter)`, run it with `bisect` and `newton` on $\cos x - x$ and $x^3 - 2x - 5$, and tabulate iteration counts. From the error sequence of $x^2 - 2$ estimate the order $p$ with $p \approx \log(e_{n+1}/e_n)/\log(e_n/e_{n-1})$. Which should you use if `f` costs a minute per evaluation, and why?

**B4. Remainders in real use.** Write `to_clock(seconds)` returning `(days, hours, minutes, secs)` using `divmod`, such that `0 <= hours < 24`, `0 <= minutes < 60`, `0 <= secs < 60` for *any* integer input, with `days` negative when `seconds` is. Check it on 10,000 seeded random integers in $[-10^6, 10^6]$ against `datetime.timedelta(seconds=s)` (whose `.days` and `.seconds` are normalised the same way). Explain which part of the notes' floor rule makes this work, and what a truncating `%` would have returned for `-1`.

**B5. A series with a stopping rule.** (a) With a `while` loop, sum $e = \sum_{k\ge0} 1/k!$ until a term is below `1e-12`; report the number of terms and the error. (b) Do the same for $\sum_{k\ge1} 1/k^2 = \pi^2/6$ with stopping term `1e-8`; report the number of terms and the error, and explain why the error is so much larger than `1e-8`. (c) Give a stopping rule that guarantees an error below `1e-6` for part (b), using the bound $\sum_{k>N} 1/k^2 < 1/N$, and check it. (d) Sum the series for $e^{-20}$ the same way as (a) and compare with `math.exp(-20)`; find a rearrangement that is accurate.

**B6. A report with computed column widths.** Given `data = [("copper", 8.96, 0.02), ("iron", 7.874, 0.004), ("aluminium", 2.70, 0.01), ("lead", 11.34, 0.05)]` (name, density in g/cm³, uncertainty), produce an aligned table whose name column is as wide as the longest name (hint: a format spec may contain a nested field, `f"{name:<{w}}"`), with density and uncertainty right-aligned with three decimals, a header, a rule line, and a final row for the mean density. Assert the exact text.

**B7. Extend `project.py`.** Add a cell after Part 5 that measures the error of `naive_sum` on `[0.1] * n` against `math.fsum` for $n = 10^3, 10^4, \ldots, 10^7$, prints a table, and estimates by least squares (written by hand, no libraries) the slope of $\log_{10}$ error against $\log_{10} n$. Assert what you measured, with a margin, and explain the slope using the naive bound $n\,u\sum|x_i|$ of the notes (note that $\sum|x_i|$ itself grows with $n$).

## C. Challenge (do at least one)

**C1. Walking the grid.** Implement `next_up(x)`, the next double above `x`, using only integer arithmetic on the bit pattern of `x` (the idea: for positive doubles, adding 1 to the 64-bit integer gives the next double). Handle `-0.0`, negative numbers, subnormals, the largest double, `inf` and `nan`. Verify against `math.nextafter(x, math.inf)` on 100,000 random bit patterns plus the special values, and explain in a paragraph why the integer trick works (hint: look at where the exponent and fraction fields sit, and what happens when $F$ overflows).

**C2. Design: `x or default`.** The function `configure(timeout, retries, name)` in the solutions has the three defaults `30`, `3` and `"anon"` applied with `x or default`. (a) For each parameter say which explicit values `0`, `""` or `[]` a caller might legitimately pass and what the idiom does to them. (b) Rewrite the function so that only a missing argument receives a default. (c) Write 100-150 words arguing for a rule about when the `or` idiom is acceptable in code you would review, with one case where you would allow it.

**C3. Does summation order matter?** Sum $\sum_{k=1}^{10^6} 1/k^2$ with six methods: forward loop, backward loop, `sorted(terms)` ascending, Kahan, a recursive pairwise sum (split in halves, base case of 8 terms), and `math.fsum` (the reference). Print the errors against `fsum` in a table, rank the methods, and explain the ranking using the grid-spacing picture. Then repeat for the alternating series $\sum_{k=1}^{10^6} (-1)^{k+1}/k$ and say what changes.
