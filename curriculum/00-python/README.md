# 00 · Python

Ten modules, about 87 hours, taking you from "I can make Python do things" to "I can predict what Python will do, write code I trust, and make it fast." Every other course here ends each module with a Python project; this course is where the Python itself gets learned properly, and it goes first.

| # | Module | Hours | Needs | Take it before |
|---|---|---|---|---|
| PY-1 | [Objects, names, numbers, and control flow](01-objects-names-and-control-flow/notes.md) | 8 | — | anything else; week 1 alongside F1 |
| PY-2 | [Collections, iteration, and the data model](02-collections-iteration-data-model/notes.md) | 9 | PY-1 | LA-2's project |
| PY-3 | [Functions as values, closures and decorators, recursion, modules](03-functions-modules-stdlib/notes.md) | 8 | PY-2 | LA-3 |
| PY-4 | [Classes and the object model; dual-number automatic differentiation](04-classes-and-the-object-model/notes.md) | 10 | PY-3 | LA-4 |
| PY-5 | [Errors, testing, debugging, type hints](05-errors-testing-debugging/notes.md) | 8 | PY-4 | LA-5 |
| PY-6 | [NumPy: arrays, memory, broadcasting, views vs copies](06-numpy/notes.md) | 10 | PY-5 | LA-5 at the latest |
| PY-7 | [Matplotlib and scientific communication](07-matplotlib-and-communication/notes.md) | 7 | PY-6 | CM-1 |
| PY-8 | [SciPy and SymPy: numerical and symbolic toolkits](08-scipy-and-sympy/notes.md) | 9 | PY-7 | CM-3, ST-1 |
| PY-9 | [pandas and working with data](09-pandas-and-data/notes.md) | 8 | PY-6, PY-7 | ST-1 |
| PY-10 | [Performance, the runtime, and writing real scientific software](10-performance-and-real-software/notes.md) | 10 | PY-8 | ML-7, SM-6 |

## What "really good at Python" means here

Three things, in order. First, you can **predict** what a piece of code does before running it: names vs values, mutability, scope, views vs copies, broadcasting. That is PY-1 to PY-6, and it is the part most self-taught programmers skip. Second, you write code you can **trust**: tests, exceptions that say what went wrong, a debugger instead of print statements, type hints. That is PY-5 and it is used in every module after. Third, you can make it **fast and real**: profile, vectorize, parallelize, package. That is PY-10.

Each module has the same structure as the math modules, with one change: where a math module derives a theorem, a Python module makes you build the thing yourself first (your own `enumerate`, your own decorator, your own matrix multiply, your own dual numbers for automatic differentiation) and only then hands you the library version. Every ▶ Try it first asks you to predict output before you run.

## Order

PY-1 through PY-5 in the first three weeks, in parallel with F1 and LA-1/LA-2 (they are light on math and heavy on practice; two tracks, as always). PY-6 and PY-7 before any physics project. PY-8 before CM-3. PY-9 before ST-1. PY-10 whenever your projects start feeling slow, and certainly before ML-7.

## Honest difficulty

None of this is conceptually hard the way LA-8 or QM-4 is. It is hard the way a sport is: the knowledge is small and the fluency comes only from volume. Do every A exercise. Type the code; do not paste it. When a ▶ Try it first asks you to predict, write the prediction down before running — the gap between what you predicted and what happened is the entire lesson.
