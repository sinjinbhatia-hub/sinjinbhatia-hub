"""Command line for numtools: ``python -m numtools {bisect,newton} EXPR ...``."""
from __future__ import annotations

import argparse
import math
import sys
from collections.abc import Callable, Sequence

from . import __version__
from .roots import bisect, newton

NAMES = {name: getattr(math, name) for name in
         ("sin", "cos", "tan", "exp", "log", "sqrt", "atan", "sinh", "cosh", "tanh", "pi", "e")}


def make_function(expr: str) -> Callable[[float], float]:
    """Compile ``expr`` once and return f(x).  eval with empty builtins is NOT a security sandbox."""
    code = compile(expr, "<expr>", "eval")
    env = {"__builtins__": {}, **NAMES}

    def f(x: float) -> float:
        return eval(code, env, {"x": x})

    return f


def central_difference(f: Callable[[float], float]) -> Callable[[float], float]:
    """Numerical f' with step 1e-6 * max(1, |x|); good to about 1e-10 for smooth f."""

    def df(x: float) -> float:
        h = 1e-6 * max(1.0, abs(x))
        return (f(x + h) - f(x - h)) / (2 * h)

    return df


def run_bisect(args: argparse.Namespace) -> tuple[float, int]:
    return bisect(make_function(args.expr), args.a, args.b, tol=args.tol)


def run_newton(args: argparse.Namespace) -> tuple[float, int]:
    f = make_function(args.expr)
    return newton(f, central_difference(f), args.x0, tol=args.tol)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m numtools", description="Find a root of EXPR (an expression in x) with the PY-1 solvers.")
    parser.add_argument("--version", action="version", version=f"numtools {__version__}")
    sub = parser.add_subparsers(dest="command", required=True, metavar="{bisect,newton}")

    def common(p: argparse.ArgumentParser) -> None:
        p.add_argument("expr", help='expression in x, for example "cos(x)-x"')
        p.add_argument("--tol", type=float, default=1e-12, help="stopping tolerance (default: %(default)s)")
        p.add_argument("-v", "--verbose", action="store_true", help="also print the iteration count and f(root)")

    pb = sub.add_parser("bisect", help="bisection on a bracket [a, b]")
    common(pb)
    pb.add_argument("a", type=float, help="left end of the bracket")
    pb.add_argument("b", type=float, help="right end of the bracket")
    pb.set_defaults(func=run_bisect)

    pn = sub.add_parser("newton", help="Newton's method from x0 (derivative by central differences)")
    common(pn)
    pn.add_argument("x0", type=float, help="starting point")
    pn.set_defaults(func=run_newton)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CLI; returns the exit status.  ``argv`` defaults to ``sys.argv[1:]``."""
    args = build_parser().parse_args(argv)
    try:
        root, iterations = args.func(args)
    except (ValueError, RuntimeError, ZeroDivisionError, SyntaxError, NameError) as err:
        print(f"error: {err}", file=sys.stderr)
        return 2
    print(f"root = {root:.12g}")
    if args.verbose:
        print(f"iterations = {iterations}")
        print(f"f(root) = {make_function(args.expr)(root):.3e}")
    return 0
