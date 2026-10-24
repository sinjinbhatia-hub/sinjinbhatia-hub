"""Function builders from PY-3 section 3.2: compose and pipeline."""
from __future__ import annotations

from collections.abc import Callable
from typing import Any


def compose(f: Callable[..., Any], g: Callable[..., Any]) -> Callable[..., Any]:
    """Return h with h(*a, **k) == f(g(*a, **k)): g runs first."""

    def h(*args: Any, **kwargs: Any) -> Any:
        return f(g(*args, **kwargs))

    return h


def pipeline(*funcs: Callable[..., Any]) -> Callable[..., Any]:
    """Return a function that applies funcs left to right; pipeline() is the identity."""

    def run(*args: Any, **kwargs: Any) -> Any:
        if not funcs:
            return args[0]
        value = funcs[0](*args, **kwargs)
        for func in funcs[1:]:
            value = func(value)
        return value

    return run
