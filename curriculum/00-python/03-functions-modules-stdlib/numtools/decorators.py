"""Decorators from PY-3 section 3.5-3.6: timer, memoize, retry(n)."""
from __future__ import annotations

import functools
import time
from collections.abc import Callable
from typing import Any


def timer(func: Callable[..., Any]) -> Callable[..., Any]:
    """Record the wall time of every call in ``wrapper.timings`` (seconds, via perf_counter)."""

    @functools.wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        start = time.perf_counter()
        try:
            return func(*args, **kwargs)
        finally:                                  # also record calls that raise
            wrapper.timings.append(time.perf_counter() - start)

    wrapper.timings = []
    return wrapper


def memoize(func: Callable[..., Any]) -> Callable[..., Any]:
    """Cache results in ``wrapper.cache``; every argument must be hashable."""
    cache: dict[Any, Any] = {}

    @functools.wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        key = (args, tuple(sorted(kwargs.items())))
        if key not in cache:
            cache[key] = func(*args, **kwargs)
        return cache[key]

    wrapper.cache = cache
    return wrapper


def retry(n: int, exceptions: tuple[type[BaseException], ...] = (Exception,), delay: float = 0.0):
    """Decorator factory: call the function up to ``n`` times in total, re-raising the last error."""
    if n < 1:
        raise ValueError("n must be at least 1")

    def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            for attempt in range(1, n + 1):
                try:
                    return func(*args, **kwargs)
                except exceptions:
                    if attempt == n:
                        raise
                    if delay:
                        time.sleep(delay)

        return wrapper

    return decorator
