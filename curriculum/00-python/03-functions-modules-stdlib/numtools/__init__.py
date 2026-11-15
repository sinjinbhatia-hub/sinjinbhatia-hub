"""numtools: the PY-1 root finders, the PY-2 sparse polynomial, and PY-3 function tools."""

__version__ = "0.1"

from .decorators import memoize, retry, timer
from .functional import compose, pipeline
from .polynomials import poly_add, poly_deriv, poly_eval, poly_mul
from .roots import bisect, newton, secant

__all__ = ["bisect", "newton", "secant", "poly_add", "poly_mul", "poly_eval", "poly_deriv",
           "timer", "memoize", "retry", "compose", "pipeline"]
