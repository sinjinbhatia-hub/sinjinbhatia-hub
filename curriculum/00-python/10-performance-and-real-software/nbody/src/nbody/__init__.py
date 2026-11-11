"""nbody: a direct-summation gravitational N-body code, the PY-10 capstone.

Public API: the force functions in nbody.forces, the integrators in nbody.integrate, and the
command line in nbody.cli (`python -m nbody --help`).
"""
from nbody.forces import NUMBA_AVAILABLE, ForceWorkspace, accel_loop, accel_numpy
from nbody.integrate import energy, leapfrog, momentum, rk4, run_ensemble, simulate, uniform_ball

__version__ = "0.1.0"
__all__ = ["NUMBA_AVAILABLE", "ForceWorkspace", "accel_loop", "accel_numpy", "energy", "leapfrog",
           "momentum", "rk4", "run_ensemble", "simulate", "uniform_ball", "__version__"]
