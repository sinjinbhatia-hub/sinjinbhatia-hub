"""Initial conditions, conserved quantities and two integrators for eq. (3.1) (PY-10 notes, section 3).

Units: G = 1, total mass 1. State: positions and velocities as (N, 3) float arrays.
"""
from __future__ import annotations

from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from multiprocessing.context import BaseContext

import numpy as np

from nbody.forces import ForceWorkspace, accel_numpy

Accel = Callable[[np.ndarray, np.ndarray, float], np.ndarray]


def uniform_ball(n: int, rng: np.random.Generator, virial: float = 0.5,
                 dtype: type = np.float64) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """n equal masses uniform in the unit ball, Gaussian velocities scaled so 2K/|U| = virial.

    The centre of mass is put at rest at the origin. Uses the unsoftened potential of a uniform
    ball, U = -3 G M^2 / (5 R), to set the velocity scale. Returns (pos, vel, mass).
    """
    u = rng.normal(size=(n, 3))
    pos = u / np.linalg.norm(u, axis=1, keepdims=True) * rng.random(n)[:, np.newaxis] ** (1 / 3)
    vel = rng.normal(size=(n, 3))
    mass = np.full(n, 1.0 / n)
    pos -= mass @ pos
    vel -= mass @ vel
    kinetic = 0.5 * mass @ (vel * vel).sum(axis=1)
    vel *= np.sqrt(virial * 0.5 * (3 / 5) / kinetic)       # K = virial * |U| / 2
    return pos.astype(dtype), vel.astype(dtype), mass.astype(dtype)


def energy(pos: np.ndarray, vel: np.ndarray, mass: np.ndarray, eps: float = 0.0) -> float:
    """Total energy K + U, U = -G sum_{i<j} m_i m_j / sqrt(r_ij^2 + eps^2), always in float64."""
    pos, vel, mass = (np.asarray(a, dtype=np.float64) for a in (pos, vel, mass))
    kinetic = 0.5 * mass @ (vel * vel).sum(axis=1)
    d = pos[np.newaxis, :, :] - pos[:, np.newaxis, :]
    r = np.sqrt((d * d).sum(axis=-1) + eps**2)
    iu = np.triu_indices(len(mass), k=1)
    potential = -np.sum(mass[iu[0]] * mass[iu[1]] / r[iu])
    return float(kinetic + potential)


def momentum(vel: np.ndarray, mass: np.ndarray) -> np.ndarray:
    """Total momentum P = sum_i m_i v_i, shape (3,)."""
    return np.asarray(mass, np.float64) @ np.asarray(vel, np.float64)


def leapfrog(pos: np.ndarray, vel: np.ndarray, mass: np.ndarray, dt: float, steps: int,
             eps: float = 0.0, accel: Accel = accel_numpy,
             record_every: int = 0) -> tuple[np.ndarray, np.ndarray, list[float]]:
    """Kick-drift-kick leapfrog (velocity Verlet), a symplectic map (CM-8 section 3.11).

    v += a dt/2; x += v dt; v += a(x) dt/2, with one force evaluation per step. Returns the final
    (pos, vel) as new arrays and, if record_every > 0, the energy every record_every steps.
    """
    pos, vel = pos.copy(), vel.copy()
    half = pos.dtype.type(0.5 * dt)
    step = pos.dtype.type(dt)
    energies = [energy(pos, vel, mass, eps)] if record_every else []
    acc = accel(pos, mass, eps)
    for s in range(1, steps + 1):
        vel += half * acc
        pos += step * vel
        acc = accel(pos, mass, eps)
        vel += half * acc
        if record_every and s % record_every == 0:
            energies.append(energy(pos, vel, mass, eps))
    return pos, vel, energies


def rk4(pos: np.ndarray, vel: np.ndarray, mass: np.ndarray, dt: float, steps: int,
        eps: float = 0.0, accel: Accel = accel_numpy,
        record_every: int = 0) -> tuple[np.ndarray, np.ndarray, list[float]]:
    """Classical fourth-order Runge-Kutta on (x, v): four force evaluations per step, not symplectic."""
    pos, vel = pos.copy(), vel.copy()
    energies = [energy(pos, vel, mass, eps)] if record_every else []
    for s in range(1, steps + 1):
        k1x, k1v = vel, accel(pos, mass, eps)
        k2x, k2v = vel + 0.5 * dt * k1v, accel(pos + 0.5 * dt * k1x, mass, eps)
        k3x, k3v = vel + 0.5 * dt * k2v, accel(pos + 0.5 * dt * k2x, mass, eps)
        k4x, k4v = vel + dt * k3v, accel(pos + dt * k3x, mass, eps)
        pos = pos + dt / 6 * (k1x + 2 * k2x + 2 * k3x + k4x)
        vel = vel + dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)
        if record_every and s % record_every == 0:
            energies.append(energy(pos, vel, mass, eps))
    return pos, vel, energies


def simulate(seed: int, n: int = 200, steps: int = 200, dt: float = 1e-3, eps: float = 0.05,
             dtype: str = "float64") -> dict:
    """One ensemble member: build uniform_ball(n) from `seed`, run leapfrog, summarize.

    A module-level function on purpose: ProcessPoolExecutor pickles it by name.
    """
    rng = np.random.default_rng(seed)
    pos, vel, mass = uniform_ball(n, rng, dtype=np.dtype(dtype).type)
    e0 = energy(pos, vel, mass, eps)
    p0 = momentum(vel, mass)
    work = ForceWorkspace(n, dtype=np.dtype(dtype).type)
    pos, vel, _ = leapfrog(pos, vel, mass, dt, steps, eps, accel=work)
    return {"seed": seed, "n": n, "steps": steps, "dtype": dtype,
            "energy_start": e0, "energy_end": energy(pos, vel, mass, eps),
            "rel_energy_change": abs(energy(pos, vel, mass, eps) - e0) / abs(e0),
            "momentum_change": float(np.abs(momentum(vel, mass) - p0).max()),
            "half_mass_radius": float(np.median(np.linalg.norm(np.asarray(pos, np.float64), axis=1)))}


def run_ensemble(seeds: list[int], workers: int = 1, chunksize: int = 1,
                 mp_context: BaseContext | None = None, **kwargs) -> list[dict]:
    """simulate(seed, **kwargs) for every seed, serially (workers = 1) or on a pool of processes.

    Call it under `if __name__ == "__main__":` (notes section 3.6): with the "spawn" and
    "forkserver" start methods every worker re-imports the main module.
    """
    job = partial(simulate, **kwargs)                      # a partial of a module-level function pickles
    if workers <= 1:
        return [job(s) for s in seeds]
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp_context) as pool:
        return list(pool.map(job, seeds, chunksize=chunksize))
