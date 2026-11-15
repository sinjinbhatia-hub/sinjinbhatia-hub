"""The first version of the program, kept as a reference: pure Python, lists of [x, y, z] lists.

This is what notes Build 1 times and profiles. Everything here is O(N) per step except
accel_loop and energy_loop, which are O(N^2).
"""
from __future__ import annotations

import math

from nbody.forces import accel_loop


def kick(vel: list[list[float]], acc: list[list[float]], h: float) -> None:
    """v <- v + h a, in place."""
    for v, a in zip(vel, acc):
        v[0] += h * a[0]
        v[1] += h * a[1]
        v[2] += h * a[2]


def drift(pos: list[list[float]], vel: list[list[float]], h: float) -> None:
    """x <- x + h v, in place."""
    for p, v in zip(pos, vel):
        p[0] += h * v[0]
        p[1] += h * v[1]
        p[2] += h * v[2]


def energy_loop(pos: list[list[float]], vel: list[list[float]], mass: list[float], eps: float) -> float:
    """K + U with U = -sum_{i<j} m_i m_j / sqrt(r_ij^2 + eps^2), by a double loop."""
    e = sum(0.5 * m * (v[0] ** 2 + v[1] ** 2 + v[2] ** 2) for m, v in zip(mass, vel))
    n = len(pos)
    for i in range(n):
        for j in range(i + 1, n):
            dx, dy, dz = (pos[j][k] - pos[i][k] for k in range(3))
            e -= mass[i] * mass[j] / math.sqrt(dx * dx + dy * dy + dz * dz + eps * eps)
    return e


def run_python(pos: list[list[float]], vel: list[list[float]], mass: list[float], dt: float,
               steps: int, eps: float) -> tuple[float, float]:
    """Leapfrog for `steps` steps, mutating pos and vel; returns (energy before, energy after)."""
    e0 = energy_loop(pos, vel, mass, eps)
    acc = accel_loop(pos, mass, eps)
    for _ in range(steps):
        kick(vel, acc, dt / 2)
        drift(pos, vel, dt)
        acc = accel_loop(pos, mass, eps)
        kick(vel, acc, dt / 2)
    return e0, energy_loop(pos, vel, mass, eps)
