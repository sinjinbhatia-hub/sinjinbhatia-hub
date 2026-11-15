"""Tests for nbody.forces: every fast version agrees with the loop, and the physics holds."""
import numpy as np
import pytest

from nbody.forces import (NUMBA_AVAILABLE, ForceWorkspace, accel_broadcast, accel_loop, accel_numpy,
                          accel_short_range)


def cluster(n: int, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    return rng.normal(size=(n, 3)), rng.uniform(0.5, 1.5, n) / n


def rel_err(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))


@pytest.mark.parametrize("n", [2, 7, 50])
@pytest.mark.parametrize("eps", [0.0, 0.05])
def test_numpy_versions_match_loop(n: int, eps: float) -> None:
    pos, mass = cluster(n)
    ref = np.array(accel_loop(pos.tolist(), mass.tolist(), eps))
    assert rel_err(accel_broadcast(pos, mass, eps), ref) < 1e-12
    assert rel_err(accel_numpy(pos, mass, eps), ref) < 1e-12
    assert rel_err(ForceWorkspace(n)(pos, mass, eps), ref) < 1e-12


def test_two_bodies_inverse_square() -> None:
    pos = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    mass = np.array([3.0, 0.5])
    acc = accel_numpy(pos, mass)
    assert np.allclose(acc, [[0.5 / 4, 0, 0], [-3.0 / 4, 0, 0]], rtol=1e-15, atol=0)


def test_newtons_third_law() -> None:
    pos, mass = cluster(100, seed=3)
    total_force = mass @ accel_numpy(pos, mass, 0.05)
    assert np.abs(total_force).max() < 1e-14


def test_workspace_result_is_overwritten() -> None:
    """The documented aliasing: keep a copy if you need the old value."""
    pos, mass = cluster(10)
    work = ForceWorkspace(10)
    first = work(pos, mass)
    kept = first.copy()
    work(pos + 1.0, mass * 2)
    assert np.shares_memory(first, work.acc)
    assert not np.allclose(first, kept)


def test_short_range_tree_matches_brute() -> None:
    rng = np.random.default_rng(5)
    pos, mass = rng.uniform(0, 10, size=(400, 3)), np.full(400, 1 / 400)
    tree = accel_short_range(pos, mass, rc=1.0, eps=0.01, method="tree")
    brute = accel_short_range(pos, mass, rc=1.0, eps=0.01, method="brute")
    assert rel_err(tree, brute) < 1e-12


@pytest.mark.skipif(not NUMBA_AVAILABLE, reason="numba not installed")
def test_numba_matches_loop() -> None:
    from nbody.forces import accel_numba
    pos, mass = cluster(30)
    ref = np.array(accel_loop(pos.tolist(), mass.tolist(), 0.05))
    assert rel_err(accel_numba(pos, mass, 0.05), ref) < 1e-12
