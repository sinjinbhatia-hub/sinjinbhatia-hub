"""Tests for nbody.integrate: conservation laws and reproducibility."""
import numpy as np

from nbody.integrate import energy, leapfrog, momentum, simulate, uniform_ball


def test_uniform_ball_is_centred_and_virialized() -> None:
    pos, vel, mass = uniform_ball(500, np.random.default_rng(0))
    assert np.abs(mass @ pos).max() < 1e-15 and np.abs(mass @ vel).max() < 1e-15
    kinetic = 0.5 * mass @ (vel**2).sum(axis=1)
    assert np.isclose(kinetic, 0.5 * 0.5 * 0.6)


def test_leapfrog_conserves_momentum_and_energy() -> None:
    pos, vel, mass = uniform_ball(64, np.random.default_rng(1))
    p0, e0 = momentum(vel, mass), energy(pos, vel, mass, 0.05)
    pos, vel, _ = leapfrog(pos, vel, mass, 1e-3, 300, eps=0.05)
    assert np.abs(momentum(vel, mass) - p0).max() < 1e-10
    assert abs(energy(pos, vel, mass, 0.05) - e0) / abs(e0) < 1e-5


def test_same_seed_same_result() -> None:
    assert simulate(7, n=20, steps=10) == simulate(7, n=20, steps=10)
