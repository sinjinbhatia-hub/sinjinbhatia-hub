"""Command line for nbody: `python -m nbody --n 100 --steps 500 --seed 1` (PY-10 notes, Build 8)."""
from __future__ import annotations

import argparse
import dataclasses
import json
import logging
import platform
import subprocess
import time
import tomllib
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import nbody
from nbody.integrate import run_ensemble

log = logging.getLogger("nbody")


@dataclass
class RunConfig:
    """Everything that determines a run. Saved next to the results so the run can be repeated."""
    n: int = 100
    steps: int = 500
    dt: float = 1e-3
    eps: float = 0.05
    seed: int = 1
    ensemble: int = 1                  # number of initial conditions: seeds seed, seed+1, ...
    workers: int = 1
    dtype: str = "float64"

    @classmethod
    def from_toml(cls, path: Path) -> RunConfig:
        """Read a [run] table; unknown keys are an error, not a silent typo."""
        table = tomllib.loads(Path(path).read_text())["run"]
        unknown = set(table) - {f.name for f in dataclasses.fields(cls)}
        if unknown:
            raise ValueError(f"unknown config keys: {sorted(unknown)}")
        return cls(**table)


def git_commit() -> str:
    """The commit hash of the working tree this code was run from, or 'unknown'."""
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True,
                             cwd=Path(__file__).parent, timeout=5)
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"
    return out.stdout.strip() if out.returncode == 0 else "unknown"


def provenance() -> dict:
    """Versions and commit: what someone needs to reproduce a result file."""
    return {"nbody": nbody.__version__, "python": platform.python_version(), "numpy": np.__version__,
            "platform": platform.platform(), "git_commit": git_commit()}


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="nbody", description="Direct-summation N-body run with leapfrog.")
    p.add_argument("--config", type=Path, help="TOML file with a [run] table; flags override it")
    p.add_argument("--n", type=int, help="number of bodies (default 100)")
    p.add_argument("--steps", type=int, help="number of leapfrog steps (default 500)")
    p.add_argument("--dt", type=float, help="time step (default 1e-3)")
    p.add_argument("--eps", type=float, help="softening length (default 0.05)")
    p.add_argument("--seed", type=int, help="seed of the first initial condition (default 1)")
    p.add_argument("--ensemble", type=int, help="number of initial conditions (default 1)")
    p.add_argument("--workers", type=int, help="processes for the ensemble (default 1)")
    p.add_argument("--dtype", choices=["float64", "float32"], help="floating type (default float64)")
    p.add_argument("--out", type=Path, help="write results and provenance to this JSON file")
    p.add_argument("-v", "--verbose", action="store_true", help="log at DEBUG level")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")
    cfg = RunConfig.from_toml(args.config) if args.config else RunConfig()
    overrides = {f.name: getattr(args, f.name) for f in dataclasses.fields(RunConfig)
                 if getattr(args, f.name) is not None}
    cfg = dataclasses.replace(cfg, **overrides)
    log.debug("config %s", cfg)
    seeds = list(range(cfg.seed, cfg.seed + cfg.ensemble))
    t0 = time.perf_counter()
    results = run_ensemble(seeds, workers=cfg.workers, n=cfg.n, steps=cfg.steps, dt=cfg.dt,
                           eps=cfg.eps, dtype=cfg.dtype)
    wall = time.perf_counter() - t0
    for r in results:
        log.info("seed %d: |dE/E| = %.2e, |dP| = %.1e, half-mass radius %.4f",
                 r["seed"], r["rel_energy_change"], r["momentum_change"], r["half_mass_radius"])
    log.info("%d run(s), N = %d, %d steps, %d worker(s): %.2f s", len(results), cfg.n, cfg.steps,
             cfg.workers, wall)
    if args.out:
        record = {"config": dataclasses.asdict(cfg), "provenance": provenance(),
                  "wall_seconds": wall, "results": results}
        args.out.write_text(json.dumps(record, indent=2))
        log.info("wrote %s", args.out)
    return 0
