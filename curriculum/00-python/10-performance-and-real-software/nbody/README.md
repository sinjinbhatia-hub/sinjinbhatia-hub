# nbody

Direct-summation gravitational N-body simulation with a leapfrog integrator, written as the
capstone of PY-10. Units: G = 1, total mass 1; Plummer softening `eps`.

## Install

    python -m venv .venv && source .venv/bin/activate
    pip install -e ".[test]"          # add ",fast" for the optional Numba kernel

## Run

    python -m nbody --n 100 --steps 500 --seed 1
    nbody --n 200 --steps 300 --ensemble 8 --workers 4 --out run.json
    nbody --help

`--out` writes the configuration, the results, and the versions and git commit that produced them.

## Test

    pytest -q

## Layout

    src/nbody/forces.py      accelerations: loop, NumPy, in-place workspace, Numba, k-d tree
    src/nbody/integrate.py   initial conditions, energy, momentum, leapfrog, RK4, ensembles
    src/nbody/cli.py         argparse front end, TOML config, provenance
    tests/                   pytest suite

## Reproduce a run

Every `--out` file records the configuration, the package and library versions, and the git
commit. To repeat a run: `git checkout <commit>`, install, and pass the recorded configuration.
