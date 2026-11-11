"""`python -m nbody ...` runs the command line. The guard matters: worker processes import this module."""
from nbody.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
