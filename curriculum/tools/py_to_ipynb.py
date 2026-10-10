"""Convert every `project.py` under the curriculum into a Jupyter notebook.

Usage:  python3 curriculum/tools/py_to_ipynb.py [root]

Cells are delimited by `# %%` (code) and `# %% [markdown]` (markdown, with the
leading `# ` stripped from each line). The notebook is written next to the
script as `project.ipynb`. Notebooks are written without outputs; run them
yourself.
"""
from __future__ import annotations
import sys
from pathlib import Path
import nbformat
from nbformat.v4 import new_notebook, new_code_cell, new_markdown_cell


def split_cells(src: str):
    cells, kind, buf = [], "code", []

    def flush():
        text = "\n".join(buf).strip("\n")
        if text.strip():
            cells.append((kind, text))

    for line in src.splitlines():
        stripped = line.strip()
        if stripped.startswith("# %%"):
            flush()
            buf = []
            kind = "markdown" if "[markdown]" in stripped else "code"
            continue
        if kind == "markdown":
            if line.startswith("# "):
                line = line[2:]
            elif line.strip() == "#":
                line = ""
        buf.append(line)
    flush()
    return cells


def convert(py: Path) -> Path:
    cells = split_cells(py.read_text(encoding="utf-8"))
    nb = new_notebook()
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    for kind, text in cells:
        nb.cells.append(new_markdown_cell(text) if kind == "markdown" else new_code_cell(text))
    out = py.with_suffix(".ipynb")
    nbformat.write(nb, out)
    return out


if __name__ == "__main__":
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
    n = 0
    for py in sorted(root.rglob("project.py")):
        convert(py)
        n += 1
    print(f"converted {n} scripts")
