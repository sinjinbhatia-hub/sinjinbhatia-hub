"""plot_helpers.py: the house style and the helpers used by every figure in PY-7.

Built step by step in the notes (section 3.7, Build 5). Copy this file next to the
project of any later module; figures then land in a ``figures/`` folder beside it.

Public API
----------
HOUSE_STYLE, OKABE_ITO, use_house_style()
    The course rcParams and the colour-blind-safe cycle (Okabe and Ito).
savefig(fig, name, ...)
    Save under ``figures/``, creating the folder; returns the path.
loglog_fit(x, y)
    Least-squares line in log10-log10 space: slope, intercept, slope standard error.
add_slope_guide(ax, order, x0, y0, decades)
    A slope triangle on log-log axes, labelled with the order.
plot_convergence(ax, h, err, order, ...)
    Data, fitted line, reference-slope triangle and slope annotation in one call.
pi_ticks(ax, axis, base)
    Tick labels as multiples of pi on an axis.
label_end(ax, line, text)
    Label a line directly at its right end instead of with a legend.
rgba_at(ax, x, y)
    The rendered colour at a data point (for the one pixel test, section 4.5).
"""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cycler
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MultipleLocator, NullFormatter

__all__ = [
    "FIGDIR",
    "HOUSE_STYLE",
    "OKABE_ITO",
    "add_slope_guide",
    "label_end",
    "loglog_fit",
    "pi_ticks",
    "plot_convergence",
    "rgba_at",
    "savefig",
    "use_house_style",
]

#: Okabe and Ito's eight colours, ordered so the first four are the most distinct.
OKABE_ITO: list[str] = [
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#009E73",  # bluish green
    "#CC79A7",  # reddish purple
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#F0E442",  # yellow
    "#000000",  # black
]

#: The course style. Apply with ``use_house_style()`` or scope it with ``plt.rc_context(HOUSE_STYLE)``.
HOUSE_STYLE: dict[str, object] = {
    "figure.figsize": (6.4, 4.0),
    "figure.constrained_layout.use": True,
    "savefig.dpi": 150,
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.6,
    "axes.prop_cycle": cycler(color=OKABE_ITO),
    "lines.linewidth": 1.6,
    "lines.markersize": 5,
    "legend.frameon": False,
    "image.cmap": "viridis",
    "image.interpolation": "nearest",
}

FIGDIR = Path(__file__).resolve().parent / "figures"


def use_house_style() -> None:
    """Apply ``HOUSE_STYLE`` to the global rcParams (use ``plt.rc_context`` to scope it)."""
    plt.rcParams.update(HOUSE_STYLE)


def savefig(
    fig: Figure,
    name: str,
    *,
    outdir: Path | str | None = None,
    fmt: str = "png",
    dpi: int = 150,
    transparent: bool = False,
    tight: bool = True,
    close: bool = True,
) -> Path:
    """Save ``fig`` as ``<outdir>/<name>.<fmt>`` (default ``figures/`` beside this file).

    Creates the folder, trims the margins (``bbox_inches="tight"``; ``tight=False`` keeps the
    exact canvas), closes the figure unless ``close=False``, and returns the path.
    """
    folder = Path(outdir) if outdir is not None else FIGDIR
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"{name}.{fmt}"
    fig.savefig(path, dpi=dpi, bbox_inches="tight" if tight else None, transparent=transparent)
    if close:
        plt.close(fig)
    return path


def loglog_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    r"""Fit :math:`\log_{10} y = a \log_{10} x + b` by least squares.

    Returns ``(a, b, se_a)``: the slope (the order of convergence when ``y`` is an error
    and ``x`` a step size), the intercept (:math:`y \approx 10^b x^a`), and the standard
    error of the slope from the residual scatter.
    """
    lx, ly = np.log10(np.asarray(x, float)), np.log10(np.asarray(y, float))
    (a, b), cov = np.polyfit(lx, ly, 1, cov=True)
    return float(a), float(b), float(np.sqrt(cov[0, 0]))


def add_slope_guide(
    ax: Axes,
    order: float,
    x0: float,
    y0: float,
    decades: float = 0.5,
    *,
    color: str = "0.3",
) -> None:
    """Draw a slope triangle whose lower-left corner is ``(x0, y0)`` on log-log axes.

    The horizontal leg spans ``decades`` decades of x and is labelled ``1``; the vertical
    leg rises (or falls, for negative ``order``) ``order * decades`` decades and is
    labelled with ``order``. A line of slope ``order`` is parallel to the hypotenuse.
    """
    if ax.get_xscale() != "log" or ax.get_yscale() != "log":
        raise ValueError("add_slope_guide needs log-log axes (use ax.loglog or set both scales)")
    x1 = x0 * 10**decades
    y1 = y0 * 10 ** (order * decades)
    ax.plot([x0, x1, x1, x0], [y0, y0, y1, y0], "-", color=color, lw=0.9, solid_joinstyle="miter")
    below = order >= 0  # put the "1" on the outside of the triangle
    ax.annotate(
        "1",
        xy=(np.sqrt(x0 * x1), y0),
        xytext=(0, -3 if below else 3),
        textcoords="offset points",
        ha="center",
        va="top" if below else "bottom",
        color=color,
        fontsize=8,
    )
    ax.annotate(
        f"{order:g}",
        xy=(x1, np.sqrt(y0 * y1)),
        xytext=(4, 0),
        textcoords="offset points",
        ha="left",
        va="center",
        color=color,
        fontsize=8,
    )


def plot_convergence(
    ax: Axes,
    h: np.ndarray,
    err: np.ndarray,
    order: float,
    *,
    label: str | None = None,
    color: str | None = None,
    guide: bool = True,
    annotate: bool = True,
) -> float:
    r"""Plot ``err`` against ``h`` on log-log axes and return the fitted order.

    Draws the data, the least-squares line :math:`\log_{10} e = a\log_{10} h + b`, a slope
    triangle for the *expected* ``order`` placed half a decade below the fitted line, and
    the fitted slope as text. (The slope's standard error from :func:`loglog_fit` is left
    off the label on purpose: it measures scatter about a line, not whether the model
    ``C h^p`` is right, and for convergence data the second-order term biases the slope
    by more than that scatter.) Raises ``ValueError`` for non-positive data, which a log
    axis would drop without a word.
    """
    h, err = np.asarray(h, float), np.asarray(err, float)
    if h.shape != err.shape or h.ndim != 1:
        raise ValueError(f"h and err must be 1-D arrays of one length, got {h.shape} and {err.shape}")
    if np.any(h <= 0) or np.any(err <= 0):
        raise ValueError("log axes need h > 0 and err > 0 (an error of exactly 0 cannot be drawn)")
    slope, icpt, se = loglog_fit(h, err)
    (line,) = ax.loglog(h, err, "o", ms=4, color=color, label=label)
    c = line.get_color()
    hh = np.array([h.min(), h.max()])
    ax.loglog(hh, 10**icpt * hh**slope, "-", lw=0.8, color=c, alpha=0.8)
    if guide:
        x0 = h.max() * 10 ** (-0.15 - 0.5)
        y0 = 10**icpt * x0**slope * 10**-0.5
        add_slope_guide(ax, order, x0, y0)
    if annotate:  # upper left of the line: the free side of a rising line
        k = np.argsort(h)[len(h) // 2]
        ax.annotate(
            f"fitted slope {slope:.2f}",
            xy=(h[k], 10**icpt * h[k] ** slope),
            xytext=(-8, 8),
            textcoords="offset points",
            ha="right",
            color=c,
            fontsize=8,
        )
    ax.xaxis.set_minor_formatter(NullFormatter())  # crowded minor labels on a short log axis
    ax.set(xlabel="step size $h$", ylabel="absolute error $|I_h - I|$")
    return slope


def pi_ticks(ax: Axes, axis: str = "x", base: float = np.pi / 2) -> None:
    r"""Label major ticks of ``axis`` ('x' or 'y') every ``base`` as multiples of :math:`\pi`.

    Labels are mathtext such as ``$\frac{\pi}{2}$`` and ``$-\pi$``; no LaTeX install needed.
    """

    def fmt(value: float, pos: int | None = None) -> str:
        frac = Fraction(value / np.pi).limit_denominator(12)
        n, d = frac.numerator, frac.denominator
        if n == 0:
            return "$0$"
        sign = "-" if n < 0 else ""
        top = ("" if abs(n) == 1 else str(abs(n))) + r"\pi"
        return f"${sign}{top}$" if d == 1 else rf"${sign}\frac{{{top}}}{{{d}}}$"

    target = ax.xaxis if axis == "x" else ax.yaxis
    target.set_major_locator(MultipleLocator(base))
    target.set_major_formatter(FuncFormatter(fmt))


def label_end(ax: Axes, line: Line2D, text: str, dx: float = 4.0) -> None:
    """Label ``line`` in its own colour just right of its last point (direct labelling)."""
    x, y = line.get_xydata()[-1]
    ax.annotate(
        text,
        xy=(x, y),
        xytext=(dx, 0),
        textcoords="offset points",
        color=line.get_color(),
        va="center",
        fontsize=9,
    )


def rgba_at(ax: Axes, x: float | np.ndarray, y: float | np.ndarray) -> np.ndarray:
    """Rendered RGBA (floats in [0, 1]) of the pixel under data point(s) ``(x, y)``.

    ``x`` and ``y`` may be scalars (returns shape ``(4,)``) or equal-length arrays (shape
    ``(n, 4)``); the figure is drawn once per call, so pass arrays when testing many
    points. Display coordinates start at the bottom left, the pixel buffer at the top
    left, hence the ``H - 1 - floor(py)``.
    """
    fig = ax.figure
    fig.canvas.draw()
    scalar = np.ndim(x) == 0 and np.ndim(y) == 0
    pts = np.column_stack(np.broadcast_arrays(np.atleast_1d(x), np.atleast_1d(y)))
    disp = ax.transData.transform(pts)
    buf = np.asarray(fig.canvas.buffer_rgba())
    out = buf[buf.shape[0] - 1 - np.floor(disp[:, 1]).astype(int), np.floor(disp[:, 0]).astype(int)] / 255.0
    return out[0] if scalar else out
