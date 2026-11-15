"""Softened gravitational accelerations, computed several ways (PY-10 notes, section 3).

Every function here computes, for bodies i = 0..N-1 with positions x_i and masses m_i,

    a_i = G * sum_{j != i} m_j (x_j - x_i) / (|x_j - x_i|^2 + eps^2)^(3/2)       (notes eq. 3.1)

with G = 1 and Plummer softening eps (eps = 0 is exact Newtonian gravity). They differ only in
how the double sum is evaluated, which is the point of the module: same physics, different cost.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.spatial import cKDTree

G = 1.0

try:                                   # Numba is optional: `pip install "nbody[fast]"`
    from numba import njit
    NUMBA_AVAILABLE = True
except ImportError:                    # the package works without it; accel_numba is then absent
    NUMBA_AVAILABLE = False


def accel_loop(pos: list[list[float]], mass: list[float], eps: float = 0.0) -> list[list[float]]:
    """Eq. (3.1) as a pure-Python double loop over lists: the baseline (notes Build 1)."""
    n = len(pos)
    eps2 = eps * eps
    acc = []
    for i in range(n):
        xi, yi, zi = pos[i]
        ax = ay = az = 0.0
        for j in range(n):
            if j == i:
                continue
            dx = pos[j][0] - xi
            dy = pos[j][1] - yi
            dz = pos[j][2] - zi
            r2 = dx * dx + dy * dy + dz * dz + eps2
            w = mass[j] / (r2 * math.sqrt(r2))
            ax += w * dx
            ay += w * dy
            az += w * dz
        acc.append([G * ax, G * ay, G * az])
    return acc


def accel_broadcast(pos: np.ndarray, mass: np.ndarray, eps: float = 0.0) -> np.ndarray:
    """Eq. (3.1) by broadcasting in the (N, N, 3) layout of PY-6 Build 4 (notes Build 2a)."""
    d = pos[np.newaxis, :, :] - pos[:, np.newaxis, :]      # (N, N, 3): d[i, j] = x_j - x_i
    r2 = (d * d).sum(axis=-1) + eps**2                     # (N, N)
    np.fill_diagonal(r2, 1.0)                              # i = j: d is zero, so any finite value works
    w = mass / (r2 * np.sqrt(r2))                          # (N, N): m_j / r_ij^3, mass broadcast along j
    return G * (w[:, :, np.newaxis] * d).sum(axis=1)       # (N, 3)


def accel_numpy(pos: np.ndarray, mass: np.ndarray, eps: float = 0.0) -> np.ndarray:
    """Eq. (3.1) by broadcasting, coordinate axis first so every pass is over contiguous memory (Build 2b)."""
    x = np.ascontiguousarray(pos.T)                        # (3, N)
    d = x[:, np.newaxis, :] - x[:, :, np.newaxis]          # (3, N, N): d[k, i, j] = x_kj - x_ki
    r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] + eps**2  # (N, N)
    np.fill_diagonal(r2, 1.0)
    w = mass / (r2 * np.sqrt(r2))                          # (N, N)
    return G * (d * w).sum(axis=2).T                       # (3, N) -> (N, 3)


class ForceWorkspace:
    """Eq. (3.1) with every intermediate preallocated and written in place (Build 3).

    Calling the workspace returns a view of its own buffer, which the next call overwrites.
    Copy the result if you keep it (notes section 8, Exercise A3).
    """

    def __init__(self, n: int, dtype: type = np.float64) -> None:
        self.dtype = np.dtype(dtype)
        self.x = np.empty((3, n), self.dtype)
        self.d = np.empty((3, n, n), self.dtype)
        self.r2 = np.empty((n, n), self.dtype)
        self.tmp = np.empty((n, n), self.dtype)
        self.acc = np.empty((3, n), self.dtype)

    def __call__(self, pos: np.ndarray, mass: np.ndarray, eps: float = 0.0) -> np.ndarray:
        x, d, r2, tmp = self.x, self.d, self.r2, self.tmp
        x[...] = pos.T
        np.subtract(x[:, np.newaxis, :], x[:, :, np.newaxis], out=d)
        np.multiply(d[0], d[0], out=r2)
        np.multiply(d[1], d[1], out=tmp)
        r2 += tmp
        np.multiply(d[2], d[2], out=tmp)
        r2 += tmp
        r2 += self.dtype.type(eps * eps)
        np.fill_diagonal(r2, 1.0)
        np.sqrt(r2, out=tmp)
        tmp *= r2                                          # r^3
        np.divide(mass, tmp, out=tmp)                      # m_j / r_ij^3
        np.multiply(d, tmp, out=d)                         # d is scratch from here on
        np.sum(d, axis=2, out=self.acc)
        if G != 1.0:
            self.acc *= G
        return self.acc.T


if NUMBA_AVAILABLE:
    @njit(cache=False)
    def accel_numba(pos: np.ndarray, mass: np.ndarray, eps: float = 0.0) -> np.ndarray:
        """Eq. (3.1) as the Build 1 double loop over arrays, compiled by Numba (Build 6)."""
        n = pos.shape[0]
        eps2 = eps * eps
        acc = np.zeros((n, 3))
        for i in range(n):
            ax = 0.0
            ay = 0.0
            az = 0.0
            for j in range(n):
                if j == i:
                    continue
                dx = pos[j, 0] - pos[i, 0]
                dy = pos[j, 1] - pos[i, 1]
                dz = pos[j, 2] - pos[i, 2]
                r2 = dx * dx + dy * dy + dz * dz + eps2
                w = mass[j] / (r2 * math.sqrt(r2))
                ax += w * dx
                ay += w * dy
                az += w * dz
            acc[i, 0] = G * ax
            acc[i, 1] = G * ay
            acc[i, 2] = G * az
        return acc


def _pair_accel(pos: np.ndarray, mass: np.ndarray, i: np.ndarray, j: np.ndarray, eps: float) -> np.ndarray:
    """Accumulate eq. (3.1) over an explicit list of unordered pairs (i < j), both directions."""
    d = pos[j] - pos[i]                                    # (P, 3)
    r2 = np.einsum("pk,pk->p", d, d) + eps**2
    inv_r3 = 1.0 / (r2 * np.sqrt(r2))
    n = len(mass)
    acc = np.empty((n, 3))
    for k in range(3):
        f = d[:, k] * inv_r3
        acc[:, k] = np.bincount(i, weights=mass[j] * f, minlength=n) - np.bincount(j, weights=mass[i] * f, minlength=n)
    return G * acc


def pairs_brute(pos: np.ndarray, rc: float, block: int = 256) -> tuple[np.ndarray, np.ndarray]:
    """All pairs i < j with |x_i - x_j| <= rc by checking every pair, in row blocks: O(N^2) work."""
    n = len(pos)
    x = np.ascontiguousarray(pos.T)                        # (3, N)
    ii, jj = [], []
    for start in range(0, n, block):
        stop = min(start + block, n)
        r2 = np.zeros((stop - start, n))
        for k in range(3):
            r2 += (x[k, start:stop, np.newaxis] - x[k, np.newaxis, :]) ** 2
        bi, bj = np.nonzero(r2 <= rc * rc)
        bi += start
        keep = bi < bj
        ii.append(bi[keep])
        jj.append(bj[keep])
    return np.concatenate(ii), np.concatenate(jj)


def pairs_tree(pos: np.ndarray, rc: float) -> tuple[np.ndarray, np.ndarray]:
    """All pairs i < j with |x_i - x_j| <= rc from a k-d tree: O(N log N) at fixed density (Build 7)."""
    p = cKDTree(pos).query_pairs(rc, output_type="ndarray")   # (P, 2), each row i < j
    return p[:, 0], p[:, 1]


def accel_short_range(pos: np.ndarray, mass: np.ndarray, rc: float, eps: float = 0.0,
                      method: str = "tree") -> np.ndarray:
    """Eq. (3.1) truncated to pairs closer than rc, neighbours found by `method` ("tree" or "brute")."""
    i, j = pairs_tree(pos, rc) if method == "tree" else pairs_brute(pos, rc)
    return _pair_accel(pos, mass, i, j, eps)
