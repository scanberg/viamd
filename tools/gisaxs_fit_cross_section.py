#!/usr/bin/env python3
"""
gisaxs_fit_cross_section.py -- fit a fibril cross-section template for the VIAMD scattering component (GISAXS / GISANS).

Input: an electron density of one fibril slice (Gaussian cube file) and, optionally, the positions of the
coarse grained beads of that slice (in the same coordinates, in slice order as in the MD topology).

The density is projected along the fibril axis onto the cross-section plane. A set of 2D Gaussians is then fitted
to the projected density in Fourier space, for |q| <= q_max (the q-range of the GISAXS calculation), which is the
only range where the template has to be right. Two models:

  --mode beads   one Gaussian per bead, at the bead position, with a free weight per bead and a free width per
                 bead name (e.g. CC, IC, OC). Requires --beads.
  --mode grid    Gaussians on a hexagonal grid (--spacing, sigma = spacing/2) with free weights. Represents any shape.

Output: a text template for VIAMD ('File' cross-section):
    anchor <u> <v>                    one per bead, slice order (if --beads is given)
    gauss  <u> <v> <weight> <sigma>   weights normalized to 1
All lengths in Ångström. The template frame is the cross-section plane with its origin at the density centroid.

Bead file format: one bead per line, 'name x y z' in Ångström (cube coordinates, i.e. bohr * 0.529177).

Examples
  python gisaxs_fit_cross_section.py slice.cube --axis z --beads slice_beads.txt --mode beads -o cellulose_7bead.txt
  python gisaxs_fit_cross_section.py slice.cube --axis z --mode grid --spacing 6 -o cellulose_grid.txt
"""
from __future__ import annotations

import argparse
import sys

import numpy as np

BOHR = 0.529177210903


def read_cube(path):
    with open(path) as f:
        f.readline(); f.readline()
        t = f.readline().split()
        natoms = int(t[0]); origin = np.array(list(map(float, t[1:4])))
        axes, n = [], []
        for _ in range(3):
            t = f.readline().split()
            n.append(int(t[0])); axes.append(list(map(float, t[1:4])))
        axes = np.array(axes)
        # Negative counts mean Å in some writers; positive means bohr (standard)
        unit = BOHR if n[0] > 0 else 1.0
        n = [abs(v) for v in n]
        atoms = []
        for _ in range(abs(natoms)):
            t = f.readline().split()
            atoms.append((int(t[0]), float(t[1]), *map(float, t[2:5])))
        if natoms < 0:
            f.readline()   # orbital line
        data = np.array(f.read().split(), dtype=float)
    data = data[: n[0] * n[1] * n[2]].reshape(n)
    origin = origin * unit
    axes = axes * unit
    if unit == BOHR:
        data = data / BOHR**3          # e/bohr^3 -> e/Å^3
    atoms = np.array([[a[0], a[2] * unit, a[3] * unit, a[4] * unit] for a in atoms]) if atoms else np.zeros((0, 4))
    return origin, axes, data, atoms


def read_beads(path):
    names, pos = [], []
    with open(path) as f:
        for line in f:
            line = line.split("#")[0].split()
            if len(line) < 4:
                continue
            names.append(line[0]); pos.append(list(map(float, line[1:4])))
    return names, np.array(pos)


def project(origin, axes, data, axis):
    """Projects the density along cube axis 'axis' (0,1,2). Returns plane coords (M,2), areal density (M,), pixel area,
    plane basis (2,3) and the slice thickness along the axis."""
    for i in range(3):
        for j in range(i + 1, 3):
            if abs(np.dot(axes[i], axes[j])) > 1e-6 * np.linalg.norm(axes[i]) * np.linalg.norm(axes[j]):
                sys.exit("non-orthogonal cube axes are not supported")
    others = [i for i in range(3) if i != axis]
    dl = np.linalg.norm(axes[axis])
    proj = data.sum(axis=axis) * dl                       # e / Å^2
    a0, a1 = axes[others[0]], axes[others[1]]
    e0, e1 = a0 / np.linalg.norm(a0), a1 / np.linalg.norm(a1)
    n0, n1 = proj.shape
    i0, i1 = np.meshgrid(np.arange(n0), np.arange(n1), indexing="ij")
    u = (i0 * np.linalg.norm(a0)).ravel()
    v = (i1 * np.linalg.norm(a1)).ravel()
    # Absolute position of the plane origin projected onto the plane
    ou, ov = np.dot(origin, e0), np.dot(origin, e1)
    area = np.linalg.norm(a0) * np.linalg.norm(a1)
    return np.column_stack([u + ou, v + ov]), proj.ravel(), area, np.array([e0, e1]), data.shape[axis] * dl


def q_points(q_max, n_r=40, n_phi=72):
    qr = (np.arange(n_r) + 0.5) / n_r * q_max
    ph = np.arange(n_phi) / n_phi * 2 * np.pi
    Q, P = np.meshgrid(qr, ph, indexing="ij")
    q = np.column_stack([(Q * np.cos(P)).ravel(), (Q * np.sin(P)).ravel()])
    w = Q.ravel()                                         # 2D area element ~ |q| dq dphi
    return q, w


def fourier(pos, weight, q, chunk=4096):
    F = np.zeros(len(q), dtype=complex)
    for s in range(0, len(pos), chunk):
        p = pos[s:s + chunk]; wt = weight[s:s + chunk]
        F += np.exp(-1j * (q @ p.T)) @ wt
    return F


def model_matrix(centers, sigma, q):
    """Columns: F of unit weight Gaussians at centers with widths sigma."""
    q2 = (q ** 2).sum(axis=1)
    return np.exp(-0.5 * q2[:, None] * sigma[None, :] ** 2) * np.exp(-1j * (q @ centers.T))


def solve_weights(A, F, w):
    sw = np.sqrt(w)[:, None]
    Ar = np.vstack([(A * sw).real, (A * sw).imag])
    Fr = np.concatenate([(F * sw[:, 0]).real, (F * sw[:, 0]).imag])
    x, *_ = np.linalg.lstsq(Ar, Fr, rcond=None)
    # Simple non-negativity: clip and refit the active set
    for _ in range(10):
        neg = x < 0
        if not neg.any():
            break
        keep = ~neg
        x = np.zeros_like(x)
        xk, *_ = np.linalg.lstsq(Ar[:, keep], Fr, rcond=None)
        x[keep] = xk
    res = np.linalg.norm(Ar @ x - Fr) / max(np.linalg.norm(Fr), 1e-30)
    return x, res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cube")
    ap.add_argument("--axis", choices=("x", "y", "z"), default="z", help="cube axis along the fibril (default z)")
    ap.add_argument("--beads", help="bead positions: 'name x y z' per line (Å), slice order")
    ap.add_argument("--mode", choices=("beads", "grid"), default=None)
    ap.add_argument("--q-max", type=float, default=2.5, help="nm^-1 (default 2.5)")
    ap.add_argument("--spacing", type=float, default=6.0, help="grid mode spacing, Å")
    ap.add_argument("--threshold", type=float, default=0.02, help="grid mode: region where the density exceeds this fraction of the maximum")
    ap.add_argument("-o", "--output", default="cross_section_template.txt")
    a = ap.parse_args(argv)

    origin, axes, data, atoms = read_cube(a.cube)
    axis = "xyz".index(a.axis)
    pos, rho, area, basis, thickness = project(origin, axes, data, axis)
    mask = rho > 0
    pos, rho = pos[mask], rho[mask]
    electrons = rho.sum() * area
    centroid = (pos * rho[:, None]).sum(axis=0) / rho.sum()
    pos = pos - centroid
    print(f"slice: {electrons:.1f} electrons, thickness {thickness:.2f} Å -> {electrons / thickness:.2f} e/Å along the fibril")
    if len(atoms):
        z_sum = atoms[:, 0].sum()
        print(f"sum of atomic numbers in the cube: {z_sum:.0f} (a large difference means valence-only / pseudopotential density)")

    q_max = a.q_max * 0.1
    q, qw = q_points(q_max)
    F = fourier(pos, rho * area, q)

    anchors = None
    names = None
    if a.beads:
        names, bp = read_beads(a.beads)
        uv = bp @ basis.T - centroid
        anchors = uv
    mode = a.mode or ("beads" if anchors is not None else "grid")

    if mode == "beads":
        if anchors is None:
            sys.exit("--mode beads requires --beads")
        types = sorted(set(names))
        tidx = np.array([types.index(n) for n in names])
        sig = {t: 6.0 for t in types}
        best = None
        # Coordinate descent over the width per bead type, weights by linear least squares
        for sweep in range(6):
            for t in types:
                cands = sig[t] * np.array([0.5, 0.7, 0.85, 1.0, 1.2, 1.45, 2.0]) if sweep == 0 else \
                        sig[t] * np.array([0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.25])
                for s in cands:
                    trial = dict(sig); trial[t] = max(float(s), 0.3)
                    sv = np.array([trial[types[i]] for i in tidx])
                    x, res = solve_weights(model_matrix(anchors, sv, q), F, qw)
                    if best is None or res < best[0]:
                        best = (res, dict(trial), x, sv)
                sig = dict(best[1])
        res, sig, x, sv = best
        centers = anchors
        print("fitted widths: " + ", ".join(f"{t} {sig[t]:.2f} Å" for t in types))
    else:
        thr = a.threshold * rho.max()
        reach = np.sqrt(((pos[rho > thr]) ** 2).sum(axis=1)).max() + a.spacing
        dy = a.spacing * np.sqrt(3) / 2
        pts = []
        for j in range(-int(reach / dy) - 1, int(reach / dy) + 2):
            off = 0.5 * a.spacing if j % 2 else 0.0
            for i in range(-int(reach / a.spacing) - 2, int(reach / a.spacing) + 3):
                p = np.array([i * a.spacing + off, j * dy])
                if np.linalg.norm(p) <= reach:
                    pts.append(p)
        centers = np.array(pts)
        # Keep grid points near density
        d2 = ((centers[:, None, :] - pos[None, rho > thr, :]) ** 2).sum(axis=2).min(axis=1)
        centers = centers[d2 <= a.spacing ** 2]
        sv = np.full(len(centers), 0.5 * a.spacing)
        x, res = solve_weights(model_matrix(centers, sv, q), F, qw)

    keep = x > 1e-6 * x.max()
    print(f"{keep.sum()} Gaussians, relative Fourier residual (|q| <= {a.q_max} nm^-1): {res:.3e}")
    wsum = x[keep].sum()
    print(f"fitted electrons {wsum:.1f} of {electrons:.1f}")

    # Rotationally averaged comparison
    qr = np.linspace(0.02, q_max, 60)
    ph = np.linspace(0, 2 * np.pi, 72, endpoint=False)
    qq = np.array([[r * np.cos(p), r * np.sin(p)] for r in qr for p in ph])
    Fr = fourier(pos, rho * area, qq)
    Fm = model_matrix(centers[keep], sv[keep], qq) @ x[keep]
    Ir = (np.abs(Fr) ** 2).reshape(len(qr), len(ph)).mean(axis=1)
    Im = (np.abs(Fm) ** 2).reshape(len(qr), len(ph)).mean(axis=1)
    print("  q [nm^-1]   I_ref        I_model      ratio")
    for i in range(0, len(qr), 6):
        print(f"  {qr[i] * 10:8.3f}   {Ir[i]:.4e}   {Im[i]:.4e}   {Im[i] / max(Ir[i], 1e-300):.3f}")

    with open(a.output, "w") as f:
        f.write(f"# Fibril cross-section template, fitted by gisaxs_fit_cross_section.py from {a.cube}\n")
        f.write(f"# mode {mode}, q_max {a.q_max} nm^-1, residual {res:.3e}, {electrons:.1f} e per slice of {thickness:.2f} Å\n")
        if anchors is not None:
            for p, n in zip(anchors, names):
                f.write(f"anchor {p[0]:.4f} {p[1]:.4f}   # {n}\n")
        for c, w, s in zip(centers[keep], x[keep] / wsum, sv[keep]):
            f.write(f"gauss {c[0]:.4f} {c[1]:.4f} {w:.6f} {s:.4f}\n")
    print(f"wrote {a.output}")


if __name__ == "__main__":
    main()
