#!/usr/bin/env python3
"""Plot the orbital plane (z = 0) of a colliding-winds run.

One PNG per input file, four panels: log density with the outlines of the AMR
blocks in the plane, log pressure, the mass fraction of wind 1 (tracer 1) and
v_x. The inputs are normally the z = 0 slices both par files write
(slicedir 3, slice_type = 'dat'):

    uv sync --extra analysis                     # numpy, matplotlib
    python plot_cwb.py                           # every output/*_d3_*.dat
    python plot_cwb.py output/cwb_quick_d3_x+0.00D+00_n0007.dat --outdir frames
    python plot_cwb.py --movie output/cwb_quick_orbit_d3_*.dat   # one mp4

A full 3D restart file (output/cwb_quick0002.dat) works too; its plane just above
z = 0 is cut out. That reads the whole domain upsampled to the finest level,
2**(3*(levmax-1)) * domain cells * nw doubles, so it only suits small runs.

What the files hold:
  - 3D restart files are always *conserved* variables.
  - Slices follow saveprim, which both par files set: *primitive* variables
    under the conserved names (m1..m3 the velocity, e the pressure). Pass
    --conserved for a run with saveprim = F.
  - The tracer is conserved (rho times the wind-1 fraction) in either case.
gamma is taken from the file header.

--movie animates all inputs into one file (mp4 through ffmpeg, which the
analysis extra's imageio-ffmpeg provides; a GIF otherwise). Its colour scales
are fixed over the whole sequence, from percentiles so that a few extreme
cells at the stars do not set them, and coarser frames are repeated up to the
finest resolution present so the image size does not jump when the finest
level first reaches the plane.
"""
import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.collections import LineCollection

agile_dir = os.environ.get(
    "AGILE_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../.."))
sys.path.insert(0, os.path.join(agile_dir, "tools", "python"))
from amrvac_pytools.datfiles.reading import datfile_utilities as dfu  # noqa: E402

# field, log scale, colour map, title
PANELS = [
    ("rho", True, "viridis", r"$\log_{10}\rho$"),
    ("p", True, "magma", r"$\log_{10}p$"),
    ("trc1", False, "coolwarm", "wind-1 fraction"),
    ("v1", False, "RdBu_r", r"$v_x$"),
]


def load(filename):
    with open(filename, "rb") as f:
        header = dfu.get_header(f)
        try:
            w = dfu.get_upsampled_data(f)
        except ValueError:
            # get_upsampled_data flags missing blocks by NaN, so a solution
            # that went NaN reports itself as "not covering the domain"
            sys.exit(f"{filename}: data contains NaN (or blocks are missing)")
        levels, indices, _ = dfu.get_tree_info(f)
    return header, w, levels, indices


def to_plane(header, w):
    """(x, y, nw) data of the plane, its z, and whether it came from a slice."""
    if header["ndim"] == 2:
        return w, 0.0, True
    xmin, xmax = np.asarray(header["xmin"]), np.asarray(header["xmax"])
    nz = w.shape[2]
    k = int(np.clip(np.floor(-xmin[2] / (xmax[2] - xmin[2]) * nz), 0, nz - 1))
    z = xmin[2] + (k + 0.5) * (xmax[2] - xmin[2]) / nz
    return w[:, :, k, :], z, False


def plane_fields(header, w2, primitive):
    """Primitive fields on the plane, transposed to (y, x) for imshow."""
    names = header["w_names"]
    get = lambda v: w2[:, :, names.index(v)].T  # noqa: E731
    rho = get("rho")
    if primitive:
        v1 = get("m1")
        p = get("e")
    else:
        mom = [get(m) for m in ("m1", "m2", "m3")]
        v1 = mom[0] / rho
        p = (header["gamma"] - 1.0) * (get("e") - 0.5 * sum(m**2 for m in mom) / rho)
    return {"rho": rho, "p": p, "trc1": get("trc1") / rho, "v1": v1}


def block_outlines(header, levels, indices, z):
    """(x, y) outlines of the leaf blocks in the plane."""
    xmin, xmax = np.asarray(header["xmin"]), np.asarray(header["xmax"])
    width1 = (xmax - xmin) * np.asarray(header["block_nx"]) / np.asarray(header["domain_nx"])
    lines = []
    for lvl, ix in zip(levels, indices):
        wd = width1 / 2 ** (lvl - 1)
        lo = xmin + (ix - 1) * wd
        hi = lo + wd
        if header["ndim"] == 3 and not lo[2] <= z < hi[2]:
            continue
        lines.append(([lo[0], hi[0], hi[0], lo[0], lo[0]],
                      [lo[1], lo[1], hi[1], hi[1], lo[1]]))
    return lines


def plot(filename, out, conserved, dpi, blocks):
    header, w, levels, indices = load(filename)
    w2, z, is_slice = to_plane(header, w)
    fields = plane_fields(header, w2, primitive=is_slice and not conserved)
    xmin, xmax = header["xmin"], header["xmax"]
    extent = [xmin[0], xmax[0], xmin[1], xmax[1]]

    fig, axes = plt.subplots(1, len(PANELS), figsize=(4.5 * len(PANELS), 4.6),
                             layout="constrained")
    for ax, (v, log, cmap, title) in zip(axes, PANELS):
        q = fields[v]
        if log:
            q = np.log10(np.maximum(q, 1e-30))
        im = ax.imshow(q, origin="lower", extent=extent, cmap=cmap)
        ax.set_title(title)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.colorbar(im, ax=ax, shrink=0.8)
    if blocks:
        for xs, ys in block_outlines(header, levels, indices, z):
            axes[0].plot(xs, ys, "w-", lw=0.3)
    fig.suptitle(f"t = {header['time']:.3f},  z = {z:.3g}")
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    print(out)


def movie(files, output, conserved, dpi, blocks, fps):
    frames = {v: [] for v, *_ in PANELS}
    outlines, times = [], []
    for fn in files:
        header, w, levels, indices = load(fn)
        w2, z, is_slice = to_plane(header, w)
        fields = plane_fields(header, w2, primitive=is_slice and not conserved)
        for v, log, *_ in PANELS:
            q = fields[v]
            frames[v].append((np.log10(np.maximum(q, 1e-30)) if log else q).astype(np.float32))
        outlines.append([list(zip(xs, ys)) for xs, ys in block_outlines(header, levels, indices, z)])
        times.append(header["time"])
    for v in frames:
        ny = max(f.shape[0] for f in frames[v])
        nx = max(f.shape[1] for f in frames[v])
        frames[v] = [np.repeat(np.repeat(f, ny // f.shape[0], axis=0), nx // f.shape[1], axis=1)
                     for f in frames[v]]
    xmin, xmax = header["xmin"], header["xmax"]
    extent = [xmin[0], xmax[0], xmin[1], xmax[1]]

    fig, axes = plt.subplots(1, len(PANELS), figsize=(4.5 * len(PANELS), 4.6),
                             layout="constrained")
    images = []
    for ax, (v, log, cmap, title) in zip(axes, PANELS):
        stack = np.stack(frames[v])
        if v == "trc1":
            vmin, vmax = 0.0, 1.0
        elif v == "v1":
            lim = float(np.percentile(np.abs(stack), 99.9)) or 1.0
            vmin, vmax = -lim, lim
        else:
            vmin, vmax = np.percentile(stack, [0.01, 99.99])
        im = ax.imshow(frames[v][0], origin="lower", extent=extent, cmap=cmap,
                       vmin=vmin, vmax=vmax, interpolation="nearest")
        ax.set_title(title)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.colorbar(im, ax=ax, shrink=0.8)
        images.append(im)
    grid = None
    if blocks:
        grid = LineCollection(outlines[0], colors="w", linewidths=0.3, alpha=0.6)
        axes[0].add_collection(grid)
    suptitle = fig.suptitle("")

    def draw(n):
        for im, (v, *_) in zip(images, PANELS):
            im.set_data(frames[v][n])
        if grid is not None:
            grid.set_segments(outlines[n])
        suptitle.set_text(f"t = {times[n]:.2f},  z = 0")
        return images

    if not animation.FFMpegWriter.isAvailable():
        try:
            import imageio_ffmpeg
            matplotlib.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
        except ImportError:
            pass
    if animation.FFMpegWriter.isAvailable():
        writer, ext = animation.FFMpegWriter(fps=fps), ".mp4"
    else:
        writer, ext = animation.PillowWriter(fps=fps), ".gif"
    output = os.path.splitext(output)[0] + ext
    animation.FuncAnimation(fig, draw, frames=len(files), blit=False).save(
        output, writer=writer, dpi=dpi)
    print(f"wrote {output} ({len(files)} frames, t = {times[0]:.2f} .. {times[-1]:.2f})")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("files", nargs="*",
                        help="z = 0 slices or 3D restart files (default: output/*_d3_*.dat)")
    parser.add_argument("--outdir", help="directory for the PNGs (default: next to each input)")
    parser.add_argument("--conserved", action="store_true",
                        help="slices hold conserved variables (run with saveprim = F)")
    parser.add_argument("--dpi", type=int, default=90)
    parser.add_argument("--no-blocks", dest="blocks", action="store_false",
                        help="do not outline the AMR blocks")
    parser.add_argument("--movie", nargs="?", const="output/cwb.mp4", metavar="FILE",
                        help="animate all inputs into FILE (default output/cwb.mp4)")
    parser.add_argument("--fps", type=int, default=20)
    args = parser.parse_args()

    files = args.files or sorted(glob.glob(os.path.join("output", "*_d3_*.dat")))
    if not files:
        sys.exit("no input files (and no output/*_d3_*.dat slices)")
    if args.movie:
        movie(files, args.movie, args.conserved, args.dpi, args.blocks, args.fps)
        return
    if args.outdir:
        os.makedirs(args.outdir, exist_ok=True)
    for fn in files:
        name = os.path.splitext(os.path.basename(fn))[0] + ".png"
        out = os.path.join(args.outdir or os.path.dirname(fn) or ".", name)
        plot(fn, out, args.conserved, args.dpi, args.blocks)


if __name__ == "__main__":
    main()
