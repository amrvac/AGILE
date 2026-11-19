#!/usr/bin/env python3
"""Multi-panel movie of the toroidal-field relativistic jet from its slices.

Reads the 2D slices agile.par writes (slicedir 1 and 2 are the jet axis
planes x = 0 and y = 0, slicedir 3 the cross-section at z = 40) when run with
slice_type = 'dat', which a second par file given after agile.par can set
(&filelist slice_type = 'dat' /):

    uv sync --extra analysis                 # numpy, matplotlib, imageio-ffmpeg
    make arch=nvidia OPENACC=1 NOGPUDIRECT=1
    mpirun -np 1 ./agile -i agile.par local.par
    python movie.py --blocks                 # output/jet_d1.mp4 (.gif without ffmpeg)

agile.par sets saveprim = T, so the slices hold *primitive* variables under the
conserved names in the header: 'e' is the gas pressure, 'm1..m3' the spatial
four-velocity u^i = lfac v^i. LABELS below gives them their real names.
Coarser AMR blocks are upsampled piecewise-constant onto the finest grid
(datfile_utilities.get_upsampled_data), as in tests/mhd/Orszag-Tang_3D/movie.py.
--blocks outlines the AMR blocks the slice cuts through, from the level and
index each leaf block carries in the file.
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

# header name -> (label, log scale, diverging colour map)
LABELS = {
    "rho": (r"$\log_{10}\rho$", True, False),
    "e": (r"$\log_{10}p$", True, False),
    "m1": (r"$u_x$", False, True),
    "m2": (r"$u_y$", False, True),
    "m3": (r"$u_z$", False, True),
    "b1": (r"$B_x$", False, True),
    "b2": (r"$B_y$", False, True),
    "b3": (r"$B_z$", False, True),
    "psi": (r"$\psi$", False, True),
    "trc1": ("jet tracer", False, False),
    "lfac": (r"Lorentz factor $\Gamma$", False, False),
}
# the in-plane axes of each slice direction
AXES = {1: ("y", "z"), 2: ("x", "z"), 3: ("x", "y")}


def load_slice(filename):
    """Return (header, w, outlines): w of shape (n1, n2, nw) on the finest
    level's grid, outlines one closed (5, 2) polyline per leaf block."""
    with open(filename, "rb") as f:
        header = dfu.get_header(f)
        w = dfu.get_upsampled_data(f)
        block_lvls, block_ixs, _ = dfu.get_tree_info(f)
    return header, w, block_outlines(header, block_lvls, block_ixs)


def block_outlines(header, block_lvls, block_ixs):
    """Rectangle of each leaf block in physical coordinates. A level-l block
    spans (xmax - xmin) * block_nx / domain_nx / 2**(l-1), and its 1-based
    index counts blocks of that size from xmin."""
    xmin, xmax = np.asarray(header["xmin"]), np.asarray(header["xmax"])
    width1 = (xmax - xmin) * np.asarray(header["block_nx"]) / np.asarray(header["domain_nx"])
    outlines = []
    for lvl, ix in zip(block_lvls, block_ixs):
        w = width1 / 2 ** (lvl - 1)
        (x0, y0), (x1, y1) = xmin + (ix - 1) * w, xmin + ix * w
        outlines.append([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)])
    return outlines


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dir", default="output", help="slice directory")
    parser.add_argument("--base", default="jet_lfac3", help="base_filename of the run")
    parser.add_argument("--slicedir", type=int, default=1, choices=(1, 2, 3))
    parser.add_argument("--vars", nargs="+", default=["rho", "e", "lfac", "b1"],
                        help="header variable names, one panel each")
    parser.add_argument("--fps", type=int, default=8)
    parser.add_argument("--dpi", type=int, default=120)
    parser.add_argument("--output", help="default: <dir>/jet_d<slicedir>.mp4 (.gif without ffmpeg)")
    parser.add_argument("--frames-dir", help="also save every frame as a PNG here")
    parser.add_argument("--blocks", action="store_true", help="outline the AMR blocks")
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(
        args.dir, f"{args.base}_d{args.slicedir}_*_n[0-9][0-9][0-9][0-9].dat")))
    if not files:
        sys.exit(f"no d{args.slicedir} slices for {args.base} in {args.dir}")

    # 70 frames of 240 x 480 x 12 doubles is under 0.5 GB; loading up front
    # gives each panel one colour scale over the whole movie.
    frames = {v: [] for v in args.vars}
    times = []
    outlines = []
    for fn in files:
        header, w, block_lines = load_slice(fn)
        outlines.append(block_lines)
        for v in args.vars:
            q = w[:, :, header["w_names"].index(v)].T  # (axis 2, axis 1) for imshow
            frames[v].append(np.log10(np.maximum(q, 1e-30)) if LABELS.get(v, ("", False))[1] else q)
        times.append(header["time"])
    # get_upsampled_data works to the finest level present in *that* slice,
    # so a plane the refinement reaches only later (z = 40 before the head
    # arrives) changes resolution mid-movie: repeat cells up to the finest.
    for v in args.vars:
        n1 = max(f.shape[0] for f in frames[v])
        n2 = max(f.shape[1] for f in frames[v])
        frames[v] = [np.repeat(np.repeat(f, n1 // f.shape[0], axis=0), n2 // f.shape[1], axis=1)
                     for f in frames[v]]
    extent = [header["xmin"][0], header["xmax"][0], header["xmin"][1], header["xmax"][1]]
    xlabel, ylabel = AXES[args.slicedir]

    width, height = extent[1] - extent[0], extent[3] - extent[2]
    panel_w = 3.0 if height > width else 4.0
    fig, axes = plt.subplots(1, len(args.vars),
                             figsize=(panel_w * len(args.vars) + 0.5, panel_w * height / width + 1.0),
                             layout="constrained", squeeze=False)
    fig.get_layout_engine().set(w_pad=0.1, h_pad=0.1)
    axes = axes[0]
    images = []
    grids = []
    for ax, v in zip(axes, args.vars):
        label, _, diverging = LABELS.get(v, (v, False, False))
        stack = np.stack(frames[v])
        if diverging:
            # a few extreme cells at the shocks would wash out the rest
            lim = np.percentile(np.abs(stack), 99.9) or 1.0
            vmin, vmax, cmap = -lim, lim, "RdBu_r"
        else:
            # percentiles, not min/max: a handful of cells at the density and
            # pressure floors (srmhd_small_*) would otherwise set the scale
            vmin, vmax = np.percentile(stack, [0.01, 99.99])
            cmap = "magma" if v != "trc1" else "viridis"
        im = ax.imshow(frames[v][0], origin="lower", extent=extent, vmin=vmin, vmax=vmax,
                       cmap=cmap, interpolation="nearest")
        fig.colorbar(im, ax=ax, location="bottom", shrink=0.9, label=label)
        ax.set_xlabel(xlabel)
        ax.set_aspect("equal")
        images.append(im)
        if args.blocks:
            # white on magma/viridis, black on the white-centred RdBu_r
            grid = LineCollection(outlines[0], colors="k" if diverging else "w",
                                  linewidths=0.3, alpha=0.5)
            ax.add_collection(grid)
            grids.append(grid)
    axes[0].set_ylabel(ylabel)
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)
    suptitle = fig.suptitle("")
    plane = {1: "x = 0", 2: "y = 0", 3: "z = 40"}[args.slicedir]

    def draw(n):
        for im, v in zip(images, args.vars):
            im.set_data(frames[v][n])
        for grid in grids:
            grid.set_segments(outlines[n])
        suptitle.set_text(f"SRMHD toroidal jet, {plane}   t = {times[n]:.1f}")
        return images + grids + [suptitle]

    if args.frames_dir:
        os.makedirs(args.frames_dir, exist_ok=True)
        for n in range(len(files)):
            draw(n)
            fig.savefig(os.path.join(args.frames_dir, f"jet_d{args.slicedir}_{n:04d}.png"), dpi=args.dpi)

    if not animation.FFMpegWriter.isAvailable():
        # no ffmpeg on the PATH: use the binary imageio-ffmpeg ships, if installed
        try:
            import imageio_ffmpeg
            matplotlib.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
        except ImportError:
            pass
    if animation.FFMpegWriter.isAvailable():
        writer, ext = animation.FFMpegWriter(fps=args.fps), ".mp4"
    else:
        writer, ext = animation.PillowWriter(fps=args.fps), ".gif"
    output = args.output or os.path.join(args.dir, f"jet_d{args.slicedir}{ext}")
    anim = animation.FuncAnimation(fig, draw, frames=len(files), blit=False)
    anim.save(output, writer=writer, dpi=args.dpi)
    print(f"wrote {output} ({len(files)} frames, t = {times[0]:.1f} .. {times[-1]:.1f})")


if __name__ == "__main__":
    main()
