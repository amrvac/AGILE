#!/usr/bin/env python3
"""Multi-panel movie of the Orszag-Tang vortex, one panel per Riemann solver.

Reads the constant-z slices that movie.par writes (slice_type = 'dat', i.e.
2D AMRVAC snapshots of the conserved variables, no ghost cells) with
tools/python's datfile_utilities, and animates them side by side on one
colour scale.

    uv sync --extra analysis                 # numpy, matplotlib, scipy
    make arch=gnu
    for s in tvdlf hll hllc; do
        mpirun -np 8 ./agile -i movie.par movie_$s.par
    done
    python movie.py                          # movie/ot_rho.mp4, or .gif without ffmpeg
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

agile_dir = os.environ.get(
    "AGILE_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../.."))
sys.path.insert(0, os.path.join(agile_dir, "tools", "python"))
from amrvac_pytools.datfiles.reading import datfile_utilities as dfu  # noqa: E402


def load_slice(filename):
    """Return (header, w) for a 2D slice on the finest level's grid.

    w has shape (nx, ny, nw); coarser AMR blocks are upsampled
    piecewise-constant (see datfile_utilities.get_upsampled_data).
    """
    with open(filename, "rb") as f:
        header = dfu.get_header(f)
        w = dfu.get_upsampled_data(f)
    return header, w


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dir", default="movie", help="slice directory")
    parser.add_argument("--solvers", nargs="+", default=["tvdlf", "hll", "hllc"])
    parser.add_argument("--var", default="rho", help="conserved variable name")
    parser.add_argument("--limiter", default="koren", help="for the title only")
    parser.add_argument("--fps", type=int, default=15)
    parser.add_argument("--dpi", type=int, default=120)
    parser.add_argument("--cmap", default="viridis")
    parser.add_argument("--output", help="default: <dir>/ot_<var>.mp4 (.gif without ffmpeg)")
    parser.add_argument("--frames-dir", help="also save every frame as a PNG here")
    args = parser.parse_args()

    files = {}
    for s in args.solvers:
        files[s] = sorted(glob.glob(os.path.join(args.dir, f"{s}_d3_*_n[0-9][0-9][0-9][0-9].dat")))
        if not files[s]:
            sys.exit(f"no slices for {s} in {args.dir}")
    nframes = min(len(f) for f in files.values())
    if any(len(f) != nframes for f in files.values()):
        print("[warning] solvers have different frame counts; using the first",
              nframes, {s: len(f) for s, f in files.items()})

    # Load everything up front: 3 x ~100 frames of 256^2 doubles is small,
    # and it gives one colour scale over the whole movie.
    frames = {s: [] for s in args.solvers}
    times = []
    for n in range(nframes):
        t_n = None
        for s in args.solvers:
            header, w = load_slice(files[s][n])
            iw = header["w_names"].index(args.var)
            frames[s].append(w[:, :, iw].T)  # (y, x) for imshow
            if t_n is None:
                t_n = header["time"]
            elif not np.isclose(header["time"], t_n, rtol=0, atol=1e-8):
                print(f"[warning] frame {n}: t = {header['time']:.4f} for {s}, {t_n:.4f} for {args.solvers[0]}")
        times.append(t_n)
    extent = [header["xmin"][0], header["xmax"][0], header["xmin"][1], header["xmax"][1]]
    vmin = min(np.min(f) for fs in frames.values() for f in fs)
    vmax = max(np.max(f) for fs in frames.values() for f in fs)

    fig, axes = plt.subplots(1, len(args.solvers), figsize=(4.2 * len(args.solvers) + 1, 4.6),
                             layout="constrained", squeeze=False)
    # the default 3 pt edge pad leaves the y label clipped in the saved frames
    fig.get_layout_engine().set(w_pad=0.1, h_pad=0.1)
    axes = axes[0]
    images = []
    for ax, s in zip(axes, args.solvers):
        images.append(ax.imshow(frames[s][0], origin="lower", extent=extent,
                                vmin=vmin, vmax=vmax, cmap=args.cmap,
                                interpolation="nearest"))
        ax.set_title(s)
        ax.set_xlabel("x")
        ax.set_aspect("equal")
    axes[0].set_ylabel("y")
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)
    fig.colorbar(images[0], ax=axes, shrink=0.85, label=args.var)
    suptitle = fig.suptitle("")

    def draw(n):
        for im, s in zip(images, args.solvers):
            im.set_data(frames[s][n])
        suptitle.set_text(f"Orszag-Tang, {args.limiter} limiter, "
                          f"{frames[args.solvers[0]][n].shape[1]}$^2$   t = {times[n]:.2f}")
        return images + [suptitle]

    if args.frames_dir:
        os.makedirs(args.frames_dir, exist_ok=True)
        for n in range(nframes):
            draw(n)
            fig.savefig(os.path.join(args.frames_dir, f"ot_{args.var}_{n:04d}.png"), dpi=args.dpi)

    if animation.FFMpegWriter.isAvailable():
        writer, ext = animation.FFMpegWriter(fps=args.fps), ".mp4"
    else:
        writer, ext = animation.PillowWriter(fps=args.fps), ".gif"
    output = args.output or os.path.join(args.dir, f"ot_{args.var}{ext}")
    anim = animation.FuncAnimation(fig, draw, frames=nframes, blit=False)
    anim.save(output, writer=writer, dpi=args.dpi)
    print(f"wrote {output} ({nframes} frames, t = {times[0]:.2f} .. {times[-1]:.2f})")


if __name__ == "__main__":
    main()
