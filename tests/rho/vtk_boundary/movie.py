#!/usr/bin/env python3
"""Movie of the data-driven inflow runs of this test: a smiling sun from
inlet.vtk carried across the box, either alone (agile.par, bc_data on the
x = 0 face, base 'data') or with a UFO added by specialbound_usr
(icarus.par, bc_icarus on the x = 1 face, base 'icarus').

Left: the image face on, following the pulse: movie.par writes y-z planes
every 1/16 in x, and each frame blends the two either side of the pulse
centre (the data pulse peaks at t = 0.125 on the inlet, and the flow is at
unit speed: x = 1.125 - t for icarus, t - 0.125 for data), so the image is shown as it is carried and
smeared by the scheme. Right: the side view, the x-z plane y = 0.72 through
the UFO, with the AMR blocks and the face-on plane marked. All are 2D slices
written by movie.par, given after the run's own par file:

    uv sync --extra analysis                 # numpy, matplotlib, imageio-ffmpeg
    make arch=nvidia OPENACC=1 NOGPUDIRECT=1 # or a CPU build (slower)
    mkdir -p output && cd output && ln -sf ../inlet.vtk .
    mpirun -np 1 ../agile -i ../icarus.par ../movie.par
    mpirun -np 1 ../agile -i ../agile.par ../movie.par    # sun only
    cd .. && python movie.py                 # output/icarus.mp4 (.gif without ffmpeg)
    python movie.py --base data              # output/data.mp4

movie.par refines to level 4 (256^3 effective); on a GTX 1660 Ti the run
takes a few minutes. Coarser blocks are upsampled piecewise-constant
(datfile_utilities.get_upsampled_data), so the AMR is visible as blockiness
where the mesh is coarse.
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

BG = "#0b0d17"
FG = "#e8e6f0"


def read_slice(filename):
    """Return (time, rho on the finest level of this slice, block outlines,
    extent)."""
    with open(filename, "rb") as f:
        header = dfu.get_header(f)
        w = dfu.get_upsampled_data(f)
        lvls, ixs, _ = dfu.get_tree_info(f)
    extent = [header["xmin"][0], header["xmax"][0], header["xmin"][1], header["xmax"][1]]
    rho = w[:, :, header["w_names"].index("rho")].T  # (axis 2, axis 1) for imshow
    return header["time"], rho, block_outlines(header, lvls, ixs), extent


def to_shape(q, shape):
    """Repeat cells up to shape (a slice upsamples only to its own finest
    level)."""
    return np.repeat(np.repeat(q, shape[0] // q.shape[0], axis=0),
                     shape[1] // q.shape[1], axis=1)


def slice_files(directory, base, slicedir):
    """{coordinate: [files in time order]} for one slice direction."""
    files = glob.glob(os.path.join(
        directory, f"{base}_d{slicedir}_*_n[0-9][0-9][0-9][0-9].dat"))
    if not files:
        sys.exit(f"no d{slicedir} slices for {base} in {directory}")
    planes = {}
    for fn in files:
        # e.g. icarus_d1_x+0.53D+00_n0042.dat
        token = os.path.basename(fn).split("_")[-2][1:]
        planes.setdefault(float(token.replace("D", "E")), []).append(fn)
    return {x: sorted(fns) for x, fns in sorted(planes.items())}


def block_outlines(header, block_lvls, block_ixs):
    """Rectangle of each leaf block in the slice plane (see
    tests/srmhd/jet_toroidal/movie.py)."""
    xmin, xmax = np.asarray(header["xmin"]), np.asarray(header["xmax"])
    width1 = (xmax - xmin) * np.asarray(header["block_nx"]) / np.asarray(header["domain_nx"])
    lines = []
    for lvl, ix in zip(block_lvls, block_ixs):
        w = width1 / 2 ** (lvl - 1)
        (x0, y0), (x1, y1) = xmin + (ix - 1) * w, xmin + ix * w
        lines.append([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)])
    return lines


def style(ax, title, xlabel, ylabel):
    ax.set_facecolor(BG)
    ax.set_title(title, color=FG, fontsize=11, pad=8)
    ax.set_xlabel(xlabel, color=FG)
    ax.set_ylabel(ylabel, color=FG)
    ax.tick_params(colors=FG, labelsize=8)
    for spine in ax.spines.values():
        spine.set_color("#3a3d52")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dir", default="output", help="slice directory")
    parser.add_argument("--base", default="icarus", help="base_filename of the run")
    parser.add_argument("--fps", type=int, default=24)
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument("--cmap", default="inferno")
    parser.add_argument("--output", help="default: <dir>/<base>.mp4 (.gif without ffmpeg)")
    parser.add_argument("--still", type=int, metavar="I",
                        help="save frame I as a PNG instead of the movie")
    args = parser.parse_args()

    (side_files,) = slice_files(args.dir, args.base, 2).values()
    face_files = slice_files(args.dir, args.base, 1)
    # flow direction and labels: icarus enters at x = 1 going -x, data at x = 0
    if args.base.startswith("icarus"):
        sign, x0, bc_name, what = -1.0, 1.125, "bc_icarus", "VTK sun + user UFO"
    else:
        sign, x0, bc_name, what = 1.0, -0.125, "bc_data", "VTK sun"
    xs = np.array(list(face_files))
    n = min([len(side_files)] + [len(f) for f in face_files.values()])

    side, side_blocks, times = [], [], []
    for fn in side_files[:n]:
        t, rho, lines, ext_side = read_slice(fn)
        side.append(rho); side_blocks.append(lines); times.append(t)
    shape = max(q.shape for q in side)
    side = [to_shape(q, shape) for q in side]
    vmin = 1.0
    vmax = np.percentile(np.stack(side), 99.95)

    def face_at(i):
        """The face-on image at the pulse centre, blended between planes."""
        xc = float(np.clip(x0 + sign * times[i], xs[0], xs[-1]))
        k = min(np.searchsorted(xs, xc), len(xs) - 1)
        k0, k1 = max(k - 1, 0), k
        a = 0.0 if k0 == k1 else (xc - xs[k0]) / (xs[k1] - xs[k0])
        q0 = read_slice(face_files[xs[k0]][i])[1]
        q1 = read_slice(face_files[xs[k1]][i])[1]
        fshape = max(q0.shape, q1.shape, shape)
        return xc, (1 - a) * to_shape(q0, fshape) + a * to_shape(q1, fshape)

    fig, (ax_f, ax_s) = plt.subplots(1, 2, figsize=(10.5, 5.6), layout="constrained")
    fig.patch.set_facecolor(BG)
    style(ax_f, "", "y", "z")
    style(ax_s, "side view:  y = 0.72, with the AMR blocks", "x", "z")

    kw = dict(origin="lower", cmap=args.cmap, vmin=vmin, vmax=vmax,
              interpolation="nearest")
    xc, q = face_at(0)
    im_f = ax_f.imshow(q, extent=[0, 1, 0, 1], **kw)
    im_s = ax_s.imshow(side[0], extent=ext_side, **kw)
    grid = LineCollection(side_blocks[0], colors="white", linewidths=0.4, alpha=0.3)
    ax_s.add_collection(grid)

    # the face-on plane, and where the data come in
    (cut,) = ax_s.plot([xc, xc], [0.12, 1.0], color="#7fd4ff", lw=1.0, ls=(0, (4, 3)), alpha=0.9)
    xin = 0.0 if sign > 0 else 1.0
    ax_s.annotate("", xy=(xin + sign * 0.20, 0.05), xytext=(xin + sign * 0.03, 0.05),
                  arrowprops=dict(arrowstyle="-|>", color=FG, lw=1.2))
    ax_s.text(xin + sign * 0.015, 0.08, f"inlet ({bc_name})", color=FG, fontsize=8,
              ha="left" if sign > 0 else "right")

    cbar = fig.colorbar(im_s, ax=[ax_f, ax_s], shrink=0.82, pad=0.02)
    cbar.set_label(r"density $\rho$", color=FG)
    cbar.ax.tick_params(colors=FG, labelsize=8)
    cbar.outline.set_edgecolor("#3a3d52")

    title = fig.suptitle("", color=FG, fontsize=13)

    def draw(i):
        xc, q = face_at(i)
        im_f.set_data(q)
        ax_f.set_title(f"face on, riding with the pulse:  x = {xc:.2f}",
                       color="#7fd4ff", fontsize=11, pad=8)
        cut.set_xdata([xc, xc])
        im_s.set_data(side[i])
        grid.set_segments(side_blocks[i])
        title.set_text(f"data-driven inflow: {what}      t = {times[i]:.2f}")
        return im_f, im_s, grid, cut, title

    if args.still is not None:
        draw(args.still)
        out = args.output or os.path.join(args.dir, f"{args.base}_{args.still:04d}.png")
        fig.savefig(out, dpi=args.dpi, facecolor=BG)
        print(out)
        return

    anim = animation.FuncAnimation(fig, draw, frames=n, blit=False)
    out = args.output
    if out is None or out.endswith(".mp4"):
        if animation.FFMpegWriter.isAvailable() or _use_imageio_ffmpeg():
            out = out or os.path.join(args.dir, f"{args.base}.mp4")
            writer = animation.FFMpegWriter(fps=args.fps, bitrate=6000,
                                            extra_args=["-pix_fmt", "yuv420p"])
        else:
            out = os.path.join(args.dir, f"{args.base}.gif")
            writer = animation.PillowWriter(fps=args.fps)
    else:
        writer = animation.PillowWriter(fps=args.fps)
    anim.save(out, writer=writer, dpi=args.dpi, savefig_kwargs={"facecolor": BG})
    print(out)


def _use_imageio_ffmpeg():
    """Point matplotlib at imageio-ffmpeg's bundled binary, if installed."""
    try:
        import imageio_ffmpeg
    except ImportError:
        return False
    matplotlib.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
    return True


if __name__ == "__main__":
    main()
