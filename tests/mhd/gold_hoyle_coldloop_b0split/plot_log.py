"""Plot the diagnostics of one or more Gold-Hoyle cold-loop 'special' logs (with Tmin, Tmax).

Usage: python plot_log.py out.png label1=path/to/a.log [label2=path/to/b.log ...]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

cols = "it t dt Ekin Emag vmax dB_L2 drho_max sigmaJ jxb_max divb_mean divb_max Tmin Tmax".split()
panels = [("vmax", "max |v| / v_A", "log"), ("Ekin", "kinetic energy", "log"),
          ("dB_L2", "L2 error of B (relative)", "log"), ("Emag", "magnetic energy", "linear"),
          ("sigmaJ", "sigma_J = <|jxB|>/<|j||B|>", "log"), ("jxb_max", "max |j x B|", "log"),
          ("drho_max", "max |rho - rho0|/rho0", "log"), ("Tminmax", "Tmin, Tmax / T0", "linear")]

out = sys.argv[1]
runs = [a.split("=", 1) for a in sys.argv[2:]]
fig, axes = plt.subplots(2, 4, figsize=(16, 7))
for label, path in runs:
    d = np.loadtxt(path)
    dd = {c: d[:, i] for i, c in enumerate(cols)}
    for ax, (c, title, scale) in zip(axes.flat, panels):
        if c == "Tminmax":
            t0 = dd["Tmax"][0]
            ax.plot(dd["t"], dd["Tmin"] / t0, label=label + " min")
            ax.plot(dd["t"], dd["Tmax"] / t0, "--", label=label + " max")
            ax.set_title(title); ax.set_xlabel("t / tau_A"); ax.legend(fontsize=8)
            continue
        y = dd[c]
        if scale == "log":
            y = np.where(y > 0, y, np.nan)
        ax.plot(dd["t"], y, label=label)
        ax.set_title(title)
        ax.set_yscale(scale)
        ax.set_xlabel("t / tau_A")
axes.flat[0].legend()
fig.tight_layout()
fig.savefig(out, dpi=110)
print("wrote", out)
