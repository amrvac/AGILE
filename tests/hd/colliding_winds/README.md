# Colliding winds in a binary

An example setup, not part of the test suite: two stars on a circular orbit,
each blowing an isotropic, supersonic (Mach 20) wind, whose collision forms
a hot shocked layer that the orbit winds into a spiral. The physics and the
parameters (`&usr_list`) are documented at the top of `mod_usr.fpp`.

It is the worked example of the GPU internal boundary (`usr_internal_bc`,
`internalboundary = T` in `&boundlist`): the hook holds each star's source
region at its wind state as the stars move, and floors the temperature
everywhere else. The source regions are kept at the finest level by
`usr_refine_grid`; tracer 1 marks the gas of wind 1.

| file | |
|---|---|
| `agile.par` | build reference and production run: 16 x 16 x 8 box, two orbits |
| `quick.par` | smoke test on the same build: 4 x 4 x 2 box, t = 2 |
| `plot_cwb.py` | z = 0 plane: one PNG per slice, or an animation with `--movie` |

```bash
make arch=nvidia OPENACC=1 NOGPUDIRECT=1        # or arch=gnu for the CPU
mpirun -np 1 ./agile -i quick.par               # ~80 s on a GTX 1660 Ti
python plot_cwb.py                              # PNGs of output/*_d3_*.dat
```

One orbit of the small box, with a slice every 0.05 (241 frames), as a movie.
A second par file overrides the first; note that its `base_filename` is
*appended* to the first one's, so it gives a suffix:

```bash
cat > orbit.par <<'PAR'
 &filelist  base_filename = '_orbit' /
 &savelist  dtsave_slice = 0.05d0 /
 &stoplist  time_max = 12.0d0 /
PAR
mpirun -np 1 ./agile -i quick.par orbit.par     # ~16 min on a GTX 1660 Ti
python plot_cwb.py --movie output/cwb_orbit.mp4 output/cwb_quick_orbit_d3_*.dat
```

`plot_cwb.py` needs the analysis extra (`uv sync --extra analysis`), whose
`imageio-ffmpeg` supplies ffmpeg for the mp4.
