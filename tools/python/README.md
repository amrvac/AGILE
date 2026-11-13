# amrvac_pytools
## Installation

This python package is regularly updated, so we recommend installing it with "develop" arguments so that it does not need to be reinstalled with each update. This can be achieved in 3 distinct ways
> *nb:* all methods require `cd amrvac_pytools-project` is executed beforehand

- **recommended:** using `conda`
```bash
conda install conda-build
conda develop .
```
- using `pip`
```bash
pip install -e .
```
- or simply vanilla python:
```bash
python setup.py develop
```

## usage

```python
import amrvac_pytools
```

### load data
```python
df = amrvac_pytools.load_vtkfile("path/to/vtk/file")
# or
df = amrvac_pytools.load_datfile("path/to/dat/file")
```

### raw arrays from a .dat file (uniform or AMR, 1D/2D/3D, also slices)
```python
from amrvac_pytools.datfiles.reading import datfile_utilities as dfu

with open("path/to/dat/file", "rb") as f:
    header = dfu.get_header(f)          # time, w_names, xmin, xmax, levmax, ...
    w = dfu.get_upsampled_data(f)       # shape (*domain_nx * 2**(levmax-1), nw)
rho = w[..., header["w_names"].index("rho")]
```
`get_upsampled_data` puts coarser AMR blocks on the finest level's grid,
piecewise-constant; `get_uniform_data` assumes a single level and misplaces
multi-level data. `.dat` files hold conserved variables.

The same upsampling is available from `load_datfile`, instead of its default
(interpolating, slower, cached to `regridded_files/`) regridding of AMR data:
```python
ds = amrvac_pytools.load_datfile("path/to/dat/file")
data = ds.load_all_data(upsample=True)   # dict of fields, e.g. data["rho"]
```

### writeme