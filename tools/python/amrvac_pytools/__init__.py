# Each reader is optional, so that e.g. the .dat tools work without vtk
# installed; import the submodule directly to see the underlying ImportError.
try:
    from .datfiles.reading.amrvac_reader import load_datfile
except ImportError:
    pass
try:
    from .vtkfiles.read import load_vtkfile
except ImportError:
    pass
