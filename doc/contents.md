title: Code Overview

# Introduction for new users {: #intro_new_users }
* [Getting Started](getting_started.html) How to install AGILE and run your first test problem.
* [Acknowledgments](acknowledgments.html) Information on collaboration and
financial support.
* [FAQ](faq.html) Frequently asked questions.
* [Examples](examples.html) An overview of documented example simulations.
* [Contributing to AGILE and its documentation](contributing.html)
* [Documentation](documentation.html) How the documentation works.

# General {: #general }
* [Equations](equations.html) The equations and parameters in physics modules.
* [Coordinate systems](coordinate_systems.html) Cartesian, spherical,
  cylindrical and logarithmic-radius geometries, the geometric source terms
  each one adds, and the polar axis.
* [Parameters](par.html) Description of all parameters in "amrvac.par" parameter file.
* [Auxiliary variables](mpiamrvac_nw.html) Description of the intended use
  for _nw, nwflux, nwaux, nwextra, nwauxio_ parameters.
* [Command line](commandline.html) Help on command-line parameters.

# Discretization methods and AMR strategy {: #discretization }
* [Discretization](discretization.html) The equation and its discretization, the
basic variables, the structure of the grid, boundaries, ghost cells.
* [AMR aspects](amrstructure.html) Some essential info on global parameters and
the data structures for the block-tree AMR.

# Additional Physics {: #special_sources }
* [Radiative cooling](radiative_cooling.html) Description of adding radiative cooling.

# IO and data analysis {: #io }
* [File format](fileformat.html) Description of the format of AGILE data file (.dat).
* [Visualisation and analysis in Python with yt](yt_usage.html) without data conversion.
* [Converting data files for visualization](convert.html) Brief notes on how to
convert to Tecplot (.plt), and VTK (.vtu) data files.
* [Slices](slices.html) How to output hypersurfaces (slices) for restart or
visualization.
* [Line of sight views](collapsed.html) How to output collapsed views for
visualisation and analysis (e.g. column densities).
* [Analysis routine](analysis.html) Using the run-time analysis routine.
* [3D Printing](print3D.html) A note on how to generate 3D printed results.
