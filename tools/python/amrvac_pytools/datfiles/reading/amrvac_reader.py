"""
Main class to load in MPI-AMRVAC native .dat files.

@author: Niels Claes    (niels.claes@kuleuven.be)
         Clément Robert (clement.robert@oca.eu)
         Adrian Kelly
         Olaf Willocx
Last edit: 1 August 2026
"""

import sys, os
import numpy as np
import copy

# Allow this file to run directly, for convenience. See main function at the
# bottom.
if __name__ == "__main__":
    sys.path.insert(0, os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "..", "..")))

from amrvac_pytools.datfiles.reading import datfile_utilities
from amrvac_pytools.datfiles.processing import process_data
from amrvac_pytools.datfiles.physics import physical_constants


class load_datfile():
    # Following methods provide easy access to functionalities in other scripts
    def amrplot(self, var, **kwargs):
        from amrvac_pytools.datfiles.plotting import amrvac_plotter
        return amrvac_plotter.amrplot(self, var, **kwargs)
    def rgplot(self, data, **kwargs):
        from amrvac_pytools.datfiles.plotting import amrvac_plotter
        return amrvac_plotter.rgplot(self, data, **kwargs)
    def halpha(self, **kwargs):
        from amrvac_pytools.datfiles.views import synthetic
        return synthetic.h_alpha(self, **kwargs)
    def faraday(self, **kwargs):
        from amrvac_pytools.datfiles.views import synthetic
        return synthetic.faraday(self, **kwargs)
    @staticmethod
    def show():
        import matplotlib.pyplot as plt
        plt.show()

    def __init__(self, filename):
        try:
            file = open(filename, "rb")
        except:
            raise FileNotFoundError("Unable to open '{}', is this an MPI-AMRVAC .dat file?".format(filename))

        self._filename = filename
        # Trim filename, assumes format .../.../filenameXXXX.datfiles
        for c in range(len(filename) - 1, -1, -1):
            if filename[c] == '/' or filename[c] == '\\':
                self._filename = filename[c+1:]
                break
        print(">> Reading {}".format(filename))

        # setup basic attributes
        self.file = file
        self.header = datfile_utilities.get_header(file)
        self.data_dict = None
        self._uniform = self._data_is_uniform()
        self._regriddir = "regridded_files"

        # load field information, including primitive variables
        self.known_fields = self._get_known_fields()

        # check if dataset is isothermal
        self.isothermal = self._data_is_isothermal()
        self.adiab_constant = 1

        # load blocktree information
        self.block_lvls, self.block_ixs, self.block_offsets, self.block_nghost, self.block_field_bytes \
          = datfile_utilities.get_tree_info(file)
        self.block_shape = np.append(self.header["block_nx"], self.header["nw"])
        self.domain_width = self.header["xmax"] - self.header["xmin"]

        # setup units
        self.units = physical_constants.units(self.header)

    def get_info(self):
        """
        Prints the current snapshot info (obtained from the header) to the console.
        """
        if self.isothermal:
            info_isothermal = ' (isothermal) '
        else:
            info_isothermal = ''
        print("")
        print(f"[INFO] Current file      : {self._filename}")
        print(f"[INFO] Datfile version   : {self.header["datfile_version"]}")
        print(f"[INFO] Current time      : {format(self.header["time"])}")
        print(f"[INFO] Current iteration : {self.header["it"]}")
        print(f"[INFO] Physics type      : {self.header["physics_type"] + info_isothermal}")
        print(f"[INFO] Boundaries        : {self.header['xmin']} -> {self.header['xmax']}")
        print(f"[INFO] Max AMR level     : {self.header["levmax"]}")
        print(f"[INFO] Uniform mesh      : {self._uniform}")
        print(f"[INFO] Block size        : {self.header["block_nx"]}")
        print(f"[INFO] Number of blocks  : {len(self.block_offsets)}")
        print(f"[INFO] Real kind         : {self.header.get("size_real", 8)}")
        if self.header.get("compression", "none") != "none":
            print(f"[INFO] Compression       : {self.header["compression"]}")
            match self.header["compression"]:
                case "zfp":
                    print(f"[INFO] ZFP precision     : {self.header["zfp_precision"]}")
        print("-" * 40)
        print(f"Known variables: {self.known_fields}")

    def get_bounds(self):
        """
        Get the bounds of the calculation domain.
        :return: A list containing the bounds in the form [[x1, x2], [y1, y2], [z1, z2]].
                 Depending on the dimension the list size is either 1, 2 or 3.
        """
        bounds = []
        for i in range(len(self.header["xmin"])):
            bounds.append([self.header["xmin"][i], self.header["xmax"][i]])
        return bounds

    def get_coordinate_arrays(self):
        """
        Calculates the domain discretisation.
        :return: A list containing the coordinate axes in the form [x, y, z] where
                 eg. 'x' is a NumPy array representing the domain discretization along the x-axis.
        """
        self._check_datadict_exists()

        w_names = self.header["w_names"]
        xmax = self.header["xmax"]
        xmin = self.header["xmin"]
        nx = np.asarray(self.data_dict[w_names[0]].shape)

        coordinate_arrays = []
        for i in range(len(nx)):
            coordinate_arrays.append(np.linspace(xmin[i], xmax[i], nx[i]))

        return np.asarray(coordinate_arrays, dtype=object)

    def load_all_data(self, nbprocs=None, regriddir=None, upsample=False):
        """
        Loads in all the data to RAM, this can take up quite some space for huge datasets.
        Data will be regridded if this is not already done.
        :param upsample: if True, non-uniform data is put on the finest level's grid by repeating each coarser
                         cell (piecewise constant, see datfile_utilities.get_upsampled_data) instead of being
                         interpolated by the regridding; fast, and nothing is written to regriddir.
                         Has no effect on uniform data, nor once the data has been loaded.
        """
        if self.data_dict is not None:
            return self.data_dict
        if self._uniform:
            data = datfile_utilities.get_uniform_data(self)
        elif upsample:
            data = datfile_utilities.get_upsampled_data(self.file)
        else:
            data = self._regrid_data(nbprocs, regriddir)
        # create dictionary containing the data
        self.data_dict = process_data.create_data_dict(data, self.header, self.adiab_constant)
        return self.data_dict

    def get_time(self):
        """
        Gets the current snapshot time.
        :return: Current time in the simulation (in AMRVAC units).
        """
        return self.header["time"]

    def get_extrema(self, var, print_result=False):
        """
        Returns the maximum and minimum value of the requested variable
        :param var: variable, see ds.known_fields which are available
        :param print_result: if True, prints the extrema to the console
        :return: min, max of the given variable
        """
        if var not in self.known_fields:
            raise KeyError("variable not known: {}".format(var))
        varmax = -1e99
        varmin = 1e99
        size_real = self.header.get("size_real", 8)
        for ileaf, offset in enumerate(self.block_offsets):
            nghost = self.block_nghost[ileaf]
            shape = datfile_utilities.block_shape_with_ghost(self.header, nghost)
            block = datfile_utilities.get_single_block_data(self.file, offset, shape, size_real)
            block = datfile_utilities.strip_ghost(block, self.header, nghost)
            block = process_data.create_data_dict(block, self.header)
            varmax = np.maximum(varmax, np.max(block[var]))
            varmin = np.minimum(varmin, np.min(block[var]))
        if print_result:
            print(">> minimum {}: {:2.3e}   |   maximum {}: {:2.3e}".format(var, varmin, var, varmax))
        return varmin, varmax

    def set_adiabatic_constant(self, adiab_cte):
        if not self.isothermal:
            print("[INFO] dataset is not isothermal, 'e' is present. Setting the adiabatic constant has no effect.")
        self.adiab_constant = adiab_cte




    # 'protected' methods
    def _data_is_uniform(self):
        """
        Checks if the data is uniformely refined.
        :return: True if grid is uniform, False otherwise.
        """
        refined_nx = 2 ** (self.header['levmax'] - 1) * self.header['domain_nx']
        nleafs_uniform = np.prod(refined_nx / self.header['block_nx'])
        if not self.header["nleafs"] == nleafs_uniform:
            return False
        return True

    def _data_is_isothermal(self):
        """
        Checks if the data is isothermal
        :return: True if hd or mhd but no energy present
        """
        if (not 'e' in self.known_fields) and \
           (self.header['physics_type'] == 'hd' or self.header['physics_type'] == 'mhd'):
            return True
        else:
            return False

    def _check_regrid_directory(self, regriddir):
        """
        Checks if the specified save directory is present. Defaults to 'regridded_files', if this folder is not
        present it is created. If 'regriddir' is explicitly specified and the folder is not present, throws an error
        and terminates the program (to prevent creation of unwanted folders).
        :param regriddir: The directory to which to save the regridded files (and load them back in).
        """
        if regriddir is None:
            regriddir = self._regriddir
            if not os.path.isdir(regriddir):
                os.mkdir(regriddir)
                print("[INFO] Created directory: {}".format(regriddir))
        else:
            if not os.path.isdir(regriddir):
                sys.exit("Specified directory does not exist: {}".format(regriddir))
            self._regriddir = regriddir

    def _regrid_data(self, nbprocs, regriddir):
        """
        Regrids non-uniform data, which is saved as a Numpy file in the folder 'regridded_files'
        (created if not present). The regrid directory is searched first if the data is already present.
        If so it is loaded, if not regridding is started.
        :param nbprocs: Number of processors to do the regridding. Defaults to max number - 2
        :param regriddir: The directory to which the regridded file is saved. This directory is searched in order to
                          try and load the file if already present.
        """
        from amrvac_pytools.datfiles.processing import regridding
        self._check_regrid_directory(regriddir)
        try:
            data = self._load_regridded_data()
        except FileNotFoundError:
            data = regridding.regrid_amr_data(self, nbprocs)
            self._save_regridded_data(data)
        return data

    def _save_regridded_data(self, data):
        """
        Saves the regridded data as a NumPy file.
        :param data: The raw, regridded data.
        """
        filename_regridded = self._regriddir + "/" + self._filename[:-4] + "_regridded.npy"
        np.save(filename_regridded, data)
        print("[INFO] Regridded data saved to {}".format(filename_regridded))

    def _load_regridded_data(self):
        """
        Loads the regridded data file back in, if present.
        :return: The data file if it is present.
        :exception: FileNotFoundError if the file is not present. The dataset will be regridded in this case.
        """
        filename_regridded = self._regriddir + "/" + self._filename[:-4] + "_regridded.npy"

        if os.path.isfile(filename_regridded):
            try:
                data = np.load(filename_regridded)
                print("[INFO] Regridded file found and loaded ({}).".format(filename_regridded))
                return data
            except OSError:
                print("[ERROR] File found but failed to load: {}".format(filename_regridded))
                sys.exit(1)
        else:
            raise FileNotFoundError

    def _check_datadict_exists(self):
        """
        Checks if the data dictionary is loaded.
        :raises AttributeError  if self.data_dict is None
        """
        if self.data_dict is None:
            print("[INFO] Dataset must be loaded to do this, call load_all_data() first.")
            raise AttributeError

    def _get_known_fields(self):
        """
        Retrieves the known fields for a given dataset.
        :return: Regular list with the known fields as strings
        """
        fields = copy.deepcopy(self.header['w_names'])
        if self.header['physics_type'] == 'hd' or self.header['physics_type'] == 'mhd':
            if self.header['ndim'] == 1:
                fields += ['v1', 'p']
            elif self.header['ndim'] == 2:
                fields += ['v1', 'v2', 'p']
            else:
                fields += ['v1', 'v2', 'v3', 'p']
            if 'e' in fields:
                fields += ['T']
        return fields


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Read .dat file info: the header, and summary of given blocks.")
    parser.add_argument("datfile",
        help=".dat file path.")
    parser.add_argument("-l", "--leaves", type=int, nargs="+", default=[],
        help="Leaf indices to summarize, e.g. -l 0 1")
    parser.add_argument("-g", "--keep-ghost", action="store_true",
        help="Keep ghost cells in the leaf summaries")
    args = parser.parse_args()

    ds = load_datfile(args.datfile)
    ds.get_info()

    ndim = ds.header["ndim"]
    for ileaf in args.leaves:
        if not 0 <= ileaf < ds.header["nleafs"]:
            print("[ERROR] leaf index {} out of range 0..{}".format(ileaf, ds.header["nleafs"]-1))
            continue
        nghost = ds.block_nghost[ileaf]
        block = datfile_utilities.read_block(ds, ileaf, keep_ghost=args.keep_ghost)
        print("-" * 40)
        print(f"leaf {ileaf} | level {ds.block_lvls[ileaf]} | spatial index {ds.block_ixs[ileaf]}")
        print(f"  byte offset {ds.block_offsets[ileaf]} | nghost lo {nghost[:ndim]} hi {nghost[ndim:]} | shape {block.shape}")
        for iw, name in enumerate(ds.header["w_names"]):
            w = block[..., iw]
            print("  {:<16s} min {: .4e}   med {: .4e}   max {: .4e}".format(
                name, np.min(w), np.median(w), np.max(w)))
