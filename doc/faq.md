title: Frequently Asked Questions

# List of questions {: #faq-list }
Below you can find a list of frequently asked questions and answers.

# Is there a mailinglist? {: #faq-mailinglist }
There is a mailing list for the related AMRVAC code that you can send questions to: <mailto:amrvacsers@ls.kuleuven.be>

You can join the mailing-list by [subscribing](https://ls.kuleuven.be/cgi-bin/wa?SUBED1=AMRVACUSERS&A=1), so that you will be informed about important changes or (bug)fixes. You can also search the [mailing list archive](https://ls.kuleuven.be/cgi-bin/wa?A0=AMRVACUSERS).

If you find issues or have feature requests specific to AGILE, its best to do this by raising a [github issue](https://github.com/amrvac/AGILE/issues).   

# Is there a way I can define my own parameters somewhere in mod_usr.t and configure them through `agile.par` ? {: #faq-own-parameters }
Indeed, there is a quick and time-saving way to read your own parameters without having to give an explicit value in the usr file and recompile each time. Instead, add this in your usr file : 

1. at the end of the usr_init subroutine, add this line : call params_read(par_files)
2. and just after the usr_init subroutine, define the params_read subroutine you just called. For instance :
  
```fortran
    !> Read parameters from a file
    subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /my_list/ my_parameter_1, my_parameter_2 ! where you tell the code to read your own parameters

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, my_list, end=111)
111    close(unitpar)
    end do

  end subroutine params_read
```

with `my_parameter_1` and `my_parameter_2` to be defined at the very beginning of the user file, just before the “contains” statement. Doing so, you can use them anywhere in the usr file and they will have the value you defined in your par file adding the following lines :

```fortran
&my_list
	my_parameter_1    = 0.2d0
	my_parameter_2    = 1.d0
/
```
# I'm looking for a way to translate some .blk files (using `convert_type= 'onegrid'`) back into standard output `.dat` that I can use to restart a simulation. My goal is to be able to change my data mid-simulation through another program. {: #faq-blk-to-dat }
The idea of the `onegrid` option for conversion is that a hierarchical block AMR grid (as stored in the `.dat` files) can be saved to an equivalent uniform grid representation at a user-chosen level (combine `level_io` with `onegrid`). You then can use any software you like to handle uniform grid data. Converting back to the .dat format is then impossible. However, you can write a user-routine to read in the uniform data, and then use it to restart a similar uniform-grid (no AMR, just domain decomposed) simulation. The `.dat` files are all that is needed to do restarts.

If you want to change the data during the simulation, in principle you do not need another program to do that. For that purpose, we provide the generic subroutines

1. usr_process_grid

2. usr_process_global

or immediately after the advance step:

3. usr_process_adv_grid

4. usr_process_adv_global

or just once, at restart (and after reading in a .dat file), then use

5. usr_before_main_loop

which allows you to modify things during runtime (of course, it should make physical sense). See their interface in _mod_usr_methods.t_.

# Is there a way to save (and visualize) the data in the boundary ghost cells? {: #faq-show-ghost-cells }
The logical variable `save_physical_boundary` can be set to true, which enforces the `.dat` files to also contain the ghost cell information for physical boundaries (i.e., those beyond the left or right edge of the domain, in each direction). You can use this file (like any other `.dat` file, to restart, and this helps if you want to use saved boundary info in your boundary value handling. However, all our present conversion options (like e.g. to `.vtu` files) do not store this extra info, and you can therefore not use them for visualizing the ghost cell info. For that, you will need to handle the `.dat` files directly, e.g. using python.

# I get an error: mpi.mod was created by a different version of GNU Fortran {: #faq-mpi-mod }
This means that your MPI library was compiled with a different version of
`gfortran` (and `GCC`) than you are using now. If you run on a cluster, contact
their support and ask them to fix it.

If you encounter this problem on your own machine, you can try to:

1. Change compilers: Use the version of gfortran that MPI was compiled with
2. Reinstall your MPI library, or install a different one (e.g. MPICH or OpenMPI), which is then hopefully compiled with the right version of GCC/gfortran.
3. Install an operating system with working package management, e.g. Debian ;)
