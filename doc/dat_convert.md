title: Dump variables to dat files

There are two ways to add variables to the dat files:

## Extra variables
This method adds extra variables in the dat file (which is also used for restart).  

This method  consists in:
1. Assign an index in the variable list w in `usr_init` subroutine (in `mod_usr.t`) using `var_set_extravar`
2. Implement one of the user methods which are called every timestep: `usr_process_grid`, `usr_modify_output`,...
 and set there  the extra variables (and properly assign the pointer to the implemented subroutine  in `usr_init`).

.
## New dat files
With this method an extra dat file is generated with new variables. (see `io/mod_convert.fpp`). 
1. Set `convert_type` to `"dat_generic_mpi"`:

```fortran
&filelist
      ...
      convert_type = "dat_generic_mpi"
      ...
/
```

2. Add a convert method in `usr_init` subroutine (`use mod_convert`):

```fortran
    call add_convert_method2(dump_vars, 4, "jx jy jz sxr", "_aux_") 
```
or

```fortran
    call add_convert_method(dump_vars, 4, (/" jx", " jy", " jz", "sxr"/), "_aux_")
```

The first way is preferred as because Fortran requires the strings in an array to have the same length and extra spaces might be necessary for this.
The last argument `_aux_` is the suffix added to `base_filename`  in order to generate  the name of the new dat file.

3. Implement `dump_vars` function: 

```fortran
        function dump_vars(ixI^L, ixO^L, w, x, nwc) result(wnew)
          use mod_global_parameters
          use mod_thermal_emission
          integer, intent(in)             :: ixI^L,ixO^L, nwc
          double precision, intent(in)    :: w(ixO^S, 1:nw)
          double precision, intent(in)    :: x(ixO^S,1:ndim)
          double precision    :: wnew(ixO^S, 1:nwc)
          double precision :: current(ixI^S,7-2*ndir:3)
          integer :: idirmin,idir
          double precision    :: tmp(ixI^S)
          call get_current(w,ixI^L,ixO^L,idirmin,current)
          wnew(ixO^S,1) =  current(ixO^S,1)
          wnew(ixO^S,2) =  current(ixO^S,2)
          wnew(ixO^S,3) =  current(ixO^S,3)
          call get_SXR(ixI^L,ixO^L,w,x,te_fl_c,tmp,9,12)
          wnew(ixO^S,4) =  tmp(ixO^S)
        end function dump_vars
```

### Examples of use in the main code

1. For the mhd model for dumping full variables by setting in the parameter file:

```fortran
&mhd_list
    mhd_dump_full_vars=.true.
/  
```

A new file is created with hardcoded suffix `new`.
