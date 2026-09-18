arch := llvm

compile = mpif90
f90_flags += -ffree-form -fimplicit-none -cpp

ifdef OPENMP
$(info Enabling OpenMP)
enabled += OPENMP
f90_flags += -fopenmp -fopenmp-version=52 --offload-arch=native
ifdef NOGPUDIRECT
$(info Disabling direct GPU-GPU copies)
enabled += NOGPUDIRECT
f90_flags += -DNOGPUDIRECT
endif
endif

link_flags += $(f90_flags)

