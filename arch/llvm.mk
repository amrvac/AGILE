arch := llvm

compile = flang
f90_flags += -ffree-form -fimplicit-none -cpp $(shell mpifort --showme:compile)

ifdef OPENMP
$(info Enabling OpenMP)
enabled += OPENMP
f90_flags += -fopenmp
endif

link_flags += $(f90_flags) $(shell mpifort --showme:link)

