include $(agile)/arch/common.mk

arch := llvm

compile = mpif90
f90_flags += -ffree-form -fimplicit-none -cpp

ifdef OPENMP
$(info Enabling OpenMP)
enabled += OPENMP
GPU_ARCH ?= native
f90_flags += -fopenmp -fopenmp-version=52 --offload-arch=$(GPU_ARCH)
endif

link_flags += $(f90_flags)

