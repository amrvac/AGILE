include $(agile)/arch/common.mk

arch := nvidia

compile = mpif90
f90_flags += -cpp -Mfree

ifdef NVTX
$(info Enabling NVTX)
enabled += NVTX
f90_flags += -DNVTX
link_flags += -lnvToolsExt
endif

ifdef OPENACC
$(info Enabling OpenACC)
enabled += OPENACC
GPU_ARCH ?= ccnative
f90_flags += -Wall -acc=gpu -gpu=$(GPU_ARCH)
ifdef DEBUG
f90_flags += -gpu=debug -Mvect=levels:0 -Mnoinline
else
f90_flags += -Mvect=levels:5 -Minline
endif
endif

ifdef OPENMP
$(info Enabling OpenMP)
enabled += OPENMP
GPU_ARCH ?= ccnative
f90_flags += -Wall -mp=gpu -gpu=$(GPU_ARCH)
ifdef DEBUG
f90_flags += -gpu=debug -Mvect=levels:0 -Mnoinline
else
f90_flags += -Mvect=levels:5 -Minline
endif
endif

ifdef INFO
f90_flags += -Minfo=all
endif

ifdef DEBUG
$(info Enable debugging symbols)
enabled += DEBUG
f90_flags += -Wall -Mbounds -g -O1 -traceback
else
f90_flags += -g -O3 -fast
endif

link_flags += $(f90_flags)
