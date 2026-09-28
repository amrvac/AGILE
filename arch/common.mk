ifdef OPENACC
ifdef OPENMP
$(error Both OpenACC and OpenMP are enabled, this is not supported)
endif
endif

ifdef NOGPUDIRECT
enable_nogpudirect = 1
ifndef OPENACC
ifndef OPENMP
$(warning Ignoring NOGPUDIRECT. No effect unless a GPU offload backend is enabled)
enable_nogpudirect = 0
endif
endif
ifeq ($(enable_nogpudirect),1)
$(info Disabling direct GPU-GPU copies)
f90_flags += -DNOGPUDIRECT
enabled += NOGPUDIRECT
endif
endif

ifdef USE_MPIWRAPPERS
$(info Enabling MPI wrappers)
f90_flags += -DUSE_MPIWRAPPERS
enabled += USE_MPIWRAPPERS
endif

