# dep_files := $(f90_files:$(build_dir)/f90/%.f90=$(build_dir)/dep/%.dep)

# TODO: is there alternative to  fortdepend without dependencies?
ifeq (, $(shell which fortdepend))
$(error "fortdepend not found. Check the readme, or pip install fortdepend.")
endif

ifdef USE_MPIWRAPPERS
  fortdepend_flags += -DUSE_MPIWRAPPERS
endif
ifdef OPENACC
  fortdepend_flags += -D_OPENACC
endif

# make sure that config is read first
ifdef CONFIG_READ
$(build_dir)/dependencies.mk: $(f90_files) $(build_dir)/f90/agile.h | $(build_dir)
	@echo "Regenerating depencies"
	@fortdepend $(fortdepend_flags) -f $(f90_files) -i mpi openacc -b $(build_dir)/obj -w -o $@

# Precompiling and dependency tracking is not needed if we're cleaning.
ifndef disable_precompile
include $(build_dir)/dependencies.mk
endif
endif

