# Makefile for automatically testing AGILE
#
# 1. You can run the hd tests with: make hd. The same applies to other physics
# modules
#
# 2. You can run them in parallel with for example: -j 4
#
# 3. You can hide the Entering/leaving directory messages with: -s

HD_DIRS := hd/KH3D hd/RT3D hd/KH3D-SMR hd/KH3D_tracer hd/sphere_advection hd/TRML hd/cloud_crushing \
           hd/spherical_uniform_flow hd/spherical_blast
MHD_DIRS := mhd/blast_wave_Cartesian_3D mhd/Orszag-Tang_3D mhd/spherical_uniform_flow mhd/spherical_blast
FFHD_DIRS := ffhd/TI3D ffhd/spherical_uniform_flow ffhd/spherical_blast
SRHD_DIRS := srhd/shockTube srhd/spherical_uniform_flow srhd/spherical_blast

TEST_DIRS := $(HD_DIRS) $(MHD_DIRS) $(FFHD_DIRS) $(SRHD_DIRS)

# Directories with altered names (useful for cleaning)
CLEAN_TEST := $(TEST_DIRS:%=clean-%)

.PHONY: $(TEST_DIRS) $(CLEAN_TEST)

all: hd mhd ffhd srhd
hd: $(HD_DIRS)
mhd: $(MHD_DIRS)
ffhd: $(FFHD_DIRS)
srhd: $(SRHD_DIRS)

$(TEST_DIRS):
	@$(MAKE) -C $@ -f test.make

clean:	$(CLEAN_TEST)

$(CLEAN_TEST):
	$(MAKE) -C $(@:clean-%=%) -f test.make clean
