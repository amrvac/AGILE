# ZFP-compressed snapshots through the parallel writer (datfile_parallel_io,
# write_snapshot_parallel) must be byte-identical to the default writer's,
# with independent and with collective MPI-IO. The blast_amr case and the
# snapshot overlays are tests/hd/spherical's (variant a: single precision, 2
# ghost layers; b: double precision, physical-boundary ghost cells), with
# compression on. Each rank's adaptive precision search runs over its blocks
# in the same order on both paths, so the streams must match.
#
# Needs a ZFP build (ZFP_DIR, see doc/par.md#par_snapshots); CI has none, so
# this directory is in no suite and reports SKIPPED without ZFP_DIR.

ifdef ZFP_DIR
TESTS := snap_zfp.log
else
TESTS := skip
endif

include ../../test_rules.make

.PHONY: skip
skip:
	@echo -e "SKIPPED snap_zfp.log (ZFP_DIR not set)"

SPH := ../spherical
SNAP_VARIANTS := a b
SNAP_PATHS := serial parallel collective
snap_zfp.log: agile $(SPH)/blast_amr.par snap_zfp.par force
	@$(RM) $@ blast_amr_snap*
	@$(foreach v,$(SNAP_VARIANTS),$(foreach p,$(SNAP_PATHS),\
		$(LAUNCHER) ./agile -i $(SPH)/blast_amr.par $(SPH)/snap_$(v).par \
		snap_zfp.par $(SPH)/snap_$(p).par > run.log &&)) true
	@if $(foreach v,$(SNAP_VARIANTS),$(foreach i,0000 0001,\
		$(foreach p,parallel collective,$(FILE_CMP) \
		blast_amr_snap$(v)_serial_snap$(i).dat \
		blast_amr_snap$(v)_$(p)_snap$(i).dat &&))) true ; \
	then echo -e "$(_green)PASSED$(_reset) $@" ; touch $@ ; \
	else echo -e "$(_red)** FAILED **$(_reset) $@" ; \
	fi
