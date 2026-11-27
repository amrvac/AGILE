TESTS := uflow.log blast.log blast_amr.log blast_amr_fixgrid.log pario.log

include ../../test_rules.make

# Generate dependency rules for the tests. They all share one build of
# ./agile, which is the point of keeping them in one directory: they agree on
# every compile-time parameter and differ only in their par file.
uflow.log: uflow.par
blast.log: blast.par
blast_amr.log: blast_amr.par
blast_amr_fixgrid.log: blast_amr_fixgrid.par

# The parallel restart-file path (datfile_parallel_io, mod_datfile_parallel)
# must write the same bytes as the default one. blast_amr gives an AMR mesh,
# and the overlays save physical-boundary ghost cells, so records differ in
# size. There is no reference log: the two runs' files are compared directly.
PARIO_DATS := 0000 0001
pario.log: agile blast_amr.par pario_serial.par pario_parallel.par force
	@$(RM) $@ blast_amr_pario_*.dat
	@$(LAUNCHER) ./agile -i blast_amr.par pario_serial.par > run.log
	@$(LAUNCHER) ./agile -i blast_amr.par pario_parallel.par > run.log
	@if $(foreach i,$(PARIO_DATS),cmp blast_amr_pario_serial$(i).dat \
		blast_amr_pario_parallel$(i).dat &&) true ; \
	then echo -e "$(_green)PASSED$(_reset) $@" ; touch $@ ; \
	else echo -e "$(_red)** FAILED **$(_reset) $@" ; \
	fi
