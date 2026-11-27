TESTS := test.log test-fixgrid.log hll.log restart.log

include ../../test_rules.make

# Generate dependency rules for the tests
test.log: test.par
test-fixgrid.log: test-fixgrid.par
hll.log: hll.par

# Restart round trip through the parallel restart-file path
# (datfile_parallel_io, mod_datfile_parallel). The run's log must match
# test.log's reference: the IO path does not change the solution. Resuming
# from its t = 0.05 file with either reader must then reproduce its t = 0.1
# file byte for byte, which also checks that the default path reads a file
# the parallel one wrote and writes the same bytes.
restart.log: $(LOG_CMP) agile test.par restart.par resume_parallel.par resume_serial.par force
	@$(RM) $@ test_restart*.dat
	@$(LAUNCHER) ./agile -i test.par restart.par > run.log
	@$(LAUNCHER) ./agile -i test.par restart.par resume_parallel.par \
		-if test_restart0000.dat > run.log
	@$(LAUNCHER) ./agile -i test.par restart.par resume_serial.par \
		-if test_restart0000.dat > run.log
	@if $(LOG_CMP) 1.0e-5 1.0e-8 test_restart.log correct_output/test.log && \
		cmp test_restart_parallel0001.dat test_restart0001.dat && \
		cmp test_restart_serial0001.dat test_restart0001.dat ; \
	then echo -e "$(_green)PASSED$(_reset) $@" ; cp test_restart.log $@ ; \
	else echo -e "$(_red)** FAILED **$(_reset) $@" ; \
	fi
