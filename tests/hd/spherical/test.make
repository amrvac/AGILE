TESTS := uflow.log blast.log blast_amr.log blast_amr_fixgrid.log

include ../../test_rules.make

# Generate dependency rules for the tests. All four share one build of
# ./agile, which is the point of keeping them in one directory: they agree on
# every compile-time parameter and differ only in their par file.
uflow.log: uflow.par
blast.log: blast.par
blast_amr.log: blast_amr.par
blast_amr_fixgrid.log: blast_amr_fixgrid.par
