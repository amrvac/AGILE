TESTS := uflow.log blast.log blast_amr.log

include ../../test_rules.make

# Generate dependency rules for the tests. All three share one build of
# ./agile: they agree on every compile-time parameter and differ only in their
# par file.
uflow.log: uflow.par
blast.log: blast.par
blast_amr.log: blast_amr.par
