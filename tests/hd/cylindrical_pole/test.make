TESTS := uflow.log uflow_amr.log blast_amr.log

include ../../test_rules.make

# Generate dependency rules for the tests. All three share one build of
# ./agile, which is the point of keeping them in one directory: they agree on
# every compile-time parameter and differ only in their par file.
uflow.log: uflow.par
uflow_amr.log: uflow_amr.par
blast_amr.log: blast_amr.par
