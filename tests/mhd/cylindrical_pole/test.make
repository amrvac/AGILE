TESTS := uflow.log uflow_hllc.log uflow_hlld.log

include ../../test_rules.make

# Generate dependency rules for the tests
uflow.log: uflow.par
uflow_hllc.log: uflow_hllc.par
uflow_hlld.log: uflow_hlld.par
