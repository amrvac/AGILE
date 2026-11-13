
TESTS := test.log test-fixgrid.log

include ../../test_rules.make

# Generate dependency rules for the tests
test.log: test.par
test-fixgrid.log: test-fixgrid.par
