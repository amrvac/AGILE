
TESTS := test.log test-lowbeta.log

include ../../test_rules.make

# Generate dependency rules for the tests
test.log: test.par
test-lowbeta.log: test-lowbeta.par
