TESTS := test.log test-lowbeta.log test-hllc.log

include ../../test_rules.make

# Generate dependency rules for the tests
test.log: test.par
test-lowbeta.log: test-lowbeta.par
test-hllc.log: test-hllc.par
