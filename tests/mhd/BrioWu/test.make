TESTS := tvdlf.log hll.log hllc.log hlld.log

include ../../test_rules.make

# Generate dependency rules for the tests
tvdlf.log: tvdlf.par
hll.log: hll.par
hllc.log: hllc.par
hlld.log: hlld.par
