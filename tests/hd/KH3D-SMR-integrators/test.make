# Every explicit time integrator, and its coefficient variants, on the AMR
# case of hd/KH3D-SMR: the static refinement puts refluxing at every stage
# under test.  pc.par sets no time_stepper, so it also checks the default.
TESTS := ssprk3.log pc.log ssprk4.log fe.log ssprk2.log rk2alfa.log ssp32.log \
         ssp42.log ssprk5.log rk4.log rk3bt.log rk3nys.log

include ../../test_rules.make

# Generate dependency rules for the tests
ssprk3.log: ssprk3.par
pc.log: pc.par
ssprk4.log: ssprk4.par
fe.log: fe.par
ssprk2.log: ssprk2.par
rk2alfa.log: rk2alfa.par
ssp32.log: ssp32.par
ssp42.log: ssp42.par
ssprk5.log: ssprk5.par
rk4.log: rk4.par
rk3bt.log: rk3bt.par
rk3nys.log: rk3nys.par
