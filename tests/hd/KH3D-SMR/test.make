TESTS := test.log pc.log ssprk4.log

include ../../test_rules.make

# Generate dependency rules for the tests.  test.par is ssprk3; pc.par sets no
# time_stepper, so it also checks that the default is Predictor_Corrector.
# Every other time integrator runs in hd/KH3D-SMR-integrators (extended suite).
test.log: test.par
pc.log: pc.par
ssprk4.log: ssprk4.par
