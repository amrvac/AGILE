TESTS := data.log icarus.log firstrow.log

include ../../test_rules.make

data.log: agile.par
icarus.log: icarus.par
firstrow.log: firstrow.par
