TESTS := data.log icarus.log

include ../../test_rules.make

data.log: agile.par
icarus.log: icarus.par
