TESTS := mixed.log primitive.log

include ../../test_rules.make

# Both runs share one build of ./agile; they differ only in their par file.
mixed.log: agile.par
primitive.log: primitive.par
