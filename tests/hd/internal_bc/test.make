TESTS := fixed.log orbit.log

include ../../test_rules.make

# Both runs share one build of ./agile; they differ only in their par file.
fixed.log: fixed.par
orbit.log: orbit.par
