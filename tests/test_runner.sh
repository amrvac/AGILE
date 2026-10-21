# this is part of the continuous integration pipeline
# call syntax:
# [SUITE=basic|extended|all] bash test_runner.sh [<name>]
#
# Run tests in $1 (default: all) and check if any fails
# return appropriate exit code (0 if and only if all tests pass)
# SUITE is read by tests/Makefile from the environment (default: basic)

if [ -z "$AGILE_DIR" ]; then
    echo "Error: AGILE_DIR is not set"
    exit 1
fi

source $AGILE_DIR/.venv/bin/activate

target=${1:-all}
echo "Running test target '$target' with SUITE=${SUITE:-basic}"

make $target | tee $target.out; test ${PIPESTATUS[0]} -eq 0
cond1=$? # first condition: build correctly with make

if [ $cond1 != 0 ] ; then
    echo "BUILD ERROR: make $target failed to build the executable"
    exit 1
fi

if [ $(uname) == Darwin ] ; then
    ! ggrep "**FAILED**" $target.out
else
    ! grep "**FAILED**" $target.out
fi
cond2=$? # second condition: no error detected while running tests

if [ $cond2 != 0 ] ; then
    echo "TEST ERROR: at least one $target test failed"
    exit 1
fi
