# The following files need to be removed before running anything else:
#
# - config.mk: stores config generated from agile.par, or another parameter
#   file if given as CONFIG=... on the command line. If this isn't removed,
#   make cannot detect a change in $CONFIG. The solution is to always rebuild
#   config.mk, but since config.mk is also included, its easiest to just
#   remove it.
#
# - agile: this is the executable symlink and suffers the same problem as
#   config.mk. Make doesn't register the change when we're simply "switching
#   config branches". Again, the solution is to remove it before building.
#   Side effect: agile is also removed on failed builds, which is a good thing.

ifeq ($(MAKE_RESTARTS),)
$(shell rm -f config.mk agile)
endif
