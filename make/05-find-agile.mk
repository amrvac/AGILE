$(info AGILE Path: $(AGILE_DIR))

agile := $(AGILE_DIR)
build := $(shell pwd)/build

ifeq ($(findstring help, $(MAKECMDGOALS)), help)
disable_precompile := 1
endif

ifeq ($(findstring clean, $(MAKECMDGOALS)), clean)
disable_precompile := 1
endif

.PHONY: help

help:
	less $(agile)/make/README.md
