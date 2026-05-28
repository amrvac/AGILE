ifdef DEBUG
fypp_flags += -DDEBUG
endif

ifeq (, $(shell which fypp))
$(error "fypp not found. Check the readme, or pip install fypp.")
endif

fypp_flags += -n

fypp_flags += -M $(agile)/make

$(info Fypp flags: $(fypp_flags))

source_files := $(shell find $(agile)/src -name '*.fpp') mod_usr.fpp
f90_files := $(patsubst %.fpp, $(build_dir)/f90/%.f90, \
	     $(notdir $(source_files)))

.PRECIOUS: $(f90_files)

define fypp_rule
-include $(patsubst %.fpp, $(build_dir)/f90/%.d, $(notdir $(1)))

$(patsubst %.fpp, $(build_dir)/f90/%.f90, $(notdir $(1))): $(1)
	@mkdir -p $$(@D)
	@echo -e "Precompiling $(_magenta)$$(<:$(agile)/src/%=%)$(_reset)"
	@python $(agile)/make/fypp-deps.py $$< $$@
	@fypp $(fypp_flags) $$<  $$@
endef

$(foreach f, $(source_files), $(eval $(call fypp_rule, $(f))))

