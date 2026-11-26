# Optional external libraries enabled through config flags.

ifneq ($(filter -DCOMPRESS_ZFP%,$(fypp_flags)),)

ifndef ZFP_DIR
$(error ZFP_DIR is not set; export ZFP_DIR=?/zfp, built with make BUILD_ZFORP=1)
endif

$(info Enabling ZFP compression, dir: $(ZFP_DIR))
f90_flags += -I$(ZFP_DIR)/modules
link_libs += $(ZFP_DIR)/lib/libzFORp.a $(ZFP_DIR)/lib/libzfp.a

endif
