# TODO: link to agile.a for later linkage with mod_usr.o

agile: $(build_dir)/obj/agile
	@rm -f agile
	@ln -s $< $@

$(build_dir)/obj/agile:
	@echo -e "Linking $(_green)$(notdir $@)$(_reset)"
	@$(link) $(link_flags) $^ -o $@

