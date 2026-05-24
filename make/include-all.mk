# Includes all numbered make files in order

include $(shell ls $(AGILE_DIR)/make/[0-9]*.mk)

