# library name
lib.name = cicmtools

# Sources:
ambicube~.class.sources := src/ambicube~.c
ambipan~.class.sources := src/ambipan~.c
vbapan~.class.sources := src/vbapan~.c

# extra files
datafiles = \
$(wildcard help/*.pd) \
license.txt \
README.md

# include Makefile.pdlibbuilder from submodule directory 'pd-lib-builder'
PDLIBBUILDER_DIR=pd-lib-builder/
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
