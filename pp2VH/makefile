# List of all program directories
SUBDIRS :=  dipole_standalone/ LO_gg2H/ virtual_standalone/ PK_standalone/

# Default target: Build all programs
.PHONY: all $(SUBDIRS) clean

all: $(SUBDIRS)

# Run `make` inside each subdirectory
$(SUBDIRS):
	$(MAKE) -C $@

# Clean all subdirectories
clean:
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir clean; done


