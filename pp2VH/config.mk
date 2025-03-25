#This is a common config file for the makefile in each PROG
# Modify accordingly [ Last modified by Niraj ]

# Compiler and flags
FC = gfortran
FFLAGS = -O -w
LDFLAGS = -lc -L$(CUBA_DIR) -lcuba
export LHFLAGS = $(shell lhapdf-config --ldflags)

# OpenLoops paths (modify for different machines)
OPENLOOPS_DIR =/home/niraj/1TB-Disc/Workspace-IITG/Packages/Install/OpenLoops
OPENLOOPS_INCLUDE = $(OPENLOOPS_DIR)/lib_src/openloops/mod
OPENLOOPS_LIB = $(OPENLOOPS_DIR)/lib

# Cuba paths
CUBA_DIR = /home/niraj/1TB-Disc/Workspace-IITG/Packages/Install/Cuba-4.2.2 
#... currently we are not using cuba path is not mandatory. You can leave empty space.

# Directories
OBJDIR = objects
SRCDIR = ..
LOCALSRCDIR = .
INCLUDEDIR = ../include/ 
