# Compiler and flags
FC = gfortran
FFLAGS = -O -w
LDFLAGS = -lc -L/home/niraj/1TB-Disc/Workspace-IITG/Packages/Install/Cuba-4.2.2 -lcuba
export LHFLAGS = $(shell lhapdf-config --ldflags)

# OpenLoops paths
OPENLOOPS_DIR = /home/niraj/1TB-Disc/Workspace-IITG/Packages/Install/OpenLoops
OPENLOOPS_INCLUDE = $(OPENLOOPS_DIR)/examples/../lib_src/openloops/mod
OPENLOOPS_LIB = $(OPENLOOPS_DIR)/examples/../lib

# Directories
OBJDIR = objects
SRCDIR = ..
LOCALSRCDIR = .

# Program and object files
PROG = runLO

OBJECTS = $(OBJDIR)/filehandle.o \
          $(OBJDIR)/mat.amp.o    \
          $(OBJDIR)/printframe.o \
          $(OBJDIR)/misc.o       \
          $(OBJDIR)/lum.o        \
          $(OBJDIR)/cone.o       \
          $(OBJDIR)/pdf_lha.o    \
          $(OBJDIR)/vsup.o       \
          $(OBJDIR)/integrand.o  \
          $(OBJDIR)/ol_interface.o  \
          $(OBJDIR)/phasespace.o

# Add OpenLoops library to the linking process
OPENLOOPS_FLAGS = -Wl,-rpath=$(OPENLOOPS_LIB) -I$(OPENLOOPS_INCLUDE) -L$(OPENLOOPS_LIB) -lopenloops

# Main target: Build executable from object files
$(PROG): $(OBJDIR) $(OBJECTS) $(LOCALSRCDIR)/main.f
	$(FC) $(LOCALSRCDIR)/main.f $(FFLAGS) -o $(PROG) $(OBJECTS) $(LIB) $(LHFLAGS) $(LDFLAGS) $(OPENLOOPS_FLAGS)

# Create the object directory if not present.
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Ensure the object directory exists and compile .f files into .o files
$(OBJDIR)/%.o: $(SRCDIR)/%.f $(LOCALSRCDIR)/%.f | $(OBJDIR)
	$(FC) $(FFLAGS) -I$(OPENLOOPS_INCLUDE) -c $< -o $@

# Pattern rule for all Fortran files in $(SRCDIR)
$(OBJDIR)/%.o: $(SRCDIR)/%.f | $(OBJDIR)
	$(FC) $(FFLAGS) -I$(OPENLOOPS_INCLUDE) -c $< -o $@

# For the specific files in $(SRCDIR)/slicing_files
$(OBJDIR)/%.o: $(SRCDIR)/slicing_files/%.f | $(OBJDIR)
	$(FC) $(FFLAGS) -I$(OPENLOOPS_INCLUDE) -c $< -o $@

# For files in the current directory
$(OBJDIR)/%.o: %.f | $(OBJDIR)
	$(FC) $(FFLAGS) -I$(OPENLOOPS_INCLUDE) -c $< -o $@

# Clean target
clean:
	rm -f $(OBJDIR)/*.o core *~ $(PROG)

