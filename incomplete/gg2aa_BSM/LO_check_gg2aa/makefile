# Compiler and flags
FC = gfortran
FFLAGS = -O -w 
LDFLAGS = -lc -L/home/niraj/1TB-Disc/Workspace-IITG/Packages/Install/Cuba-4.2.2 -lcuba
export LHFLAGS = $(shell lhapdf-config --ldflags)

# Directories
OBJDIR = objects
SRCDIR = ..
LOCALSRCDIR = .

# Program and object files
PROG = runLO

OBJECTS = $(OBJDIR)/filehandle.o \
          $(OBJDIR)/cuba.o       \
          $(OBJDIR)/mat.amp.o    \
          $(OBJDIR)/cuts.lo2.o   \
          $(OBJDIR)/cone.o       \
          $(OBJDIR)/printframe.o \
          $(OBJDIR)/misc.o       \
          $(OBJDIR)/lum.o        \
          $(OBJDIR)/pdf_lha.o    \
          $(OBJDIR)/vsup.o       \
          $(OBJDIR)/integrand.o  \
          $(OBJDIR)/phasespace.o \
          $(OBJDIR)/couplings.o

# Main target: Build executable from object files
$(PROG): $(OBJDIR) $(OBJECTS) $(LOCALSRCDIR)/main.f
	$(FC) $(LOCALSRCDIR)/main.f $(FFLAGS) -o $(PROG) $(OBJECTS) $(LIB) $(LHFLAGS) $(LDFLAGS)

# Create the object directory if not present.
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Ensure the object directory exists and compile .f files into .o files
$(OBJDIR)/%.o: $(SRCDIR)/%.f $(LOCALSRCDIR)/%.f| $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Pattern rule for all Fortran files in $(SRCDIR)
$(OBJDIR)/%.o: $(SRCDIR)/%.f | $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@

# For the specific files in $(SRCDIR)/slicing_files
$(OBJDIR)/%.o: $(SRCDIR)/slicing_files/%.f | $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@

# For files in the current directory 
$(OBJDIR)/%.o: %.f | $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Clean target
clean:
	rm -f $(OBJDIR)/*.o core *~ $(PROG)
