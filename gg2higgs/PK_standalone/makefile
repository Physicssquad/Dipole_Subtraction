 include ../config.mk
# Program and object files
PROG = runPK

OBJECTS = $(OBJDIR)/filehandle.o \
          $(OBJDIR)/mat.amp.o    \
          $(OBJDIR)/printframe.o \
          $(OBJDIR)/misc.o       \
          $(OBJDIR)/lum.o        \
          $(OBJDIR)/pdf_lha.o    \
          $(OBJDIR)/vsup.o       \
          $(OBJDIR)/integrand.o  \
          $(OBJDIR)/ol_interface.o  \
          $(OBJDIR)/phasespace.o \
          $(OBJDIR)/PK_new.o 

# Add OpenLoops library to the linking process
OPENLOOPS_FLAGS = -Wl,-rpath=$(OPENLOOPS_LIB) -I$(OPENLOOPS_INCLUDE) -L$(OPENLOOPS_LIB) -lopenloops

# Main target: Build executable from object files
$(PROG): $(OBJDIR) $(OBJECTS) $(LOCALSRCDIR)/main.f
	$(FC) $(LOCALSRCDIR)/main.f $(FFLAGS) -o $(PROG) $(OBJECTS) $(LIB) $(LHFLAGS) $(LDFLAGS) $(OPENLOOPS_FLAGS)

# Create the object directory if not present.
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Ensure the object directory exists and compile .f files into .o files
$(OBJDIR)/%.o: $(SRCDIR)/%.f $(LOCALSRCDIR)/%.f| $(OBJDIR)
	$(FC) $(FFLAGS) -I$(OPENLOOPS_INCLUDE) -c $< -o $@

# Pattern rule for all Fortran files in $(SRCDIR)
$(OBJDIR)/%.o: $(SRCDIR)/%.f | $(OBJDIR)
	$(FC) $(FFLAGS) -I$(OPENLOOPS_INCLUDE) -c $< -o $@

# For files in the current directory 
$(OBJDIR)/%.o: %.f | $(OBJDIR)
	$(FC) $(FFLAGS) -I$(OPENLOOPS_INCLUDE) -c $< -o $@

# Clean target
clean:
	rm -f $(OBJDIR)/*.o core *~ $(PROG)
