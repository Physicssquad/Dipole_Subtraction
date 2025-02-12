#FC = gfortran
#FFLAGS   = -O
#LIBS = 
###LIB = -lm 
#LIBS = 
#LIB = 
#OBJECTS = phasespace.o vsup.o integrand.o misc.o check-dipole.o
#PROG = run.x
#$(PROG): $(OBJECTS) main.f  
#	$(FC) main.f $(FFLAGS) -o $(PROG) $(OBJECTS) $(LIB) $(LIBS) 
#clean: 
#	rm -f *.o core *~ *.x
# Compiler and flags
FC = gfortran
FFLAGS = -O
LDFLAGS =
export LHFLAGS = $(shell lhapdf-config --ldflags)

# Directories
OBJDIR = objects

# Program and object files
PROG = run.x

OBJECTS = $(OBJDIR)/phasespace.o \
          $(OBJDIR)/vsup.o       \
          $(OBJDIR)/integrand.o  \
          $(OBJDIR)/misc.o       \
          $(OBJDIR)/check-dipole.o \
          $(OBJDIR)/lum.o        \
          $(OBJDIR)/pdf_lha.o

# Main target: Build executable from object files
$(PROG): $(OBJDIR) $(OBJECTS) main.f
	$(FC) main.f $(FFLAGS) -o $(PROG) $(OBJECTS) $(LHFLAGS) $(LDFLAGS)

# Create the object directory if not present.
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Compile .f files into .o files, ensuring they go into $(OBJDIR)
$(OBJDIR)/%.o: %.f | $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Clean target
clean:
	rm -f $(OBJDIR)/*.o core *~ $(PROG)

