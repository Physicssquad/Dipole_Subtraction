# Compiler and flags
FC = gfortran
FFLAGS = -O -w 
LDFLAGS = -lc -L/home/niraj/1TB-Disc/Workspace-IITG/Packages/Install/Cuba-4.2.2 -lcuba
export LHFLAGS = $(shell lhapdf-config --ldflags)

# Directories
OBJDIR = objects
SRCDIR = ../
LOCALSRCDIR = .


# OpenLoops paths
OPENLOOPS_DIR = /home/niraj/1TB-Disc/Workspace-IITG/Packages/Install/OpenLoops
OPENLOOPS_INCLUDE = $(OPENLOOPS_DIR)/lib_src/openloops/mod
OPENLOOPS_LIB = $(OPENLOOPS_DIR)/lib


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


#FC = gfortran
#FFLAGS   = -O 
#LDFLAGS = -lc #-fopenmp 
#MYFLAG = #-L./home/niraj/1TB-Disc/Workspace-IITG/Packages/Source/mylibrary -lfile
##export LHFLAGS = $(shell lhapdf-config --ldflags)
#LHFLAGS = $(shell lhapdf-config --ldflags)
#PROG    = runPK
#OBJECTS4  = cuts.lo2.o mat.amp.o phasespace.o misc.o printframe.o PK_new.o pdf_lha.o lum.o integrand.o vsup.o filehandle.o cone.o couplings.o 
#	
#$(PROG): $(OBJECTS4) main.f
#	$(FC) main.f $(FFLAGS) -o $(PROG) $(OBJECTS4) $(LDFLAGS) $(LHFLAGS) $(MYFLAG) 
#lum.o: ../lum.f
#	$(FC) $(FFLAGS) -c ../lum.f -o lum.o
#misc.o: ../misc.f
#	$(FC) $(FFLAGS) -c ../misc.f -o misc.o
#vsup.o: ../vsup.f
#	$(FC) $(FFLAGS) -c ../vsup.f -o vsup.o
#pdf_lha.o: ../pdf_lha.f
#	$(FC) $(FFLAGS) -c ../pdf_lha.f -o pdf_lha.o
#cuba.o: ../cuba.f
#	$(FC) $(FFLAGS) -c ../cuba.f -o cuba.o
#phasespace.o: ../phasespace.f
#	$(FC) $(FFLAGS) -c ../phasespace.f -o phasespace.o
#printframe.o: ../printframe.f
#	$(FC) $(FFLAGS) -c ../printframe.f -o printframe.o
#filehandle.o: ../filehandle.f
#	$(FC) $(FFLAGS) -c ../filehandle.f -o filehandle.o
#mat.amp.o: ../mat.amp.f
#	$(FC) $(FFLAGS) -c ../mat.amp.f -o mat.amp.o
#
#
#couplings.o: ../slicing_files/couplings.f
#	$(FC) $(FFLAGS) -c ../slicing_files/couplings.f -o couplings.o
#cone.o: ../slicing_files/cone.f
#	$(FC) $(FFLAGS) -c ../slicing_files/cone.f -o cone.o
#cuts.lo2.o: ../slicing_files/cuts.lo2.f
#	$(FC) $(FFLAGS) -c ../slicing_files/cuts.lo2.f -o cuts.lo2.o
#
#	
##	rm -f *.o
#clean:
#	rm -f $(OBJECTS4)
#	rm -f $(PROG).o $(PROG) *.o 
#
#distclean:
#	make clean
#	rm -f  $(PROG)
#
#$(PROG).o :  $(PROG).f param_card.dat
#cmatrix.o: *.inc


