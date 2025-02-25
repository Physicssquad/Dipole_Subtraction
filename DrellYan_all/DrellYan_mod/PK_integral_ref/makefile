#FC = gfortran#g77
##FC = f77#gfortran#g77
#FFLAGS   = -O -w 
#LDFLAGS = -lc 
#MYFLAG = #-L./home/niraj/1TB-Disc/Workspace-IITG/Packages/Source/mylibrary -lfile
##export LHFLAGS = $(shell lhapdf-config --ldflags)
#LHFLAGS = $(shell lhapdf-config --ldflags)
#PROG    = runPK
#OBJECTS4  =filehandle.o misc.o phasespace.o printframe.o PK_new.o pdf_lha.o lum.o integrand.o phasespacePK.o vsup.o
#	
#$(PROG): $(OBJECTS4) main.f
#	$(FC) main.f $(FFLAGS) -o $(PROG) $(OBJECTS4) $(LDFLAGS) $(LHFLAGS) $(MYFLAG) 
#lum.o: ../lum.f
#	$(FC) $(FFLAGS) -c ../lum.f -o lum.o
#vsup.o: ../vsup.f
#	$(FC) $(FFLAGS) -c ../vsup.f -o vsup.o
#pdf_lha.o: ../pdf_lha.f
#	$(FC) $(FFLAGS) -c ../pdf_lha.f -o pdf_lha.o
#cuba.o: ../cuba.f
#	$(FC) $(FFLAGS) -c ../cuba.f -o cuba.o
#alphaS.o: ../alphaS.f
#	$(FC) $(FFLAGS) -c ../alphaS.f -o alphaS.o
#phasespace.o: ../phasespace.f
#	$(FC) $(FFLAGS) -c ../phasespace.f -o phasespace.o
#misc.o: ../misc.f
#	$(FC) $(FFLAGS) -c ../misc.f -o misc.o
#printframe.o: ../printframe.f
#	$(FC) $(FFLAGS) -c ../printframe.f -o printframe.o
#filehandle.o: ../filehandle.f
#	$(FC) $(FFLAGS) -c ../filehandle.f -o filehandle.o
#	
##	rm -f *.o
#clean:
#	rm -f $(OBJECTS4)
#	rm -f $(PROG).o $(PROG) 
#
#distclean:
#	make clean
#	rm -f  $(PROG)
#
#$(PROG).o :  $(PROG).f param_card.dat
#cmatrix.o: *.inc
#
#OBJECTS4  =filehandle.o misc.o phasespace.o printframe.o PK_new.o pdf_lha.o lum.o integrand.o phasespacePK.o vsup.o
#	
# Compiler and flags
FC = gfortran
FFLAGS = -O -w
LDFLAGS = -lc -L/home/niraj/1TB-Disc/Workspace-IITG/Packages/Install/Cuba-4.2.2 -lcuba 
export LHFLAGS = $(shell lhapdf-config --ldflags)

# OpenLoops paths
OPENLOOPS_DIR = /home/niraj/1TB-Disc/Workspace-IITG/Packages/Install/OpenLoops
OPENLOOPS_INCLUDE = $(OPENLOOPS_DIR)/lib_src/openloops/mod
OPENLOOPS_LIB = $(OPENLOOPS_DIR)/lib

# Directories
OBJDIR = objects
SRCDIR = ..
LOCALSRCDIR = .

# Program and object files
PROG = runPK

OBJECTS = $(OBJDIR)/mat.amp.o    \
          $(OBJDIR)/printframe.o \
          $(OBJDIR)/filehandle.o \
          $(OBJDIR)/misc.o       \
          $(OBJDIR)/lum.o        \
          $(OBJDIR)/pdf_lha.o    \
          $(OBJDIR)/vsup.o       \
          $(OBJDIR)/integrand.o  \
          $(OBJDIR)/phasespace.o \
          $(OBJDIR)/phasespacePK.o \
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
$(OBJDIR)/%.o: $(SRCDIR)/%.f | $(OBJDIR)
	$(FC) $(FFLAGS) -I$(OPENLOOPS_INCLUDE) -c $< -o $@

# For files in the current directory
$(OBJDIR)/%.o: %.f | $(OBJDIR)
	$(FC) $(FFLAGS) -I$(OPENLOOPS_INCLUDE) -c $< -o $@

# Clean target
clean:
	rm -f $(OBJDIR)/*.o core *~ $(PROG)
