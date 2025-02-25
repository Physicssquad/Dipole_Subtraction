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
PROG = runDipole

OBJECTS = $(OBJDIR)/filehandle.o \
          $(OBJDIR)/cuba.o       \
          $(OBJDIR)/mat.amp.o    \
          $(OBJDIR)/cuts.nlo3.o  \
          $(OBJDIR)/cone.o       \
          $(OBJDIR)/printframe.o \
          $(OBJDIR)/misc.o       \
          $(OBJDIR)/lum.o        \
          $(OBJDIR)/pdf_lha.o    \
          $(OBJDIR)/vsup.o       \
          $(OBJDIR)/integrand.o  \
          $(OBJDIR)/phasespace.o \
          $(OBJDIR)/dipole.o	 \
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

#**********[ OLD MAKEFILE SCRIPT ]***********#
#FC = gfortran
#FFLAGS   = -O
##LDFLAGS = #-L/home/niraj/1TB-Disc/Workspace-IITG/Packages/Install/autodipole/AutoDipole_V1.2.3/process/Proc_gg_uuxg -ldipole -ldhelas3 -lg2c
#LDFLAGS = -L/home/niraj/1TB-Disc/Workspace-IITG/Packages/Install/Cuba-4.2.2 -lcuba
#export LHFLAGS=$(shell lhapdf-config --ldflags)
#LIB = -lm 
#
##OBJECTS = Iterm.o matrixLO1.o couplings.o cmatrix.o PK.o dipole.o lum.o pdf_lha.o vsup.o mat.amp.o integrand.o misc.o phasespace.o
#OBJECTS = cuba.o printframe.o dipole.o lum.o pdf_lha.o vsup.o mat.amp.o integrand.o misc.o phasespace.o
#
## objects from slicing files
#OBJECTS2 = filehandle.o cone.o couplings.o cuts.nlo3.o
#
#PROG = runDipole
#$(PROG): $(OBJECTS) $(OBJECTS2) main.f  
#	$(FC) main.f $(FFLAGS) -o $(PROG) $(OBJECTS) $(OBJECTS2) $(LIB) $(LIBS) $(LHFLAGS) $(LDFLAGS) 
##	rm -f *o
#
#lum.o: ../lum.f
#	$(FC) $(FFLAGS) -c ../lum.f -o lum.o
#
#vsup.o: ../vsup.f
#	$(FC) $(FFLAGS) -c ../vsup.f -o vsup.o
#
#pdf_lha.o: ../pdf_lha.f
#	$(FC) $(FFLAGS) -c ../pdf_lha.f -o pdf_lha.o
#
#mat.amp.o: ../mat.amp.f
#	$(FC) $(FFLAGS) -c ../mat.amp.f -o mat.amp.o
#
#misc.o: ../misc.f
#	$(FC) $(FFLAGS) -c ../misc.f -o misc.o
#
#filehandle.o: ../filehandle.f
#	$(FC) $(FFLAGS) -c ../filehandle.f -o filehandle.o
#
#phasespace.o: ../phasespace.f
#	$(FC) $(FFLAGS) -c ../phasespace.f -o phasespace.o
#
#printframe.o: ../printframe.f
#	$(FC) $(FFLAGS) -c ../printframe.f -o printframe.o
#
#couplings.o: ../slicing_files/couplings.f
#	$(FC) $(FFLAGS) -c ../slicing_files/couplings.f -o couplings.o
#
#cuts.nlo3.o: ../slicing_files/cuts.nlo3.f
#	$(FC) $(FFLAGS) -c ../slicing_files/cuts.nlo3.f -o cuts.nlo3.o
#
#cone.o: ../slicing_files/cone.f
#	$(FC) $(FFLAGS) -c ../slicing_files/cone.f -o cone.o
#
#clean: 
#	rm -f *.o core *~ *.x $(PROG)
#$(PROG): 
#$(PROG).o: $(PROG).f param_card.dat

