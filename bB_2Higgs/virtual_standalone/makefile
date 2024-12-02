FC = gfortran
FFLAGS   = -O
LDFLAGS = -lc
LIBS = 
export LHFLAGS=$(shell lhapdf-config --ldflags)
LIB = -lm 
#OBJECTS = misc.o Iterm.o matrixLO1.o couplings.o cmatrix.o cuba.o lum.o pdf_lha.o vsup.o integrand.o phasespace.o
OBJECTS =filehandle.o couplings.o cuts.nlo3.o mat.amp.o printframe.o misc.o lum.o pdf_lha.o vsup.o integrand.o phasespace.o cone.o
PROG = runVir
$(PROG): $(OBJECTS) main.f  
	$(FC) main.f $(FFLAGS) -o $(PROG) $(OBJECTS) $(LIB) $(LIBS) $(LHFLAGS) $(LDFLAGS) 
lum.o: ../lum.f
	$(FC) $(FFLAGS) -c ../lum.f -o lum.o
vsup.o: ../vsup.f
	$(FC) $(FFLAGS) -c ../vsup.f -o vsup.o
pdf_lha.o: ../pdf_lha.f
	$(FC) $(FFLAGS) -c ../pdf_lha.f -o pdf_lha.o
misc.o: ../misc.f
	$(FC) $(FFLAGS) -c ../misc.f -o misc.o
phasespace.o: ../phasespace.f
	$(FC) $(FFLAGS) -c ../phasespace.f -o phasespace.o
printframe.o: ../printframe.f
	$(FC) $(FFLAGS) -c ../printframe.f -o printframe.o
cone.o: ../slicing_files/cone.f
	$(FC) $(FFLAGS) -c ../slicing_files/cone.f -o cone.o
mat.amp.o: ../mat.amp.f
	$(FC) $(FFLAGS) -c ../mat.amp.f -o mat.amp.o
couplings.o: ../slicing_files/couplings.f
	$(FC) $(FFLAGS) -c ../slicing_files/couplings.f -o couplings.o
filehandle.o: ../filehandle.f
	$(FC) $(FFLAGS) -c ../filehandle.f -o filehandle.o

cuts.nlo3.o: ../slicing_files/cuts.nlo3.f
	$(FC) $(FFLAGS) -c ../slicing_files/cuts.nlo3.f -o cuts.nlo3.o


clean: 
	rm -f *.o core *~ *.x $(PROG)
$(PROG): 
$(PROG).o: $(PROG).f param_card.dat
