FC = gfortran
FFLAGS   = -O
LDFLAGS = -lc -L./../ -lfile
LIBS = 
export LHFLAGS=$(shell lhapdf-config --ldflags)
LIB = -lm 
#OBJECTS = misc.o Iterm.o matrixLO1.o couplings.o cmatrix.o cuba.o lum.o pdf_lha.o vsup.o integrand.o phasespace.o
OBJECTS =printframe.o misc.o lum.o pdf_lha.o vsup.o integrand.o phasespace.o
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
clean: 
	rm -f *.o core *~ *.x $(PROG)
$(PROG): 
$(PROG).o: $(PROG).f param_card.dat
