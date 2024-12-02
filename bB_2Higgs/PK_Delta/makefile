FC = gfortran
FFLAGS   = -O 
LDFLAGS = -lc 
LHFLAGS = $(shell lhapdf-config --ldflags)
PROG    = runPK
OBJECTS4  = couplings.o cone.o cuts.lo2.o filehandle.o mat.amp.o phasespace.o misc.o printframe.o PK_new.o pdf_lha.o lum.o integrand.o vsup.o
	
$(PROG): $(OBJECTS4) main.f
	$(FC) main.f $(FFLAGS) -o $(PROG) $(OBJECTS4) $(LDFLAGS) $(LHFLAGS) $(MYFLAG) 
lum.o: ../lum.f
	$(FC) $(FFLAGS) -c ../lum.f -o lum.o
vsup.o: ../vsup.f
	$(FC) $(FFLAGS) -c ../vsup.f -o vsup.o
pdf_lha.o: ../pdf_lha.f
	$(FC) $(FFLAGS) -c ../pdf_lha.f -o pdf_lha.o
phasespace.o: ../phasespace.f
	$(FC) $(FFLAGS) -c ../phasespace.f -o phasespace.o
misc.o: ../misc.f
	$(FC) $(FFLAGS) -c ../misc.f -o misc.o
printframe.o: ../printframe.f
	$(FC) $(FFLAGS) -c ../printframe.f -o printframe.o
mat.amp.o: ../mat.amp.f
	$(FC) $(FFLAGS) -c ../mat.amp.f -o mat.amp.o
filehandle.o: ../filehandle.f
	$(FC) $(FFLAGS) -c ../filehandle.f -o filehandle.o



couplings.o: ../slicing_files/couplings.f
	$(FC) $(FFLAGS) -c ../slicing_files/couplings.f -o couplings.o
cone.o: ../slicing_files/cone.f
	$(FC) $(FFLAGS) -c ../slicing_files/cone.f -o cone.o
cuts.lo2.o: ../slicing_files/cuts.lo2.f
	$(FC) $(FFLAGS) -c ../slicing_files/cuts.lo2.f -o cuts.lo2.o


	
#	rm -f *.o
clean:
	rm -f $(OBJECTS4)
	rm -f $(PROG).o $(PROG) 

distclean:
	make clean
	rm -f  $(PROG)

$(PROG).o :  $(PROG).f param_card.dat
cmatrix.o: *.inc

