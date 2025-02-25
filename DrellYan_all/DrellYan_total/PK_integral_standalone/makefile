FC = gfortran#g77
#FC = f77#gfortran#g77
FFLAGS   = -O 
LDFLAGS = -lc 
MYFLAG = #-L./home/niraj/1TB-Disc/Workspace-IITG/Packages/Source/mylibrary -lfile
#export LHFLAGS = $(shell lhapdf-config --ldflags)
LHFLAGS = $(shell lhapdf-config --ldflags)
PROG    = runPK
OBJECTS4  =filehandle.o misc.o phasespace.o printframe.o PK_new.o pdf_lha.o lum.o integrand.o phasespacePK.o vsup.o
	
$(PROG): $(OBJECTS4) main.f
	$(FC) main.f $(FFLAGS) -o $(PROG) $(OBJECTS4) $(LDFLAGS) $(LHFLAGS) $(MYFLAG) 
lum.o: ../lum.f
	$(FC) $(FFLAGS) -c ../lum.f -o lum.o
vsup.o: ../vsup.f
	$(FC) $(FFLAGS) -c ../vsup.f -o vsup.o
pdf_lha.o: ../pdf_lha.f
	$(FC) $(FFLAGS) -c ../pdf_lha.f -o pdf_lha.o
cuba.o: ../cuba.f
	$(FC) $(FFLAGS) -c ../cuba.f -o cuba.o
alphaS.o: ../alphaS.f
	$(FC) $(FFLAGS) -c ../alphaS.f -o alphaS.o
phasespace.o: ../phasespace.f
	$(FC) $(FFLAGS) -c ../phasespace.f -o phasespace.o
misc.o: ../misc.f
	$(FC) $(FFLAGS) -c ../misc.f -o misc.o
printframe.o: ../printframe.f
	$(FC) $(FFLAGS) -c ../printframe.f -o printframe.o
filehandle.o: ../filehandle.f
	$(FC) $(FFLAGS) -c ../filehandle.f -o filehandle.o
	
#	rm -f *.o
clean:
	rm -f $(OBJECTS4)
	rm -f $(PROG).o $(PROG) 

distclean:
	make clean
	rm -f  $(PROG)

$(PROG).o :  $(PROG).f param_card.dat
cmatrix.o: *.inc

