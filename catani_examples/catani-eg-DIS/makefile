FC = gfortran
FFLAGS   = -O
LIBS = 
##LIB = -lm 
LIBS = 
LIB = 
OBJECTS = phasespace.o vsup.o integrand.o misc.o check-dipole.o
PROG = run.x
$(PROG): $(OBJECTS) main.f  
	$(FC) main.f $(FFLAGS) -o $(PROG) $(OBJECTS) $(LIB) $(LIBS) 
clean: 
	rm -f *.o core *~ *.x
