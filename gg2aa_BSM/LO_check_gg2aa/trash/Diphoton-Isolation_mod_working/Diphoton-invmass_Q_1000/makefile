#F77=fort77
#F77=mpxlf
#F77=mpif77
#F77=mpif77.openmpi
F77=gfortran

#FFLAGS= -c -g -Wall
#FFLAGS= -c -O3
FFLAGS= -c 
#FFLAGS= -c -O3 -qstrict

#F77=pgf77
#FFLAGS= -c -O2 -Munroll -tp p6 -pc 64

LINK=$(F77)
export LFLAGS=$(shell lhapdf-config --ldflags)

all: cteq.x

parton_new.o: parton_new.f 
	$(F77) $(FFLAGS) parton_new.f 
dilog.o: dilog.f
	$(F77) $(FFLAGS) dilog.f

event.lo2.o: event.lo2.f
	$(F77) $(FFLAGS) event.lo2.f

sm.gg.amp.o: sm.gg.amp.f
	$(F77) $(FFLAGS) sm.gg.amp.f
bsm.gg.intf.amp.o: bsm.gg.intf.amp.f
	$(F77) $(FFLAGS) bsm.gg.intf.amp.f

event.nlo2.o: event.nlo2.f
	$(F77) $(FFLAGS) event.nlo2.f

event.nlo3.o: event.nlo3.f
	$(F77) $(FFLAGS) event.nlo3.f

main.o: main.f
	$(F77) $(FFLAGS) main.f
misc.o: misc.f
	$(F77) $(FFLAGS) misc.f
pdf.o: pdf.f
	$(F77) $(FFLAGS) pdf.f

sig.lo2.o: sig.lo2.f
	$(F77) $(FFLAGS) sig.lo2.f

sig.nlo2.qqb.o: sig.nlo2.qqb.f
	$(F77) $(FFLAGS) sig.nlo2.qqb.f
sig.nlo2.gg.o: sig.nlo2.gg.f
	$(F77) $(FFLAGS) sig.nlo2.gg.f
sig.nlo2.qg.o: sig.nlo2.qg.f
	$(F77) $(FFLAGS) sig.nlo2.qg.f

sig.nlo3.o: sig.nlo3.f
	$(F77) $(FFLAGS) sig.nlo3.f

matrix.lo.o:matrix.lo.f  
	$(F77) $(FFLAGS) matrix.lo.f
vrt.col.sft.o:vrt.col.sft.f  
	$(F77) $(FFLAGS) vrt.col.sft.f
diphoton7.o: diphoton7.f
	$(F77) $(FFLAGS) diphoton7.f
lum.o:lum.f  
	$(F77) $(FFLAGS) lum.f
vsup.o: vsup.f
	$(F77) $(FFLAGS) vsup.f

pdft.Pgg.o: pdft.Pgg.f
	$(F77) $(FFLAGS) pdft.Pgg.f
pdft.Pqq.o: pdft.Pqq.f
	$(F77) $(FFLAGS) pdft.Pqq.f
pdft.Pqg.o: pdft.Pqg.f
	$(F77) $(FFLAGS) pdft.Pqg.f
pdft.Pgq.o: pdft.Pgq.f
	$(F77) $(FFLAGS) pdft.Pgq.f

kin.var.lo2.o: kin.var.lo2.f
	$(F77) $(FFLAGS) kin.var.lo2.f
kin.var.nlo3.o: kin.var.nlo3.f
	$(F77) $(FFLAGS) kin.var.nlo3.f
cuts.lo2.o: cuts.lo2.f
	$(F77) $(FFLAGS) cuts.lo2.f
cuts.nlo2.o: cuts.nlo2.f
	$(F77) $(FFLAGS) cuts.nlo2.f
cuts.nlo3.o: cuts.nlo3.f
	$(F77) $(FFLAGS) cuts.nlo3.f
couplings.o: couplings.f
	$(F77) $(FFLAGS) couplings.f
cone.o: cone.f
	$(F77) $(FFLAGS) cone.f


cteq.x: main.o vsup.o dilog.o event.lo2.o sm.gg.amp.o bsm.gg.intf.amp.o event.nlo2.o event.nlo3.o sig.lo2.o sig.nlo2.qqb.o sig.nlo2.qg.o sig.nlo2.gg.o sig.nlo3.o pdf.o parton_new.o misc.o lum.o pdft.Pqq.o pdft.Pgg.o pdft.Pqg.o pdft.Pgq.o kin.var.lo2.o kin.var.nlo3.o cuts.lo2.o cuts.nlo2.o cuts.nlo3.o matrix.lo.o vrt.col.sft.o couplings.o cone.o diphoton7.o

	$(LINK) -o cteq.x main.o vsup.o dilog.o event.lo2.o sm.gg.amp.o bsm.gg.intf.amp.o event.nlo2.o event.nlo3.o sig.lo2.o sig.nlo2.qqb.o sig.nlo2.qg.o sig.nlo2.gg.o sig.nlo3.o pdf.o parton_new.o misc.o lum.o pdft.Pqq.o pdft.Pgg.o pdft.Pqg.o pdft.Pgq.o kin.var.lo2.o kin.var.nlo3.o cuts.lo2.o cuts.nlo2.o cuts.nlo3.o matrix.lo.o vrt.col.sft.o couplings.o cone.o diphoton7.o $(LFLAGS)

clean: 
	rm -f *.o core *~ *.x
