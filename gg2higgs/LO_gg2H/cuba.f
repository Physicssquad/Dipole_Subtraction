
c      subroutine cubacheck(integral(1),error(1))
      subroutine cubacheck(answer,sd)
      implicit none

c     Common Parameters      
      integer ndim, ncomp, nvec, last,flags, seed, mineval, maxeval
      real*8 epsrel, epsabs, userdata, answer, sd
      parameter (ndim = 2)
      parameter (ncomp = 1)
      parameter (userdata = 0)
      parameter (nvec = 1)
      parameter (epsrel = 1D-1)
      parameter (epsabs = 1D-18)
      parameter (last = 4)
      parameter (flags = 1)
      parameter (seed = 0)
      parameter (mineval = 1000000)
      parameter (maxeval = 500000000)

c     Suave Specific Parameters

      integer nstart, nincrease, nbatch, gridno
      integer*8 spin
      character*(*) statefile 
      parameter (nstart = 1000)
      parameter (nincrease = 1000)
      parameter (nbatch = 1000)
      parameter (gridno = 0)
      parameter (statefile = " ")
      parameter (spin = -1)

c     Divonne Specific Parameters

      integer nnew, nmin
      real*8 flatness
      parameter (nnew = 1000)
      parameter (nmin = 2)
      parameter (flatness = 25D0)

c     Cuhre Specific Parameters

      integer key1, key2, key3, maxpass
      real*8 border, maxchisq, mindeviation
      integer ngiven, ldxgiven, nextra
      parameter (key1 = 47)
      parameter (key2 = 1)
      parameter (key3 = 1)
      parameter (maxpass = 5)
      parameter (border = 0D0)
      parameter (maxchisq = 10D0)
      parameter (mindeviation = .25D0)
      parameter (ngiven = 0)
      parameter (ldxgiven = ndim)
      parameter (nextra = 0)
      integer key
      parameter (key = 0)
        
      external integrand
      real*8 integral(ncomp), error(ncomp), prob(ncomp)
      integer verbose, nregions, neval, fail
      character*16 env
      parameter (verbose = 0)

      integer c,I
       call cubacores(5,1000)
c      write(*,*)"Check Integral using:"
c      print *, "----------------------------------------"
c      write(*,*)"     VEGAS     ~~~~~> 1"
c      write(*,*)"     SUAVE     ~~~~~> 2"
c      write(*,*)"    DIVONNE    ~~~~~> 3"
c      write(*,*)"     CUHRE     ~~~~~> 4"
c      write(*,*)" Run all at Once ~~~> 5"
c      print *, "----------------------------------------"
c      read(*,*)i
        i=1
c      write(*,*)"READ VERBOSE:"
c      read(*,*)flags
      print *, " "
      write(*,*)" Initializing Integration using selected method "
      print *, " "

      IF (i .eq. 1 .or. i .eq. 5) then
      print *, "-------------------- Vegas test --------------------"
      CALL vegas(ndim, ncomp, integrand, userdata, nvec,
     .      epsrel, epsabs,flags, seed, mineval, maxeval,
     .      nstart, nincrease, nbatch, gridno, statefile, spin,
     .      neval, fail, integral, error, prob)
        write(*,*)"neval      = ", neval
        write(*,*)"fail       = ", fail
        write(*,*)"Integral   = ",integral(1),"+-",error(1)
c         print '(F20.12," +- ",F20.12,"   p = ",F8.3)',
c     &    (integral(c), error(c), prob(c), c = 1, ncomp)

c        write(*,*)"with chisq = ",prob(1)
        print *, " "

      ENDIF
      IF (i .eq. 2 .or. i .eq. 5) then
        print *, "-------------------- Suave test --------------------"
       call suave(ndim, ncomp, integrand, userdata, nvec,
     .      epsrel, epsabs,flags,seed, mineval, maxeval,
     .      nnew,nmin,flatness,statefile,spin,
     .      nregions, neval, fail, integral, error, prob)
        write(*,*)"nregions   = ", nregions
        write(*,*)"neval      = ", neval
        write(*,*)"fail       = ", fail
        write(*,*)"Integral   = ",integral(1),"+-",error(1)
c        write(*,*)"prob       = ",prob(1)
        print *, " "

      ENDIF
      IF (i .eq. 3 .or. i .eq. 5) then
       print *, "------------------- Divonne test -------------------"
       call divonne(ndim, ncomp, integrand, userdata, nvec,
     .    epsrel, epsabs, verbose, seed,
     .    mineval, maxeval, key1, key2, key3, maxpass,
     .    border, maxchisq, mindeviation,
     .    ngiven, ldxgiven, 0, nextra, 0,
     .    statefile, spin,
     .    nregions, neval, fail, integral, error, prob)

        write(*,*)"nregions   = ", nregions
        write(*,*)"neval      = ", neval
        write(*,*)"fail       = ", fail
        write(*,*)"Integral   = ",integral(1),"+-",error(1)
c        write(*,*)"prob       = ",prob(1)
        print *, " "
      ENDIF
      IF (i .eq. 4 .or. i .eq. 5) then
       print *, "-------------------- Cuhre test --------------------"
       call cuhre(ndim, ncomp, integrand, userdata, nvec,
     .    epsrel, epsabs, verbose + last,
     .    mineval, maxeval, key,
     .    statefile, spin,
     .    nregions, neval, fail, integral, error, prob)
        write(*,*)"nregions   = ", nregions
        write(*,*)"neval      = ", neval
        write(*,*)"fail       = ", fail
        write(*,*)"Integral   = ",integral(1),"+-",error(1)
c        write(*,*)"prob       = ",prob(1)
        print *, " "
      ENDIF
	answer = integral(1)
	    sd = error(1)
      end
