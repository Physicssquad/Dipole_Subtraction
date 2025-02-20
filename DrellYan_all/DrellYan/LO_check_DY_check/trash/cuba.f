      subroutine cubacheck
      implicit none
c      double precision epsrel,epsabs

c     Common Parameters      
      integer ndim, ncomp, nvec, last,flags, seed, mineval, maxeval
      integer n4,i2
      real*8 epsrel, epsabs, userdata
      parameter (ndim = 6)
      parameter (ncomp = 1)
      parameter (userdata = 0)
      parameter (nvec = 1)
      parameter (last = 4)
      parameter (flags = 5)
      parameter (seed = 0)
      parameter (mineval = 0)
      parameter (maxeval = 500000000)
      common/countc/n4

c     Suave Specific Parameters

      integer nstart, nincrease, nbatch, gridno
      integer*8 spin
      character*(*) statefile 
      parameter (nstart = 1000)
      parameter (nincrease = 700)
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

      integer c,I,cut(1:2)
      common/jetcut/cut
      epsrel = 1D-6
      epsabs = 1D-25
      write(*,*)"Check Integral using:"
      print *, "----------------------------------------"
      write(*,*)"     VEGAS     ~~~~~> 1"
      write(*,*)"     SUAVE     ~~~~~> 2"
      write(*,*)"    DIVONNE    ~~~~~> 3"
      write(*,*)"     CUHRE     ~~~~~> 4"
      write(*,*)" Run all at Once ~~~> 5"
      print *, "----------------------------------------"
      read(*,*)i
c        i=1
c      write(*,*)"READ VERBOSE:"
c      read(*,*)flags
c      print*,"Enter the precision."
c      print*,"Current Value of epsrel:",epsrel
c      print*,"Current Value of epsabs:",epsabs
c      print*," "
c      print*,"Press 1 to go with it else press 2 to modify"
c      read*,i2
c      if ( i2 .eq. 2 ) then
c              print*,"Enter epsrel:"
c              read*,epsrel
c              print*,"Enter epsabs:"
c              read*,epsabs
c              print*," "
c              print*,"Values are "
c              print*,"epsrel:",epsrel
c              print*,"epsabs:",epsabs
c              call sleep(1)
c      endif
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
c        write(*,*)'n4 = ',n4
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
      ELSEIF (i .eq. 4 .or. i .eq. 5) then
      ENDIF
      end
