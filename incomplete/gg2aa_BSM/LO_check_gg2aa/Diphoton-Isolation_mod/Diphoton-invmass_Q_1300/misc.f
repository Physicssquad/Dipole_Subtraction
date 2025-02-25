      function pqg0(z)
      implicit double precision (a-h,o-z)
      omz=1d0-z
      pqg0=0.5d0*(z*z+omz*omz)
      return
      end

      function pqg1(z)
      implicit double precision (a-h,o-z)
      pqg1=z*(z-1d0)
      return
      end

      function pgg0(z)
      implicit double precision (a-h,o-z)
      parameter (n=3)
      omz=1d0-z
      pgg0=2*n*(z/omz+omz/z+z*omz)
      return
      end

      function pgg1(z)
      implicit double precision (a-h,o-z)
      pgg1=0d0
      return
      end

      function pqq0(z)
      implicit double precision (a-h,o-z)
      parameter (n=3)
      parameter (v=n*n-1)
      pqq0=0.5d0*v*(1d0+z*z)/(n*(1d0-z))
      return
      end
      function pqq1(z)
      implicit double precision (a-h,o-z)
      parameter (n=3)
      parameter (v=n*n-1)
      pqq1=0.5d0*v*(z-1d0)/n
      return
      end

      function pgq0(z)
      implicit double precision (a-h,o-z)
      pgq0=pqq0(1d0-z)
      return
      end

      function pgq1(z)
      implicit double precision (a-h,o-z)
      pgq1=pqq1(1d0-z)
      return
      end

      function gfun(y,a)
      implicit double precision (a-h,o-z)
      tmp=dlog(a*(1-y)/y)
      gfun=-0.5*tmp*tmp
      return
      end

      function ffun(y,dc)
      implicit double precision (a-h,o-z)
      omy=1-y
      omdc=1-dc
      ffun=dilog((1-dc/omy)/omdc)+dlog(1-dc/omy)*dlog(dc*y/omy/omdc)
      return
      end


cc. cone size
c      function rjet(pa,pb)
c      implicit double precision (a-h,o-z)
c      dimension pa(0:3),pb(0:3)
c      parameter (pi=3.14159265358979d0)
c      ya=rapid(pa)
c      yb=rapid(pb)
c      phia=datan2(pa(2),pa(1))
c      phib=datan2(pb(2),pb(1))
c      delt=dabs(phia-phib)
c      if(delt.gt.pi)delt=delt-2*pi
c      rjet=dsqrt((ya-yb)**2+delt*delt)
c      return
c      end
c. rapidity
c      function rapid(p)
c      implicit double precision (a-h,o-z)
c      dimension p(0:3)
c      rapid=0.5d0*dlog( (p(0)+p(3)) / (p(0)-p(3)) )
c      return
c      end
c. pt
c      function eperp(p)
c      implicit double precision (a-h,o-z)
c      dimension p(0:3)
c      eperp=dsqrt(p(1)*p(1)+p(2)*p(2))
c      return
c      end
c. function to reconstruct azimuthal angle
c. n=0 returns   0<phi<2pi
c. n=1 returns -pi<phi<+pi
c      function phi(n,p)
c      implicit double precision (a-h,o-z)
c      dimension p(0:3)
c      parameter (pi=3.141592653589793d0)
c      pt=eperp(p)
c      if(p(2).ge.0.d0)then
c         phi=dacos(p(1)/pt)
c      else
c         phi=dacos(-p(1)/pt)+pi
c      endif
c      if(n.eq.1)phi=phi-2*pi*idnint(0.5d0*phi/pi)
c      return
c      end
c. dot product

      function dot(p,q)
      implicit double precision (a-h,o-z)
      dimension p(0:3),q(0:3)
      dot=p(0)*q(0)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)
      return
      end

      subroutine vwrite(p)
      implicit double precision (a-h,o-z)
      dimension p(0:3)
      write(*,*) ' '
      write(*,*) 'p(0) = ',p(0)
      write(*,*) 'p(1) = ',p(1)
      write(*,*) 'p(2) = ',p(2)
      write(*,*) 'p(3) = ',p(3)
      write(*,*) 'p.p  = ',dot(p,p)
      return
      end

      subroutine jwrite(et,eta,phi)
      implicit double precision (a-h,o-z)
      write(*,*) ' '
      p0=et*dcosh(eta)
      px=et*dcos(phi)
      py=et*dsin(phi)
      pz=et*dsinh(eta)
      write(*,*) 'p(0) = ',p0
      write(*,*) 'p(1) = ',px
      write(*,*) 'p(2) = ',py
      write(*,*) 'p(3) = ',pz
      write(*,*) 'p.p  = ',p0*p0-px*px-py*py-pz*pz
      return
      end

c. concatenate str1 and str2 into str
      subroutine strcat(str1,str2,str)
      character *(*) str1,str2,str
      l1=istrl(str1)
      l2=istrl(str2)
      l =len(str)
      if(l.lt.l1+l2) then
          write(*,*) 'error: l1+l2>l in strcat'
          write(*,*) 'l1=',l1,' str1=',str1
          write(*,*) 'l2=',l2,' str2=',str2
          write(*,*) 'l=',l
          stop
      endif
      if(l1.ne.0) str(1:l1)=str1(1:l1)
      if(l2.ne.0) str(l1+1:l1+l2)=str2(1:l2)
      if(l1+l2+1.le.l) str(l1+l2+1:l)= ' '
      end
c. returns the position of the last non-blank character in string
      function istrl(string)
      character * (*) string
      i = len(string)
      dowhile(i.gt.0.and.string(i:i).eq.' ')
         i=i-1
      enddo
      istrl = i
      end

c	For calculating the the bin size \Delta M (in GeV) of the first
c	RS resonance measured in TeV for ATLAS detector.
c----------------------------------------------------------------------
c       M should be converted into TeV, hence the name aMT
c       Delta M will be in GeV.
c----------------------------------------------------------------------
        subroutine bin_size(aM, Delta_M) 
        implicit double precision(a-h,o-z)
        aMT = aM/1000d0
        Delta_M = 24d0*dsqrt(0.625d0*aMT + aMT**2d0 + 0.0056d0)
c        write(*,321)aM, Delta_M
c 321    format(f10.4,f6.1)
        return
        end
c=========================================================================
c subroutine "class" does the classification of the processors
c with label "iirank" in order to assign  them each to one subprocess.
c This classification depends both on the order (LO/NLO) and on 
c the model SM/BSM.
c-------------------------------------------------------------------------
c At NLO in SM, the processors with iirank=1 calculate both qg(3-body)
c and gg(2-body).  At NLO in BSM, the processors with iirank=4 calculate
c both qg INTF(3-body) and gg INTF (2-body).  This way we can minimize 
c the number of processors and hence the time of computation.
c=========================================================================

        subroutine class(iirank, io_min, io_max, is_min, is_max)
        implicit double precision(a-h, o-z) 
        common/model/model
        common/max_order/iorder

        if (iorder.eq.1) then

        if (model.eq.0) then

        if (iirank.eq.0) then
        io_min = 1
        io_max = 1              ! select only SM
        is_min = 1              ! qqb
        is_max = 1
        elseif (iirank.eq.1) then
        io_min = 1
        io_max = 1              ! select only SM
        is_min = 2              ! qg
        is_max = 3              ! gg
        endif

        elseif (model.eq.1 .or. model.eq.2 .or. model.eq.3) then

        if (iirank.eq.0) then
        io_min = 2                              ! select BSM direct
        io_max = 2
        is_min = 1                              ! select qqb
        is_max = 1                              
        elseif (iirank.eq.1) then
        io_min = 2                              ! select BSM direct
        io_max = 2
        is_min = 2                              ! select qg
        is_max = 2                              
        elseif (iirank.eq.2) then
        io_min = 2                              ! select BSM direct
        io_max = 2
        is_min = 3                              ! select gg
        is_max = 3                              
        elseif (iirank.eq.3) then
        io_min = 3                              ! select SM*BSM interference
        io_max = 3
        is_min = 1                              ! select qqb intf
        is_max = 1                              
        elseif (iirank.eq.4) then
        io_min = 3                              ! select SM*BSM interference
        io_max = 3
        is_min = 2                              ! select qg intf
        is_max = 3                              ! select gg intf
        endif

        endif           ! model endif

c       Assigning only one subprocess to one processor (machine)

        elseif (iorder.eq.0) then

        if (iirank.eq.0) then

        if (model.eq.0) then
        io_min = 1
        io_max = 1
        elseif (model.eq.1 .or. model.eq.2 .or. model.eq.3) then
        io_min = 2
        io_max = 3
        endif

        is_min = 1
        is_max = 3

        endif               ! endif iirank


        endif               ! endif iorder

        return
        end




