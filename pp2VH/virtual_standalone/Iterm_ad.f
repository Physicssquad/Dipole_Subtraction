c Generated by AutoDipole
c Kouhei Hasegawa, Sven Moch, and Peter Uwer, 2009
c Filename:Iterm.f
c Subroutine: Iterm evaluates all I-terms
c
c input
c   real p      :a phase space point
c output
c   real coef   :coefficient of all I-terms
c   real SumI   :sum of all coefficients
c 
c Subroutine to calculate I terms 
c LO and Virtual process:{{g, p[1]}, {g, p[2]}} --> {{u, p[3]}, {ubar, p[4]}}
c I(i) = coef(i,-2)*eps^(-2) + coef(i,-1)*eps^(-1) + coef(i,0)  
c SumI[-2,-1,0] = Sum_{i=1}^{12}coef[-2,-1,0] 

      subroutine Iterm(p_ex,id_LO,coef,SumI,AL) 
      use openloops
      implicit none 

      integer i,j 
      double precision p(0:3,1:4),p_ex(0:3,1:3),coef(12,-2:0),SumI(-2:0)
      double precision Pi,rtwo,Eul,AL,CF,CA,TR,mt,mb,dot,mu,Nf,log2
      double precision CLV(12),q(0:3,1:4),dilog,dot2
      double precision s12,s13,s14,s23,s24,s34,Born_gg2h_ 
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),answer
      integer id_LO
 
      common /MASS/ mt,mb  
c      common /usedalpha/ AL  
    
      Pi=3.141592653589793238D0 
      Eul=0.5772156649015328606065120d0 
      rtwo=dsqrt(2.d0) 
      log2=Log(2.d0) 
      CF=4.D0/3.D0 
      CA=3.D0 
      TR=0.5D0 
      mu=125d0/2d0
      Nf=1D0 

!      p4(0) = p1(0) 
!      p4(1) = p1(1) 
!      p4(2) = p1(2) 
!      p4(3) = p1(3) 

!      s12=2.d0*dot(p1,p2) 
!      s13=2.d0*dot(p1,p3) 
!      s14=2.d0*dot(p1,p4) 
!      s23=2.d0*dot(p2,p3) 
!      s24=2.d0*dot(p2,p4) 
!      s34=2.d0*dot(p3,p4) 
  
       s12 = 2d0*dot2(p_ex,1,2,1)
       s13 = 2d0*dot2(p_ex,1,3,1)
       s23 = 2d0*dot2(p_ex,2,3,1)
c       s12 = 2d0*dot2(p_ex,1,3,1)

       s14 = 2d0*dot2(p_ex,1,2,1)
       s24 = 2d0*dot2(p_ex,1,2,1)
       s34 = 2d0*dot2(p_ex,1,2,1)
       call p2d_to_p1d_3(p_ex,p1,p2,p3)
      call p1d_to_p2d_3(p1,p2,p3,p_ex) 
c      Born = Born_gg2h_(0,AL,p1,p2,p3)
      answer = Born_gg2h_(0,AL,p1,p2,p3)
c      call evaluate_tree(id_LO,p_ex,answer)
      do i=1,12
       CLV(i) = -answer 
      enddo

      coef(1,-2)=
     - -0.5*(AL*CLV(1))/Pi
c  	print*,coef(1,-2),CLV(1),answer

c	stop
      coef(1,-1)=
     -         (AL*CLV(1)*(-3 + 2*Eul - 2*Log(4*Pi) - 
     -      2*Log(mu**2/s34)))/(4.*Pi)
   
      coef(1,0)=
     -         (AL*CLV(1)*(-60 - 6*Eul**2 + 7*Pi**2 - 18*Log(4*Pi) - 
     -      6*Log(4*Pi)**2 + 6*Eul*(3 + 2*Log(4*Pi)) + 
     -      6*(-3 + 2*Eul - 2*Log(4*Pi))*Log(mu**2/s34) - 
     -      6*Log(mu**2/s34)**2))/(24.*Pi)
   
      coef(2,-2)=
     - -0.5*(AL*CLV(2))/Pi
   
      coef(2,-1)=
     -         (AL*CLV(2)*(-3 + 2*Eul - 2*Log(4*Pi) - 
     -      2*Log(mu**2/s34)))/(4.*Pi)
   
      coef(2,0)=
     -         (AL*CLV(2)*(-60 - 6*Eul**2 + 7*Pi**2 - 18*Log(4*Pi) - 
     -      6*Log(4*Pi)**2 + 6*Eul*(3 + 2*Log(4*Pi)) + 
     -      6*(-3 + 2*Eul - 2*Log(4*Pi))*Log(mu**2/s34) - 
     -      6*Log(mu**2/s34)**2))/(24.*Pi)
   
      coef(3,-2)=
     - -0.5*(AL*CLV(3))/Pi
   
      coef(3,-1)=
     -         (AL*CLV(3)*(-3 + 2*Eul - 2*Log(4*Pi) - 
     -      2*Log(mu**2/s13)))/(4.*Pi)
   
      coef(3,0)=
     -         (AL*CLV(3)*(-60 - 6*Eul**2 + 7*Pi**2 - 18*Log(4*Pi) - 
     -      6*Log(4*Pi)**2 + 6*Eul*(3 + 2*Log(4*Pi)) + 
     -      6*(-3 + 2*Eul - 2*Log(4*Pi))*Log(mu**2/s13) - 
     -      6*Log(mu**2/s13)**2))/(24.*Pi)
   
      coef(4,-2)=
     - -0.5*(AL*CLV(4))/Pi
   
      coef(4,-1)=
     -         (AL*CLV(4)*(-3 + 2*Eul - 2*Log(4*Pi) - 
     -      2*Log(mu**2/s23)))/(4.*Pi)
   
      coef(4,0)=
     -         (AL*CLV(4)*(-60 - 6*Eul**2 + 7*Pi**2 - 18*Log(4*Pi) - 
     -      6*Log(4*Pi)**2 + 6*Eul*(3 + 2*Log(4*Pi)) + 
     -      6*(-3 + 2*Eul - 2*Log(4*Pi))*Log(mu**2/s23) - 
     -      6*Log(mu**2/s23)**2))/(24.*Pi)
   
      coef(5,-2)=
     - -0.5*(AL*CLV(5))/Pi
   
      coef(5,-1)=
     -         (AL*CLV(5)*(-3 + 2*Eul - 2*Log(4*Pi) - 
     -      2*Log(mu**2/s14)))/(4.*Pi)
   
      coef(5,0)=
     -         (AL*CLV(5)*(-60 - 6*Eul**2 + 7*Pi**2 - 18*Log(4*Pi) - 
     -      6*Log(4*Pi)**2 + 6*Eul*(3 + 2*Log(4*Pi)) + 
     -      6*(-3 + 2*Eul - 2*Log(4*Pi))*Log(mu**2/s14) - 
     -      6*Log(mu**2/s14)**2))/(24.*Pi)
   
      coef(6,-2)=
     - -0.5*(AL*CLV(6))/Pi
   
      coef(6,-1)=
     -         (AL*CLV(6)*(-3 + 2*Eul - 2*Log(4*Pi) - 
     -      2*Log(mu**2/s24)))/(4.*Pi)
   
      coef(6,0)=
     -         (AL*CLV(6)*(-60 - 6*Eul**2 + 7*Pi**2 - 18*Log(4*Pi) - 
     -      6*Log(4*Pi)**2 + 6*Eul*(3 + 2*Log(4*Pi)) + 
     -      6*(-3 + 2*Eul - 2*Log(4*Pi))*Log(mu**2/s24) - 
     -      6*Log(mu**2/s24)**2))/(24.*Pi)
   
      coef(7,-2)=
     - -0.5*(AL*CLV(7))/Pi
   
      coef(7,-1)=
     -         (AL*CLV(7)*(4*Nf*TR + CA*(-11 + 6*Eul - 6*Log(4*Pi)) - 
     -      6*CA*Log(mu**2/s13)))/(12.*CA*Pi)
   
      coef(7,0)=
     -         (AL*CLV(7)*(8*Nf*TR*(8 - 3*Eul + 3*Log(4*Pi)) + 
     -      CA*(-200 + 66*Eul - 18*Eul**2 + 21*Pi**2 - 
     -         66*Log(4*Pi) + 36*Eul*Log(4*Pi) - 18*Log(4*Pi)**2)
     -        + 6*(4*Nf*TR + CA*(-11 + 6*Eul - 6*Log(4*Pi)))*
     -       Log(mu**2/s13) - 18*CA*Log(mu**2/s13)**2))/
     -  (72.*CA*Pi)
   
      coef(8,-2)=
     - -0.5*(AL*CLV(8))/Pi
   
      coef(8,-1)=
     -         (AL*CLV(8)*(4*Nf*TR + CA*(-11 + 6*Eul - 6*Log(4*Pi)) - 
     -      6*CA*Log(mu**2/s14)))/(12.*CA*Pi)
   
      coef(8,0)=
     -         (AL*CLV(8)*(8*Nf*TR*(8 - 3*Eul + 3*Log(4*Pi)) + 
     -      CA*(-200 + 66*Eul - 18*Eul**2 + 21*Pi**2 - 
     -         66*Log(4*Pi) + 36*Eul*Log(4*Pi) - 18*Log(4*Pi)**2)
     -        + 6*(4*Nf*TR + CA*(-11 + 6*Eul - 6*Log(4*Pi)))*
     -       Log(mu**2/s14) - 18*CA*Log(mu**2/s14)**2))/
     -  (72.*CA*Pi)
   
      coef(9,-2)=
     - -0.5*(AL*CLV(9))/Pi
   
      coef(9,-1)=
     -         (AL*CLV(9)*(4*Nf*TR + CA*(-11 + 6*Eul - 6*Log(4*Pi)) - 
     -      6*CA*Log(mu**2/s23)))/(12.*CA*Pi)
   
      coef(9,0)=
     -         (AL*CLV(9)*(8*Nf*TR*(8 - 3*Eul + 3*Log(4*Pi)) + 
     -      CA*(-200 + 66*Eul - 18*Eul**2 + 21*Pi**2 - 
     -         66*Log(4*Pi) + 36*Eul*Log(4*Pi) - 18*Log(4*Pi)**2)
     -        + 6*(4*Nf*TR + CA*(-11 + 6*Eul - 6*Log(4*Pi)))*
     -       Log(mu**2/s23) - 18*CA*Log(mu**2/s23)**2))/
     -  (72.*CA*Pi)
   
      coef(10,-2)=
     - -0.5*(AL*CLV(10))/Pi
   
      coef(10,-1)=
     -         (AL*CLV(10)*(4*Nf*TR + CA*(-11 + 6*Eul - 6*Log(4*Pi)) - 
     -      6*CA*Log(mu**2/s24)))/(12.*CA*Pi)
   
      coef(10,0)=
     -         (AL*CLV(10)*(8*Nf*TR*(8 - 3*Eul + 3*Log(4*Pi)) + 
     -      CA*(-200 + 66*Eul - 18*Eul**2 + 21*Pi**2 - 
     -         66*Log(4*Pi) + 36*Eul*Log(4*Pi) - 18*Log(4*Pi)**2)
     -        + 6*(4*Nf*TR + CA*(-11 + 6*Eul - 6*Log(4*Pi)))*
     -       Log(mu**2/s24) - 18*CA*Log(mu**2/s24)**2))/
     -  (72.*CA*Pi)
   
      coef(11,-2)=
     - -0.5*(AL*CLV(11))/Pi
   
      coef(11,-1)=
     -         (AL*CLV(11)*(4*Nf*TR + CA*(-11 + 6*Eul - 6*Log(4*Pi)) - 
     -      6*CA*Log(mu**2/s12)))/(12.*CA*Pi)
   
      coef(11,0)=
     -         (AL*CLV(11)*(8*Nf*TR*(8 - 3*Eul + 3*Log(4*Pi)) + 
     -      CA*(-200 + 66*Eul - 18*Eul**2 + 21*Pi**2 - 
     -         66*Log(4*Pi) + 36*Eul*Log(4*Pi) - 18*Log(4*Pi)**2)
     -        + 6*(4*Nf*TR + CA*(-11 + 6*Eul - 6*Log(4*Pi)))*
     -       Log(mu**2/s12) - 18*CA*Log(mu**2/s12)**2))/
     -  (72.*CA*Pi)
   
      coef(12,-2)=
     - -0.5*(AL*CLV(12))/Pi
   
      coef(12,-1)=
     -         (AL*CLV(12)*(4*Nf*TR + CA*(-11 + 6*Eul - 6*Log(4*Pi)) - 
     -      6*CA*Log(mu**2/s12)))/(12.*CA*Pi)
   
      coef(12,0)=
     -         (AL*CLV(12)*(8*Nf*TR*(8 - 3*Eul + 3*Log(4*Pi)) + 
     -      CA*(-200 + 66*Eul - 18*Eul**2 + 21*Pi**2 - 
     -         66*Log(4*Pi) + 36*Eul*Log(4*Pi) - 18*Log(4*Pi)**2)
     -        + 6*(4*Nf*TR + CA*(-11 + 6*Eul - 6*Log(4*Pi)))*
     -       Log(mu**2/s12) - 18*CA*Log(mu**2/s12)**2))/
     -  (72.*CA*Pi)

c        do i=1,12
c        do j=0,2
c	print*,"I:",i,"eps:",j,coef(i,-j)
c        enddo
c	print*
c        enddo
   
      SumI(-2) = 0.d0 
      SumI(-1) = 0.d0 
      SumI(-0) = 0.d0 
 
      do i=1,12 
      SumI(-2) = SumI(-2) + coef(i,-2) 
      SumI(-1) = SumI(-1) + coef(i,-1) 
      SumI(0) = SumI(0) + coef(i,0) 
      enddo 
 
 
      return 
      end 

      subroutine mexchange(q,p,i,j)
      implicit none
      integer i,j,k,l
      double precision p(0:3,1:4),q(0:3,1:4)
      do k=0,3
       do l=3,4
        q(k,l) = p(k,l)
       enddo
       q(k,i) = p(k,j)
       q(k,j) = p(k,i)
      enddo
      return
      end
 
      subroutine mmap(q,p)
      implicit none
      integer i
      double precision  p(0:3),q(0:3)
      do i=0,3
       q(i) = p(i)
      enddo
      return
      end
 
 
