
c-------------------------------------------c
c-     2 -> 2 matrix element squared       -c
c-------------------------------------------c      

c     SM AND BSM
      subroutine setsigma(s,t,u,sigma)
      implicit double precision (a-h,o-z)
      dimension sigma(10)
      parameter (n=3)
      parameter (pi=3.14159265358979d0)
      common/model/model

      rp34=dsqrt(s)
      ALEM=1.0D0/128.0D0
      e=DSQRT(4.0D0*PI*ALEM)
      if (model.eq.1 .or. model.eq.2 .or. model.eq.3) then
      call coupfact(model,rp34,AK2D,AK2DINTF)
      endif

!qqb sm
      ! digamma
      ! qqb->gamma+gamma

      sigma(1)=(t/u + u/t)
      sigma(1)=sigma(1)*e**4/n      ! symm. factor 1/2 due to 2 identical 
                                    ! photons is included here


!gg sm
      ! digamma
      ! gg->gamma+gamma

      sigma(2)=smsigma_gg(s,t,u)
c      write(*,*)'Sigma_gg=',sigma(2)

      ! digamma
      ! qqb->gamma+gamma   BSM

      sigma(3)=u*t*(u*u + t*t)
      sigma(3)=sigma(3)*ak2D**2/16.d0/n
                                ! symm. factor 1/2 due to 2 identical 
                                ! photons is included here
    
     
      ! gg->gamma+gamma   BSM

      sigma(4)=u**4 + t**4
      sigma(4)=sigma(4)*ak2D**2/16.d0/(n*n-1.0d0) 
                                ! symm. factor 1/2 due to 2 identical 
                                ! photons is included here
c      sigma(4)=sigma(4)+bsmsigma_gg_intf(s,t,u)
      ! digamma :end
    
      ! qqb->gamma+gamma BSM interference with sm

      sigma(5)=-(u*u + t*t)/4.0d0
      sigma(5)=sigma(5)*2.0d0*AK2DINTF/n
                                ! symm. factor 1/2 due to 2 identical 
                                ! photons is included here
                                ! inteference has got factor 2

      ! gg->gamma+gamma BSM interference with sm
      sigma(6)=bsmsigma_gg_intf(s,t,u)
c      write(*,*)'bsm_intf_gg_matrix=',sigma(6)

      return
      end


