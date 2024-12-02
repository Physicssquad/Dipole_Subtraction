      subroutine pdft_Pqg(x,xmu,nt,nf,deltas,deltac,s12,fin,fout)
      implicit double precision (a-h,o-z)
      dimension fin(-6:6),fout(-6:6),ftmp(-6:6)
      dimension fvec(1)
      parameter (n=3)
      parameter (cf=4d0/3d0)
      common/corr/icorr

      omx=1-x
      xjac=omx
      call brm48(fvec,1)

      y=x+xjac*fvec(1)
      atmp=deltac*s12/(xmu*xmu)
      xlt=dlog(atmp)
      omy=1-y
      xly=xlt+dlog(omy/y)

      fout(0)=0.d0
      do i=1,6
         fout(+i) = 0.d0
         fout(-i) = 0.d0         
      enddo 
 
      ptqg=xly*pqg0(y)-pqg1(y)

      call pdf(x/y,xmu,ftmp)
      fsum=0
      do i=1,nf
         fsum=fsum+ftmp(-i)+ftmp(+i)
      enddo      

      if((x.le.1.d0).and.(y.le.1.d0))then         
         do i=1,nf
            fout(+i) = xjac*ftmp(0)*ptqg/y
            fout(-i) = fout(+i)
         enddo                  
      endif
      
      return
      end

 
