      subroutine pdft_Pqq(x,xmu,nt,nf,deltas,deltac,s12,fin,fout)
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

      call pdf(x/y,xmu,ftmp)
      fsum=0
      do i=1,nf
         fsum=fsum+ftmp(-i)+ftmp(+i)
      enddo

      oms=1-deltas

      if((x.le.oms).and.(y.le.oms))then
         f34=gfun(oms,atmp)-gfun(x,atmp)
         if(icorr.eq.1)then
            f34=f34-(ffun(oms,deltac)-ffun(x,deltac))
            xly=xly-dlog(1-deltac/omy)
         endif
         
         f34=f34/(omx-deltas)
         ptqq=xly*pqq0(y)-pqq1(y)
         
         do i=1,nf
            j=+i
            fout(j) = 
     .           2*cf*fin(j)*f34 
     .           + (ftmp(j)*ptqq-fin(j)*2*cf*xly/omy)/y 
            fout(j) = fout(j)*xjac
            j=-i
            fout(j) = 
     .           2*cf*fin(j)*f34 
     .           + (ftmp(j)*ptqq-fin(j)*2*cf*xly/omy)/y 
            fout(j) = fout(j)*xjac
         enddo         
      endif

      return
      end

    
