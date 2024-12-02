      subroutine pdft_Pgg(x,xmu,nt,nf,deltas,deltac,s12,fin,fout)
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
      
       
c     call pdf(nt,x/y,xmu,ftmp)
      call pdf(x/y,xmu,ftmp)
      
      oms=1-deltas
      IF((x.le.oms).and.(y.le.oms))THEN         
         f34=gfun(oms,atmp)-gfun(x,atmp)
         if(icorr.eq.1)then
            f34=f34-(ffun(oms,deltac)-ffun(x,deltac))
            xly=xly-dlog(1-deltac/omy)
         endif

         ptgg=xly*pgg0(y)-pgg1(y)

         f34=f34/(omx-deltas)
         
         fout(0) = 2*n*fin(0)*f34
     &        + (ftmp(0)*ptgg-fin(0)*2*n*xly/omy)/y
         
         fout(0) = fout(0)*xjac     
         
      ENDIF

      return
      end

