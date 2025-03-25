      subroutine vegcall(IT,TI,TSI,AVGI,SD,CHI2A,WGT,SWGT)
      implicit double precision (a-h,o-z)
      integer IT,i,nbin
      double precision TI,TSI,AVGI,SD,CHI2A,CALLS,SCALLS
      double precision intres,intav,sigav,SWGT,WGT,bin2m,bin2w
      dimension dfun(100,10),dfus(100,10),dfsq(100,10),dsdi(100,10)
      common/distfun/dfun,dfsq
      common/distfus/dfus,dsdi
      common/callnu/CALLS
      common/vbin/binc,binw,nbin
      common/vegparam/npt1,its1
      common/dimen/nd
      common/isub/io,is
c... here the function is having multiple binnings whuch are being computed inside
c... a do loop.

      do j=1,nbin
      do k=is,is
       dfsq(j,k) = dfsq(j,k) - dfun(j,k)**2*(CALLS)**(-1)
      enddo
      enddo

      if (IT .eq. 1) then
      do j=1,nbin
      do k=is,is
        dfus(j,k) = dfun(j,k)*WGT
        dsdi(j,k) = dfsq(j,k)**(-1)
      enddo
      enddo

      else

      do j=1,nbin
      do k=is,is
        dfus(j,k) = dfus(j,k)+dfun(j,k)*WGT
        dsdi(j,k) = dsdi(j,k) + dfsq(j,k)**(-1)
      enddo
      enddo

      endif

      do j=1,nbin
      do k=is,is
      dfun(j,k) = 0.d0
      dfsq(j,k) = 0.d0
      enddo
      enddo

      if (IT .eq. its1) then

      do j=1,nbin
      do k=is,is
       dfus(j,k) = dfus(j,k)*SWGT**(-1)
      enddo
      enddo

      endif 
      return
      end
