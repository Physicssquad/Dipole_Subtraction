      program intPK
      implicit double precision (a-h,o-z)
      dimension c(1:2)
      character*50 name,mode
      character*100 command,run_tag,dir_path,filename,filename1
      common/energy/s
      external flo2_PK
      common/leg_choice/leg
      common/usedalpha/AL,ge   
      common/distribution/xq
      dimension PK(1:50),err(1:50)
      

      !input data card
      open(unit=10,file='../run.vegas.dat',status='unknown')
      do i=1,6
      read (10,*)
      enddo
      read (10,*) pt1           ! vegas points     
      read (10,*) its1          ! vegas iterations 
      npt1 = pt1
      close(10)
      

      open(unit=10,file='../param_card.dat',status='unknown')    
      read (10,*) ge          ! [ 1/Alpha_ew ]
      close(10)

      open(unit=15,file='../run.machine.dat',status='unknown')
      read (15,*) mid           ! machine id Tevatron:0 LHC:1
      read (15,*) ecm           ! ecm
      read (15,*) name          ! lhapdf set
      read (15,*) it_max        ! it_max no of q for distribution
      read (15,*) xq            ! initialise xq value
      read (15,*) xincr         ! increment in Gev from xq 
      read (15,*) run_tag
      read (15,*) iprint        ! to print data in file
      close(15)

      open(unit=20,file='../output_files.dat',status='unknown')
      read (20,*)
      read (20,*) 
      read (20,*) filename
      close(20)



        xq_initial = xq
      call initpdfsetbyname(name)
      Call initPDF(0)
       s=ecm*ecm

        print*,"  "
        print*,"____________________________________"
        print*,"    Calculating P and K terms       "
        print*,"____________________________________"
        print*,"````````````````````````````````````"

c ~~~~~~~~~~~Writing in a file to store~~~~~~~~~~~~c        
      if (iprint .eq. 1) call output(run_tag,filename)            

      call printframe1(pt1,its1)   ! Prints Vegas points
        do j=1,it_max
      call printframe2(xq)

          call brm48i(40,0,0) ! initialize random number generator
          call vsup(4,npt1,its1,flo2_PK,ai_lo2,sd,chi2)

          PK(j) = ai_lo2
          err(j)=sd

      mode = "[P+K]"
      call printframe3(mode,ai_lo2,sd,chi2)   ! Prints Summary of integration

           xq=xq + xincr
        enddo

        xq = xq_initial
      call printframe4(mode)
        
        do j=1,it_max
          write(*,'(i7,3e27.15,3e27.15)')
     .             int(xq),PK(j),err(j)
          xq = xq + xincr
        enddo
c~~~~~~~~~~~~~~~~~~~~~~~~~ * END * ~~~~~~~~~~~~~~~~~~~~~~~~~~c       

c ~~~~~~~~~~~Writing in a file to compare~~~~~~~~~~~~c        
        if (iprint .eq. 0 ) goto 123
       open(unit=21,file='../summary/'//trim(run_tag)//
     .   '/'//trim(filename),status='unknown')
c     .   '/'//trim(filename),status='unknown', access='append')
         xq = xq_initial
         do i=1,it_max
          write(21,*)xq,PK(i),err(i)
          xq = xq + xincr
         enddo
         close(21)
123        continue
       end

      double precision function dot(p1,p2)
      implicit none
      double precision p1(0:3),p2(0:3)
      dot=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      end
