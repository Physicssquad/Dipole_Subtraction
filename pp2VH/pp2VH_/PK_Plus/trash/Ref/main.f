      program intPK
      implicit double precision (a-h,o-z)
      dimension c(1:2)
      character*50 name
      character*100 command,run_tag,dir_path,filename,filename1
      common/energy/s
      external flo2_PK
      common/leg_choice/leg
      common/usedalpha/AL,ge   
      common/distribution/xq
      dimension PK(1:50),err(1:50)
      

c      include 'coupl.inc'
c      include 'nexternal.inc'
c      call setpara('param_card.dat',.true.)


      !input data card
      open(unit=10,file='../run.vegas.dat',status='unknown')
      read (10,*)
      read (10,*)
      read (10,*)
      read (10,*)
      read (10,*)
      read (10,*)
      read (10,*) pt1          ! vegas points     LO 2 body
      read (10,*) its1          ! vegas iterations LO 2 body
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


        xq_initial = xq
      call initpdfsetbyname(name)
      Call initPDF(0)
       s=ecm*ecm

        print*,"  "
        print*,"____________________________________"
        print*,"    Calculating P and K terms       "
        print*,"____________________________________"
        print*,"````````````````````````````````````"

        filename = 'PK.dat'
c ~~~~~~~~~~~Writing in a file to store~~~~~~~~~~~~c        
      if (iprint .eq. 1) call output(run_tag,filename)            

        do j=1,it_max
          print*," "
      write(*,*) achar(27)//'[1;32m' // "For xq = ",int(xq) ,achar(27) 
     .   //'[0m'
          print*," "

          call brm48i(40,0,0) ! initialize random number generator
          call vsup(4,npt1,its1,flo2_PK,ai_lo2,sd,chi2)

          PK(j) = ai_lo2
          err(j)=sd
            print*,"  "
            print*,"  "
            write(*,*)achar(27)//'[1;32m'//"Integral [P+K]=", 
     .      ai_lo2,achar(27) //'[0m', "+-",err(j)
            write(*,*)"with chisq    =",chi2
            print*," "
            print*," "

           xq=xq + xincr
        enddo
        xq = xq_initial
       write(*,*)achar(27)//'[1;32m'//"   xq"," ","   Integral PK ",
     .  "                    error",
     . achar(27)//'[0m'
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
