        program compare
        implicit double precision (a-h,o-z)
        dimension xqLO(1:50),xintLO(1:50),xqvir(1:50),xintvir(1:50),
     .    xqreal(1:50),xintreal(1:50),xqPK(1:50),xintPK(1:50)
     .                                ,xqPKterm2(1:50),xintPKterm2(1:50)
     .                                ,xqLO_ch(1:50),xintLO_ch(1:50)   
     .                                ,xqch(1:50),xintch(1:50)   
        character*100 run_tag
        character*50 name

        
      open(unit=15,file='../../run.machine.dat',status='unknown')
      read (15,*) mid                   ! machine id Tevatron:0 LHC:1
      read (15,*) ecm                   ! ecm
      read (15,*) name                  !lhapdf set
      read (15,*) it_max                !it_max no of q for distribution
      read (15,*) xq_initial            ! initialise xq value
      read (15,*) step_size             ! size in the multiplle of loop variable 
      read (15,*) run_tag               ! name of run directory to save output
      close(15)

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c checking existance of dir and file run.machine.dat
c  Setting up environment      
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       inlo=0
       call system("test -d ../"// trim(run_tag) //
     . " && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr1
       close(13)
       call system("rm command.txt")
       if(ierr1 .eq. 0) then 
               print*,"Directory not found :" //run_tag
               print*,"Currently stored data are in:"
               call system("cd ../ && ls")
               goto 123 
       endif
       call system("test -f ../"// trim(run_tag) //"/run.machine.dat 
     . && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr
       close(13)
       call system("rm command.txt")
       if(ierr .eq. 0 ) then
       call system("cp ../../run.machine.dat ../"// trim(run_tag))
       endif
        
      open(unit=15,file='../'//trim(run_tag)//'/run.machine.dat'
     . ,status='unknown')
      read (15,*) mid                   ! machine id Tevatron:0 LHC:1
      read (15,*) ecm                   ! ecm
      read (15,*) name                  !lhapdf set
      read (15,*) it_max                !it_max no of q for distribution
      read (15,*) xq_initial            ! initialise xq value
      read (15,*) step_size             ! size in the multiplle of loop variable
      read (15,*) run_tag               ! name of run directory to save output
      close(15)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      print*," "
      print*,"Reading Data from directory: /summary/"//trim(run_tag)
      print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
      print*,"            ecm:", int(ecm),"[GeV]"                   ! ecm
      print*,"     LHApdfname:   ", name                  !lhapdf set
      print*,"         it_max:" ,int(it_max)                !it_max no of q for distribution
      print*,"initial Q value:", int(xq_initial),"[GeV]"            ! initialise xq value
      print*,"     step size :", int(step_size)             ! size in the multiplle of loop variable
      print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
      print*," "
      call sleep(1)

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c                [LO contribution]      
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       call system("test -f ../"// trim(run_tag) //"/LO.dat 
     . && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr
       close(13)
       call system("rm command.txt")
       if(ierr .eq. 1) then
        open(unit=17,file='../'//trim(run_tag)//'/LO.dat',
     .     status='unknown')
        do i=1,it_max
        read(17,*) xqLO(i),xintLO(i)
        enddo
        close(17)
        print*,"/LO.dat"
       write(*,*)achar(27)//'[1;32m'//"   xq"," ","        Integral_LO",
     . achar(27)//'[0m'
        do i=1,it_max
        write(*,'(i7,3e27.15)')int(xqLO(i)),xintLO(i)
        enddo
c        call sleep(1)
      endif

c~~~~~~~~~~~~~~~~~[ Chinmoy LO ]        

             call system("test -f ../"// trim(run_tag) //"/LO2.dat
     . && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr
       close(13)
       call system("rm command.txt")
       if(ierr .eq. 1) then
        open(unit=17,file='../'//trim(run_tag)//'/LO2.dat',
     .     status='unknown')
        do i=1,it_max
        read(17,*) xqLO_ch(i),xintLO_ch(i)
        enddo
        close(17)
        print*,"/LO2.dat"
       write(*,*)achar(27)//'[1;32m'//"   xq"," ","        Chinmoy_LO",
     . achar(27)//'[0m'
        do i=1,it_max
        write(*,'(i7,3e27.15)')int(xqLO_ch(i)),xintLO_ch(i)
        enddo
c        call sleep(1)
      endif

      
c~~~~~~~~~~~~~~~~~[ ratio ]        

       call system("test -f ../"// trim(run_tag) //"/LO2.dat
     . && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr
       close(13)
       call system("rm command.txt")
       if(ierr .eq. 1) then
        open(unit=17,file='../'//trim(run_tag)//'/LO2.dat',
     .     status='unknown')
        do i=1,it_max
        read(17,*) xqLO_ch(i),xintLO_ch(i)
        enddo
        close(17)
        print*,"/LO2.dat"
       write(*,*)achar(27)//'[1;32m'//"   xq          LO my/Chinmoy_LO",
     . achar(27)//'[0m'
        do i=1,it_max
        write(*,'(i7,3e27.15)')int(xqLO_ch(i)),xintLO(i)/xintLO_ch(i)
        enddo
c        call sleep(1)
        stop
      endif


c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c                [Virtual contribution]
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       call system("test -f ../"// trim(run_tag) //"/virtual.dat 
     . && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr
       close(13)
       call system("rm command.txt")
       if(ierr .eq. 1) then
               inlo=inlo+1
        open(unit=17,file='../'//trim(run_tag)//'/virtual.dat',
     .     status='unknown')
        do i=1,it_max
        read(17,*) xqVir(i),xintVir(i)
c        read(17,*) xqVir(i),xintVir(i)
        enddo
        close(17)
        print*," "
        print*,"/virtual.dat"
       write(*,*)achar(27)//'[1;32m'//"   xq           Integral_VIR",
     . achar(27)//'[0m'
        do i=1,it_max
          write(*,'(i7,3e27.15)')int(xqVir(i)),xintVir(i)
        enddo
        call sleep(1)
        endif

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c                  ["Real - Dipole"]
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       call system("test -f ../"// trim(run_tag) //"/real.dat 
     . && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr
       close(13)
       call system("rm command.txt")
       if(ierr .eq. 1) then
               inlo=inlo+1
        open(unit=17,file='../'//trim(run_tag)//'/real.dat',
     .     status='unknown')
        do i=1,it_max
        read(17,*) xqreal(i),xintreal(i)
c        read(17,*) xqreal(i),xintreal(i)
        enddo
        close(17)

        print*," "
        print*,"/real.dat"
       write(*,*)achar(27)//'[1;32m'//"   xq           Integral Real",
     . achar(27)//'[0m'
        do i=1,it_max
          write(*,'(i7,3e27.15)')int(xqreal(i)),xintreal(i)
        enddo
        call sleep(1)
        endif
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c                          PK terms
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       call system("test -f ../"// trim(run_tag) //"/PK.dat 
     . && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr
       close(13)
       call system("rm command.txt")
       if(ierr .eq. 1) then
               inlo=inlo+1
        open(unit=17,file='../'//trim(run_tag)//'/PK.dat',
     .     status='unknown')
        do i=1,it_max
        read(17,*) xqPK(i),xintPK(i)
c        read(17,*) xqPK(i),xintPK(i)
        enddo
        close(17)
        print*," "
        print*,"/PK.dat"
       write(*,*)achar(27)//'[1;32m'//"   xq","          Integral PK",
     . achar(27)//'[0m'
        do i=1,it_max
          write(*,'(i7,3e27.15)')int(xqPK(i)),xintPK(i)
        enddo
        call sleep(1)
        endif

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c       Total sig_NLO will be
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       if (inlo .eq. 3) then
       print*," "
       write(*,*)achar(27)//'[1;32m'//"   xq  ",
     . "     sigma NLO dipole", 
     . achar(27)//'[0m'
        xq = xq_initial
        do i=1,it_max
          write(*,'(i7,3e27.15)')int(xq),
     .  xintPK(i)+xintvir(i)+xintreal(i)!+xintLO(i))
c     .  xintvir(i)+xintreal(i))
c          write(*,'(i7,3e27.15)')int(xqPKterm2(i)),
c     .  xintvir(i)+xintreal(i)+xintLO(i)
        xq = xq + step_size
        enddo
        call sleep(1)
        endif

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c       "Sigma_Chinmoy" and Ratio with Dipole
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        inloch = 0
       call system("test -f ../"// trim(run_tag) //"/smqqb.nlo.out
     . && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr
       close(13)
       call system("rm command.txt")
       if(ierr .eq. 1) then
               inloch = 1
        open(unit=17,file='../'//trim(run_tag)//'/smqqb.nlo.out',
     .     status='unknown')
        do i=1,it_max
        read(17,*) xqch(i),xintch(i)
        enddo
        close(17)
       print*," "
       write(*,*)achar(27)//'[1;32m'//"   xq  ",
     ." sigma NLO chinmoy", 
     . achar(27)//'[0m'
        do i=1,it_max
          write(*,'(i7,3e27.15)')int(xqch(i)),xintch(i)
        enddo
        call sleep(1)
        endif
c~~~~~~~~~~~~~~~~~[ ratio ]        

        if (inlo + inloch .eq. 4) then
       call system("test -f ../"// trim(run_tag) //"/smqqb.nlo.out
     . && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr
       close(13)
       call system("rm command.txt")
       if(ierr .eq. 1) then
        print*," "
        write(*,*)achar(27)//'[1;32m'//"   xq  ",
     .  "  chinmoy/dipole sigma_NLO",
     . achar(27)//'[0m'
        xq = xq_initial
        do i=1,it_max
          write(*,'(i7,3e27.15)')int(xq),
     .  xintch(i)/(xintPK(i)+xintvir(i)+xintreal(i))
c          write(*,'(i7,3e27.15)')int(xqPKterm2(i)),
c     .  xintvir(i)+xintreal(i)+xintLO(i)
        xq = xq + step_size
        enddo
        endif
        endif
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
123        continue
        end
