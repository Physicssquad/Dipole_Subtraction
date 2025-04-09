        program compare
        implicit double precision (a-h,o-z)
        dimension xqLO(1:50),xintLO(1:50),xqvir(1:50),xintvir(1:50),
     .    xqreal(1:50),xintreal(1:50),xqPK(1:50),xintPK(1:50)
     .                                ,xqPKterm2(1:50),xintPKterm2(1:50)
     .                                ,xqLO_ch(1:50),xintLO_ch(1:50)   
     .                                ,xqch(1:50),xintch(1:50)   
     .                                ,xqLO_ch2(1:50),xintLO_ch2(1:50)  
     .                                ,err1(1:50),err2(1:50),err3(1:50) 
        character*100 run_tag
        character*50 name,firstfile,secondfile,thirdfile,fourthfile
     .               ,decision
     . ,real_dipole,virtual,PK,LO,ref

        
      open(unit=15,file='../../run.machine.dat',status='unknown')
      read (15,*) mid                   ! machine id Tevatron:0 LHC:1
      read (15,*) ecm                   ! ecm
      read (15,*) name                  !lhapdf set
      read (15,*) it_max                !it_max no of q for distribution
      read (15,*) xq_initial            ! initialise xq value
      read (15,*) step_size             ! size in the multiplle of loop variable 
      read (15,*) run_tag               ! name of run directory to save output
      close(15)

      open(unit=20,file='../../output_files.dat',status='unknown')
      read (20,*) real_dipole
      read (20,*) virtual 
      read (20,*) PK
      read (20,*) LO
      read (20,*) ref
      read (20,*) firstfile
      read (20,*) secondfile
      read (20,*) thirdfile
      close(20)
      fourthfile = PK


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

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     Ready to start comparing files    
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       print*,"Available files in "//trim(run_tag)//"  are"
       call system ("cd ../"//trim(run_tag)//" && ls -ltr")

       print*,"These file names are supplied by ../output_files.dat"
       print*,firstfile
       print*,secondfile
       print*,thirdfile
       print*,"will be combined to the output filename"
       print*,fourthfile
       print*," "
       print*,"Shall I combine ?       [y/n]"
       read*,decision
       if (decision .eq. 'n') then
       print*,"Enter the first file name"
       read*,firstfile
       print*,"Enter the second file name"
       read*,secondfile
       print*,"Enter the third file name"
       read*,thirdfile
       print*,"Enter the name of file to combine all"
       read*,fourthfile
       endif


c~~~~~~~~~~~~~~~~~[ first  file ]        

      call system("test -f ../"//trim(run_tag)//"/"//trim(firstfile)//
     . " && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr1
       close(13)
       call system("rm command.txt")
       if(ierr1 .eq. 1) then
        open(unit=17,file='../'//trim(run_tag)//'/'//trim(firstfile),
     .     status='unknown')
        do i=1,it_max
        read(17,*) xqLO(i),xintLO(i),err1(i)
        enddo
        close(17)
        print*,"/"//trim(firstfile)
       write(*,*)achar(27)//'[1;32m'//"   xq        Integral firstfile",
     . achar(27)//'[0m'
        do i=1,it_max
        write(*,'(i7,3e27.15)')int(xqLO(i)),xintLO(i),err1(i)    
        enddo
      endif

c~~~~~~~~~~~~~~~~~[ second file ]        

      call system("test -f ../"//trim(run_tag)//"/"//trim(secondfile)//
     . " && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr2
       close(13)
       call system("rm command.txt")
       if(ierr2 .eq. 1) then
        open(unit=17,file='../'//trim(run_tag)//'/'//trim(secondfile)
     .     ,status='unknown')
        do i=1,it_max
        read(17,*) xqLO_ch(i),xintLO_ch(i),err2(i)
        enddo
        close(17)
        print*,"/"//trim(secondfile)
       write(*,*)achar(27)//'[1;32m'//"   xq       integral secondfile",
     . achar(27)//'[0m'
        do i=1,it_max
        write(*,'(i7,3e27.15)')int(xqLO_ch(i)),xintLO_ch(i),err2(i)
        enddo
      endif
c~~~~~~~~~~~~~~~~~[ third file ]        

      call system("test -f ../"//trim(run_tag)//"/"//trim(thirdfile)//
     . " && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr3
       close(13)
       call system("rm command.txt")
       if(ierr3 .eq. 1) then
        open(unit=17,file='../'//trim(run_tag)//'/'//trim(thirdfile)
     .     ,status='unknown')
        do i=1,it_max
        read(17,*) xqLO_ch2(i),xintLO_ch2(i),err3(i)
        enddo
        close(17)
        print*,"/"//trim(thirdfile)
       write(*,*)achar(27)//'[1;32m'//"   xq       integral thirdfile",
     . achar(27)//'[0m'
        do i=1,it_max
        write(*,'(i7,3e27.15)')int(xqLO_ch2(i)),xintLO_ch2(i),err3(i)
        enddo
      endif


      
c~~~~~~~~~~~~~~~~~[ combining all ]

       if(ierr1 + ierr2 +ierr3 .eq. 3) then
               print*," "
        print*,"   /"//trim(firstfile)
        print*,"+  /"//trim(secondfile)
        print*,"+  /"//trim(thirdfile)
        print*,"=  /"//trim(fourthfile)
               print*," "


c       write(*,*)achar(27)//'[1;32m'//"   xq         first / second",
c       write(*,*)achar(27)//'[1;32m'//"   xq
c     . [first-second]/first*100",
c     . achar(27)//'[0m'
        open(unit=17,file='../'//trim(run_tag)//'/'//trim(fourthfile)
     .     ,status='unknown')
        do i=1,it_max
c      write(*,'(i7,3f10.6)')int(xqLO_ch(i)),dabs(xintLO(i)-xintLO_ch(i))
      write(*,'(i7,3e27.15)')int(xqLO_ch(i))
     .         , xintLO(i)+xintLO_ch(i)+xintLO_ch2(i)
     .         ,err1(i)+err2(i)+err3(i)

      write(17,'(i7,3e27.15)')int(xqLO_ch(i))
     .         ,xintLO(i)+xintLO_ch(i)+xintLO_ch2(i)
     .         ,err1(i)+err2(i)+err3(i)
        enddo
      write(17,*)"This data is combined from "
      write(17,*)trim(firstfile)//" + "//trim(secondfile)//" + "
     . //trim(thirdfile)//" = "//trim(fourthfile)
        close(17)
        else
            print*,"                  ⚠️"
           Print*,"Opps one or more filenames are not correct"
           print*,"or maybe one of the file doesn't exist "
           print*,"Take a closer look 👁️ ^👁️"
            print*,"                  "
           stop
      endif
        stop
123        continue
        end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
