        program compare
        implicit double precision (a-h,o-z)
        dimension xqLO(1:50),xintLO(1:50),xqvir(1:50),xintvir(1:50),
     .    xqreal(1:50),xintreal(1:50),xqPK(1:50),xintPK(1:50)
     .                                ,xqPKterm2(1:50),xintPKterm2(1:50)
     .                                ,xqLO_ch(1:50),xintLO_ch(1:50)   
     .                                ,xqch(1:50),xintch(1:50)   
        character*100 run_tag
        character*50 name,firstfile,secondfile

        
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

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     Ready to start comparing files    
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c       print*,"Available files in "//trim(run_tag)//"  are"
c       call system ("cd ../"//trim(run_tag)//" && ls -ltr")
c       print*,"Enter the first file name"
c       read*,firstfile
c       print*,"Enter the second file name"
c       read*,secondfile

	firstfile = 'LO_all.dat'
	secondfile = 'LO_ref.dat'

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
        read(17,*) xqLO(i),xintLO(i)
        enddo
        close(17)
        print*,"/"//trim(firstfile)
       write(*,*)achar(27)//'[1;32m'//"   xq        Integral firstfile",
     . achar(27)//'[0m'
        do i=1,it_max
        write(*,'(i7,3e27.15)')int(xqLO(i)),xintLO(i)
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
        read(17,*) xqLO_ch(i),xintLO_ch(i)
        enddo
        close(17)
        print*,"/"//trim(secondfile)
       write(*,*)achar(27)//'[1;32m'//"   xq       integral secondfile",
     . achar(27)//'[0m'
        do i=1,it_max
        write(*,'(i7,3e27.15)')int(xqLO_ch(i)),xintLO_ch(i)
        enddo
      endif

      
c~~~~~~~~~~~~~~~~~[ ratio ]        

       if(ierr1 + ierr2 .eq. 2) then
        print*,"/"//trim(firstfile)//" and  /"//trim(secondfile)
c       write(*,*)achar(27)//'[1;32m'//"   xq         first / second",
       write(*,*)achar(27)//'[1;32m'//"   xq
     . [first second]",
     . achar(27)//'[0m'
        do i=1,it_max
c      write(*,'(i7,3f10.6)')int(xqLO_ch(i)),xintLO(i)+xintLO_ch(i)
c      write(*,'(i7,3f10.6)')int(xqLO_ch(i)),xintLO(i)/xintLO_ch(i)
      write(*,'(i7,3e27.15)')int(xqLO_ch(i)),xintLO_ch(i)/xintLO(i)
c      write(*,'(i7,3e27.15)')int(xqLO_ch(i)),
c     .          dabs(xintLO(i)-xintLO_ch(i))
c     .     /xintLO_ch(i)*100d0
        enddo
	print*," "
       write(*,*)achar(27)//'[1;32m'//"   xq
     . [second first]",
     . achar(27)//'[0m'
        do i=1,it_max
c      write(*,'(i7,3f10.6)')int(xqLO_ch(i)),xintLO(i)+xintLO_ch(i)
c      write(*,'(i7,3e27.15)')int(xqLO_ch(i)),xintLO(i)/xintLO_ch(i)
      write(*,'(i7,3f10.6)')int(xqLO_ch(i)),xintLO(i)/xintLO_ch(i)
c      write(*,'(i7,3e27.15)')int(xqLO_ch(i)),
c     .          dabs(xintLO(i)-xintLO_ch(i))
c     .     /xintLO_ch(i)*100d0
        enddo

      endif
        stop
123        continue
        end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
