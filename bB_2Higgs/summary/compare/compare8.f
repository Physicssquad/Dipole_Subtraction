        program compare
        implicit double precision (a-h,o-z)
        dimension xqLO(1:50),xintLO(1:50),xqvir(1:50),xintvir(1:50),
     .    xqreal(1:50),xintreal(1:50),xqPK(1:50),xintPK(1:50)
     .                                ,xqPKterm2(1:50),xintPKterm2(1:50)
     .                                ,xqch(1:50),xintch(1:50)
     .   ,xLO_err(1:50),xVir_err(1:50),xreal_err(1:50),xPK_err(1:50)
        character*100 run_tag
        character*50 name,yes
        character*100 message
        character*100 real_dipole,virtual,PK,ref,LO

        
c       Leading Order
      open(unit=15,file='../../run.machine.dat',status='unknown')
      read (15,*) mid                   ! machine id Tevatron:0 LHC:1
      read (15,*) ecm                   ! ecm
      read (15,*) name                  !lhapdf set
      read (15,*) it_max                !it_max no of q for distribution
      read (15,*) xq_initial            ! initialise xq value
      read (15,*) step_size             ! size in the multiplle of loop variable 
      read (15,*) run_tag               ! name of run directory to save output
      close(15)

c      open(unit=20,file='../../output_files2.dat',status='unknown')
      open(unit=20,file='../../output_files.dat',status='unknown')
      read (20,*) real_dipole
      read (20,*) virtual 
      read (20,*) PK
      read (20,*) LO
      read (20,*) ref
      close(20)


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
       print*,"Enter short text about this run for future references."
       print*,"Run Discription:"
       read(*, '(A)', advance='yes')message
       call system("echo " //message// " >> ../" //trim(run_tag)//
     .    "/run.machine.dat")
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
      read (15,*) iprint
      read (15, '(A)', iostat=ierr) message  ! this identifies the message for the process
      close(15)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


      print*," "
      print*,"Reading Data from directory: /summary/"//trim(run_tag)
      print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
      call system("cat ../"//trim(run_tag)//"/run.machine.dat")
c      print*,"            ecm:", int(ecm),"[GeV]"                   ! ecm
c      print*,"     LHApdfname:   ", name                  !lhapdf set
c      print*,"         it_max:" ,int(it_max)                !it_max no of q for distribution
c      print*,"initial Q value:", int(xq_initial),"[GeV]"            ! initialise xq value
c      print*,"     step size :", int(step_size)             ! size in the multiplle of loop variable
c      print*,"Run Description:   ",message
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
        read(17,*) xqLO(i),xintLO(i),xLO_err(i)
        enddo
        close(17)
        print*,"/LO.dat"
       write(*,*)achar(27)//'[1;32m'//"   xq"," ","        Integral_LO",
     .     "                   error",
     . achar(27)//'[0m'
        do i=1,it_max
        write(*,'(i7,3e27.15,3e27.15)')int(xqLO(i)),xintLO(i),xLO_err(i)
        enddo
      endif

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c                [Virtual contribution]
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       call system("test -f ../"// trim(run_tag) //"/"//trim(virtual)//
     . " && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr
       close(13)
       call system("rm command.txt")
       if(ierr .eq. 1) then
               inlo=inlo+1
        open(unit=17,file='../'//trim(run_tag)//'/'//trim(virtual)
     .     ,status='unknown')
        do i=1,it_max
        read(17,*) xqVir(i),xintVir(i),xVir_err(i)
c        read(17,*) xqVir(i),xintVir(i)
        enddo
        close(17)
        print*," "
        print*,"/"//trim(virtual)
       write(*,*)achar(27)//'[1;32m'//"   xq           Integral_VIR",
     .     "                   error",
     . achar(27)//'[0m'
        do i=1,it_max
           write(*,'(i7,3e27.15,3e27.15)')int(xqVir(i)),
     .     xintVir(i),xVir_err(i)
        enddo
        endif

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c                  ["Real - Dipole"]
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       call system("test -f ../"// trim(run_tag) //"/"
     .   //trim(real_dipole)//
     . " && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr
       close(13)
       call system("rm command.txt")
       if(ierr .eq. 1) then
               inlo=inlo+1
        open(unit=17,file='../'//trim(run_tag)//'/'
     .   //trim(real_dipole)
     .     ,status='unknown')
        do i=1,it_max
        read(17,*) xqreal(i),xintreal(i),xreal_err(i)
c        read(17,*) xqreal(i),xintreal(i)
        enddo
        close(17)

        print*," "
        print*,"/"//trim(real_dipole)
       write(*,*)achar(27)//'[1;32m'//"   xq           Integral Real",
     .  "                  error",
     . achar(27)//'[0m'
        do i=1,it_max
          write(*,'(i7,3e27.15,3e27.15)')int(xqreal(i)),xintreal(i)
     .       ,xreal_err(i)
        enddo
        endif
9899    continue
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c                          PK terms
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       call system("test -f ../"// trim(run_tag) //"/"//trim(PK)//
     . " && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr
       close(13)
       call system("rm command.txt")
       if(ierr .eq. 1) then
               inlo=inlo+1
        open(unit=17,file='../'//trim(run_tag)//'/'//trim(PK),
     .     status='unknown')
        do i=1,it_max
        read(17,*) xqPK(i),xintPK(i),xPK_err(i)
c        read(17,*) xqPK(i),xintPK(i)
        enddo
        close(17)
        print*," "
        print*,"/"//trim(PK)
       write(*,*)achar(27)//'[1;32m'//"   xq","          Integral PK",
     .  "                     error",
     . achar(27)//'[0m'
        do i=1,it_max
        write(*,'(i7,3e27.15,3e27.15)')int(xqPK(i)),xintPK(i),xPK_err(i)
        enddo
        print*," "
        print*,"Do you want to see individual data of PK ? 
     .          [y/n]"
        read*,yes
        if (yes .eq. 'y' .or. yes .eq. 'Y') then 
       
        call compare6  

        endif
        
        print*," "

        elseif( ierr .eq. 0 ) then 

        print*," /"//trim(PK)//" file doesn't exist ‼️  "
        print*,"Do you want to combine now ?                      [y/n]"
        read*,yes
        if (yes .eq. 'y' .or. yes .eq. 'Y') then 
                call compare5
                goto 9899 
        else
                stop
        endif 

        endif

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c       Total sig_NLO will be
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       if (inlo .eq. 3) then
       print*," "
       write(*,*)achar(27)//'[1;32m'//"   xq  ",
     . "       sigma NLO dipole","                 error",
     . achar(27)//'[0m'
        call sleep(1)
        xq = xq_initial
        do i=1,it_max
          write(*,'(i7,3e27.15,3e27.15)')int(xq),
     .  xintPK(i)+xintvir(i)+xintreal(i)
     .   ,xVir_err(i)+xreal_err(i)+xPK_err(i)
c     .  xintvir(i)+xintreal(i))
c          write(*,'(i7,3e27.15)')int(xqPKterm2(i)),
c     .  xintvir(i)+xintreal(i)+xintLO(i)
        xq = xq + step_size
        enddo
        endif

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c       "Sigma_Chinmoy" and Ratio with Dipole
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        inloch = 0
       call system("test -f ../"// trim(run_tag) //"/"//trim(ref)//
     . " && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr
       close(13)
       call system("rm command.txt")
       if(ierr .eq. 1) then
               inloch = 1
        open(unit=17,file='../'//trim(run_tag)//'/'//trim(ref)
     .     ,status='unknown')
        do i=1,it_max
        read(17,*) xqch(i),xintch(i)
        enddo
        close(17)
       print*," "
       print*,"For comparison with other data "
       print*," "
       print*,"/"//trim(ref)
       write(*,*)achar(27)//'[1;32m'//"   xq  ",
     ."     sigma slicing code", 
     . achar(27)//'[0m'
        do i=1,it_max
          write(*,'(i7,3e27.15)')int(xqch(i)),xintch(i)
        enddo
        endif
c~~~~~~~~~~~~~~~~~[ ratio ]        

        if (inlo + inloch .eq. 4) then
       print*," "
        write(*,*)achar(27)//'[1;32m'//"   xq  ",
     .  "  slicing/dipole sigma_NLO, and dipole/slicing",
     . achar(27)//'[0m'
        xq = xq_initial
        do i=1,it_max
          write(*,'(i7,3e27.15)')int(xq),
     .  xintch(i)/(xintPK(i)+xintvir(i)+xintreal(i)),
     .  (xintPK(i)+xintvir(i)+xintreal(i))/xintch(i)
c     .  (xintPK(i)+xintvir(i)+xintreal(i))/xintch(i)
c          write(*,'(i7,3e27.15)')int(xqPKterm2(i)),
c     .  xintvir(i)+xintreal(i)+xintLO(i)
        xq = xq + step_size
        enddo
        endif
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
123        continue
        end

        subroutine compare5
        implicit double precision (a-h,o-z)
        dimension xqLO(1:50),xintLO(1:50),xqvir(1:50),xintvir(1:50),
     .    xqreal(1:50),xintreal(1:50),xqPK(1:50),xintPK(1:50)
     .                                ,xqPKterm2(1:50),xintPKterm2(1:50)
     .                                ,xqLO_ch(1:50),xintLO_ch(1:50)   
     .                                ,xqch(1:50),xintch(1:50)   
     .                                ,xqLO_ch2(1:50),xintLO_ch2(1:50)  
     .                                ,err1(1:50),err2(1:50),err3(1:50) 
        character*100 run_tag,location,fpath
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
      print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
      print*,"Reading Data from directory: /summary/"//trim(run_tag)
      print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
c    print*,"            ecm:", int(ecm),"[GeV]"                   ! ecm
c    print*,"     LHApdfname:   ", name                  !lhapdf set
c    print*,"         it_max:" ,int(it_max)                !it_max no of q for distribution
c    print*,"initial Q value:", int(xq_initial),"[GeV]"            ! initialise xq value
c    print*,"     step size :", int(step_size)             ! size in the multiplle of loop variable
c    print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
c    print*," "

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     Ready to start combining files Specific for PK Only   
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       print*,"Available files in "//trim(run_tag)//"  are"
       call system ("cd ../"//trim(run_tag)//" && ls -ltr")

       print*," "
       print*," "
       print*,"Following file names are supplied by ../output_files.dat"
       print*,firstfile
       print*,secondfile
       print*,thirdfile
       print*,"will be combined to the output filename"
       print*,fourthfile
       print*," "
       print*,"Shall I combine ?       [y/n]"
       read*,decision
       if (decision .eq. 'n') then
       print*,"Do you want to Enter filenames manually?           [y/n]"
       read*,decision
       if ( decision .eq. 'y') then
       print*,"Enter the first file name"
       read*,firstfile
       print*,"Enter the second file name"
       read*,secondfile
       print*,"Enter the third file name"
       read*,thirdfile
       print*,"Enter the name of file to combine all"
       read*,fourthfile
       else
               stop
       endif
       endif


c~~~~~~~~~~~~~~~~~[ first  file ]        

       location =  '../../PK_Isolated/summary'
      call system("test -f " //trim(location)//"/"//trim(run_tag)//"/"
     . //trim(firstfile)//
     . " && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr1
       close(13)
       call system("rm command.txt")
       if(ierr1 .eq. 1) then
      fpath=  trim(location)//'/'//trim(run_tag)//'/'//trim(firstfile)
        open(unit=17,file=fpath,status='unknown')
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

      call system("test -f " //trim(location)//"/"//trim(run_tag)//"/"
     . //trim(secondfile)//
     . " && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr2
       close(13)
       call system("rm command.txt")
       if(ierr2 .eq. 1) then
      fpath=  trim(location)//'/'//trim(run_tag)//'/'//trim(secondfile)
        open(unit=17,file=fpath,status='unknown')
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

      call system("test -f " //trim(location)//"/"//trim(run_tag)//"/"
     . //trim(thirdfile)//
     . " && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr3
       close(13)
       call system("rm command.txt")
       if(ierr3 .eq. 1) then
      fpath=  trim(location)//'/'//trim(run_tag)//'/'//trim(firstfile)
        open(unit=17,file=fpath,status='unknown')
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


       open(unit=17,file='../'//trim(run_tag)//'/'//trim(fourthfile)
     .     ,status='unknown')
        do i=1,it_max

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
           print*,"or maybe some of the files doesn't exist "
           print*,"Take a closer look 👁️ ^👁️"
            print*,"                  "
           stop
      endif
123        continue
        end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        subroutine compare6
        implicit double precision (a-h,o-z)
        dimension xqLO(1:50),xintLO(1:50),xqvir(1:50),xintvir(1:50),
     .    xqreal(1:50),xintreal(1:50),xqPK(1:50),xintPK(1:50)
     .                                ,xqPKterm2(1:50),xintPKterm2(1:50)
     .                                ,xqLO_ch(1:50),xintLO_ch(1:50)   
     .                                ,xqch(1:50),xintch(1:50)   
     .                                ,xqLO_ch2(1:50),xintLO_ch2(1:50)  
     .                                ,err1(1:50),err2(1:50),err3(1:50) 
        character*100 run_tag,location,fpath
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
c
       print*," "
       print*," "
       print*,"These file names are supplied by ../output_files.dat"
       print*,firstfile
       print*,secondfile
       print*,thirdfile
c       print*,"will be combined to the output filename"
c       print*,fourthfile
       print*," "
c       print*,"Shall I combine ?       [y/n]"
c       read*,decision
c       if (decision .eq. 'n') then
       print*,"Do you want to Enter filenames manually?           [y/n]"
       read*,decision
       if ( decision .eq. 'y') then
       print*,"Enter the first file name"
       read*,firstfile
       print*,"Enter the second file name"
       read*,secondfile
       print*,"Enter the third file name"
       read*,thirdfile
c       print*,"Enter the name of file to combine all"
c       read*,fourthfile
c       else
c               stop
c       endif
       endif
c

c~~~~~~~~~~~~~~~~~[ first  file ]        

       location =  '../../PK_Isolated/summary'
      call system("test -f " //trim(location)//"/"//trim(run_tag)//"/"
     . //trim(firstfile)//
     . " && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr1
       close(13)
       call system("rm command.txt")
       if(ierr1 .eq. 1) then
      fpath=  trim(location)//'/'//trim(run_tag)//'/'//trim(firstfile)
        open(unit=17,file=fpath,status='unknown')
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

      call system("test -f " //trim(location)//"/"//trim(run_tag)//"/"
     . //trim(secondfile)//
c      call system("test -f ../"//trim(run_tag)//"/"//trim(secondfile)//
     . " && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr2
       close(13)
       call system("rm command.txt")
       if(ierr2 .eq. 1) then
      fpath=  trim(location)//'/'//trim(run_tag)//'/'//trim(secondfile)
        open(unit=17,file=fpath,status='unknown')
c        open(unit=17,file=//trim(location)//'/'//trim(run_tag)//'/'
c     . //trim(secondfile)
c        open(unit=17,file='../'//trim(run_tag)//'/'//trim(secondfile)
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

      call system("test -f " //trim(location)//"/"//trim(run_tag)//"/"
     . //trim(thirdfile)//
c      call system("test -f ../"//trim(run_tag)//"/"//trim(thirdfile)//
     . " && echo 1 > command.txt || echo 0  > command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr3
       close(13)
       call system("rm command.txt")
       if(ierr3 .eq. 1) then
      fpath=  trim(location)//'/'//trim(run_tag)//'/'//trim(firstfile)
        open(unit=17,file=fpath,status='unknown')
c        open(unit=17,file='../'//trim(run_tag)//'/'//trim(thirdfile)
c        open(unit=17,file=//trim(location)//'/'//trim(run_tag)//'/'
c     . //trim(thirdfile)
c     .     ,status='unknown')
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
c
       if(ierr1 + ierr2 +ierr3 .lt. 3) then
               print*," "
           print*,"                  ⚠️"
           Print*,"Opps one or more filenames are not correct"

           stop
      endif
123        continue
        end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
