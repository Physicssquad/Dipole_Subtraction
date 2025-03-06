      program handle_output
       implicit double precision (a-h,o-z)
       character*50 name,filename,run_tag,filename_tmp,mode,tag
     .  ,filename1,filename2,filename3,filename4
     .  ,filename5,filename6,filename7,filename8,input_machine
     .  ,input_files
       character*100 filename_tmp1

       call get_command_argument(1,mode)
       call get_command_argument(2,tag)
       call get_command_argument(3,input_machine)
       call get_command_argument(4,input_files)

       open(unit=15,file='../../trash/broken/'//trim(input_machine)
     .     ,status='unknown')
       read (15,*) mid           ! machine id Tevatron:0 LHC:1
       read (15,*) ecm           ! ecm
       read (15,*) name          ! lhapdf set
       read (15,*) it_max        ! it_max no of q for distribution
       read (15,*) xq            ! initialise xq value
       read (15,*) xincr         ! increment in Gev from xq 
       read (15,*) run_tag
       read (15,*) iprint        ! to print data in file
       close(15)

c       open(unit=20,file='../../output_files.dat',status='unknown')
       open(unit=20,file='../../trash/broken/'//trim(input_files)
     .     ,status='unknown')
       read (20,*) filename1
       read (20,*) filename2
       read (20,*) filename3
       read (20,*) filename4
       read (20,*) filename5
       read (20,*) filename6
       read (20,*) filename7
       read (20,*) filename8
       close(20)

       if (tag .eq. 'dipole')     filename = filename1
       if (tag .eq. 'virtual')    filename = filename2
       if (tag .eq. 'PK')         filename = filename3
       if (tag .eq. 'LO')         filename = filename4
       if (tag .eq. 'PK_Plus')    filename = filename6
       if (tag .eq. 'PK_Regular') filename = filename7
       if (tag .eq. 'PK_Delta')   filename = filename8

       if (iprint .eq. 1 ) then

       if (tag .eq. 'PK_Plus' .or. tag .eq. 'PK_Regular'
     .  .or. tag .eq. 'PK_Delta') then    
       call system("cd ../../trash/broken && cat "//trim(mode)// " >>
     . ../../PK_Isolated/summary/"//trim(run_tag)//"/"//trim(filename))
c     .  ../../summary/"//trim(run_tag)//"/"//trim(filename)) !//
c     .  " && rm -f "//trim(mode))
        filename_tmp1 = "../../PK_Isolated/summary/"
     .    //trim(run_tag)//"/"//trim(filename)

       call system(
     . 'sed -i -r "s/\\x1B\\[[0-9;]*[a-zA-Z]//g" '//trim(filename_tmp1))
       else
       call system("cd ../../trash/broken && cat "//trim(mode)// " >>
     .  ../../summary/"//trim(run_tag)//"/"//trim(filename)) !//
c     .  " && rm -f "//trim(mode))
        filename_tmp = "../"//trim(run_tag)//"/"//trim(filename)

       call system(
     . 'sed -i -r "s/\\x1B\\[[0-9;]*[a-zA-Z]//g" '//trim(filename_tmp))
       endif
       endif
       
       end

