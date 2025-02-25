c       subroutine printframe1(pt1,its1)
c       subroutine printframe2(xq)
c       subroutine printframe3(name,xintegral,error,chi2)
c       subroutine printframe4(name)

c ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        
       subroutine printframe0(name)
       implicit double precision(a-h,o-z)
       character*50 name

        print*,"  "
        print*,"  "
        print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print*,"    Calculating "//trim(name)
        print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print*,"  "
       end
c ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        
       subroutine printframe1(pt1,its1)
       implicit double precision(a-h,o-z)
          print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
          write(*,'(a,es25.0)')"using vegas points:",pt1
          write(*,*)"        iteration:         ",its1
          print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       return
       end
c ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        
        subroutine printframe2(xq)
        implicit double precision(a-h,o-z)
          print*," "
      write(*,*) achar(27)//'[1;32m' // "For xq = ",int(xq) ,achar(27) 
     .   //'[0m'
          print*," "

        end
c ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        
        subroutine printframe3(name,xintegral,error,chi2)
        implicit double precision(a-h,o-z)
        character*50 name

            print*,"  "
            print*,"  "
            write(*,*)achar(27)//'[1;32m'//"Integral "//trim(name)// 
     .      ":",xintegral,achar(27) //'[0m', "+-",error
            write(*,*)"with chisq    =",chi2
            print*," "
            print*," "
        end
c ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        
       subroutine printframe4(name)
       character*50 name
       
            print*," "
       write(*,*)achar(27)//'[1;32m'//"   xq         Integral "
     .  //trim(name),
     .  "                 error",
     . achar(27)//'[0m'
       end
cc ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        
       subroutine printframe5(e_cut,t_cut)
       implicit double precision(a-h,o-z)
          print*,"                                        "
          print*,"~~~~~~~~~~~~~[Technical Cuts Used]~~~~~~~~~~~"
          write(*,'(a,es25.0)')"  Cut for Energy e5 :",e_cut
          write(*,'(a,es25.0)')"  Collinear cut sij :",t_cut
          print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       return
       end
c ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        
