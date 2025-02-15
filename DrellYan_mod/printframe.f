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
          write(*,'(A,ES25.0)')"Using Vegas points:",pt1
          write(*,*)"        Iteration:         ",its1
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
c ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine printframe8(p1,p2,p3,p4)
      implicit double precision(a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
        print*,"p1:",p1
        print*,"p2:",p2
        print*,"p3:",p3
        print*,"p4:",p4
        print*

      end

