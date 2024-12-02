c       subroutine printframe1(pt1,its1)
c       subroutine printframe2(xq)
c       subroutine printframe3(name,xintegral,error,chi2)
c       subroutine printframe4(name)
c       subroutine printframe6(ecm,xmur,xmuf,pdf_name,amH)
c       subroutine printframe5(e_cut,t_cut)

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
cc ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        
       subroutine printframe6(ecm,xmur,xmuf,pdf_name,amH)
       implicit double precision(a-h,o-z)
c       character*50 pdf_name
       character(len=*) pdf_name
          print*,"                                        "
          print*,"~~~~~~~~~~~~~[Parameters  Used]~~~~~~~~~~~"
          write(*,*)"Centre of Mass Energy: ",int(ecm)
          write(*,*)"                  mur: ",xmur
          write(*,*)"                  muf: ",xmuf
          write(*,*)"          LHAPDF used: ",pdf_name
          write(*,*)"           Higgs mass: ",amH
          print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       return
       end
c ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        
