      module kinparams
 
      use inprms_ps
 
      integer, parameter :: iqprec=0
      integer, parameter :: xp=kind(1.d0)
      integer, parameter :: nexp=4
      integer, parameter :: nfsp=2
      integer, parameter :: nkins=1
      integer, parameter :: nkinst=1
      integer, parameter :: nkinag=1
      integer, parameter :: nkin=1
      integer, parameter :: igind=4
      real(xp) :: m(nexp)
      character(10) :: expc(nexp)
      character(1) :: expl(nexp)
 
      integer :: iwkin(nkin)
      real(kind(1.d0)) :: aw(nkin), xmc(0:nkin)
      real(kind(1.d0)) :: aw_w(nkin,3)
 
      contains
       subroutine exp_masses
       m=(/mg,mg,mh,mg/)
       end subroutine exp_masses
       subroutine exp_names
        expc(1)='g'
        expl(1)='g'
        expc(2)='g'
        expl(2)='g'
        expc(3)='h'
        expl(3)='h'
        expc(4)='g'
        expl(4)='g'
       end subroutine exp_names
 
      end module kinparams
