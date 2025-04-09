module constants_mod
  use iso_fortran_env, only: real64
  implicit none

  ! Make everything private by default
  private

  ! Public constants
  public :: PI, CA, CF, NF, TR, hbarc2, EulerGamma, mass, decayW
  public :: init_mass

  ! Define the constants
  real(real64), parameter :: PI = 3.141592653589793238462643d0
  real(real64), parameter :: CA = 3.0d0
  real(real64), parameter :: CF = 4.0d0 / 3.0d0
  real(real64), parameter :: NF = 5.0d0
  real(real64), parameter :: TR = 0.5d0
  real(real64), parameter :: hbarc2 = 0.3894d12 ! in fb
  real(real64), parameter :: EulerGamma = 0.5772156649015328606065120d0
  real(real64), dimension(-25:25) :: mass
  real(real64), dimension(-25:25) :: decayW 
 
    contains
   
     subroutine init_mass 
     implicit none
        mass   = 0.0_real64
        decayW = 0.0_real64
        ! Fermions - Quarks
        mass(1)   = 0.0_real64       ! down (d)
        mass(1)   = 0.0_real64       ! up (u)
        mass(3)   = 0.0_real64       ! strange (s)
        mass(4)   = 0.0_real64       ! charm (c)
        mass(5)   = 0.0_real64       ! bottom (b)
        mass(6)   = 173.0_real64     ! top (t)
        
        ! Fermions - Leptons
        mass(11)  = 0.0_real64       ! electron (e-)
        mass(13)  = 0.0_real64       ! muon (μ-)
        mass(15)  = 1.77686d0        ! tau (τ-)
        
        ! Neutrinos (set to zero, or tiny if needed)
        mass(12)  = 0d0              ! electron neutrino (νe)
        mass(14)  = 0d0              ! muon neutrino (νμ)
        mass(16)  = 0d0              ! tau neutrino (ντ)
        
        ! Anti-fermions: same mass as their particles
        mass(-1)  = mass(1)          ! anti-down
        mass(-2)  = mass(2)          ! anti-up
        mass(-3)  = mass(3)          ! anti-strange
        mass(-4)  = mass(4)          ! anti-charm
        mass(-5)  = mass(5)          ! anti-bottom
        mass(-6)  = mass(6)          ! anti-top
        
        mass(-11) = mass(11)         ! positron
        mass(-13) = mass(13)         ! anti-muon
        mass(-15) = mass(15)         ! anti-tau
        mass(-12) = mass(12)
        mass(-14) = mass(14)
        mass(-16) = mass(16)
        
        ! Gauge Bosons
        mass(21)  = 0d0              ! gluon (g)
        mass(22)  = 0d0              ! photon (γ)
        mass(23)  = 9.118760e+01     ! Z boson
        mass(24)  = 80.3850115091    ! W+ boson
        mass(-24) = mass(24)         ! W- boson
        
        ! Higgs
        mass(25)  = 125.10d0    ! Higgs boson (h)
        
        ! Optional: ghost and Goldstone bosons, if needed
            ! mass(9000001) = ...d0  ! ghost or BSM
      
        ! DECAY WIDTHS
         decayW(23) = 2.495200e+00     ! Z-boson
         decayW(24) = 2.085400e+00     ! W-boson
      end subroutine init_mass
end module constants_mod

