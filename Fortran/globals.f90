module globals
  implicit none
  integer, parameter :: dp=kind(0.d0) ! double precision
  !
  ! Useful constants:
  !

  ! Indices into the state-variable vector:
  integer, parameter :: idxN = 1
  integer, parameter :: idxDOC = 2
  integer, parameter :: idxSi = 3

  ! Type of spectra:
  integer, parameter :: typeGeneralist = 1
  integer, parameter :: typeGeneralist_csp = 2
  integer, parameter :: typeDiatom = 3
  integer, parameter :: typeDiatom_simple = 4  
  integer, parameter :: typeCopepod = 10

  ! Useful mathematical constants:
  real(dp), parameter :: onethird = 1.d0/3.d0
  real(dp), parameter :: twothirds = 2.d0/3.d0
  real(dp), parameter :: threequarters = 3.d0/4.d0
  real(dp), parameter :: pi = 4*ATAN(1.d0)

  type typeRates
     real(dp), dimension(:), allocatable:: flvl, JF, F, JEnc
     real(dp), dimension(:), allocatable:: JN, JDOC, JSi, JL
     real(dp), dimension(:), allocatable:: JNtot, JLreal, JCtot, Jtot
     real(dp), dimension(:), allocatable:: JCloss_feeding, JCloss_photouptake, JNlossLiebig, JClossLiebig
     real(dp), dimension(:), allocatable:: JNloss, JCloss, JSiloss
     real(dp), dimension(:), allocatable:: mortpred, mortHTL
     real(dp), dimension(:), allocatable:: g, mortStarve, mort ! Multicellular rates
     real(dp), dimension(:), allocatable:: dudt
  end type typeRates

  integer:: nGrid ! Total number of grid points incl.  points for nutrients
end module globals
