module globals
  implicit none
  integer, parameter :: dp=kind(0.d0) ! double precision
  !
  ! Useful constants:
  !
  integer, parameter :: idxN = 1
  integer, parameter :: idxDOC = 2
  integer, parameter :: idxB = 3
  integer, parameter :: typeGeneralist = 1
  integer, parameter :: typeCopepod = 10
  real(dp), parameter :: onethird = 1.d0/3.d0
  real(dp), parameter :: twothirds = 2.d0/3.d0
  real(dp), parameter :: threequarters = 3.d0/4.d0

  type typeRates
     real(dp), dimension(:), allocatable:: flvl, JF, F, JEnc
     real(dp), dimension(:), allocatable:: JN, JDOC, JL
     real(dp), dimension(:), allocatable:: JNtot, JLreal, JCtot, Jtot
     real(dp), dimension(:), allocatable:: JCloss_feeding, JCloss_photouptake, JNlossLiebig, JClossLiebig
     real(dp), dimension(:), allocatable:: JNloss, JCloss
     real(dp), dimension(:), allocatable:: mortpred, mortHTL
     real(dp), dimension(:), allocatable:: g, mortStarve, mort ! Multicellular rates
     real(dp), dimension(:), allocatable:: dudt
  end type typeRates

  integer:: nGrid ! Total number of grid points incl. two points for N and DOC

!contains
  
!!$  subroutine initGlobals(nnGrid)
!!$    integer, intent(in):: nnGrid
!!$
!!$    nGrid = nnGrid+2
!!$
!!$    allocate(m(nGrid))
!!$    allocate(beta(nGrid))
!!$    allocate(sigma(nGrid))
!!$    allocate(AF(nGrid))
!!$    allocate(JFmax(nGrid))
!!$    allocate(epsilonF(nGrid))
!!$  end subroutine initGlobals

end module globals
