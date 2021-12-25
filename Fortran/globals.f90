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

  ! Types of spectra:
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

  ! Small number to avoid divisions by zero
  real(dp), parameter :: eps = 1d-200

  ! type typeRates
  !    real(dp), dimension(:), allocatable:: flvl, JF, F, JEnc
  !    real(dp), dimension(:), allocatable:: JN, JDOC, JSi, JL
  !    real(dp), dimension(:), allocatable:: JNtot, JLreal, JCtot, Jtot
  !    real(dp), dimension(:), allocatable:: JCloss_feeding, JCloss_photouptake, JNlossLiebig, JClossLiebig
  !    real(dp), dimension(:), allocatable:: JNloss, JCloss, JSiloss
  !    real(dp), dimension(:), allocatable:: mortpred, mortHTL
  !    real(dp), dimension(:), allocatable:: g, mortStarve, mort ! Multicellular rates
  !    real(dp), dimension(:), allocatable:: dudt
  ! end type typeRates

  ! Temperature Q10 corrections (for Q10=2 and Q10=1.5)
  real(dp) :: fTemp2, fTemp15
  real(dp), parameter:: Tref = 10. ! Reference temperature

  contains

  ! -----------------------------------------------
  ! Temperature Q10 function
  ! -----------------------------------------------
  function fTemp(Q10, T) result(f)
    real(dp), intent(in):: Q10, T
    real(dp):: f

    f = Q10**((T-Tref)/10.)
  end function fTemp

  ! -----------------------------------------------
  ! Update the temperature corrections only if T has changed
  ! -----------------------------------------------
  subroutine updateTemperature(T)
    real(dp), intent(in) :: T
    real(dp), save :: Told = -1000.

    if (T .ne. Told) then
      Told = T
      fTemp2 = fTemp(2.d0, T)
      fTemp15 = fTemp(1.5d0, T)
    end if
  end subroutine updateTemperature

end module globals
