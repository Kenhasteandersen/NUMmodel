module globals
  implicit none
  integer, parameter :: dp=kind(0.d0) ! double precision
  
  ! Useful mathematical constants:
  real(dp), parameter :: onethird = 1.d0/3.d0
  real(dp), parameter :: twothirds = 2.d0/3.d0
  real(dp), parameter :: threequarters = 3.d0/4.d0
  real(dp), parameter :: pi = 4*ATAN(1.d0)

  ! Small number to avoid divisions by zero
  real(dp), parameter :: eps = 1d-200

  ! Temperature scalings parameters:
  real(dp) :: fTemp2, fTemp15 ! Temperature Q10 corrections (for Q10=2 and Q10=1.5)
  ! Temperature scalings parameters:
  real(dp) :: fTemp2, fTemp15 ! Temperature Q10 corrections (for Q10=2 and Q10=1.5)
  real(dp), parameter:: Tref = 10. ! Reference temperature
  
  character(len=16) :: inputfile='../input/input.h'
  
 

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
