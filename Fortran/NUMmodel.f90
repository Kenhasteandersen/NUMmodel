! Calculates the Mandelbrot set
module NUMmodel
  implicit none
  integer, parameter :: dp=kind(0.d0) ! double precision

    type parameters
     integer:: n
     real(dp), dimension(:), allocatable:: m(:)
  end type parameters

contains

  function setparameters(n) result(p)
    integer, intent(in):: n
    type(parameters):: p

    p%n = n
    allocate(p%m(n))
  end function setparameters
    
  function calcrates(u) result(dudt)
    real(dp), intent(in):: u(:)
    real(dp):: dudt(size(u))

    dudt = 0.01*u

  end function calcrates
  
end module NUMmodel


