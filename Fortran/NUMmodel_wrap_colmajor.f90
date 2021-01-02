module NUMmodel_wrap
  use iso_c_binding, only: c_double, c_int
  use NUMmodel, only:  parametersGeneralistsOnly, parametersGeneralistsCopepod, calcderivatives, rates

  implicit none

contains

  subroutine f_parametersGeneralistsOnly() bind(c)
    call parametersGeneralistsOnly()
  end subroutine f_parametersGeneralistsOnly

  subroutine f_parametersGeneralistsCopepod() bind(c)
    call parametersGeneralistsCopepod()
  end subroutine f_parametersGeneralistsCopepod  
  
  subroutine f_calcDerivatives(nGrid, u, L, dt, dudt) bind(c)
    integer(c_int), intent(in), value:: nGrid
    real(c_double), intent(in), value:: L, dt
    real(c_double), intent(in):: u(nGrid)
    real(c_double), intent(out):: dudt(nGrid)
    
    call calcDerivatives(u, L, dt)
    dudt = rates%dudt
  end subroutine f_calcDerivatives
  
end module NUMmodel_wrap


