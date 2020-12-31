module NUMmodel_wrap
  use iso_c_binding, only: c_double, c_int
  use NUMmodel, only:  parametersGeneralistsOnly, calcderivatives, rates

  implicit none

contains

  subroutine f_parametersGeneralistsOnly() bind(c)
    call parametersGeneralistsOnly()
  end subroutine f_parametersGeneralistsOnly
  
  subroutine f_calcDerivatives(nGrid, u, L, dt, dudt) bind(c)
    integer(c_int), intent(in), value:: nGrid
    real(c_double), intent(in), value:: L, dt
    real(c_double), intent(in):: u(nGrid)
    real(c_double), intent(out):: dudt(nGrid)
    
    call calcDerivatives(u, L, dt)
    dudt = rates%dudt
  end subroutine f_calcDerivatives
!!$    call setParameters(n, m, rhoCN, epsilonL, epsilonF, AN, AL, AF, Jmax, JFmax, Jresp, &
!!$         JlossPassive, theta, mort, mort2, mortHTL, remin, remin2, cLeakage)
!!$  end subroutine f_setParameters
!!$    
!!$  subroutine f_calcRates(T, L, n, u,  gammaN, gammaDOC, dudt) bind(c)
!!$    integer(c_int), intent(in), value:: n
!!$    real(c_double), intent(in), value:: T, L, gammaN, gammaDOC
!!$    real(c_double), intent(in):: u(n)
!!$    real(c_double), intent(out):: dudt(n)
!!$    
!!$    dudt = calcrates(T,L,u,gammaN,gammaDOC)
!!$  end subroutine f_calcRates
  
end module NUMmodel_wrap


