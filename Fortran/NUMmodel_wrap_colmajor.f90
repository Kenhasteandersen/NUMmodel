module NUMmodel_wrap
  use iso_c_binding, only: c_double, c_int
  use NUMmodel, only:  setupGeneralistsOnly, setupGeneralistsCopepod, &
       setupGeneric, calcderivatives, rates, simulateChemostatEuler

  implicit none

contains

  subroutine f_setupGeneralistsOnly() bind(c)
    call setupGeneralistsOnly()
  end subroutine f_setupGeneralistsOnly

  subroutine f_setupGeneralistsCopepod() bind(c)
    call setupGeneralistsCopepod()
  end subroutine f_setupGeneralistsCopepod

  subroutine f_setupGeneric(nCopepods, mAdult) bind(c)
    integer(c_int), intent(in), value:: nCopepods
    real(c_double), intent(in):: mAdult(nCopepods)

    call setupGeneric(mAdult)
  end subroutine f_setupGeneric
  
  subroutine f_calcDerivatives(nGrid, u, L, dt, dudt) bind(c)
    integer(c_int), intent(in), value:: nGrid
    real(c_double), intent(in), value:: L, dt
    real(c_double), intent(in):: u(nGrid)
    real(c_double), intent(out):: dudt(nGrid)
    
    call calcDerivatives(u, L, dt)
    dudt = rates%dudt
  end subroutine f_calcDerivatives

  subroutine f_simulateChemostatEuler(nGrid, u, L, Ndeep, diff, tEnd, dt) bind(c)
    integer(c_int), intent(in), value:: nGrid
    real(c_double), intent(inout):: u(nGrid)
    real(c_double), intent(in), value:: L, Ndeep, diff, tEnd, dt
    
    call simulateChemostatEuler(u, L, Ndeep, diff, tEnd, dt)
  end subroutine f_simulateChemostatEuler
end module NUMmodel_wrap


