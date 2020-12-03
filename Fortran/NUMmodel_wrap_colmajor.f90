module NUMmodel_wrap
  ! To wrap calc_num_iter for use in Julia with column-major arrays.  
  ! Still using iso_c_binding.

  use iso_c_binding, only: c_double, c_int
  use NUMmodel, only:  calcrates, setparameters

  implicit none

contains

  ! need to make a subroutine as only scalars can be returned

  subroutine f_setParameters(n, m, rhoCN, epsilonL, epsilonF, &
       AN, AL, AF, &
       Jmax, JFmax, Jresp, JlossPassive, &
       theta, &
       mort, mort2, mortHTL, remin, remin2, cLeakage) bind(c)
    integer(c_int), intent(in), value:: n
    real(c_double), intent(in):: m(n)
    real(c_double), intent(in), value:: rhoCN, epsilonL, epsilonF
    real(c_double), intent(in):: AN(n), AL(n), AF(n)
    real(c_double), intent(in):: Jmax(n), JFmax(n), Jresp(n), JlossPassive(n)
    real(c_double), intent(in):: theta(n, n)
    real(c_double), intent(in):: mort(n)
    real(c_double), intent(in), value:: mort2
    real(c_double), intent(in):: mortHTL(n)
    real(c_double), intent(in), value:: remin, remin2, cLeakage
    
    call setParameters(n, m, rhoCN, epsilonL, epsilonF, AN, AL, AF, Jmax, JFmax, Jresp, &
         JlossPassive, theta, mort, mort2, mortHTL, remin, remin2, cLeakage)
  end subroutine f_setParameters
    
  subroutine f_calcRates(T, L, n, u,  gammaN, gammaDOC, dudt) bind(c)
    integer(c_int), intent(in), value:: n
    real(c_double), intent(in), value:: T, L, gammaN, gammaDOC
    real(c_double), intent(in):: u(n)
    real(c_double), intent(out):: dudt(n)
    
    dudt = calcrates(T,L,u,gammaN,gammaDOC)
  end subroutine f_calcRates
  
end module NUMmodel_wrap


