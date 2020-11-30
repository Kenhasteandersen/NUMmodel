module NUMmodel_wrap
  ! To wrap calc_num_iter for use in Julia with column-major arrays.  
  ! Still using iso_c_binding.

  use iso_c_binding, only: c_double, c_int
  use NUMmodel, only:  calcrates, setparameters, parameters

  implicit none

contains

  ! need to make a subroutine as only scalars can be returned

  subroutine f_setparameters(n, p) bind(c)
    integer(c_int), intent(in), value:: n
    type(parameters), intent(out):: p
  end subroutine f_setparameters
    
  subroutine f_calcrates(n, u, dudt) bind(c)
    integer(c_int), intent(in), value:: n
    real(c_double), intent(in):: u(n)
    real(c_double), intent(out):: dudt(n)
    
    dudt = calcrates(u)
  end subroutine f_calcrates
  
end module NUMmodel_wrap


