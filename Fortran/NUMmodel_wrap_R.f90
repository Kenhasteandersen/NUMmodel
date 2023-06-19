
subroutine f_setupGeneralistsSimpleOnly(n, errorio, errorstr)
  use NUMmodel, only:  setupGeneralistsSimpleOnly
  use globals
  use iso_c_binding, only: c_bool, c_char
  integer, intent(in):: n
  logical(c_bool), intent(out) :: errorio
  character(c_char), dimension(*) :: errorstr
  call setupGeneralistsSimpleOnly(n, errorio, errorstr)
end subroutine f_setupGeneralistsSimpleOnly

subroutine f_setupGeneralistsOnly(n, errorio, errorstr)
  use NUMmodel, only:  setupGeneralistsOnly
  use globals
  use iso_c_binding, only: c_bool, c_char
  integer, intent(in):: n
  logical(c_bool), intent(out) :: errorio
  character(c_char), dimension(*) :: errorstr
  call setupGeneralistsOnly(n, errorio, errorstr)
end subroutine f_setupGeneralistsOnly

subroutine f_setupGeneric(nAdult, mAdult, errorio, errorstr)
  use NUMmodel, only:  setupGeneric
  use globals
  use iso_c_binding, only: c_bool, c_char
  integer, intent(in):: nAdult
  real(dp), intent(in):: mAdult(nAdult)
  logical(c_bool), intent(out) :: errorio
  character(c_char), dimension(*) :: errorstr
  call setupGeneric(mAdult, errorio, errorstr)
end subroutine 

!!$  subroutine f_setupGeneralistsCopepod()
!!$    call setupGeneralistsCopepod()
!!$  end subroutine f_setupGeneralistsCopepod
!!$
!!$  subroutine f_setupGeneric(nCopepods, mAdult)
!!$    integer, intent(in):: nCopepods
!!$    real(dp), intent(in):: mAdult(nCopepods)
!!$
!!$    call setupGeneric(mAdult)
!!$  end subroutine f_setupGeneric

!!$  subroutine test(x)
!!$    integer, intent(in):: x
!!$   write(6,*) 'test'
!!$    write(6,*) x
!!$  end subroutine test

  subroutine f_setHTL(mHTL, mortHTL, bQuadraticHTL, bDecliningHTL)
    use globals
    use NUMmodel, only: setHTL
    real(dp), intent(in):: mHTL, mortHTL
    logical, intent(in):: bQuadraticHTL, bDecliningHTL

    call setHTL(mHTL, mortHTL, bQuadraticHTL, bDecliningHTL)
  end subroutine f_setHTL

  subroutine f_calcDerivatives(u, L, T, dt, dudt)
    use globals
    use NUMmodel, only: calcDerivatives, nGrid
    real(dp), intent(in):: L, T, dt
    real(dp), intent(in):: u(nGrid)
    real(dp), intent(out):: dudt(nGrid)

    call calcDerivatives(u, L, T, dt, dudt)
  end subroutine f_calcDerivatives

! Returns the rates calculated from last call to calcDerivatives
  subroutine f_getRates(jN, jDOC, jL, jSi, jF, jFreal, f, &
    jTot, jMax, jFmaxx, jR, jRespTot, jLossPassive, &
    jNloss,jLreal, jPOM, &
    mortpred, mortHTL, mort2, mort)
    use globals
    use NUMmodel, only: getRates, nNutrients, nGrid
    real(dp), intent(out):: jN(nGrid-nNutrients), jDOC(nGrid-nNutrients), jL(nGrid-nNutrients)
    real(dp), intent(out):: jSi(nGrid-nNutrients)
    real(dp), intent(out):: jF(nGrid-nNutrients), jFreal(nGrid-nNutrients), f(nGrid-nNutrients)
    real(dp), intent(out):: jTot(nGrid-nNutrients), jMax(nGrid-nNutrients), jFmaxx(nGrid-nNutrients)
    real(dp), intent(out):: jR(nGrid-nNutrients), jRespTot(nGrid-nNutrients)
    real(dp), intent(out):: jLossPassive(nGrid-nNutrients), jNloss(nGrid-nNutrients), jLreal(nGrid-nNutrients)
    real(dp), intent(out):: jPOM(nGrid-nNutrients)
    real(dp), intent(out):: mortpred(nGrid-nNutrients), mortHTL(nGrid-nNutrients)
    real(dp), intent(out):: mort2(nGrid-nNutrients), mort(nGrid-nNutrients)

    call getRates(jN, jDOC, jL, jSi, jF, jFreal, f, &
    jTot, jMax, jFmaxx, jR, jRespTot, jLossPassive, &
    jNloss,jLreal, jPOM, &
    mortpred, mortHTL, mort2, mort)
  end subroutine f_getRates

  !subroutine f_simulateChemostatEuler(nGri, u, L, Ndeep, diff, tEnd, dt)
  !  use globals
  !  use NUMmodel, only: simulateChemostatEuler
  !  integer, intent(in):: nGri
  !  real(dp), intent(inout):: u(nGri)
  !  real(dp), intent(in):: L, Ndeep, diff, tEnd, dt

    !call simulateChemostatEuler(u, L, Ndeep, diff, tEnd, dt)
!end subroutine f_simulateChemostatEuler

  subroutine f_getFunctions(u, ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro)
    use globals
    use NUMmodel, only: getFunctions, nGrid
    real(dp), intent(in) :: u(nGrid)
    real(dp), intent(out):: ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro

    call getFunctions(u, ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro)
  end subroutine f_getFunctions
  
!  subroutine f_getBalance(Nbalance, Cbalance)
!    use globals
!    use NUMmodel, only: getBalance
!   real(dp), intent(out):: Nbalance, Cbalance

!    call getBalance(Nbalance, Cbalance)
!  end subroutine f_getBalance

   subroutine f_simulateChemostatEuler(u, L, T, nNutrients, Ndeep, diff, tEnd, dt, bLosses)
    use globals
    use NUMmodel, only: simulateChemostatEuler, nGrid

    integer, intent(in):: nNutrients
    real(dp), intent(inout):: u(nGrid)
    real(dp), intent(in):: L, T, Ndeep(nNutrients), diff, tEnd, dt
    logical(1), intent(in):: bLosses

    call simulateChemostatEuler(u, L, T, Ndeep, diff, tEnd, dt, bLosses)
  end subroutine f_simulateChemostatEuler
!!$
!!$  subroutine f_simulateEuler(nGrid, u, L, tEnd, dt)
!!$    integer, intent(in):: nGrid
!!$    real(dp), intent(inout):: u(nGrid)
!!$    real(dp), intent(in)!:: L, tEnd, dt
!!$
!!$    call simulateEuler(u, L, tEnd, dt)
!!$  end subroutine f_simulateEuler
!!$
