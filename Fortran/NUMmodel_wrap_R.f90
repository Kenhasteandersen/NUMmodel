
subroutine f_setupGeneralistsOnly(n)
  use NUMmodel, only:  setupGeneralistsOnly
  use globals
  integer, intent(in):: n
  call setupGeneralistsOnly(n)
end subroutine f_setupGeneralistsOnly

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

  subroutine f_calcDerivatives(nn, u, L, T, dt, dudt)
    use globals
    use NUMmodel, only: calcDerivatives, rates
    integer, intent(in):: nn
    real(dp), intent(in):: L, T, dt
    real(dp), intent(in):: u(nn)
    real(dp), intent(out):: dudt(nn)

    call calcDerivatives(u, L, T, dt)
    dudt = rates%dudt
  end subroutine f_calcDerivatives

! Returns the rates calculated from last call to calcDerivatives
  subroutine f_getRates(jN, jDOC, jL, jSi, jF, jFreal,&
    jTot, jMax, jFmaxx, jR, jLossPassive, &
    jNloss,jLreal, &
    mortpred, mortHTL, mort2, mort)
    use globals
    use NUMmodel, only: getRates, nNutrients
    real(dp), intent(out):: jN(nGrid-nNutrients), jDOC(nGrid-nNutrients), jL(nGrid-nNutrients)
    real(dp), intent(out):: jSi(nGrid-nNutrients)
    real(dp), intent(out):: jF(nGrid-nNutrients), jFreal(nGrid-nNutrients)
    real(dp), intent(out):: jTot(nGrid-nNutrients), jMax(nGrid-nNutrients), jFmaxx(nGrid-nNutrients)
    real(dp), intent(out):: jR(nGrid-nNutrients)
    real(dp), intent(out):: jLossPassive(nGrid-nNutrients), jNloss(nGrid-nNutrients), jLreal(nGrid-nNutrients)
    real(dp), intent(out):: mortpred(nGrid-nNutrients), mortHTL(nGrid-nNutrients)
    real(dp), intent(out):: mort2(nGrid-nNutrients), mort(nGrid-nNutrients)

    call getRates(jN, jDOC, jL, jSi, jF, jFreal,&
    jTot, jMax, jFmaxx, jR, jLossPassive, &
    jNloss,jLreal, &
    mortpred, mortHTL, mort2, mort)
  end subroutine f_getRates

  subroutine f_getFunctions(ProdGross, ProdNet,ProdHTL,eHTL,Bpico,Bnano,Bmicro)
    use globals
    use NUMmodel, only: getFunctions
   real(dp), intent(out):: ProdGross, ProdNet,ProdHTL,eHTL,Bpico,Bnano,Bmicro

    call getFunctions(ProdGross, ProdNet,ProdHTL,eHTL,Bpico,Bnano,Bmicro)
  end subroutine f_getFunctions
  

   subroutine f_simulateChemostatEuler(u, L, T, nNutrients, Ndeep, diff, tEnd, dt)
    use globals
    use NUMmodel, only: simulateChemostatEuler

    integer, intent(in):: nNutrients
    real(dp), intent(inout):: u(nGrid)
    real(dp), intent(in):: L, T, Ndeep(nNutrients), diff, tEnd, dt

    call simulateChemostatEuler(u, L, T, Ndeep, diff, tEnd, dt)
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
