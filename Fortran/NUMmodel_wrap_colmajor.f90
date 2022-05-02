module NUMmodel_wrap
  use iso_c_binding, only: c_double, c_int, c_bool
  use NUMmodel, only:  nGrid, idxB, setupGeneralistsOnly,  &
      setupGeneralistsPOM, setHTL, setmortHTL, getRates, &         
       setupGeneralistsOnly_csp, &
       setupDiatomsOnly, &
      setupDiatoms_simpleOnly, setupGeneralistsDiatoms_simple, &
       setupGeneralistsDiatoms, &
       setupGeneralistsCopepod, &
       setupGeneric, setupNUMmodel, setupGeneric_csp, &
       calcderivatives, &
       simulateChemostatEuler, simulateEuler, getFunctions, &
       getBalance

  use globals

  implicit none

contains

  subroutine f_setupGeneralistsOnly(n) bind(c)
    integer(c_int), intent(in), value:: n
    call setupGeneralistsOnly(n)
  end subroutine f_setupGeneralistsOnly

  subroutine f_setupGeneralistsPOM(n, nPOM) bind(c)
    integer(c_int), intent(in), value:: n, nPOM
    call setupGeneralistsPOM(n, nPOM)
  end subroutine f_setupGeneralistsPOM

  subroutine f_setupGeneralistsOnly_csp() bind(c)
    call setupGeneralistsOnly_csp()
  end subroutine f_setupGeneralistsOnly_csp

  subroutine f_setupDiatomsOnly(n) bind(c)
    integer(c_int), intent(in), value:: n
    call setupDiatomsOnly(n)
  end subroutine f_setupDiatomsOnly

  subroutine f_setupDiatoms_simpleOnly(n) bind(c)
    integer(c_int), intent(in), value:: n
    call setupDiatoms_simpleOnly(n)
  end subroutine f_setupDiatoms_simpleOnly

  subroutine f_setupGeneralistsDiatoms(n) bind(c)
    integer(c_int), intent(in), value:: n
    call setupGeneralistsDiatoms(n)
  end subroutine f_setupGeneralistsDiatoms

  subroutine f_setupGeneralistsDiatoms_simple(n) bind(c)
    integer(c_int), intent(in), value:: n
    call setupGeneralistsDiatoms_simple(n)
  end subroutine f_setupGeneralistsDiatoms_simple

  subroutine f_setupGeneralistsCopepod() bind(c)
    call setupGeneralistsCopepod()
  end subroutine f_setupGeneralistsCopepod

  subroutine f_setupGeneric(nCopepods, mAdult) bind(c)
    integer(c_int), intent(in), value:: nCopepods
    real(c_double), intent(in):: mAdult(nCopepods)

    call setupGeneric(mAdult)
  end subroutine f_setupGeneric

  subroutine f_setupNUMmodel(n,nCopepod,nPOM, nCopepods, mAdult) bind(c)
    integer(c_int), intent(in), value:: n,nCopepod,nPOM, nCopepods
    real(c_double), intent(in):: mAdult(nCopepods)

    call setupNUMmodel(n,nCopepod,nPOM,mAdult)
  end subroutine f_setupNUMmodel

  subroutine f_setupGeneric_csp(nCopepods, mAdult) bind(c)
    integer(c_int), intent(in), value:: nCopepods
    real(c_double), intent(in):: mAdult(nCopepods)

    call setupGeneric_csp(mAdult)
  end subroutine f_setupGeneric_csp

  subroutine f_setHTL(mHTL, mortHTL, bQuadraticHTL, bDecliningHTL) bind(c)
    real(c_double), intent(in), value:: mHTL, mortHTL
    logical, intent(in), value:: bQuadraticHTL, bDecliningHTL

    call setHTL(mHTL, mortHTL, bQuadraticHTL, bDecliningHTL)
  end subroutine f_setHTL 

  subroutine f_setMortHTL(mortHTL) bind(c)
    real(c_double), intent(in):: mortHTL(nGrid-idxB+1)

    call setMortHTL(mortHTL)
  end subroutine f_setMortHTL

  subroutine test(x) bind(c)
    integer(c_int), intent(in), value:: x
   write(6,*) 'test'
    write(6,*) x
  end subroutine test

  subroutine f_calcDerivatives(u, L, T, dt, dudt) bind(c)
    real(c_double), intent(in), value:: L, T, dt
    real(c_double), intent(in):: u(nGrid)
    real(c_double), intent(out):: dudt(nGrid)

    call calcDerivatives(u, L, T, dt, dudt)
    
  end subroutine f_calcDerivatives

  ! subroutine f_calcRates(nGrid, u, L, T, jN, jL, jF, jTot, mortHTL, mortpred, g) bind(c)
  !   use NUMmodel, only: idxB
  !   integer(c_int), intent(in), value:: nGrid
  !   real(c_double), intent(in):: u(nGrid)
  !   real(c_double), intent(in), value:: L, T
  !   real(c_double), intent(out):: jN(nGrid), jL(nGrid), jF(nGrid)
  !   real(c_double), intent(out):: jTot(nGrid), mortHTL(nGrid), mortpred(nGrid), g(nGrid)

  !   call calcDerivatives(u, L, T, 0.d0)
  !   ! jN(idxB:nGrid) = rates%JN(idxB:nGrid) / m(idxB:nGrid)
  !   ! jL(idxB:nGrid) = rates%JL(idxB:nGrid) / m(idxB:nGrid)
  !   ! jF(idxB:nGrid) = rates%JF(idxB:nGrid) / m(idxB:nGrid)
  !   ! jtot(idxB:nGrid) = rates%Jtot(idxB:nGrid) / m(idxB:nGrid)
  !   ! mortHTL(idxB:nGrid) = rates%mortHTL(idxB:nGrid)
  !   ! mortpred(idxB:nGrid) = rates%mortpred(idxB:nGrid)
  !   ! g(idxB:nGrid) = rates%g(idxB:nGrid)
  ! end subroutine f_calcRates

  subroutine f_simulateChemostatEuler(u, L, T, nNutrients, Ndeep, diff, tEnd, dt, bLosses) bind(c)
    !integer(c_int), intent(in), value:: nGrid
    real(c_double), intent(inout):: u(nGrid)
    integer(c_int), intent(in), value:: nNutrients
    real(c_double), intent(in):: Ndeep(nNutrients)
    real(c_double), intent(in), value:: L, T, diff, tEnd, dt
    logical(c_bool), intent(in), value:: bLosses

    call simulateChemostatEuler(u, L, T, Ndeep, diff, tEnd, dt, bLosses)
  end subroutine f_simulateChemostatEuler

  subroutine f_simulateEuler(u, L, T, tEnd, dt) bind(c)
    real(c_double), intent(inout):: u(nGrid)
    real(c_double), intent(in), value:: L, T, tEnd, dt

    call simulateEuler(u, L, T, tEnd, dt)
  end subroutine f_simulateEuler

  subroutine f_getFunctions(u, ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro) bind(c)
    real(c_double), intent(in):: u(nGrid)
    real(c_double), intent(out):: ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro

    call getFunctions(u, ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro)
  end subroutine f_getFunctions

  subroutine f_getBalance(u, dudt, Nbalance, Cbalance) bind(c)
     real(c_double), intent(in):: u(nGrid), dudt(nGrid)
     real(c_double), intent(out):: Nbalance, Cbalance

     call getBalance(u, dudt, Nbalance, Cbalance)
   end subroutine f_getBalance  
  
  subroutine f_getMass(m, mDelta) bind(c)
    use globals
    use NUMmodel, only: getMass

    real(c_double), intent(inout):: m(nGrid), mDelta(nGrid)
    
    call getMass(m, mDelta) 
  end subroutine f_getMass

  subroutine f_getSinking(velocity) bind(c)
    use globals
    use NUMmodel, only: getSinking, nGrid

    real(c_double), intent(inout):: velocity(nGrid)

    call getSinking(velocity)
  end subroutine f_getSinking

  subroutine f_getRates(jN, jDOC, jL, jSi, jF, jFreal,&
    jTot, jMax, jFmax, jR, jLossPassive, &
    jNloss,jLreal, &
    mortpred, mortHTL, mort2, mort) bind(c)
    use globals
    use NUMmodel, only: nNutrients, getRates
    real(dp), intent(out):: jN(nGrid-nNutrients), jDOC(nGrid-nNutrients), jL(nGrid-nNutrients)
    real(dp), intent(out):: jSi(nGrid-nNutrients)
    real(dp), intent(out):: jF(nGrid-nNutrients), jFreal(nGrid-nNutrients)
    real(dp), intent(out):: jTot(nGrid-nNutrients), jMax(nGrid-nNutrients), jFmax(nGrid-nNutrients)
    real(dp), intent(out):: jR(nGrid-nNutrients)
    real(dp), intent(out):: jLossPassive(nGrid-nNutrients), jNloss(nGrid-nNutrients), jLreal(nGrid-nNutrients)
    real(dp), intent(out):: mortpred(nGrid-nNutrients), mortHTL(nGrid-nNutrients)
    real(dp), intent(out):: mort2(nGrid-nNutrients), mort(nGrid-nNutrients)

   call getRates(jN, jDOC, jL, jSi, jF, jFreal,&
   jTot, jMax, jFmax, jR, jLossPassive, &
   jNloss,jLreal, &
   mortpred, mortHTL, mort2, mort)
  end subroutine f_getRates
  
end module NUMmodel_wrap
