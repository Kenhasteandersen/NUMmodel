module NUMmodel_wrap
  use iso_c_binding, only: c_double, c_int, c_bool, c_char, c_null_char
  use NUMmodel, only:  nGrid, idxB, nGrid, &
       setupGeneralistsSimpleOnly, setupGeneralistsSimplePOM, &
       setupGeneralistsOnly,setupGeneralistsPOM,  &
       setupDiatomsOnly, &
       setupDiatoms_simpleOnly, setupGeneralistsDiatoms_simple, &
       setupGeneralistsDiatoms, &
       setupGeneralistsSimpleCopepod, &
       setupGeneric, setupNUMmodel, setupNUMmodelSimple, setupGenDiatCope, &
       calcderivatives, &
       simulateChemostatEuler, simulateEuler, simulateEulerFunctions, getFunctions, &
       setHTL, setmortHTL, setSinking, getRates, getBalance, getLost, theta

  use globals

  implicit none

contains

  subroutine f_setupGeneralistsOnly(n, errorio, errorstr) bind(c)
    integer(c_int), intent(in), value:: n
    logical(c_bool), intent(out) :: errorio
    character(c_char), dimension(*) :: errorstr
    call setupGeneralistsOnly(n,errorio, errorstr)
  end subroutine f_setupGeneralistsOnly

  subroutine f_setupGeneralistsSimplePOM(n, nPOM, errorio, errorstr) bind(c)
    integer(c_int), intent(in), value:: n, nPOM
    logical(c_bool), intent(out) :: errorio
    character(c_char), dimension(*) :: errorstr
    call setupGeneralistsSimplePOM(n, nPOM,errorio, errorstr)
  end subroutine f_setupGeneralistsSimplePOM
  
  
  subroutine f_setupGeneralistsPOM(n, nPOM, errorio, errorstr) bind(c)
    integer(c_int), intent(in), value:: n, nPOM
    logical(c_bool), intent(out) :: errorio
    character(c_char), dimension(*) :: errorstr
    call setupGeneralistsPOM(n, nPOM,errorio, errorstr)
  end subroutine f_setupGeneralistsPOM

  subroutine f_setupGeneralistsSimpleOnly(n, errorio, errorstr) bind(c)
    integer(c_int), intent(in), value:: n
    logical(c_bool), intent(out) :: errorio
    character(c_char), dimension(*) :: errorstr
    call setupGeneralistsSimpleOnly(n,errorio, errorstr)
  end subroutine f_setupGeneralistsSimpleOnly
 
  subroutine f_setupDiatomsOnly(n, errorio, errorstr) bind(c)
    integer(c_int), intent(in), value:: n
    logical(c_bool), intent(out) :: errorio
    character(c_char), dimension(*) :: errorstr
    call setupDiatomsOnly(n,errorio, errorstr)
  end subroutine f_setupDiatomsOnly

  subroutine f_setupDiatoms_simpleOnly(n, errorio, errorstr) bind(c)
    integer(c_int), intent(in), value:: n
    logical(c_bool), intent(out) :: errorio
    character(c_char), dimension(*) :: errorstr
    call setupDiatoms_simpleOnly(n,errorio, errorstr)
  end subroutine f_setupDiatoms_simpleOnly

  subroutine f_setupGeneralistsDiatoms(n, errorio, errorstr) bind(c)
    integer(c_int), intent(in), value:: n
    logical(c_bool), intent(out) :: errorio
    character(c_char), dimension(*) :: errorstr
    call setupGeneralistsDiatoms(n,errorio, errorstr)
  end subroutine f_setupGeneralistsDiatoms

  subroutine f_setupGeneralistsDiatoms_simple(n, errorio, errorstr) bind(c)
    integer(c_int), intent(in), value:: n
    logical(c_bool), intent(out) :: errorio
    character(c_char), dimension(*) :: errorstr
    call setupGeneralistsDiatoms_simple(n,errorio, errorstr)
  end subroutine f_setupGeneralistsDiatoms_simple

  subroutine f_setupGeneralistsSimpleCopepod(errorio, errorstr) bind(c)
    logical(c_bool), intent(out) :: errorio
    character(c_char), dimension(*) :: errorstr
    call setupGeneralistsSimpleCopepod(errorio, errorstr)
  end subroutine f_setupGeneralistsSimpleCopepod

  subroutine f_setupGeneric(nCopepods, mAdult, errorio, errorstr) bind(c)
    integer(c_int), intent(in), value:: nCopepods
    real(c_double), intent(in):: mAdult(nCopepods)
    logical(c_bool), intent(out) :: errorio
    character(c_char), dimension(*) :: errorstr
    call setupGeneric(mAdult,errorio, errorstr)
  end subroutine f_setupGeneric

  subroutine f_setupNUMmodel(n,nCopepod,nPOM, nCopepodsPassive, mAdultPassive, &
     nCopepodsActive, mAdultActive, errorio, errorstr) bind(c)
    integer(c_int), intent(in), value:: n,nCopepod,nPOM, nCopepodsPassive, nCopepodsActive
    real(c_double), intent(in):: mAdultPassive(nCopepodsPassive), mAdultActive(nCopepodsActive)
    logical(c_bool), intent(out) :: errorio
    character(c_char), dimension(*) :: errorstr
    call setupNUMmodel(n,nCopepod,nPOM, mAdultPassive, mAdultActive,errorio, errorstr)
  end subroutine f_setupNUMmodel

    subroutine f_setupNUMmodelSimple(n,nCopepod,nPOM, nCopepods, mAdult, errorio, errorstr) bind(c)
    integer(c_int), intent(in), value:: n,nCopepod,nPOM, nCopepods
    real(c_double), intent(in):: mAdult(nCopepods)
    logical(c_bool), intent(out) :: errorio
    character(c_char), dimension(*) :: errorstr
    call setupNUMmodelSimple(n,nCopepod,nPOM,mAdult,errorio, errorstr)
  end subroutine f_setupNUMmodelSimple

  subroutine f_setupGenDiatCope(n,nCopepod,nPOM, nCopepods, mAdult, errorio, errorstr) bind(c)
    integer(c_int), intent(in), value:: n,nCopepod,nPOM, nCopepods
    real(c_double), intent(in):: mAdult(nCopepods)
    logical(c_bool), intent(out) :: errorio
    character(c_char), dimension(*) :: errorstr
    call setupGenDiatCope(n,nCopepod,nPOM,mAdult,errorio, errorstr)
  end subroutine f_setupGenDiatCope

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

  subroutine f_simulateEulerFunctions(u, L, T, tEnd, dt, &
    ProdGross, ProdNet,ProdHTL,prodBact,eHTL,Bpico,Bnano,Bmicro) bind(c)
    real(c_double), intent(inout):: u(nGrid) ! Initial conditions and result after integration
    real(c_double), intent(in), value:: L,T,tEnd,dt
    real(c_double), intent(out):: ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro

    call simulateEulerFunctions(u, L, T, tEnd, dt, &
      ProdGross, ProdNet,ProdHTL,prodBact,eHTL,Bpico,Bnano,Bmicro)
  end subroutine f_simulateEulerFunctions

   subroutine f_getFunctions(u, ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro) bind(c)
    real(c_double), intent(in):: u(nGrid)
    real(c_double), intent(out):: ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro

    call getFunctions(u, ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro)
  end subroutine f_getFunctions

  subroutine f_getBalance(u, dudt, Cbalance, Nbalance, Sibalance) bind(c)
     real(c_double), intent(in):: u(nGrid), dudt(nGrid)
     real(c_double), intent(out):: Cbalance, Nbalance, Sibalance

     call getBalance(u, dudt, Cbalance, Nbalance, Sibalance)
   end subroutine f_getBalance  

  subroutine f_getLost(u, Clost, Nlost, SiLost) bind(c)
    real(c_double), intent(in):: u(nGrid)
    real(c_double), intent(out):: Clost, Nlost, SiLost
    
    call getLost(u, Clost, Nlost, SiLost)
  end subroutine f_getLost
  
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

  subroutine f_setSinking(velocity) bind(c)
    use globals
    use NUMmodel, only: setSinking, nGrid

    real(c_double), intent(in):: velocity(nGrid)

    call setSinking(velocity)
  end subroutine f_setSinking
   
  subroutine f_getRates(jN, jDOC, jL, jSi, jF, jFreal, f, &
    jTot, jMax, jFmax, jR, jResptot, jLossPassive, &
    jNloss,jLreal, jPOM, &
    mortpred, mortHTL, mort2, mort) bind(c)
    use globals
    use NUMmodel, only: nNutrients, getRates
    real(c_double), intent(out):: jN(nGrid-nNutrients), jDOC(nGrid-nNutrients), jL(nGrid-nNutrients)
    real(c_double), intent(out):: jSi(nGrid-nNutrients)
    real(c_double), intent(out):: jF(nGrid-nNutrients), jFreal(nGrid-nNutrients), f(nGrid-nNutrients)
    real(c_double), intent(out):: jTot(nGrid-nNutrients), jMax(nGrid-nNutrients), jFmax(nGrid-nNutrients)
    real(c_double), intent(out):: jR(nGrid-nNutrients), jResptot(nGrid-nNutrients)
    real(c_double), intent(out):: jLossPassive(nGrid-nNutrients), jNloss(nGrid-nNutrients), jLreal(nGrid-nNutrients)
    real(c_double), intent(out):: jPOM(nGrid-nNutrients)
    real(c_double), intent(out):: mortpred(nGrid-nNutrients), mortHTL(nGrid-nNutrients)
    real(c_double), intent(out):: mort2(nGrid-nNutrients), mort(nGrid-nNutrients)

   call getRates(jN, jDOC, jL, jSi, jF, jFreal, f, &
   jTot, jMax, jFmax, jR, jResptot, jLossPassive, &
   jNloss,jLreal, jPOM, &
   mortpred, mortHTL, mort2, mort)
  end subroutine f_getRates
  
  subroutine f_getTheta(thetaMatrix) bind(c)
    real(c_double), intent (inout) :: thetaMatrix(nGrid,nGrid)

    thetaMatrix = theta
  end subroutine

end module NUMmodel_wrap
