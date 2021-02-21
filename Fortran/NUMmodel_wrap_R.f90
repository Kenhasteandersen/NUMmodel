
  subroutine f_setupGeneralistsOnly()
    use NUMmodel, only:  setupGeneralistsOnly
    use globals
    call setupGeneralistsOnly()
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
  
  subroutine f_calcDerivatives(nn, u, L, dt, dudt)
    use globals
    use NUMmodel, only: calcDerivatives, rates
    integer, intent(in):: nn
    real(dp), intent(in):: L, dt
    real(dp), intent(in):: u(nn)
    real(dp), intent(out):: dudt(nn)

    call calcDerivatives(u, L, dt)
    dudt = rates%dudt
  end subroutine f_calcDerivatives
!!$
!!$  subroutine f_calcRates(nGrid, u, L, jN, jL, jF) 
!!$    integer, intent(in):: nGrid
!!$    real(dp), intent(in):: u(nGrid)
!!$    real(dp), intent(in):: L
!!$    real(dp), intent(out):: jN(nGrid), jL(nGrid), jF(nGrid)
!!$
!!$    call calcDerivatives(u, L, 0.d0)
!!$    jN(idxB:nGrid) = rates%JN(idxB:nGrid) / m(idxB:nGrid)
!!$    jL(idxB:nGrid) = rates%JL(idxB:nGrid) / m(idxB:nGrid)
!!$    jF(idxB:nGrid) = rates%JF(idxB:nGrid) / m(idxB:nGrid)
!!$  end subroutine f_calcRates
!!$    
!!$  subroutine f_simulateChemostatEuler(nGrid, u, L, Ndeep, diff, tEnd, dt) 
!!$    integer, intent(in):: nGrid
!!$    real(dp), intent(inout):: u(nGrid)
!!$    real(dp), intent(in):: L, Ndeep, diff, tEnd, dt
!!$    
!!$    call simulateChemostatEuler(u, L, Ndeep, diff, tEnd, dt)
!!$  end subroutine f_simulateChemostatEuler
!!$  
!!$  subroutine f_simulateEuler(nGrid, u, L, tEnd, dt) 
!!$    integer, intent(in):: nGrid
!!$    real(dp), intent(inout):: u(nGrid)
!!$    real(dp), intent(in)!:: L, tEnd, dt
!!$    
!!$    call simulateEuler(u, L, tEnd, dt)
!!$  end subroutine f_simulateEuler
!!$  



