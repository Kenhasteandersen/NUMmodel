program NUMmodeltest
  use NUMmodel
  use globals
  implicit none

  real(dp), allocatable:: u0(:), u00(:), dudt(:)
  real(dp):: ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro
  integer:: i
  real(dp):: Nbalance,Cbalance, Sibalance

  !call setupGeneric( (/0.1d0, 1.0d0 /) )
  !call setHTL(0.0001d0, 1.d0, .true.)

  !call setupGeneralistsCopepod()
  !call setupGeneralistsOnly(10)
  !call setupGeneralistsSimpleOnly(10)

  !call setupGeneralistsOnly_csp()
 ! call setupGeneralistsOnly_csp()
  !call setupGeneralistsOnly_csp()
  !call setupGeneralistsOnly_csp()
  !call parametersFinalize(0.d0, .false.)
  
  !call setupGeneralistsDiatoms(10)
  !call setupGeneralistsDiatoms_simple(10)
  call setupGeneralistsOnly(10)
  !call setupGeneralistssimpleOnly(10)
  !call setupDiatoms_simpleOnly(10)
  !call setupDiatomsOnly(10)
  !call setupDiatoms_simpleOnly(10)

  allocate(u0(nGrid))
  allocate(u00(nGrid))
  allocate(dudt(nGrid))
  u00(idxN) = 50.d0
  u00(idxDOC) = 10.d0
  !u00(idxSi) = 10.d0
  do i = idxB, nGrid
     u00(i) = 1.0d0 !*(i-2)
  end do
  dudt = 0.d0

  !call simulateEuler(u00, 60.d0, 100.d0, 10.d0, 0.1d0)
  !                          ( u ,   L   ,   T  ,   Ndeep  , diff ,  tEnd  ,   dt , bLosses    )
  !call simulateChemostatEuler(u00, 100.d0, 10.d0, u00(1:2), 0.5d0, 1000.d0, 0.1d0, logical(.true.,1))
  !                      u  ,  L  ,   T  ,   dt , dudt
  call calcDerivatives(u00, 100.d0, 10.d0, 0.1d0, dudt)
  !write(*,*) u00

  ProdGross = 0
  ProdNet = 0
  ProdHTL=0
  ProdBact = 0
  eHTL=0
  Bpico=0
  Bnano=0
  Bmicro=0

  !call getFunctions(u00, ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro)
  !write(*,*) ProdGross, ProdNet,ProdHTL, ProdBact, eHTL
  !write(*,*) u00
  !call calcDerivatives(u00, 60.d0, 15.d0, 0.1d0, dudt)
  !call printRates()
  !!$  u0=u00
  !!$  call printU(u0)
  ! call calcDerivatives(u00, 150.d0, 0.1d0)

 ! write(6,*) theta(3:5, 3:5)
 ! 
 call printRates()
 !
 ! write(6,*) 'xxxx'
 ! call setupGeneric( (/0.1d0, 1.0d0 /) )
 !write(*,*) Bpico, Bnano, Bmicro
    call getBalance(u00, dudt, Nbalance,Cbalance,Sibalance)
    write(*,*) 'Nbalance:', Nbalance
    write(*,*) 'Cbalance:', Cbalance
    !write(*,*) 'Sibalance:', Sibalance

  end program NUMmodeltest
 