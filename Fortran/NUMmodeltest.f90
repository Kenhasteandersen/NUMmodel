program NUMmodeltest
  use NUMmodel
  use globals
  implicit none

  real(dp), allocatable:: u0(:), u00(:), dudt(:)
  real(dp):: ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro
  integer:: i
  real(dp):: Nbalance,Cbalance, Sibalance

  !call setupNUMmodel( (/0.1d0, 1.0d0 /) )

  !call setupGeneralistsCopepod()
  !call setupGeneralistsOnly(10)
  !call setupGeneralistsSimpleOnly(10)

  !call parametersFinalize(0.d0, .false.)
  
  !call setupGeneralistsDiatoms(10)
  !call setupGeneralistsDiatoms_simple(10)
  ! setupGeneralistsOnly(10)
  !call setupGenDiatCope(3,3,(/0.1d0, 1.0d0 /))
  !call setupGenDiatCope(3,5,1,(/0.1d0, 1.0d0 /))
   !               2 gens cop POM   mAdult     
   !call setupNUMmodel(3 , 1 , 2 ,(/0.1d0 /), (/1.d0 /))
   !call setupGenDiatCope(3 , 1 , 2 ,(/0.1d0 /))

   !              gen-diat-cop      POM      mAdult    
  !call setupGenDiatCope(3,   2,    1,    (/0.1d0, 1.d0/))

  !call setupGeneralistssimpleOnly(10)
  !call setupDiatoms_simpleOnly(10)
  !call setupDiatomsOnly(10)
  !call setupDiatoms_simpleOnly(10)
  !call setHTL(0.1d0, 0.1d0, .false., .false.)
  !call setupGeneralistsOnly(5)
  !call setupGeneralistsDiatoms(10)
  !call setupGeneralistsDiatoms_simple(10)
  !call setupGeneralistsPOM(10,5)
  call setupNUMmodel(2,2,1, (/1.d0 /), (/1.d0/) )
  !call setupNUMmodelsimple(10,10,10, (/0.1d0, 1.0d0/) )
  !call setupGeneric( (/0.1d0/) )

  !call setHTL(0.0d0, 0.1d0, .true., .true.)

  allocate(u0(nGrid))
  allocate(u00(nGrid))
  allocate(dudt(nGrid))
  u00(idxN) = 50.d0
  u00(idxDOC) = 10.d0
  u00(idxSi) = 10.d0
  do i = idxB, nGrid
     u00(i) = 10! + 0.1*(i-2)
  end do
  dudt = 0.d0

  !call simulateEuler(u00, 60.d0, 100.d0, 10.d0, 0.1d0)
  !                          ( u ,   L   ,   T  ,   Ndeep  , diff ,  tEnd  ,   dt , bLosses    )
  !call simulateChemostatEuler(u00, 100.d0, 10.d0, u00(1:2), 0.5d0, 1000.d0, 0.1d0, logical(.true.,1))
  !                      u  ,  L  ,   T  ,   dt , dudt
  
  !call simulateChemostatEuler(u00, 100.d0, 10.d0, u00(1:2), 0.1d0, 1000.d0, 0.1d0, logical(.false.,1))
  !call calcDerivatives(u00, 20.d0, 20.d0, 0.0000001d0, dudt)
  !call printRates()

  !select type (spec => group(1)%spec)
  !    type is (spectrumGeneralists)
  !      write(*,*) getNbalanceGeneralists(spec, u00(idxN), dudt(idxN), u00(idxB:nGrid), dudt(idxB:nGrid))    
  !      write(*,*) getCbalanceGeneralists(spec, u00(idxDOC), dudt(idxDOC), u00(idxB:nGrid), dudt(idxB:nGrid))    
  !end select
  
  
  call calcDerivatives(u00, 100.d0, 10.d0, 0.0d0, dudt)
  write(*,*) dudt
  !write(*,*) 'ngrid',nGrid
  !write(*,*) 'ngroups',nGroups
  !write(*,*) 'nbutrients',nNutrients

  ProdGross = 0
  ProdNet = 0
  ProdHTL=0
  ProdBact = 0
  eHTL=0
  Bpico=0
  Bnano=0
  Bmicro=0

 ! call getFunctions(u00, ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro)
  !write(*,*) ProdGross, ProdNet,ProdHTL, ProdBact, eHTL
  !write(*,*) u00
  !call calcDerivatives(u00, 60.d0, 15.d0, 0.1d0, dudt)
  !call printRates()
  !!$  u0=u00
  !!$  call printU(u0)
  ! call calcDerivatives(u00, 150.d0, 0.1d0)

 ! write(6,*) theta(3:5, 3:5)
 ! 
 !call printRates()
 !
 ! write(6,*) 'xxxx'
 ! call setupGeneric( (/0.1d0, 1.0d0 /) )
 !write(*,*) Bpico, Bnano, Bmicro
    call getBalance(u00, dudt, Nbalance,Cbalance,Sibalance)
    write(*,*) 'Nbalance:', Nbalance
    write(*,*) 'Cbalance:', Cbalance
    write(*,*) 'Sibalance:', Sibalance

!do i = 5,9
!   write(*,*) i, theta(i+3,6:9)
!end do

  end program NUMmodeltest
 