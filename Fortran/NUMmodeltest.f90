!
! This file contains test code to run the NUM model directly without interfacing through matlab 
! or R. It is a mess as it is only used for testing purposes.
!
program NUMmodeltest
  use NUMmodel
  use globals
  implicit none

  real(dp), allocatable:: u0(:), u00(:), dudt(:)
  real(dp):: ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro,mHTL
  integer:: i
  real(dp):: Nbalance,Cbalance, Sibalance
  logical(1) :: TRUE1
  logical(1):: FALSE1
  !real(dp):: myout
  character(len=20) :: errorstr
  logical(1):: errorio=.false. 

  TRUE1 = .true.
  FALSE1 = .false.
  !call setupNUMmodel( (/0.1d0, 1.0d0 /) )

  !call setupGeneralistsOnly(10)
  !call setupGeneralistsSimpleOnly(10)

  ! !call setupGenDiatCope(3 , 1 , 2 ,(/0.1d0 /), errorio, errorstr)

  !call setupGeneralistsDiatoms(6,errorio,errorstr)
  !call setupDiatomsOnly(10,errorio,errorstr)
  
  !call setupGeneralistsOnly(5,errorio,errorstr)
  !call setupGeneralistsPOM(5,1, errorio, errorstr)
  !call setupNUMmodel(10,6,1, (/.2d0, 5.d0 /), (/1.d0, 31.6d0, 1000.d0/) ,errorio,errorstr)
  !call setupGeneralistsDiatoms(10, errorio, errorstr)
  !call setupGeneric( (/1.d0 /), errorio, errorstr )
  !call setupGeneralistsDiatoms(10, errorio, errorstr)
  !call setupNUMmodel(5,5,1, (/1.d0 /), (/10.d0/),errorio,errorstr)
  call setupNUMmodelGelatinous(5,5,1, (/1.d0 /), (/10.d0/),(/5.d0/),errorio,errorstr)
  !call setHTL(0.005d0, 0.1d0, TRUE1, FALSE1, )

  if (errorio .eqv. .false.) then
    print*, 'Parameters loaded correctly'
  else
    print*, 'Error loading parameter ', errorstr
  end if


  allocate(u0(nGrid))
  allocate(u00(nGrid))
  allocate(dudt(nGrid))
  u00(idxN) = 1.d0
  u00(idxDOC) = 0.01d0
  u00(idxSi) = 0.01d0
  do i = idxB, nGrid
     u00(i) = 0.0686d0 !+ 0.1*(i-2)
  end do
  dudt = 0.d0

  !call getSinking(u00)
  !write(*,*) u00
  !u00(8:12) = 5.d0

  !call simulateEuler(u00, 60.d0, 100.d0, 10.d0, 0.1d0)
  !                          ( u ,   L   ,   T  ,   Ndeep  , diff ,  tEnd  ,   dt , bLosses    )
  !call simulateChemostatEuler(u00, 100.d0, 10.d0, u00(1:2), 0.5d0, 1000.d0, 0.1d0, logical(.true.,1))
  !                      u  ,  L  ,   T  ,   dt , dudt
  
  call simulateChemostatEuler(u00, 100.d0, 10.d0, u00(1:2), 0.1d0, 1000.d0, 0.1d0, logical(.false.,1))
  !call calcDerivatives(u00, 20.d0, 20.d0, 0.0000001d0, dudt)
  !call printRates()

  !select type (spec => group(1)%spec)
  !    type is (spectrumGeneralists)
  !      write(*,*) getNbalanceGeneralists(spec, u00(idxN), dudt(idxN), u00(idxB:nGrid), dudt(idxB:nGrid))    
  !      write(*,*) getCbalanceGeneralists(spec, u00(idxDOC), dudt(idxDOC), u00(idxB:nGrid), dudt(idxB:nGrid))    
  !end select
  
call setHTL(0.1d0, 1.d0, logical(.false.,1),logical(.false.,1),logical(.true.,1))

  call calcDerivatives(u00, 60.d0, 15.d0, 0.1d0, dudt)
  call printRates()
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
  mHTL=0.d0

 !call getFunctions(u00, ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro,mHTL)
!write(*,*) ProdGross, ProdNet,ProdHTL, ProdBact, eHTL,mHTL
  !write(*,*) dudt
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
  write(*,*) dudt
  call getBalance(u00, dudt, Cbalance,Nbalance,Sibalance)
  write(*,*) 'Cbalance:', Cbalance
  write(*,*) 'Nbalance:', Nbalance
  write(*,*) 'Sibalance:', Sibalance
   

!
  end program NUMmodeltest