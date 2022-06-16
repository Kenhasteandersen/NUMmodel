program NUMmodeltest
  use NUMmodel
  use globals
  implicit none

  real(dp), allocatable:: u0(:), u00(:), dudt(:)
  real(dp):: ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro
  integer:: i


  call setupGeneric( (/0.1d0, 1.0d0 /) )
  !call setHTL(0.0001d0, 1.d0, .true.)

  !call setupGeneralistsCopepod()
  call setHTL(0.1d0, 0.1d0, .false., .false.)
  !call setupGeneralistsOnly(5)
  !call setupGeneralistsPOM(10,5)
  !call setupNUMmodel(10,10,10, (/0.1d0, 1.0d0 /) )

  allocate(u0(nGrid))
  allocate(u00(nGrid))
  allocate(dudt(nGrid))
  u00(idxN) = 150.d0
  u00(idxDOC) = 10.d0
  !u00(idxSi) = 10.d0
  do i = idxB, nGrid
     u00(i) = 10 + 0.1*(i-2)
  end do
  !u00(17:22) = 0.d0 ! No POM
  dudt = 0.d0

 
  !write(*,*) group(1)%spec%velocity
  !write(*,*) group(2)%spec%velocity
  !write(*,*) group(3)%spec%velocity
  !write(*,*) group(4)%spec%velocity
  !call getSinking(u00)
  !write(*,*) u00
  !write(*,*) u00
  !call simulateChemostatEuler(u00, 60.d0, 100.d0, (/150.d0, 0.d0/), 0.01d0, .01d0, 0.01d0, logical(.true.,1))
  !write(*,*) u00
  
  !call simulateChemostatEuler(u00, 100.d0, 10.d0, u00(1:2), 0.1d0, 1000.d0, 0.1d0, logical(.false.,1))
  call calcDerivatives(u00, 20.d0, 20.d0, 0.0000001d0, dudt)
  !call printRates()

  !select type (spec => group(1)%spec)
  !    type is (spectrumGeneralists)
  !      write(*,*) getNbalanceGeneralists(spec, u00(idxN), dudt(idxN), u00(idxB:nGrid), dudt(idxB:nGrid))    
  !      write(*,*) getCbalanceGeneralists(spec, u00(idxDOC), dudt(idxDOC), u00(idxB:nGrid), dudt(idxB:nGrid))    
  !end select
  
  !write(*,*) 'dudt:',dudt
  !write(*,*) 'u',u00
  !call printRates()

  !write(*, '(6f10.6)') theta

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

!!$  u0=u00
!!$  call simulateChemostatEuler(u0, 100.d0, 150.d0, 0.05d0, 300.d0, 0.01d0)
!!$  call printU(u0)
 ! call calcDerivatives(u00, 150.d0, 0.1d0)

 ! write(6,*) theta(3:5, 3:5)
 !call printRates(m, rates)
 ! write(6,*) 'xxxx'
 ! call setupGeneric( (/0.1d0, 1.0d0 /) )
!  call setupGeneralistsOnly()

  !call getFunctions(ProdGross, ProdNet,ProdHTL,eHTL,Bpico,Bnano,Bmicro)
  !write(*,*) Bpico, Bnano, Bmicro
  
end program NUMmodeltest
