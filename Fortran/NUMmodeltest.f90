program NUMmodeltest
  use NUMmodel
  implicit none

  real(dp), allocatable:: u0(:), u00(:)
  !real(dp):: ProdGross, ProdNet,ProdHTL,eHTL,Bpico,Bnano,Bmicro
  integer:: i


  !call setupGeneric( (/0.1d0, 1.0d0 /) )
  !call setupGeneralistsCopepod()
  !call setupGeneralistsOnly(10)
  !call setupGeneralistsOnly_csp()
 ! call setupGeneralistsOnly_csp()
  call setupGeneralistsOnly_csp()
  !call setupGeneralistsOnly(25)

  !call setupGeneralistsDiatoms(10)
  !call setupDiatoms_simpleOnly(10)
  !call setupDiatoms_simpleOnly(10)

  allocate(u0(nGrid))
  allocate(u00(nGrid))
  u00(idxN) = 150.d0
  u00(idxDOC) = 1.d0
  u00(idxSi) = 7.d0
  do i = idxB, nGrid
     u00(i) = 1.0d0*(i-2)
  end do

  call calcDerivatives(u00, 60.d0, 10.d0, 10.d0)
  !call printRates(m,rates)
 
  !call simulateEuler(u00, 60.d0, 100.d0, 0.1d0)
  !call simulateChemostatEuler(u00, 60.d0, 20.d0, u00(1:3), 0.1d0, 0.2d0, 0.1d0)
  
  !write(*,*) u00
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
