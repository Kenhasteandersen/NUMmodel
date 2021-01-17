program NUMmodeltest
  use NUMmodel

  real(dp), allocatable:: u0(:), u00(:)
  integer:: i

  !call setupGeneric( (/0.1d0, 1.0d0 /) )
  !call setupGeneralistsCopepod()
  call setupGeneralistsOnly()
  !call setupGeneralistsOnly()
  allocate(u0(nGrid))
  allocate(u00(nGrid))
  u00(1) = 150.d0
  u00(2) = 1.d0
  do i = 3, nGrid
     u00(i) = 1.0d0
  end do
  !call calcDerivatives(u, 100.d0, 0.d0)
  u0=u00
  call simulateChemostatEuler(u0, 100.d0, 150.d0, 0.05d0, 300.d0, 0.01d0)
  call printU(u0)
  !call calcDerivatives(usave(size(usave,1),:), 100.d0, 0.1d0)
  !call printRates(m, rates)
 ! write(6,*) 'xxxx'

end program NUMmodeltest
