program NUMmodeltest
  use NUMmodel

  real(dp), allocatable:: u0(:)
  integer:: i

  !call setupGeneric( (/0.1d0, 1.0d0 /) )
  call setupGeneralistsCopepod()
  !call setupGeneralistsOnly()
  allocate(u0(nGrid))
  u0(1) = 150.d0
  u0(2) = 1.d0
  do i = 3, nGrid
     u0(i) = 1.0d0
  end do
  !call calcDerivatives(u, 100.d0, 0.d0)
  call simulateChemostatEuler(u0, 100.d0, 150.d0, 0.05d0, 30.d0, 0.01d0)
  !call calcDerivatives(usave(size(usave,1),:), 100.d0, 0.1d0)
  !call printRates(m, rates)
 ! write(6,*) 'xxxx'
  call printU(u0)

end program NUMmodeltest
