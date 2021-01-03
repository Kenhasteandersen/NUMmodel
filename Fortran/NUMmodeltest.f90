program NUMmodeltest
  use NUMmodel

  real(dp), allocatable:: u0(:), usave(:,:)
  integer:: i

  !call setupGeneric( (/0.1d0, 1.0d0 /) )
  call setupGeneralistsCopepod()
  !call setupGeneralistsOnly()
  allocate(u0(nGrid))
  u0(1) = 150.d0
  u0(2) = 1.d0
  do i = 1, 20
     u0(2+i) = i*1.0d0
  end do
  !call calcDerivatives(u, 100.d0, 0.d0)
  usave = simulateChemostatEuler(u0, 100.d0, 0.05d0, 365.d0, 0.001d0)
  !call calcDerivatives(usave(size(usave,1),:), 100.d0, 0.1d0)
  call printRates(m, rates)
 ! write(6,*) 'xxxx'
  call printU(usave(size(usave,1),:))

end program NUMmodeltest
