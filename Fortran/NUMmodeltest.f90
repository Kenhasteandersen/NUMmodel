program NUMmodeltest
  use NUMmodel

!  real(dp), allocatable:: usave(:,:)
  integer:: i

  !call setupGeneric( (/0.1d0, 1.0d0 /) )
  call setupGeneralistsCopepod()
  !call setupGeneralistsOnly()
  u0(1) = 150.d0
  u0(2) = 1.d0
  do i = 1, 20
     u0(2+i) = i*1.0d0
  end do
  call calcDerivatives(u0, 100.d0, 0.d0)
  !call simulateChemostatEuler(100.d0, 0.05d0, 3650.d0, 0.1d0, usave)
  !call calcDerivatives(usave(size(usave,1),:), 100.d0, 0.1d0)
  call printRates(m, rates)
 ! write(6,*) 'xxxx'
  !call printU(usave(size(usave,1),:))

end program NUMmodeltest
