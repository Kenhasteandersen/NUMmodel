module copepods
  use globals
  use sizespectrum
  implicit none

  real(dp), parameter:: rhoCN = 5.68
  real(dp), parameter:: epsilonFF = 0.67 ! Assimilation efficiency
  real(dp), parameter:: epsilonR = 0.25 ! Reproductive efficiency
  real(dp), parameter:: alphaF = 0.1 !  PROBABLY WRONG!
  real(dp), parameter:: q = 0.75 ! Exponent of clerance rate
  real(dp), parameter:: h = 1.37 ! Factor for maximum ingestion rate
  real(dp), parameter:: hExponent = 0.75 ! Exponent for maximum ingestions rate
  real(dp), parameter:: Kappa = 0.16 ! Factor for respiration
  real(dp), parameter:: p = 0.75 ! Exponent for respiration
  real(dp), parameter:: AdultOffspring = 100.
  real(dp), parameter:: remin = 0.0 ! fraction of mortality losses reminerilized to N and DOC
  real(dp), parameter:: remin2 = 1.d0 ! fraction of virulysis remineralized to N and DOC
  
  real(dp),  dimension(:), allocatable:: Jresp(:), gamma(:)
  
contains
  
  function initCopepod(n, ixStart, mAdult) result(this)
    type(typeSpectrum):: this
    integer, intent(in):: n, ixStart
    real(dp), intent(in):: mAdult
    real(dp):: lnDelta, mMin
    !
    ! Calc grid. Grid runs from mLower(1) = offspring size to m(n) = adult size
    !
    lnDelta = (log(mAdult)-log(mAdult/AdultOffspring)) / (n-0.5)
    mMin = exp(log(mAdult/AdultOffspring)+0.5*lnDelta);
    this = initSpectrum(n, ixStart, mMin, mAdult)

    if (.not. allocated(Jresp)) then
       allocate(Jresp(n))
       allocate(gamma(n))
    end if
    
    this%AF = alphaF*this%m**q
    this%JFmax = h*this%m**hExponent
    Jresp = Kappa*this%m**p
    !
    ! Export grid to globals:
    !
    m(this%ixStart : this%ixEnd) = this%m
    beta(this%ixStart : this%ixEnd) = 10000.
    sigma(this%ixStart : this%ixEnd) = 1.5
    !
    ! Export feeding parameters:
    !
    AF(this%ixStart : this%ixEnd) = this%AF
    JFmax(this%ixStart : this%ixEnd) = this%JFmax
    epsilonF(this%ixStart : this%ixEnd) = epsilonFF
  end function initCopepod

  subroutine calcDerivativesCopepod(this, u, rates)
    type(typeSpectrum), intent(in):: this
    real(dp), intent(in):: u(:)
    type(typeRates), intent(inout):: rates
    integer:: ix, i
    real(dp):: nu, b
    
    do ix = this%ixStart, this%ixEnd
       i = ix-this%ixStart+1
       !
       ! Growth and reproduction:
       !
       nu = epsilonFF*rates%JF(ix) - Jresp(i)
       rates%g(ix) = max(0.d0, nu)/this%m(i)
       rates%mortStarve(ix) = -min(0.d0, nu)/this%m(i)
       b = epsilonR * rates%g(this%ixEnd) ! Birth rate
       !
       ! Mortality:
       !
       rates%mort(ix) = rates%mortpred(ix) + rates%mortHTL(ix) + rates%mortStarve(ix)
         ! Flux:
         gamma(i) = (rates%g(ix)-rates%mort(ix)) / (1 - this%z(i)**(1-rates%mort(ix)/rates%g(ix)))
         rates%Jtot(ix) = nu
      end do
      !
      ! Assemble derivative:
      ! 
      rates%dudt(this%ixStart) = b*u(this%ixEnd) &
           + (rates%g(this%ixStart)-gamma(1)-rates%mort(this%ixStart))*u(this%ixStart)
      do ix = this%ixStart+1, this%ixEnd-1
         i = ix-this%ixStart+1
         rates%dudt(ix) = &
              gamma(i-1)*u(ix-1) & ! growth into group
              + (rates%g(ix)-gamma(i)-rates%mort(ix))*u(ix)  ! growth out of group
      end do
      rates%dudt(this%ixEnd) = &
           gamma(this%n-1)*u(this%ixEnd-1) & ! growth into adult group
           - rates%mort(this%ixEnd)*u(this%ixEnd); ! adult mortality
    end subroutine calcDerivativesCopepod
    
  end module copepods
  
