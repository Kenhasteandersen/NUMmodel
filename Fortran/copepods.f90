!
! Module to handle copepods defined by their adult mass.
! All parameters are for active copepods
!
module copepods
  use globals
  use spectrum
  implicit none

  private
  real(dp), parameter:: rhoCN = 5.68
  real(dp), parameter:: epsilonF = 0.67 ! Assimilation efficiency
  real(dp), parameter:: epsilonR = 0.25 ! Reproductive efficiency
  real(dp), parameter:: beta = 10000.d0
  real(dp), parameter:: sigma = 1.5d0
  real(dp), parameter:: alphaF = 0.01 !  
  real(dp), parameter:: q = 0.75 ! Exponent of clerance rate
  real(dp), parameter:: h = 1.37 ! Factor for maximum ingestion rate
  real(dp), parameter:: hExponent = 0.75 ! Exponent for maximum ingestions rate
  real(dp), parameter:: Kappa = 0.16 ! Factor for respiration
  real(dp), parameter:: p = 0.75 ! Exponent for respiration
  real(dp), parameter:: AdultOffspring = 100.
  real(dp), parameter:: remin = 0.0 ! fraction of mortality losses reminerilized to N and DOC
  real(dp), parameter:: remin2 = 1.d0 ! fraction of virulysis remineralized to N and DOC

  real(dp),  dimension(:), allocatable:: Jresp(:), gamma(:)

  public initCopepod, calcDerivativesCopepod
contains

  function initCopepod(n, ixOffset, mAdult) result(this)
    type(typeSpectrum):: this
    integer, intent(in):: n, ixOffset
    real(dp), intent(in):: mAdult
    real(dp):: lnDelta, mMin
    !
    ! Calc grid. Grid runs from mLower(1) = offspring size to m(n) = adult size
    !
    lnDelta = (log(mAdult)-log(mAdult/AdultOffspring)) / (n-0.5)
    mMin = exp(log(mAdult/AdultOffspring)+0.5*lnDelta);
    this = initSpectrum(n, ixOffset, mMin, mAdult)

    if (.not. allocated(Jresp)) then
       allocate(Jresp(n))
       allocate(gamma(n))
    end if

    this%beta = beta
    this%sigma = sigma
    this%epsilonF = epsilonF
    this%AF = alphaF*this%m**q
    this%JFmax = h*this%m**hExponent
    Jresp = Kappa*this%m**p
  end function initCopepod

  subroutine calcDerivativesCopepod(this, u, rates)
    type(typeSpectrum), intent(in):: this
    real(dp), intent(in):: u(:)
    type(typeRates), intent(inout):: rates
    integer:: ix, i
    real(dp):: nu, b

    do i = 1, this%n
       ix = i+this%ixOffset
       !
       ! Growth and reproduction:
       !
       nu = epsilonF*rates%JF(ix) - Jresp(i)
       rates%g(ix) = max(0.d0, nu)/this%m(i)
       rates%mortStarve(ix) = -min(0.d0, nu)/this%m(i)
       !
       ! Mortality:
       !
       rates%mort(ix) = rates%mortpred(ix) + rates%mortHTL(ix) + rates%mortStarve(ix)
       ! Flux:
       gamma(i) = (rates%g(ix)-rates%mort(ix)) / (1 - this%z(i)**(1-rates%mort(ix)/rates%g(ix)))
       rates%Jtot(ix) = nu
    end do
    b = epsilonR * rates%g(this%ixEnd) ! Birth rate
    !
    ! Assemble derivatives:
    !
    ! 1st stage:
    rates%dudt(this%ixStart) = b*u(this%n) &
         + (rates%g(this%ixStart)-gamma(1)-rates%mort(this%ixStart))*u(1)
    ! Middle stages:
    do i = 2, this%n-1
       ix = i+this%ixOffset
       rates%dudt(ix) = &
            gamma(i-1)*u(i-1) & ! growth into group
            + (rates%g(ix)-gamma(i)-rates%mort(ix))*u(i)  ! growth out of group
    end do
    !Adults
    rates%dudt(this%ixEnd) = &
         gamma(this%n-1)*u(this%n-1) & ! growth into adult group
         - rates%mort(this%ixEnd)*u(this%n); ! adult mortality
  end subroutine calcDerivativesCopepod
  
  end module copepods
