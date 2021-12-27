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

  type, extends(spectrumMulticellular) :: spectrumCopepod
    real(dp), allocatable :: gamma(:), g(:), mortStarve(:), mort(:)
  contains
    procedure, pass :: initCopepod
    procedure :: calcDerivativesCopepod
    procedure :: printRates => printRatesCopepod
  end type spectrumCopepod

  public spectrumCopepod, initCopepod, calcDerivativesCopepod, printRatesCopepod
contains

  subroutine initCopepod(this, n, mAdult)
    class(spectrumCopepod), intent(inout):: this
    integer, intent(in):: n
    real(dp), intent(in):: mAdult
    real(dp):: lnDelta, mMin
    !
    ! Calc grid. Grid runs from mLower(1) = offspring size to m(n) = adult size
    !
    lnDelta = (log(mAdult)-log(mAdult/AdultOffspring)) / (n-0.5)
    mMin = exp(log(mAdult/AdultOffspring)+0.5*lnDelta);
    call this%initSpectrum(n, mMin, mAdult)

    allocate(this%gamma(n))
    allocate(this%g(n))
    allocate(this%mortStarve(n))
    allocate(this%mort(n))

    this%beta = beta
    this%sigma = sigma
    this%epsilonF = epsilonF
    this%AF = alphaF*this%m**q
    this%JFmax = h*this%m**hExponent
    this%Jresp = Kappa*this%m**p
  end subroutine initCopepod

  subroutine calcDerivativesCopepod(this, u, dudt)
     class(spectrumCopepod), intent(inout):: this
    real(dp), intent(in):: u(this%n)
    real(dp), intent(inout):: dudt(this%n)
    integer:: i
    real(dp):: nu, b

    do i = 1, this%n
       !
       ! Growth and reproduction:
       !
       nu = epsilonF*this%JF(i) - this%Jresp(i)
       this%g(i) = max(0.d0, nu)/this%m(i)
       this%mortStarve(i) = -min(0.d0, nu)/this%m(i)
       !
       ! Mortality:
       !this%mortHTL(i) = this%mortHTL(i)*u(i)
       this%mort(i) = this%mortpred(i) + this%mortHTL(i)*h + this%mortStarve(i)
       ! Flux:
       if ( this%g(i) .ne. 0.) then
         this%gamma(i) = (this%g(i)-this%mort(i)) / (1 - this%z(i)**(1-this%mort(i)/this%g(i)))
       else
         this%gamma(i) = 0.d0
       end if
       this%Jtot(i) = nu
    end do
    b = epsilonR * this%g(this%n) ! Birth rate
    !
    ! Assemble derivatives:
    !
    ! 1st stage:
    dudt(1) = b*u(this%n) &
         + (this%g(1) - this%gamma(1) - this%mort(1))*u(1)
    ! Middle stages:
    do i = 2, this%n-1
       dudt(i) = &
            this%gamma(i-1)*u(i-1) & ! growth into group
            + (this%g(i) - this%gamma(i) - this%mort(i))*u(i)  ! growth out of group
    end do
    !Adults
    dudt(this%n) = &
         this%gamma(this%n-1)*u(this%n-1) & ! growth into adult group
         - this%mort(this%n)*u(this%n); ! adult mortality
  end subroutine calcDerivativesCopepod

  subroutine printRatesCopepod(this)
   class(spectrumCopepod), intent(in):: this

   write(*,*) "Copepod with ", this%n, " size classes and adult size ", this%m(this%n), "ugC."
     call this%printRatesSpectrum()

     99 format (a10, 20f10.6)
 
     write(*,99) "gamma:", this%gamma
     write(*,99) "mortStarve:", this%mortStarve
     write(*,99) "g:", this%g
  end subroutine printRatesCopepod

end module copepods
