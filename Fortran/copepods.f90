!
! Module to handle copepods defined by their adult mass.
! All parameters are for active copepods. Follows Serra-Pompei et al (2020) with
! the update to respiration in Serra-Pompei (2022).
!
module copepods
  use globals
  use spectrum
  use input
  implicit none

  private

  !real(dp) :: rhoCN
  real(dp) :: epsilonF  ! Assimilation efficiency
  real(dp) :: epsilonR  ! Reproductive efficiency
  real(dp) :: beta
  real(dp) :: sigma 
  real(dp) :: alphaF  ! Clearance rate coefficient
  real(dp) :: q  ! Exponent of clerance rate
  real(dp) :: h  ! Factor for maximum ingestion rate
  real(dp) :: hExponent  ! Exponent for maximum ingestions rate
  real(dp) :: kBasal ! 0.006 ! Factor for basal metabolism.This value represents basal
  real(dp) :: kSDA  ! Factor for SDA metabolism (Serra-Pompei 2020). This value assumes that the
  real(dp) :: p  ! Exponent for respiration
  real(dp) :: AdultOffspring 
  real(dp) :: remin! fraction of mortality losses reminerilized to N and DOC
  !real(dp), parameter:: rhoCN = 5.68
  !real(dp), parameter:: epsilonF = 0.67 ! Assimilation efficiency
  !real(dp), parameter:: epsilonR = 0.25 ! Reproductive efficiency
  !real(dp), parameter:: beta = 10000.d0
  !real(dp), parameter:: sigma = 1.5d0
  !real(dp), parameter:: alphaF = 0.011 ! Clearance rate coefficient
  !real(dp), parameter:: q = 0.75 ! Exponent of clerance rate
  !real(dp), parameter:: h = 1.37 ! Factor for maximum ingestion rate
  !real(dp), parameter:: hExponent = 0.75 ! Exponent for maximum ingestions rate
  !real(dp), parameter:: kBasal = 0.006! 0.006 ! Factor for basal metabolism. This value represents basal
                                       ! metabolism at starvation. Following Kiørboe (1985)
                                       ! the starvation metabolism is approx 0.2*0.18=0.036 times 
                                       ! the maximum metabolism (kSDA). Increased to 0.01 to avoid too long transients.
  !real(dp), parameter:: kSDA = 0.16 ! Factor for SDA metabolism (Serra-Pompei 2020). This value assumes that the
                                    ! data in Kiørboe and Hirst (2014) are for fully fed copepods.
  !real(dp), parameter:: p = 0.75 ! Exponent for respiration
  !real(dp), parameter:: AdultOffspring = 100.
  !real(dp), parameter:: remin = 0.0 ! fraction of mortality losses reminerilized to N and DOC
  !real(dp), parameter:: remin2 = 1.d0 ! fraction of virulysis remineralized to N and DOC

  type, extends(spectrumMulticellular) :: spectrumCopepod
    real(dp), allocatable :: gamma(:), g(:), mortStarve(:), mort(:), JrespFactor(:)
  contains
    procedure, pass :: initCopepod
    procedure :: calcDerivativesCopepod
    procedure :: printRates => printRatesCopepod
  end type spectrumCopepod

  public spectrumCopepod, initCopepod, calcDerivativesCopepod, printRatesCopepod
contains

  subroutine read_namelist()
    integer :: file_unit,io_err

    namelist /input_copepods / epsilonF, epsilonR, beta, sigma, alphaF,q, &
             & h, hExponent, kBasal, kSDA, p, AdultOffspring, remin

    call open_inputfile(file_unit, io_err)
    read(file_unit, nml=input_copepods, iostat=io_err)
    call close_inputfile(file_unit, io_err)

  end subroutine read_namelist

  subroutine initCopepod(this, n, mAdult)
    class(spectrumCopepod), intent(inout):: this
    integer, intent(in):: n
    real(dp), intent(in):: mAdult
    real(dp):: lnDelta, mMin

    call read_namelist()
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
    allocate(this%JrespFactor(n))

    this%beta = beta
    this%sigma = sigma
    this%epsilonF = epsilonF
    this%AF = alphaF*this%m**q
    this%JFmax = h*this%m**hExponent
    this%JrespFactor = epsilonF*this%JFmax
    this%mort2constant = 0.d0 ! No quadratic mortality
    this%mort2 = 0.d0

    this%mPOM = 3.5e-3*this%m ! Size of fecal pellets (Serra-Pompei 2022 approximated)
  end subroutine initCopepod

  subroutine calcDerivativesCopepod(this, u, dNdt, dudt)
    class(spectrumCopepod), intent(inout):: this
    real(dp), intent(in):: u(this%n)
    real(dp), intent(inout):: dNdt, dudt(this%n)
    integer:: i
    real(dp):: nu, b

    do i = 1, this%n
       !
       ! Growth and reproduction:
       !

       ! Basal and SDA respiration:
       this%Jresp(i) = this%JrespFactor(i) * kBasal * fTemp2 +  kSDA * this%JF(i)
       ! Available energy:
       nu = this%JF(i) - this%Jresp(i)
       ! Production of POM:
       this%jPOM = (1-epsilonF)*this%JF(i)/(this%m(i) * epsilonF)
       ! Available energy rate (1/day):
       this%g(i) = max(0.d0, nu)/this%m(i)
       ! Starvation:
       this%mortStarve(i) = -min(0.d0, nu)/this%m(i)
       !
       ! Mortality:
       !this%mortHTL(i) = this%mortHTL(i)*u(i)

       !this%mortHTL(i) = this%mortHTL(i) * fTemp2     

       this%mort(i) = this%mortpred(i) + this%mortStarve(i) + this%mortHTL(i)
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

    dNdt = dNdt + sum( this%Jresp*u/(this%m*rhoCN) )  &! All respiration of carbon results in a corresponding
                                    ! surplus of nutrients. This surplus (pee) is routed to nutrients
                + (1-epsilonR)*this%g(this%n)*u(this%n)/rhoCN  ! Should perhaps also go to DOC
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
