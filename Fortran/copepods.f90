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

  integer, parameter :: passive = 1
  integer, parameter :: active = 2

  private

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
  real(dp) :: AdultOffspring ! Adult-offspring mass ratio
  real(dp) :: vulnerability ! Passed to "palatability" in the parent spectrum class
  real(dp) :: DiatomsPreference

  type, extends(spectrumMulticellular) :: spectrumCopepod
    integer :: feedingmode ! Active=1; Passive=2
    real(dp), allocatable :: gamma(:), g(:), mortStarve(:), mort(:), JrespFactor(:)
    real(dp) :: DiatomsPreference ! Feeding preference on diatoms
  contains
    procedure, pass :: initCopepod
    procedure :: calcDerivativesCopepod
    procedure :: printRates => printRatesCopepod
    procedure :: getNbalance
    procedure :: getCbalance
  end type spectrumCopepod
  
  public active, passive, spectrumCopepod, initCopepod, calcDerivativesCopepod, printRatesCopepod

contains

  subroutine read_namelist(feedingmode)
    integer, intent(in) :: feedingmode
    integer :: file_unit,io_err

    namelist /input_copepods_passive / epsilonF, epsilonR, beta, sigma, alphaF,q, &
             & h, hExponent, kBasal, kSDA, AdultOffspring, vulnerability, DiatomsPreference
    namelist /input_copepods_active /  epsilonF, epsilonR, beta, sigma, alphaF,q, &
             & h, hExponent, kBasal, kSDA, AdultOffspring, vulnerability, DiatomsPreference
    
    call open_inputfile(file_unit, io_err)

    if (feedingmode .eq. active) then
      read(file_unit, nml=input_copepods_active, iostat=io_err)
    else
      if (feedingmode .eq. passive) then
        read(file_unit, nml=input_copepods_passive, iostat=io_err)
      else
        stop
      endif
    endif

    call close_inputfile(file_unit, io_err)
  end subroutine read_namelist

  subroutine initCopepod(this, feedingmode, n, mAdult)
    class(spectrumCopepod), intent(inout):: this
    integer, intent(in) :: feedingmode ! Whether the copepods is active or passive
    integer, intent(in):: n
    real(dp), intent(in):: mAdult
    integer:: i

    this%feedingmode = feedingmode
    
    call read_namelist(this%feedingmode)
    this%DiatomsPreference = DiatomsPreference
    !
    ! Calc grid. Grid runs from mLower(1) = offspring size to m(n) = adult size
    !
    call this%initMulticellular(n, mAdult/AdultOffspring, mAdult)

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
    this%palatability = vulnerability

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
       this%Jresptot(i) = this%JrespFactor(i) * kBasal * fTemp2 + kSDA * this%JF(i)
       ! Available energy:
       nu = this%JF(i) - this%Jresptot(i)
       ! Available energy rate (1/day):
       this%g(i) = max(0.d0, nu)/this%m(i)
       ! Starvation:
       this%mortStarve(i) = -min(0.d0, nu)/this%m(i)
       !
       ! Mortality:
       !
       this%mortHTL(i) = this%mortHTL(i) * fTemp2     

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
    ! Production of POM:
    this%jPOM = &
          (1-epsilonF)*this%JF/(this%m * epsilonF) & ! Unassimilated food (fecal pellets)
        + this%mortStarve                            ! Copepods dead from starvation
    this%jPOM(this%n) = this%jPOM(this%n) + (1.d0-epsilonR)*this%g(this%n) ! Lost reproductive flux
  
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

    dNdt = dNdt + sum( this%Jresptot*u/this%m )/rhoCN  ! All respiration of carbon results in a corresponding
                                                   ! surplus of nutrients. This surplus (pee) is routed to nutrients
    !
    ! Check balance: (should be zero)
    !
    !write(*,*) 'Copepod N balance:', &
    !      + sum(this%JF/this%m*u)/rhoCN &  ! Gains from feeding
    !      - sum(dudt)/rhoCN & ! Accumulation of biomass
    !      - sum( this%Jresp*u/(this%m*rhoCN) ) & ! Losses from respiration
    !      - (1-epsilonR)*this%g(this%n)*u(this%n)/rhoCN  & ! Losses from reproduction
    !      - sum(this%mort*u)/rhoCN  ! Mortality losses

  end subroutine calcDerivativesCopepod

  subroutine printRatesCopepod(this)
   class(spectrumCopepod), intent(in):: this
   character(7) :: type

   if (this%feedingmode .eq. passive) then
      type = 'Passive'
   else
      type = 'Active'
   endif

   write(*,*) type, " copepod with ", this%n, " size classes and adult size ", this%m(this%n), "ugC."
     call this%printRatesSpectrum()

     99 format (a10, 20f10.6)
 
     write(*,99) "gamma:", this%gamma
     write(*,99) "mortStarve:", this%mortStarve
     write(*,99) "g:", this%g
  end subroutine printRatesCopepod

  function getNbalance(this, u, dudt) result(Nbalance)
    class(spectrumCopepod), intent(in):: this
    real(dp):: Nbalance
    real(dp), intent(in):: u(this%n), dudt(this%n)

    Nbalance = sum(dudt)/rhoCN ! HTL losses
  end function getNbalance

  function getCbalance(this, u, dudt) result(Cbalance)
    class(spectrumCopepod), intent(in):: this
    real(dp):: Cbalance
    real(dp), intent(in):: u(this%n), dudt(this%n)

    Cbalance = sum( dudt & ! Change in standing stock of N
          + this%Jresptot*u/this%m )  ! Losses from respiration
   end function getCbalance

end module copepods
