!
! Module to handle zooplankton defined by their adult mass.
! All parameters are for active zooplanktons. Follows Serra-Pompei et al (2020) with
! the update to respiration in Serra-Pompei (2022).
!
module zooplanktons
  use globals
  use spectrum
  use read_input_module
  implicit none

  integer, parameter :: gelatinous = 1
  integer, parameter :: nongelatinous = 2
  real(dp) :: epsilonR, kBasal, kSDA
  real(dp) :: DiatomsPreference

  private

  type, extends(spectrumMulticellular) :: spectrumZooplankton
    integer :: feedingmode ! Nongelatinous=2; Gelatinous=1
    real(dp), allocatable :: gamma(:), g(:), mortStarve(:), mort(:), JrespFactor(:)
    real(dp) :: DiatomsPreference ! Feeding preference on diatoms
  contains
    procedure, pass :: initZooplankton
    procedure :: calcDerivativesZooplankton
    procedure :: printRates => printRatesZooplankton
  end type spectrumZooplankton
  
  public nongelatinous, gelatinous, spectrumZooplankton, initZooplankton, calcDerivativesZooplankton, printRatesZooplankton

contains

  subroutine initZooplankton(this, feedingmode, n, mAdult,errorio,errorstr)
    use iso_c_binding, only: c_char
    class(spectrumZooplankton), intent(inout):: this
    integer, intent(in) :: feedingmode ! Whether the zooplankton is nongelatinous or gelatinous
    integer, intent(in):: n
    real(dp), intent(in):: mAdult
    logical(1), intent(out):: errorio 
    character(c_char), dimension(*), intent(out) :: errorstr
    real(dp) :: alphaF, q, h, hExponent, AdultOffspring
    real(dp) :: vulnerability
    
    character(len=30)::this_listname
    
    this%feedingmode = feedingmode
    
    if (feedingmode .eq. nongelatinous) then
      this_listname='zooplanktons_nongelatinous'
    else if (feedingmode .eq. gelatinous) then
      this_listname='zooplanktons_gelatinous'
    else
       print*, 'no feeding mode defined'
       stop
    end if
    
    print*, 'Loading parameter for ',this_listname,' from ', inputfile, ':'

    !
    ! Calc grid. Grid runs from mLower(1) = offspring size to m(n) = adult size
    !
    call read_input(inputfile,this_listname,'AdultOffspring',AdultOffspring,errorio,errorstr)
    call this%initMulticellular(n, mAdult/AdultOffspring, mAdult)

    call read_input(inputfile,this_listname,'alphaF',alphaF,errorio,errorstr)
    call read_input(inputfile,this_listname,'q',q,errorio,errorstr)
    call read_input(inputfile,this_listname,'h',h,errorio,errorstr)
    call read_input(inputfile,this_listname,'hExponent',hExponent,errorio,errorstr)
    call read_input(inputfile,this_listname,'vulnerability',vulnerability,errorio,errorstr)
    
    call read_input(inputfile,this_listname,'epsilonR',epsilonR,errorio,errorstr)
    call read_input(inputfile,this_listname,'kBasal',kBasal,errorio,errorstr)
    call read_input(inputfile,this_listname,'kSDA',kSDA,errorio,errorstr)
    call read_input(inputfile,this_listname,'DiatomsPreference',DiatomsPreference,errorio,errorstr)
    
    call read_input(inputfile,this_listname,'epsilonF',this%epsilonF,errorio,errorstr)
    call read_input(inputfile,this_listname,'beta',this%beta,errorio,errorstr)
    call read_input(inputfile,this_listname,'sigma',this%sigma,errorio,errorstr)
    this%DiatomsPreference=DiatomsPreference
    

    allocate(this%gamma(n))
    allocate(this%g(n))
    allocate(this%mortStarve(n))
    allocate(this%mort(n))
    allocate(this%JrespFactor(n))

    this%AF = alphaF*this%m**q
    this%JFmax = h*this%m**hExponent
    this%JrespFactor = this%JFmax/this%JFmax
    this%mort2constant = 0.d0 ! No quadratic mortality
    this%mort2 = 0.d0
    this%palatability = vulnerability

    this%mPOM = 3.5e-3*this%m ! Size of fecal pellets (Serra-Pompei 2022 approximated)
  end subroutine initZooplankton

  subroutine calcDerivativesZooplankton(this, u, dNdt, dudt)
    class(spectrumZooplankton), intent(inout):: this
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
       !this%Jresptot(i) = this%Jresptot(i) - min(0.d0, nu) ! Limit respiration to the energy available
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
          (1-this%epsilonF)*this%JF/(this%m * this%epsilonF) !& ! Unassimilated food (fecal pellets)
       ! + this%mortStarve                            ! Zooplanktons dead from starvation are not counted here, because
                                                      ! the starvation is already respired
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
    !write(*,*) 'Zooplankton N balance:', &
    !      + sum(this%JF/this%m*u)/rhoCN &  ! Gains from feeding
    !      - sum(dudt)/rhoCN & ! Accumulation of biomass
    !      - sum( this%Jresp*u/(this%m*rhoCN) ) & ! Losses from respiration
    !      - (1-epsilonR)*this%g(this%n)*u(this%n)/rhoCN  & ! Losses from reproduction
    !      - sum(this%mort*u)/rhoCN  ! Mortality losses

  end subroutine calcDerivativesZooplankton

  subroutine printRatesZooplankton(this)
   class(spectrumZooplankton), intent(in):: this
   character(7) :: type

   if (this%feedingmode .eq. gelatinous) then
      type = 'Gelatinous'
   else
      type = 'Nongelatinous'
   endif

   write(*,*) type, " zooplankton with ", this%n, " size classes and adult size ", this%m(this%n), "ugC."
     call this%printRatesSpectrum()

     99 format (a10, 20f10.6)
 
     write(*,99) "gamma:", this%gamma
     write(*,99) "mortStarve:", this%mortStarve
     write(*,99) "g:", this%g
  end subroutine printRatesZooplankton

end module zooplanktons
