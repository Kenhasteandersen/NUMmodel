!
! Module to handle copepods from a specific species. These can have different species-
! specific parameters, and species-dependent temperature functions.
!
module copepodSpecies
  use globals
  use spectrum
  use read_input_module
  use copepods
  implicit none

  private

  type, extends(spectrumCopepod) :: spectrumCopepodSpecies
    real(dp) :: alphaF, q, h, hExponent, AdultOffspring, mAdult
    real(dp) :: epsilonR, kBasal, kSDA

  contains
    procedure, pass :: initCopepodSpecies
    procedure :: calcDerivativesCopepodSpecies
  end type spectrumCopepodSpecies
  
  public spectrumCopepodSpecies, initCopepodSpecies, calcDerivativesCopepodSpecies

contains

  subroutine initCopepodSpecies(this, nameSpecies, n, errorio,errorstr)
    use iso_c_binding, only: c_char
    class(spectrumCopepodSpecies), intent(inout):: this
    character(len=20), intent(in):: nameSpecies
    integer, intent(in):: n
    logical(1), intent(out):: errorio 
    character(c_char), dimension(*), intent(out) :: errorstr
    real(dp):: alphaF,q,h,hExponent
    
    print*, 'Loading parameter for ',nameSpecies,' from ', inputfile, ':'
    !
    ! Calc grid. Grid runs from mLower(1) = offspring size to m(n) = adult size
    !
    call read_input(inputfile,nameSpecies,'AdultMass',this%mAdult,errorio,errorstr)
    call read_input(inputfile,nameSpecies,'AdultOffspring',this%AdultOffspring,errorio,errorstr)
    call this%initMulticellular(n, this%mAdult/this%AdultOffspring, this%mAdult)

    call read_input(inputfile,nameSpecies,'alphaF',alphaF,errorio,errorstr)
    call read_input(inputfile,nameSpecies,'q',q,errorio,errorstr)
    call read_input(inputfile,nameSpecies,'h',h,errorio,errorstr)
    call read_input(inputfile,nameSpecies,'hExponent',hExponent,errorio,errorstr)
    call read_input(inputfile,nameSpecies,'vulnerability',this%palatability,errorio,errorstr)
    
    call read_input(inputfile,nameSpecies,'epsilonR',this%epsilonR,errorio,errorstr)
    call read_input(inputfile,nameSpecies,'kBasal',this%kBasal,errorio,errorstr)
    call read_input(inputfile,nameSpecies,'kSDA',this%kSDA,errorio,errorstr)
    call read_input(inputfile,nameSpecies,'DiatomsPreference',this%DiatomsPreference,errorio,errorstr)
    
    call read_input(inputfile,nameSpecies,'epsilonF',this%epsilonF,errorio,errorstr)
    call read_input(inputfile,nameSpecies,'beta',this%beta,errorio,errorstr)
    call read_input(inputfile,nameSpecies,'sigma',this%sigma,errorio,errorstr)

    allocate(this%gamma(n))
    allocate(this%g(n))
    allocate(this%mortStarve(n))
    allocate(this%mort(n))
    allocate(this%JrespFactor(n))

    this%AF = alphaF*this%m**q
    this%JFmax = h*this%m**hExponent
    this%JrespFactor = this%epsilonF*this%JFmax
    this%mort2constant = 0.d0 ! No quadratic mortality
    this%mort2 = 0.d0

    this%mPOM = 3.5e-3*this%m ! Size of fecal pellets (Serra-Pompei 2022 approximated)
  end subroutine initCopepodSpecies

  subroutine calcDerivativesCopepodSpecies(this, u, T, dNdt, dudt)
    class(spectrumCopepodSpecies), intent(inout):: this
    real(dp), intent(in):: u(this%n)
    real(dp), intent(in):: T
    real(dp), intent(inout):: dNdt, dudt(this%n)
    integer:: i
    real(dp):: nu, b

    do i = 1, this%n
       !
       ! Growth and reproduction:
       !

       ! Basal and SDA respiration:
       this%Jresptot(i) = this%JrespFactor(i) * this%kBasal * fTemp2 + this%kSDA * this%JF(i)
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
    b = this%epsilonR * this%g(this%n) ! Birth rate
    ! Production of POM:
    this%jPOM = &
          (1-this%epsilonF)*this%JF/(this%m * this%epsilonF) !& ! Unassimilated food (fecal pellets)
       ! + this%mortStarve                            ! Copepods dead from starvation are not counted here, because
                                                      ! the starvation is already respired
    this%jPOM(this%n) = this%jPOM(this%n) + (1.d0-this%epsilonR)*this%g(this%n) ! Lost reproductive flux
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

  end subroutine calcDerivativesCopepodSpecies

end module copepodSpecies
