!
! Module that handles the general methods for all spectrum groups.
! It defines the two abstract classes for unicellullar and multicellular plankton.
!
module spectrum
  use globals
  implicit none
  ! ------------------------------------------------
  ! Abstract class for all spectra:
  !
  type, abstract :: typeSpectrum
     integer:: type
     integer:: n  ! Number of size classes

     ! Grid:
     real(dp), dimension(:), allocatable:: m(:)  ! Geometric center mass of size-bin
     real(dp), dimension(:), allocatable:: mLower(:)  ! Smallest size in the bin
     real(dp), dimension(:), allocatable:: mDelta(:)   ! Width of the bin
     real(dp), dimension(:), allocatable:: z(:) ! Ratio btw upper and lower size of bin
     ! Feeding:
     real(dp):: palatability ! [0:1] Reduction of risk of predation
     real(dp):: beta, sigma ! Pred:prey mass ratio and width
     real(dp):: epsilonF ! Assimilation efficiency
     real(dp):: epsilonL ! Light Assimilation efficiency
     real(dp), dimension(:), allocatable:: flvl(:), AF(:), JFmax(:), JF(:), f(:)
     ! Growth:
     real(dp), dimension(:), allocatable:: Jtot, JCloss_feeding, JNlossLiebig
     real(dp), dimension(:), allocatable:: JNloss, JCloss, Jresp, Jresptot
     real(dp), dimension(:), allocatable:: mPOM, jPOM! mass and flux of POM created
     ! Mortality:
     real(dp), dimension(:), allocatable:: mortpred, mortHTL, mort2
     real(dp) :: mort2constant
     ! Sinking:
     real(dp), dimension(:), allocatable:: velocity ! sinking velocity m/day
 
     contains 

     procedure, pass :: initSpectrum
     procedure :: calcGrid
     procedure :: calcFeeding
     procedure :: printRates => printRatesSpectrum
     procedure :: printRatesSpectrum
     procedure :: getClost => getClostSpectrum
  end type typeSpectrum

  ! ------------------------------------------------
  ! Abstract class for all unicellular spectra:
  !
  type, abstract, extends(typeSpectrum) :: spectrumUnicellular
    real(dp),  dimension(:), allocatable:: nu, r
 
    real(dp),  dimension(:), allocatable:: JlossPassive
    ! Resource uptake affinities
    real(dp), dimension(:), allocatable:: AN, AL
    ! Resource uptake fluxes
    real(dp), dimension(:), allocatable:: JN, JDOC, JL
    real(dp), dimension(:), allocatable:: JNtot, JLreal, JCtot 
    real(dp), dimension(:), allocatable:: JCloss_photouptake, JClossLiebig

    real(dp), dimension(:), allocatable:: Jmax
    real(dp), dimension(:), allocatable:: JNreal, JDOCreal

    contains

    procedure, pass :: initUnicellular
    procedure :: printRatesUnicellular
    procedure :: getCLost => getClostUnicellular
    procedure :: getProdNet
    procedure :: getProdBact => getProdBactUnicellular
  end type spectrumUnicellular

  ! ------------------------------------------------
  ! Abstact class for all multicellular spectra:
  ! 
  type, abstract, extends(typeSpectrum) :: spectrumMulticellular
  contains
    procedure :: initMulticellular
    !procedure :: printRatesMulticellular
    !procedure :: getCbalance => getCbalanceMulticellular
  end type spectrumMulticellular
  ! ------------------------------------------------
  ! Type needed to make an array of spectra:
  !
  type spectrumContainer 
    class(typeSpectrum), allocatable :: spec
  end type spectrumContainer

contains

! ==========================================================================
!  Member functions for the abstract spectrum parent class.
! ==========================================================================
  subroutine initSpectrum(this, n)
    class(typeSpectrum) :: this
    integer, intent(in):: n

    this%n = n
    allocate(this%m(n))
    allocate(this%mLower(n))
    allocate(this%mDelta(n))
    allocate(this%z(n))

    allocate(this%AF(n))
    allocate(this%JFmax(n))
    allocate(this%flvl(n))
    allocate(this%JF(n))
    allocate(this%f(n))

    !allocate(this%JFreal(n))

    allocate(this%Jtot(n))
    allocate(this%Jresp(n))
    allocate(this%Jresptot(n))
    allocate(this%JCloss_feeding(n))
    allocate(this%JNlossLiebig(n))
    allocate(this%JNloss(n))
    allocate(this%JCloss(n))

    allocate(this%mPOM(n))
    allocate(this%jPOM(n))

    allocate(this%mortpred(n))
    allocate(this%mortHTL(n))
    allocate(this%mort2(n))

    allocate(this%velocity(n))
    ! Set feeding to dummy values. Relevant for non-feeding groups (diatoms)
    this%AF = 0.d0
    this%JFmax = 0.d0
    this%flvl = 0.d0
    this%JF = 0.d0
    this%f = 0.d0
    this%epsilonF = 1.d0 ! Probably overridden by the specific group, but must be >0.
    this%palatability = 1.d0 ! set to default
    this%mPOM = 0.d0
    this%jPOM = 0.d0
    this%velocity = 0.d0 ! Probably overridden by the specific group (POM at least)
    this%JResp = 0.d0
    this%Jresptot = 0.d0
    this%Jtot = 0.d0
    this%JCloss_feeding = 0.d0
    this%JNlossLiebig = 0.d0
    this%JNloss = 0.d0
    this%jCloss = 0.d0
  end subroutine initSpectrum

  !
  ! Set up a grid for the size groups
  ! The minimum size corresponds to the lower size of the 
  ! first grid cell and the maxmimum size corresponds to
  ! the upper size of the last grid cell
  !
  subroutine calcGrid(this, mMin, mMax)
    class(typeSpectrum) :: this
    real(dp), intent(in):: mMin, mMax
    integer:: i
    real(dp):: x, deltax
    
  deltax = (log(mMax)-log(mMin)) / this%n
  do i=1, this%n
     x = log(mMin) + (i-0.5)*deltax
     this%m(i) = exp(x)
     this%mLower(i) = exp(x - 0.5*deltax)
     this%mDelta(i) =  exp(x + 0.5*deltax)-this%mLower(i)
     this%z(i) = this%mLower(i)/(this%mLower(i) + this%mDelta(i))
  end do
  this%mort2constant = 0.004/log(this%m(2) / this%m(1))
end subroutine calcGrid

!function getNbalanceSpectrum(this, u, dudt) result(Nbalance)
!    real(dp):: Nbalance
!    class(typeSpectrum), intent(in):: this
!    real(dp), intent(in):: u(this%n), dudt(this%n)

!    Nbalance = sum( dudt ) / rhoCN
!  end function getNbalanceSpectrum

  subroutine getLossesSpectrum(this, u, Nloss, Closs, SiLoss)
    class(typeSpectrum), intent(in):: this
    real(dp), intent(in):: u(this%n)
    real(dp), intent(out) :: Nloss, Closs, SiLoss

    Nloss = 0.d0
    Closs = sum( this%Jresptot * u / this%m )
    SiLoss = 0.d0
  end subroutine getLossesSpectrum
  !
  ! Returns the amount of encounter and potentially assimilated food available for a group, JF
  !
  subroutine calcFeeding(this, F)
    class (typeSpectrum), intent(inout):: this
    real(dp), intent(in):: F(this%n)

    this%flvl = this%epsilonF * this%AF*F / & ! Note: adding a small number in the
      ((this%AF*F+eps) + fTemp2*this%JFmax)   ! demonominator to avoid negative values if F = JFmax = 0.
    this%JF = this%flvl * fTemp2*this%JFmax
  end subroutine calcFeeding

  !
  ! Returns the carbon that is lost from the system (by default only respiration, but 
  !  spectrumUnicellular overrides this function to also account for imports of carbon by photosynthesis)
  !
  function getClostSpectrum(this, u) result(Clost)
    class (typeSpectrum), intent(in):: this
    real(dp), intent(in):: u(this%n)
    real(dp) :: Clost

    Clost = sum( this%Jresptot*u/this%m )
  end function getClostSpectrum

  subroutine printRatesSpectrum(this)
    class (typeSpectrum), intent(in):: this

    99 format (a10, 20f10.6)
    write(*,'(a10, 20d10.3)') "m: ", this%m
    write(*,99) "jF:", this%JF / this%m
    write(*,99) "jTot:", this%Jtot / this%m
    write(*,99) "mortpred", this%mortpred
    write(*,99) "mortHTL", this%mortHTL
    write(*,99) "jPOM:", this%jPOM
  end subroutine printRatesSpectrum

  ! ==========================================================================
  !  Member functions for the abstract unicellular class:
  ! ==========================================================================

  subroutine initUnicellular(this, n, mMin, mMax)
    class(spectrumUnicellular), intent(inout) :: this
    integer, intent(in):: n
    real(dp), intent(in):: mMin, mMax

    call this%initSpectrum(n)
    call this%calcGrid(mMin, mMax)
  
    allocate(this%nu(n))
    allocate(this%r(n))

    allocate(this%JlossPassive(n))

    allocate(this%AN(n))
    allocate(this%AL(n))

    allocate(this%JN(n))
    allocate(this%JDOC(n))
    allocate(this%JL(n))
    allocate(this%JNtot(n))
    allocate(this%JLreal(n))

    allocate(this%JCtot (n))
    allocate(this%JCloss_photouptake(n))
    allocate(this%JClossLiebig(n))

    allocate(this%Jmax(n))
    allocate(this%JNreal(n))
    allocate(this%JDOCreal(n))

    this%mPOM = this%m ! Assume that POM created by dead cells are 
                       ! the same size as the cells
  end subroutine initUnicellular
  !
  ! Carbon lost from the system (with gains being negative):
  !
  function getClostUnicellular(this, u) result(Clost)
    class (spectrumUnicellular), intent(in):: this
    real(dp), intent(in):: u(this%n)
    real(dp) :: Clost

    Clost = sum( ( &
      - this%JLreal & ! Fixed carbon
      - this%JCloss_photouptake & !Fixed carbon which is later exuded
      + this%Jresptot & ! Respiration
      )/this%m * u )
  end function getClostUnicellular

  subroutine printRatesUnicellular(this)
    class (spectrumUnicellular), intent(in):: this 
    
    99 format (a10, 20f10.6)
 
    call this%printRatesSpectrum()

    write(*,'(a10, 20d10.3)') "r:", this%r
    write(*,99) "jN:", this%JN / this%m
    write(*,99) "jNreal:", this%JNreal / this%m
    write(*,99) "jL:", this%JL / this%m
    write(*,99) "jLreal:", this%JLreal / this%m
    write(*,99) "jDOC:", this%JDOCreal / this%m
    write(*,99) "jDOCreal:", this%JDOCreal / this%m
    write(*,99) "jLossPass.", this%JlossPassive / this%m
  end subroutine printRatesUnicellular

  !function getCbalanceUnicellular(this, u, dudt) result(Cbalance)
  !  real(dp):: Cbalance
  !  class(spectrumUnicellular), intent(in):: this
  !  real(dp), intent(in):: u(this%n), dudt(this%n)

  !  Cbalance = sum(dudt &
  !  - this%JLreal*u/this%m &
  !  - this%JCloss_photouptake*u/this%m &
  !  + this%Jresptot*u/this%m &
  !  )
  !end function getCbalanceUnicellular
  !
  ! Returns the net primary production calculated as the total amount of carbon fixed
  ! by photsynthesis minus the respiration. Units: mugC/day/m3
  ! (See Andersen and Visser (2023) table 5)
  !
  function getProdNet(this, u) result(ProdNet)
    real(dp):: ProdNet
    class(spectrumUnicellular), intent(in):: this
    real(dp), intent(in):: u(this%n)
    integer:: i
 
    ProdNet = 0.d0
    do i = 1, this%n
       ProdNet = ProdNet + max( 0.d0, &
                   (this%JLreal(i)-this%Jresptot(i))*u(i)/this%m(i) )
    end do
  end function getProdNet
  !
  ! Returns the net bacterial production calculated as the total amount of DOC
  ! taken up minus the respiration. Units: mugC/day/m3
  ! (See Andersen and Visser (2023) table 5)
  !
  function getProdBactUnicellular(this, u) result(ProdBact)
    real(dp):: ProdBact
    class(spectrumUnicellular), intent(in):: this
    real(dp), intent(in):: u(this%n)
    integer:: i

    ProdBact = 0.d0
    do i = 1, this%n
      ProdBact = ProdBact + max( 0.d0, &
                  (this%JDOCreal(i)-this%Jresptot(i))*u(i)/this%m(i) )
    end do
  end function getProdBactUnicellular

  ! ==========================================================================
  !  Member functions for the abstract multicellular class:
  ! ==========================================================================

  subroutine initMulticellular(this, n, mOffspring, mAdult)
    class(spectrumMulticellular), intent(inout) :: this
    integer, intent(in):: n
    real(dp), intent(in):: mOffspring, mAdult
    !integer:: i
    real(dp):: lnDelta
    
    call this%initSpectrum(n)
    !
    ! Set up a grid.
    ! The minimum size corresponds to the offspring size.
    ! The central size of the last cell is the mass of the adult.
    ! The (log) width of the last cell is the same as the other cells;
    ! this is needed to ensure that the mort2 is the same for all
    ! cells.
    !
    lnDelta = (log(mAdult)-log(mOffspring)) / (n-0.5)
    call this%calcGrid(mOffspring, exp( log(mAdult) + 0.5*lnDelta) )
 
    this%mort2constant = 0.004/log(this%m(2) / this%m(1))
  end subroutine initMulticellular
  
  !function getCbalanceMulticellular(this, u, dudt) result(Cbalance)
  !  class(spectrumMulticellular), intent(in):: this
  !  real(dp):: Cbalance
  !  real(dp), intent(in):: u(this%n), dudt(this%n)

  !  Cbalance = sum( dudt & ! Change in standing stock of N
  !        + this%Jresptot*u/this%m )  ! Losses from respiration
  !end function getCbalanceMulticellular

  !subroutine printRatesMulticellular(this)
  !  class(spectrumMulticellular), intent(in) :: this
  !end subroutine printRatesMulticellular

 end module spectrum
