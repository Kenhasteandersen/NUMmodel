!
! Module that handles the general methods for all spectrum groups.
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
     real(dp), dimension(:), allocatable:: flvl(:), AF(:), JFmax(:), JF(:), f(:)
     ! Growth:
     real(dp), dimension(:), allocatable:: Jtot, JCloss_feeding, JNlossLiebig
     real(dp), dimension(:), allocatable:: JNloss, JCloss, Jresp
     real(dp), dimension(:), allocatable:: mPOM, jPOM ! mass and flux of POM created
     ! Mortality:
     real(dp), dimension(:), allocatable:: mortpred, mortHTL, mort2
     real(dp) :: mort2constant
     ! Sinking:
     real(dp), dimension(:), allocatable:: velocity ! sinking velocity m/day
 
     contains 

     procedure, pass :: initSpectrum
     procedure :: calcFeeding
     procedure :: printRates => printRatesSpectrum
     procedure :: printRatesSpectrum
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

    real(dp), dimension(:), allocatable:: Jmax, Jresptot
    real(dp), dimension(:), allocatable:: JNreal, JDOCreal, JSireal


    contains

    procedure, pass :: initUnicellular
    procedure :: printRatesUnicellular
    procedure :: getProdNet
    procedure :: getProdBact => getProdBactUnicellular
  end type spectrumUnicellular
  ! ------------------------------------------------
  ! Abstact class for all multicellular spectra:
  ! 
  type, abstract, extends(typeSpectrum) :: spectrumMulticellular
  contains
    procedure :: initMulticellular
    procedure :: printRatesMulticellular
  end type spectrumMulticellular
  ! ------------------------------------------------
  ! Type needed to make an array of spectra:
  !
  type spectrumContainer 
    class(typeSpectrum), allocatable :: spec
  end type spectrumContainer

contains

! ==========================================================================
!  Member functions for the abstract spectrum parent class:
! ==========================================================================

  subroutine initSpectrum(this, n, mMin, mMax)
    class(typeSpectrum) :: this
    integer, intent(in):: n
    real(dp), intent(in):: mMin, mMax

    this%n = n
    allocate(this%m(n))
    allocate(this%mLower(n))
    allocate(this%mDelta(n))
    allocate(this%z(n))
    call calcGrid()

    allocate(this%AF(n))
    allocate(this%JFmax(n))
    allocate(this%flvl(n))
    allocate(this%JF(n))
    !allocate(this%JFreal(n))

    allocate(this%Jtot(n))
    allocate(this%Jresp(n))
    allocate(this%JCloss_feeding(n))
    allocate(this%JNlossLiebig(n))
    allocate(this%JNloss(n))
    allocate(this%JCloss(n))

    allocate(this%f(n))


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
    this%epsilonF = 1.d0 ! Probably overridden by the specific group, but must be >0.
    this%palatability = 1.d0 ! set to default
    this%mPOM = 0.d0
    this%jPOM = 0.d0
    this%velocity = 0.d0 ! Probably overridden by the specific group (POM at least)
    this%JResp = 0.d0
    this%Jtot = 0.d0
    this%JCloss_feeding = 0.d0
    this%JNlossLiebig = 0.d0
    this%JNloss = 0.d0
    this%jCloss = 0.d0
    this%f = 0.d0
    contains

      !
  ! Set up a grid given minimum and maximum center masses
  !
  subroutine calcGrid()
    !class (typeSpectrum), intent(inout):: this
    !real(dp), intent(in):: mMin, mMax
    real(dp):: deltax, x
    integer:: i

    deltax = (log(mMax)-log(mMin))/(this%n-1)
    do i=1, this%n
       x = log(mMin) + (i-1)*deltax
       this%m(i) = exp(x)
       this%mLower(i) = exp(x - 0.5*deltax)
       this%mDelta(i) =  exp(x + 0.5*deltax)-this%mLower(i)
       this%z(i) = this%mLower(i)/(this%mLower(i) + this%mDelta(i))
    end do
    this%mort2constant = 0.004/log(this%m(2) / this%m(1))
  end subroutine calcGrid

  end subroutine initSpectrum

  subroutine calcFeeding(this, F)
    class (typeSpectrum), intent(inout):: this
    real(dp), intent(in):: F(this%n)

    this%flvl = this%epsilonF * this%AF*F / &
      ((this%AF*F+eps) + fTemp2*this%JFmax)
    this%JF = this%flvl * fTemp2*this%JFmax
  end subroutine calcFeeding

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

    call this%initSpectrum(n, mMin, mMax)

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
    allocate(this%Jresptot(n))
    allocate(this%JNreal(n))
    allocate(this%JDOCreal(n))
    allocate(this%JSireal(n))



    this%mPOM = this%m ! Assume that POM created by dead cells are 
                       !the same size as the cells
  end subroutine initUnicellular

  subroutine printRatesUnicellular(this)
    class (spectrumUnicellular), intent(in):: this 
    
    99 format (a10, 20f10.6)
 
    call this%printRatesSpectrum()

    write(*,'(a10, 20d10.3)') "r:", this%r
    write(*,99) "jNreal:", this%JNreal / this%m
    write(*,99) "jL:", this%JL / this%m
    write(*,99) "jLreal:", this%JLreal / this%m
    write(*,99) "jDOCreal:", this%JDOCreal / this%m
    write(*,99) "jSireal:", this%JSireal / this%m
    write(*,99) "jLossPass.", this%JlossPassive / this%m
  end subroutine printRatesUnicellular

  function getProdNet(this, u) result(ProdNet)
    real(dp):: ProdNet
    class(spectrumUnicellular), intent(in):: this
    real(dp), intent(in):: u(this%n)
    integer:: i
 
    ProdNet = 0.d0
    do i = 1, this%n
       ProdNet = ProdNet + max( 0.d0, &
                   (this%JLreal(i)-ftemp2*this%Jresp(i))*u(i)/this%m(i) )
     end do
    end function getProdNet

  function getProdBactUnicellular(this, u) result(ProdBact)
    real(dp):: ProdBact
    class(spectrumUnicellular), intent(in):: this
    real(dp), intent(in):: u(this%n)

    ProdBact = 0.d0
  end function getProdBactUnicellular

  ! ==========================================================================
  !  Member functions for the abstract unicellular class:
  ! ==========================================================================

  subroutine initMulticellular(this, n, mMin, mMax)
    class(spectrumMulticellular), intent(inout) :: this
    integer, intent(in):: n
    real(dp), intent(in):: mMin, mMax
    
    call this%initSpectrum(n, mMin, mMax)
  end subroutine initMulticellular
  
  subroutine printRatesMulticellular(this)
    class(spectrumMulticellular), intent(in) :: this
  end subroutine printRatesMulticellular

 end module spectrum