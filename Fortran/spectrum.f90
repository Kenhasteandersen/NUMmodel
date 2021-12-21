!
! Module that handles the general methods for all spectrum groups.
!
module spectrum
  use globals
  implicit none

  type typeSpectrum
     integer:: type
     integer:: n  ! Number of size classes
     integer:: ixStart, ixEnd ! Start and end indices into u and rates
     integer:: ixOffset ! ixStart-1
     ! Grid:
     real(dp), dimension(:), allocatable:: m(:)  ! Geometric center mass of size-bin
     real(dp), dimension(:), allocatable:: mLower(:)  ! Smallest size in the bin
     real(dp), dimension(:), allocatable:: mDelta(:)   ! Width of the bin
     real(dp), dimension(:), allocatable:: z(:) ! Ratio btw upper and lower size of bin
     ! Feeding:
     real(dp):: palatability
     real(dp):: beta, sigma, epsilonF
     real(dp), dimension(:), allocatable:: flvl(:), AF(:), JFmax(:), JF(:)
     real(dp), dimension(:), allocatable:: mortpred
  end type typeSpectrum

  public initSpectrum, calcFeeding
contains

  function initSpectrum(type, n, ixOffset, mMin, mMax) result(this)
    type (typeSpectrum):: this
    integer, intent(in):: type, n, ixOffset
    real(dp), intent(in):: mMin, mMax

    this%type = type
    this%n = n
    this%ixOffset = ixOffset
    this%ixStart = ixOffset+1
    this%ixEnd = ixOffset+n
    allocate(this%m(n))
    allocate(this%mLower(n))
    allocate(this%mDelta(n))
    allocate(this%z(n))
    call calcGrid(this, mMin, mMax)

    allocate(this%AF(n))
    allocate(this%JFmax(n))
    allocate(this%flvl(n))
    allocate(this%JF(n))

    allocate(this%mortpred(n))
    ! Set feeding to dummy values. Relevant for non-feeding groups (diatoms)
    this%AF = 0.d0
    this%JFmax = 1.d0
    this%flvl = 0.d0
    this%JF = 0.d0
    this%palatability = 1.d0 ! set to default
  end function initSpectrum
  !
  ! Set up a grid given minimum and maximum center masses
  !
  subroutine calcGrid(this, mMin, mMax)
    type (typeSpectrum), intent(inout):: this
    real(dp), intent(in):: mMin, mMax
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
  end subroutine calcGrid

  subroutine calcFeeding(this, F)
    type (typeSpectrum), intent(inout):: this
    real(dp), intent(in):: F(this%n)

    this%flvl = this%epsilonF * this%AF*F / &
      (this%AF*F + fTemp2*this%JFmax)
    this%JF = this%flvl * fTemp2*this%JFmax
  end subroutine

 end module spectrum
