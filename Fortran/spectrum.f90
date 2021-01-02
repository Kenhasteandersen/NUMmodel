module sizespectrum
  use globals
  implicit none
  
  type typeSpectrum
     integer:: typeGroup
     integer:: ixStart, ixEnd, n
     ! Grid:
     real(dp), dimension(:), allocatable:: m(:), mLower(:), mDelta(:), z(:)
     ! Feeding:
     real(dp), dimension(:), allocatable:: AF(:), JFmax(:)
     
  end type typeSpectrum
  
  public initSpectrum
contains
  
  function initSpectrum(n, ixStart, mMin, mMax) result(this)
    type (typeSpectrum):: this
    integer, intent(in):: n, ixStart
    real(dp), intent(in):: mMin, mMax

    this%n = n
    this%ixStart = ixStart
    this%ixEnd = ixStart + n - 1
    allocate(this%m(n))
    allocate(this%mLower(n))
    allocate(this%mDelta(n))
    allocate(this%z(n))
    call calcGrid(this, mMin, mMax)

    allocate(this%AF(n))
    allocate(this%JFmax(n))
  end function initSpectrum

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


 end module sizespectrum
  
