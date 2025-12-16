module gelatinous
use globals
use spectrum
use copepods
use read_input_module
implicit none

private

! The gelatinous zooplankton type inherits everything from copepods
type, extends(spectrumCopepod) :: spectrumGelatinous
contains
  procedure, pass :: initGelatinous
  procedure :: calcDerivativesGelatinous
  procedure :: printRates => printRatesGelatinous
end type spectrumGelatinous

public spectrumGelatinous, initGelatinous, calcDerivativesGelatinous, printRatesGelatinous

contains

subroutine initGelatinous(this, n, mAdult, errorio,errorstr)
  use iso_c_binding, only: c_char
  class(spectrumGelatinous), intent(inout):: this
  integer, intent(in):: n
  real(dp), intent(in):: mAdult
  logical(1), intent(out):: errorio 
  character(c_char), dimension(*), intent(out) :: errorstr
character(len=20)::this_listname
  ! Initialize as a copepod with gelatinous feeding mode
  ! Gelatinous zooplankton are considered active feeders
  this%feedingmode = active
  this_listname = 'gelatinous'
  
  call this%readCopepodInput(this_listname, n, mAdult, errorio, errorstr)
  call read_input(inputfile,this_listname,'selectionHTL',this%selectionHTL,errorio,errorstr)

  this%mPOM = 3.5e-3*this%m ! WHAT IS THE SIZE OF FECAL PELLETS FROM GELATINOUS ZOOPLANKTON?
end subroutine initGelatinous

subroutine calcDerivativesGelatinous(this, u, dNdt, dudt)
  class(spectrumGelatinous), intent(inout):: this
  real(dp), intent(in):: u(this%n)
  real(dp), intent(inout):: dNdt, dudt(this%n)

  call calcDerivativesCopepod(this, u, dNdt, dudt)

end subroutine calcDerivativesGelatinous

subroutine printRatesGelatinous(this)
  class(spectrumGelatinous), intent(in):: this

  write(*,*) "Gelatinous zooplankton with ", this%n, " size classes and adult size ", this%m(this%n), "ugC."
  call this%printRatesSpectrum()

  99 format (a10, 20f10.6)
 
  write(*,99) "gamma:", this%gamma
  write(*,99) "mortStarve:", this%mortStarve
  write(*,99) "g:", this%g
end subroutine printRatesGelatinous

end module gelatinous 