!
! Module for handling particular organic matter (POM)
! To have other groups produce POM, do the following:
!  1) Define the masses of POM that each group produces in the vector mPOM 
!     This vector is already set to default as mPOM = m, which works for unicellular groups
!  2) Define the fluxes to POM in jPOM (note this is a rate 1/day)
!
! NOTE: THE DYNAMICS OF POM IN GENERALISTS AND COPEPODS NEEDS TO BE REVISITED
!
module POM
    use globals
    use spectrum
    use read_input_module
    implicit none
  
    private 
  
    type, extends(typeSpectrum) :: spectrumPOM
            
      real(dp):: remin

    contains
      procedure, pass :: initPOM
      procedure :: calcDerivativesPOM
      procedure :: printRates => printRatesPOM
      procedure :: getCbalance
    end type spectrumPOM
   
    public initPOM, spectrumPOM, calcDerivativesPOM, printRatesPOM
  
  contains

  subroutine initPOM(this, n, mMax,errorio,errorstr)
    use iso_c_binding, only: c_char
    class(spectrumPOM):: this
    integer, intent(in):: n
    real(dp), intent(in):: mMax
    logical(1), intent(out):: errorio 
    character(c_char), dimension(*), intent(out) :: errorstr
    integer:: file_unit,io_err
    real(dp) :: mMin
    print*, 'Loading parameter for POM from ', inputfile, ':'
    call read_input(inputfile,'POM','mMin',mMin,errorio,errorstr)
    call read_input(inputfile,'POM','remin',this%remin,errorio,errorstr)

    call this%initSpectrum(n)
    call this%calcGrid(mMin, mMax)

    this%velocity = 400*this%m**0.513 ! Copepod fecal pellets from Serra-Pompei (2022)
    this%mort2 = 0.d0 ! No virulysis of POM
  end subroutine initPOM

  subroutine calcDerivativesPOM(this, u, dNdt, dDOCdt, dudt)
    class(spectrumPOM):: this
    real(dp), intent(in):: u(this%n)
    real(dp), intent(inout) :: dNdt, dDOCdt, dudt(this%n)

    this%Jresptot = fTemp2*this%remin*this%m
    dudt = dudt - fTemp2*this%remin*u - this%mortpred*u
    dNdt = dNdt + sum(fTemp2*this%remin*u)/rhoCN
    dDOCdt = dDOCdt ! remineralized carbon is respired, so lost
  end subroutine calcDerivativesPOM

  function getCbalance(this, u, dudt) result(Cbalance)
    real(dp):: Cbalance
    class(spectrumPOM), intent(in):: this
    real(dp), intent(in):: u(this%n), dudt(this%n)

    Cbalance = sum(dudt + this%Jresptot/this%m*u)
  end function getCbalance

  subroutine printRatesPOM(this)
    class(spectrumPOM), intent(in):: this
  
    write(*,*) "POM with ", this%n, " size classes:"
    call this%printRatesSpectrum()
  end subroutine printRatesPOM
end
