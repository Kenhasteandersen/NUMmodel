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
    implicit none
  
    private 
  
    real(dp), parameter:: rhoCN = 5.68
    real(dp), parameter:: remin = 0.5d0 ! remineralisation rate (1/day)
    real(dp), parameter:: mMin = 1e-9 ! Smallest POM mass
  
    type, extends(typeSpectrum) :: spectrumPOM
            
    contains
      procedure, pass :: initPOM
      procedure :: calcDerivativesPOM
      procedure :: printRates => printRatesPOM
    end type spectrumPOM
   
    public initPOM, spectrumPOM, calcDerivativesPOM, printRatesPOM
  
  contains

  subroutine initPOM(this, n, mMax)
    class(spectrumPOM):: this
    integer, intent(in):: n
    real(dp), intent(in):: mMax
    
    call this%initSpectrum(n, mMin, mMax)

    this%velocity = 10.d0 ! Size-independent fast sinking (10 m/day)
    this%mort2 = 0.d0 ! No virulysis of POM
  end subroutine initPOM

  subroutine calcDerivativesPOM(this, u, dNdt, dDOCdt, dudt)
    class(spectrumPOM):: this
    real(dp), intent(in):: u(this%n)
    real(dp), intent(inout) :: dNdt, dDOCdt, dudt(this%n)

    dudt = dudt - remin*u - this%mortpred*u
    dNdt = dNdt + sum(remin*u/rhoCN)
    dDOCdt = dDOCdt + sum(remin*u)
  end subroutine calcDerivativesPOM

  subroutine printRatesPOM(this)
    class(spectrumPOM), intent(in):: this
  
    write(*,*) "POM with ", this%n, " size classes:"
    call this%printRatesSpectrum()
  end subroutine printRatesPOM
end
