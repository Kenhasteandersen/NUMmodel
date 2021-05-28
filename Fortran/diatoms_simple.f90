!
! Module to handle diatoms
! Maximal simplicity; based on generalists with the only difference:
!  - a fixed-sized vacuole which increases radius
!  - lack of ability of phagotrophy
!  - rely on silicate
!  - lower predation risk due to silicate shell
!
module diatoms_simple
    use globals
    use spectrum
    implicit none

    private
    !
    ! Stoichiometry:
    !
    real(dp), parameter:: rhoCN = 5.68
    real(dp), parameter:: rhoCSi = 3.4 
     !cell properties
    ! real(dp), parameter:: cm = 6d-7 ! rhoCmem 6*d.-7 ! Carbon content in cell membrane mugC(mum^-3)
    ! real(dp), parameter:: cb = 1.25d-7 ! rhoCcyt 1.25*d.-7 ! Carbon content in cell cytoplasm mugC(mum^-3)  
     real(dp), parameter:: v=0.6 ! Vacuole fraction
    !
    ! Light uptake:
    !
    real(dp), parameter:: epsilonL = 0.9 ! Light uptake efficiency
    real(dp), parameter:: alphaL = 0.206
    real(dp), parameter:: rLstar = 8.25
  !
    ! Costs
    !
    real(dp), parameter :: bL = 0 !0.08 ! cost of light harvesting mugC(mugC)^-1
    real(dp), parameter :: bN = 0 !0.45 ! cost of N uptake mugC(mugN)^-1
    real(dp), parameter :: bSi = 0 !0.45 ! cost of Si uptake mugC(mugSi)^-1
    !
    ! Dissolved nutrient uptake:
    !
    real(dp), parameter:: alphaN = 0.682 / (1-v) ! L/d/mugC/mum^2
    real(dp), parameter:: rNstar = 2 ! mum
    !
    ! Metabolism
    !
    real(dp), parameter:: cLeakage = 0.03 ! passive leakage of C and N
    real(dp), parameter:: c = 0.0015 ! Parameter for cell wall fraction of mass.
    real(dp), parameter:: delta = 0.05 ! Thickness of cell wall in mum
    real(dp), parameter:: alphaJ = 1.5 ! Constant for jmax.  per day
    real(dp), parameter:: cR = 0.1
    !
    ! Bio-geo:
    !
    real(dp), parameter:: remin = 0.0 ! fraction of mortality losses reminerilized to N and DOC
    real(dp), parameter:: remin2 = 1.d0 ! fraction of virulysis remineralized to N and DOC
    real(dp), parameter:: reminHTL = 0.d0 ! fraction of HTL mortality remineralized
    !
    ! Needed internal variables:
    !
    real(dp), dimension(:), allocatable:: AN(:), AL(:), Jmax(:), JlossPassive(:)
    real(dp), dimension(:), allocatable:: nu(:), mort(:)
    real(dp), dimension(:), allocatable:: JN(:), JL(:), Jresp(:)

    real(dp):: mort2

    public initDiatoms_simple, calcRatesDiatoms_simple, calcDerivativesDiatoms_simple
  contains
      
    function initDiatoms_simple(n, ixOffset, mMax) result(this)
      type(typeSpectrum):: this
      real(dp), intent(in):: mMax
      integer, intent(in):: n, ixOffset
      real(dp), parameter:: mMin = 3.1623d-9
      real(dp):: r(n)
      real(dp):: hs(n) ! Shell thickness mum ! TO BE FIXED!!
      real(dp), parameter:: rho = 0.57*1d6*1d-12
      real(dp) :: fl !investment in photosynthesis
  
      this = initSpectrum(typeDiatom_simple, n, ixOffset, mMin, mMax)
  
      if ( allocated(AN) ) then
         deallocate(AN)
         deallocate(AL)
         deallocate(Jmax)
         deallocate(Jresp)
         deallocate(JlossPassive)
         deallocate(nu)
         deallocate(mort)
       end if
      
      allocate(AN(n)) 
      allocate(AL(n))
      allocate(Jmax(n))
      allocate(Jresp(n))
      allocate(JlossPassive(n))
      allocate(nu(n))
      allocate(mort(n))
      
      ! Radius:
      !
      r = (threequarters/pi * this%m/rho/(1-v))**onethird  ! Andy's approximation

      write(*,*) r

      AN = alphaN * r**(-2.) / (1.+(r/rNstar)**(-2.)) * this%m
      AL = alphaL/r * (1-exp(-r/rLstar)) * this%m

      JlossPassive = cLeakage/r * this%m ! in units of C
  
      !nu = c * this%m**(-onethird)
      nu = 3*delta/r
      Jmax = alphaJ * this%m * (1.d0-nu) ! mugC/day
      Jresp = cR*alphaJ*this%m
      mort = 0*0.005*(Jmax/this%m) * this%m**(-0.25);
      mort2 = 0.0002*n
  
      this%beta = 0.d0 ! No feeding

    end function initDiatoms_simple
 
    subroutine calcRatesDiatoms_simple(this, rates, L, N, DOC, Si, gammaN, gammaDOC, gammaSi)
      type(typeSpectrum), intent(in):: this
      real(dp), intent(in):: gammaN, gammaDOC, gammaSi
      type(typeRates), intent(inout):: rates
      real(dp), intent(in):: L, N, DOC, Si
      integer:: ix, i
  
      do i = 1, this%n
         ix = i+this%ixOffset
         !
         ! Uptakes
         !
         rates%JN(ix) =   gammaN * AN(i)*N*rhoCN ! Diffusive nutrient uptake in units of C/time
         rates%JDOC(ix) = gammaDOC * AN(i)*DOC ! Diffusive DOC uptake, units of C/time
         rates%JSi(ix) = gammaSi * AN(i)*Si*rhoCSi! Diffusive Si uptake, units of C/time
         rates%JL(ix) =   epsilonL * AL(i)*L  ! Photoharvesting
              !
         ! Estimate the limiting growth rate:
         !
         rates%Jtot(ix) = min( rates%JN(ix), rates%JL(ix)-Jresp(i), rates%JSi(ix) )
         !    
         ! Account for possible carbon limitation
         !
         rates%Jtot(ix) = min( rates%Jtot(ix), &
            rates%JL(ix)-Jresp(i) - (bN/rhoCN + bSi/rhoCSi)*rates%Jtot(ix) )
         
         ! Jtot saturation
         rates%Jtot(ix)=Jmax(i)*rates%Jtot(ix) / ( rates%Jtot(ix) + Jmax(i) )
      end do
    end subroutine calcRatesDiatoms_simple
  
    subroutine calcDerivativesDiatoms_simple(this, u, rates)
      type(typeSpectrum), intent(in):: this
      type(typeRates), intent(inout):: rates
      real(dp), intent(in):: u(this%n)
      real(dp):: mortloss
      integer:: i, ix
  
      do i = 1, this%n
        ix = i+this%ixOffset
        mortloss = u(i)*(remin2*mort2*u(i) +reminHTL* rates%mortHTL(ix))
        !
        ! Update nitrogen:
        !
        rates%dudt(idxN) = rates%dudt(idxN)  &
            - ((rates%Jtot(ix))*u(i)/this%m(i) &
            + mortloss)/rhoCN
        !
        ! Update DOC:
        !
        rates%dudt(idxDOC) = 0.d0 !rates%dudt(idxDOC) &
             !+ (-rates%JDOC(ix) &
             !+ rates%JCloss(ix))*u(i)/this%m(i) &
             !+ mortloss
        !
        ! Update Si:
        !
        rates%dudt(idxSi) = rates%dudt(idxSi) &
             + ((-rates%Jtot(ix))*u(i)/this%m(i) &
             + mortloss)/rhoCSi
        !
        ! Update the diatoms:
        !
        rates%dudt(ix) = (rates%Jtot(ix)/this%m(i)  &
             - mort(i) &
             - rates%mortpred(ix) &
             - mort2*u(i) &
             - rates%mortHTL(ix))*u(i)

     end do
   end subroutine calcDerivativesDiatoms_simple
     
  end module diatoms_simple