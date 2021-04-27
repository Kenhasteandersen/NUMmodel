!
! Module to handle generalist unicellulars
!
module diatoms
    use globals
    use spectrum
    implicit none

    private
    real(dp), parameter:: rhoCN = 5.68
    real(dp), parameter:: rhoCSi = 100. ! TO BE FIXED!!!
    !
    ! Light uptake:
    !
    real(dp), parameter:: epsilonL = 0.9 ! Light uptake efficiency
    real(dp), parameter:: alphaL = 0.206
    real(dp), parameter:: rLstar = 8.25
    !
    ! Dissolved nutrient uptake:
    !
    real(dp), parameter:: alphaN = 0.682 ! L/d/mugC/mum^2
    real(dp), parameter:: rNstar = 2 ! mum
    !
    ! Metabolism
    !
    real(dp), parameter:: cLeakage = 0.03 ! passive leakage of C and N
    real(dp), parameter:: c = 0.0015 ! Parameter for cell wall fraction of mass.
    real(dp), parameter:: delta = 0.05 ! Thickness of cell wall in mum
              !The constant is increased a bit to limit the lower cell size
    real(dp), parameter:: alphaJ = 1.5 ! Constant for jmax.  per day
    real(dp), parameter:: cR = 0.1
    !
    ! Biogeo:
    !
    real(dp), parameter:: remin = 0.0 ! fraction of mortality losses reminerilized to N and DOC
    real(dp), parameter:: remin2 = 1.d0 ! fraction of virulysis remineralized to N and DOC
    real(dp), parameter:: reminHTL = 0.d0 ! fraction of HTL mortality remineralized
  
    real(dp),  dimension(:), allocatable:: AN(:), AL(:), Jmax(:),  JlossPassive(:)
    real(dp),  dimension(:), allocatable:: nu(:), mort(:)
    real(dp),  dimension(:), allocatable:: Jresp(:)
    real(dp):: mort2
  
    public initDiatoms, calcRatesDiatoms, calcDerivativesDiatoms
  contains
  
    function initDiatoms(n, ixOffset, mMax) result(this)
      type(typeSpectrum):: this
      real(dp), intent(in):: mMax
      integer, intent(in):: n, ixOffset
      real(dp), parameter:: mMin = 3.1623d-9
      real(dp):: r(n)
      real(dp), parameter:: rho = 0.57*1d6*1d-12
  
      this = initSpectrum(typeDiatom, n, ixOffset, mMin, mMax)
  
      if ( allocated(AN) ) then
         deallocate(AN)
         deallocate(AL)
         deallocate(Jmax)
         deallocate(Jresp)
         deallocate(JlossPassive)
         deallocate(nu)
         deallocate(mort)
         
         !deallocate(JN)
         !deallocate(JL)
         !deallocate(JFreal)
      end if
      
      allocate(AN(n))
      allocate(AL(n))
      allocate(Jmax(n))
      allocate(Jresp(n))
      allocate(JlossPassive(n))
      allocate(nu(n))
      allocate(mort(n))
  
      !allocate(JN(n))
      !allocate(JL(n))
      !allocate(JFreal(n))
    
      r = (3./(4.*pi)*this%m/rho)**onethird
      
      AN = alphaN*r**(-2.) / (1.+(r/rNstar)**(-2.)) * this%m
      AL = alphaL/r * (1-exp(-r/rLstar)) * this%m
      
      this%beta = 0.d0 ! No feeding

      JlossPassive = cLeakage/r * this%m ! in units of C
  
      !nu = c * this%m**(-onethird)
      nu = 3*delta/r
      Jmax = alphaJ * this%m * (1.d0-nu) ! mugC/day
      Jresp = cR*alphaJ*this%m
      mort = 0*0.005*(Jmax/this%m) * this%m**(-0.25);
      mort2 = 0.0002*n
    end function initDiatoms
 
    subroutine calcRatesDiatoms(this, rates, L, N, DOC, Si, gammaN, gammaDOC, gammaSi)
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
         rates%JSi(ix) = gammaSi * AN(i)*Si*rhoCSi ! Diffusive Si uptake, units of C/time
         rates%JL(ix) =   epsilonL * AL(i)*L  ! Photoharvesting
         ! Total nitrogen uptake:
         rates%JNtot(ix) = rates%JN(ix)+rates%JF(ix)-Jlosspassive(i) ! In units of C
         ! Total carbon uptake
         rates%JCtot(ix) = rates%JL(ix)+rates%JF(ix)+rates%JDOC(ix)-Jresp(i)-JlossPassive(i)
         ! Liebig + synthesis limitation:
         rates%Jtot(ix) = min( rates%JNtot(ix), rates%JCtot(ix), rates%JSi(ix) )
  
         ! Actual uptakes:
         rates%JCtot(ix) = &
              + rates%JL(ix)  &
              + rates%JDOC(ix)  &
              - Jresp(i)  &
              - JlossPassive(i)
         rates%JNtot(ix) = &
              + rates%JN(ix) &
              + JlossPassive(i)
         !
         ! Losses:
         !
         rates%JCloss_photouptake(ix) = (1.-epsilonL)/epsilonL * rates%JLreal(ix)
         rates%JNlossLiebig(ix) = max( 0.d0, rates%JNtot(ix)-rates%Jtot(ix))  ! In units of C
         rates%JClossLiebig(ix) = max( 0.d0, rates%JCtot(ix)-rates%Jtot(ix)) ! C losses from Liebig, not counting losses from photoharvesting
         rates%JSiloss(ix) = 0.d0 ! NEEDS TO BE FIXED

         rates%JNloss(ix) = &
              + rates%JNlossLiebig(ix)&
              + JlossPassive(i) ! In units of C
         rates%JCloss(ix) = &
              + rates%JCloss_photouptake(ix) &
              + rates%JClossLiebig(ix) &
              + JlossPassive(i)
      end do
    end subroutine calcRatesDiatoms
  
    subroutine calcDerivativesDiatoms(this, u, rates)
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
             + ((-rates%JN(ix) &
             + rates%JNloss(ix))*u(i)/this%m(i) &
             + mortloss)/rhoCN
        !
        ! Update DOC:
        !
        rates%dudt(idxDOC) = rates%dudt(idxDOC) &
             + (-rates%JDOC(ix) &
             + rates%JCloss(ix))*u(i)/this%m(i) &
             + mortloss
        !
        ! Update Si:
        !
        rates%dudt(idxSi) = rates%dudt(idxSi) &
             + ((-rates%JSi(ix) &
             + rates%JSiLoss(ix))*u(i)/this%m(i) &
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
   end subroutine calcDerivativesDiatoms
     
  end module diatoms
  