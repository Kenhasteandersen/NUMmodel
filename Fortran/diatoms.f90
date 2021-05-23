!
! Module to handle diatoms
! AMALIA, diatomsA2 with Andys Jtot assuming co-limitation
! matlab crushes with this one
!
module diatoms
    use globals
    use spectrum
    implicit none

    private
    !
    ! Stoichiometry:
    !
    real(dp), parameter:: rhoCN = 5.68
    real(dp), parameter:: rhoCSi =3.4 ! TO BE FIXED!!!
     !cell properties
     real(dp), parameter:: cm = 6d-7 ! rhoCmem 6*d.-7 ! Carbon content in cell membrane mugC(mum^-3)
     real(dp), parameter:: cb = 1.25d-7 ! rhoCcyt 1.25*d.-7 ! Carbon content in cell cytoplasm mugC(mum^-3)  
     real(dp), parameter:: v=0.6 ! Vacuole fraction
    !
    ! Light uptake:
    !
    real(dp), parameter:: epsilonL = 0.9 ! Light uptake efficiency
    real(dp), parameter:: alphaL = 0.206
    real(dp), parameter:: rLstar = 8.25
     ! Andy's eqs
    real(dp), parameter :: y=0.05 !photosynthetic quantum yield (mugCmE^-1)
    real(dp), parameter :: kappa=0.2 ! light investment coefficient (0.2 mum^-1)
    real(dp), parameter :: a=2270 !cross-sectional area of a chromophore ( m^2(mol chromophore)^-1 ) 
    real(dp), parameter :: dcr=40 !chromophore number density in cytoplasm up to~ ( mol chromophore m^-3 )
    real(dp), parameter :: cL=0.18 ! photosynthetic affinity constant (= pi*y)
    !
    ! Costs
    !
    real(dp), parameter :: bL=0.08 ! cost of light harvesting mugC(mugC)^-1
    real(dp), parameter :: bN=0.45 ! cost of N uptake mugC(mugN)^-1
    real(dp), parameter :: bSi=0.45 ! cost of Si uptake mugC(mugSi)^-1
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
    real(dp), dimension(:), allocatable:: JsN(:), JsSi(:), JsC(:)
    real(dp):: mort2
    real(dp), allocatable :: fl !investment in photosynthesis
    real(dp), allocatable :: dN, dSi, dL 
    public initDiatoms, calcRatesDiatoms, calcDerivativesDiatoms
  contains
      
    function initDiatoms(n, ixOffset, mMax) result(this)
      type(typeSpectrum):: this
      real(dp), intent(in):: mMax
      integer, intent(in):: n, ixOffset
      real(dp), parameter:: mMin = 3.1623d-9
      real(dp):: r(n)
      real(dp):: hs(n) ! Shell thickness mum ! TO BE FIXED!!
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
         deallocate(JsN)
      end if
      
      allocate(AN(n)) 
      allocate(AL(n))
      allocate(Jmax(n))
      allocate(Jresp(n))
      allocate(JlossPassive(n))
      allocate(nu(n))
      allocate(mort(n))
      allocate(JsN(n)) !! YOU FORGOT TO ALLOCATE THIS VARIABLE
      allocate(JsSi(n)) !! YOU FORGOT TO ALLOCATE THIS VARIABLE
      allocate(JsC(n)) !! YOU FORGOT TO ALLOCATE THIS VARIABLE
      
      ! Radius:
      !
      r = (3./(4.*pi*cb*(1-v))*this%m)**onethird ! Andy's approximation

      fl=4.d0/3.d0*a*dcr/kappa !investment in photoharvesting
      !
      ! Affinities:
      !
      AN = alphaN*r**(-2.) / (1.+(r/rNstar)**(-2.)) * this%m
       !AL = alphaL/r * (1-exp(-r/rLstar)) * this%m
      AL = cL*r**2.*( 1.-exp(-kappa*fl*r*(1-v)) ) *this%m
      this%beta = 0.d0 ! No feeding

      JlossPassive = cLeakage/r * this%m ! in units of C
  
      nu = 3*delta/r ! wall fraction
      hs = 0.0023*r**(0.72) !eq from ecotip report
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
      real(dp) :: f
      integer:: ix, i, jx, j
  
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
         rates%JNtot(ix) = rates%JN(ix)-Jlosspassive(i) ! In units of C
         ! Total carbon uptake
         rates%JCtot(ix) = rates%JL(ix)+rates%JDOC(ix)-Jresp(i)-JlossPassive(i)
         ! Liebig + synthesis limitation:
         !rates%Jtot(ix) = min(Jmax(i),rates%JNtot(ix), rates%JCtot(ix), rates%JSi(ix) )
         !
         ! calculate Jtot for three different cases of resourse limitation
         !
         dL=1 ! for all three cases
         ! N-limited case 
         dN=1
         dSi= rhoCN*rates%JN(ix)/ (rhoCSi*rates%JSi(ix))  
         JsN(ix)= dL*rates%JL(ix)*(1-bL)-dN*rates%JN(ix)*bN-dSi*rates%JSi(ix)*bSi-Jresp(i) 
         ! Si-limited case
         dSi=1
         dN= rhoCSi*rates%JSi(ix) / (rhoCN*rates%JN(ix))
         JsSi(ix)= dL*rates%JL(ix)*(1-bL)-dN*rates%JN(ix)*bN-dSi*rates%JSi(ix)*bSi-Jresp(i) 
         ! C-limited case
         dN=(rhoCSi*rates%JSi(ix) * (rates%JL(ix)*(1-bL)-Jresp(i) ) ) / &
         ( rhoCSi*rates%JSi(ix)*rates%JN(ix)*(bN+rhoCN) + rhoCN*rates%JN(ix) )
         dSi= ( rhoCN*rates%JN(ix)* (rates%JL(ix)*(1-bL)-Jresp(i) ) ) / &
         ( rhoCN*rates%JSi(ix)*rates%JN(ix)*(bSi+rhoCSi)+rhoCSi*rates%JSi(ix) )
         JsC(ix)= dL*rates%JL(ix)*(1-bL)-dN*rates%JN(ix)*bN-dSi*rates%JSi(ix)*bSi-Jresp(i)
         !      
         ! find which case yields the maximum Jtot
            !do j = 1, this%n
             !  jx = j+this%ixOffset
               rates%Jtot(ix)=max(JsN(ix), JsSi(ix), JsC(ix))
            !end do    
         !
         ! jtot after down-regulation for diatoms, eqs 4.31 Ecotip
         !rates%Jtot(ix)=dL*rates%JL(ix)*(1-bL)-dN*rates%JN(ix)*bN-dSi*rates%JSi(ix)*bSi-Jresp(i)   

         f=rates%Jtot(ix)/(rates%Jtot(ix)+Jmax(i))
         rates%Jtot(ix)= f*Jmax(i)
         !rates%JLreal(ix)= rates%JL(ix)-max(0.d0, &
         !min((rates%JCtot(ix)-rates%Jtot(ix)), rates%JL(ix))) 
         rates%JLreal(ix)= rates%JL(ix)-&
               min((rates%JCtot(ix)-rates%Jtot(ix)), rates%JL(ix)) 
         !write(*,*) 'JLrealD=' rates%JLreal(ix)/ this%m(i)
         ! Actual uptakes:
         rates%JCtot(ix) = &
              + rates%JLreal(ix)  &
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
  