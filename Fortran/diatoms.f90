!
! Module to handle diatoms
! AMALIA, diatomsAa with Andys Jtot assuming co-limitation
! with Down-regulation Andy's structure
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
     real(dp), parameter:: rhoCSi = 3.4 
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

     real(dp):: mort2
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
       real(dp) :: fl !investment in photosynthesis
   
       this = initSpectrum(typeDiatom, n, ixOffset, mMin, mMax)
   
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
       r = (3./(4.*pi*cb*(1-v))*this%m)**onethird ! Andy's approximation
       fl=4.d0/3.d0*a*dcr/kappa !investment in photoharvesting
       !
       ! Affinities:
       !
       AN = alphaN*r**(-2.) / (1.+(r/rNstar)**(-2.)) * this%m
        !AL = alphaL/r * (1-exp(-r/rLstar)) * this%m
       AL = cL*r**2.*( 1.-exp(-kappa*fl*r*(1-v)) )! *this%m
       !AL = alphaL/r * (1-exp(-r/rLstar)) * this%m
       this%beta = 0.d0 ! No feeding
       this%palatability = 0.5d0
 
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
       integer:: ix, i

       real(dp) :: dN,dN1,dN2,dN3, dSi, dSi1,dSi2,dSi3,dL ,Jlieb
    
   
       do i = 1, this%n
          ix = i+this%ixOffset
          !
          ! Uptakes
          !
          rates%JN(ix) =   gammaN * AN(i)*N ! Diffusive nutrient uptake in units of N/time
          rates%JDOC(ix) = gammaDOC * AN(i)*DOC ! Diffusive DOC uptake, units of C/time
          rates%JSi(ix) = gammaSi * AN(i)*Si! Diffusive Si uptake, units of Si/time
          rates%JL(ix) =   epsilonL * AL(i)*L  ! Photoharvesting
          !
          ! calculate Jtot Ecotip eqs(4.31) for three different cases of resourse limitation
          !
       
          !
          ! N-limited case 
          !
          dL = 1. 
          !dL=-1.
          !dN = max(0.,min(1.,(rhoCSi * (rates%JL(ix)*(1-bL)-Jresp(i) ) ) / &
          !( rates%JN(ix)*(rhoCN*bSi+rhoCSi*bN+rhoCN*rhoCsi ) )))
          !dSi = max(0.,min(1.,( rhoCN *(rates%JL(ix)*(1-bL)-Jresp(i) ) ) / &
          !( rates%JSi(ix)* (rhoCN*bSi+rhoCSi*bN+rhoCN*rhoCsi ))))

          !if ((dN .gt. 1.) .or. (dSi .gt. 1.)) then
               ! 
               ! Si-limited case
               !
               
           !    dL = max(0.,(Jresp(i)*rhoCN + rates%JSi(ix)*rhoCSi*rhoCN + rates%JSi(ix)*bSi*rhoCN + rates%JSi(ix)*bN*rhoCSi) &
           !    /(rhoCN*(rates%JL(ix)*(1 -bL))))
           !    dN = max(0.,min(1.,rates%JSi(ix)*rhoCSi / rates%JN(ix)*rhoCN))
           !    dSi = 1.

            !   if ((dL .gt. 1.) .or. (dN .gt. 1.)) then
                    !
                    ! C limited
                    !
                    !dL = ... stuff ... ! Only needed for DOC dynamics
             !       dN = max(0.,min(1.,(rhoCSi * (rates%JL(ix)*(1-bL)-Jresp(i) ) ) / &
             !            ( rates%JN(ix)*(rhoCN*bSi+rhoCSi*bN+rhoCN*rhoCsi ) )))
             !  end if
          !end if
          !
          !if ((dL .gt. 1) .or. (dN .gt. 1) .or. (dSi .gt. 1) )then
          !     write(*,*) rates%JN(ix)
          !     write(*,*) rates%JL(ix)
          !     write(*,*) rates%JSi(ix)
          !     write(*,*) '-----'
          ! end if  

          !
          ! Estimate the limiting growth rate:
          !
          Jlieb= min( rates%JN(ix)*rhoCN, rates%JL(ix)-Jresp(i), rates%JSi(ix)*rhoCSi )
          !    
          ! Account for possible carbon limitation
          !
          Jlieb = min( Jlieb, rates%JL(ix)-Jresp(i) - bN*Jlieb/rhoCN - bSi*Jlieb/rhoCSi)
          !
          !
          ! Liebig's Law
          !
          !Jlieb = max(0., min( rates%JL(ix)*(1-bL)-Jresp(i) - dN*bN*rates%JN(ix) - dSi*bSi*rates%JSi(ix), &
          !     dN*rates%JN(ix)*rhoCN, dSi*rates%JSi(ix)*rhoCSi ) )
          !
          !BEFORE THIS I need to find JCtot 
          !      
          ! Jtot saturation
          rates%Jtot(ix)=Jmax(i)*Jlieb / ( Jlieb + Jmax(i) )
          !
          ! Fix units of JN and JSi:
          !
          rates%JN(ix) = rates%JN(ix)*rhoCN
          rates%JSi(ix) = rates%Jsi(ix)*rhoCSi
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
         ! update rates based on Ecotip (5.7)- no JxLoss
         !
         ! Update nitrogen:
         !
         rates%dudt(idxN) = rates%dudt(idxN)  &
         - ((rates%Jtot(ix))*u(i)/this%m(i) &
         !+ rates%JNloss(ix))*u(i)/this%m(i) &
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
    end subroutine calcDerivativesDiatoms
      
   end module diatoms