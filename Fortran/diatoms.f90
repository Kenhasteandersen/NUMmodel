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
     !real(dp), parameter:: rhoCN = 5.68
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

     type, extends(spectrumUnicellular) :: spectrumDiatoms
       real(dp), dimension(:), allocatable:: JSi

     contains
       procedure, pass :: initDiatoms
       procedure :: calcRates => calcRatesDiatoms
       procedure :: calcDerivativesDiatoms
       procedure :: printRates => printRatesDiatoms
     end type spectrumDiatoms

     public spectrumDiatoms, initDiatoms, calcRatesDiatoms, calcDerivativesDiatoms, printRatesDiatoms
   contains
       
     subroutine initDiatoms(this, n, mMax)
       class(spectrumDiatoms):: this
       real(dp), intent(in):: mMax
       integer, intent(in):: n
       real(dp), parameter:: mMin = 3.1623d-9
       real(dp):: hs(n) ! Shell thickness mum ! TO BE FIXED!!
       real(dp), parameter:: rho = 0.57*1d6*1d-12
       real(dp) :: fl !investment in photosynthesis
   
       call this%initUnicellular(n, mMin, mMax)
       allocate(this%JSi(this%n))
       !
       ! Radius:
       !
       this%r = (3./(4.*pi*cb*(1-v))*this%m)**onethird ! Andy's approximation
       fl=4.d0/3.d0*a*dcr/kappa !investment in photoharvesting
       !
       ! Affinities:
       !
       this%AN = alphaN*this%r**(-2.) / (1.+(this%r/rNstar)**(-2.)) * this%m
        !AL = alphaL/r * (1-exp(-r/rLstar)) * this%m
       this%AL = cL*this%r**2.*( 1.-exp(-kappa*fl*this%r*(1-v)) )! *this%m
       !AL = alphaL/r * (1-exp(-r/rLstar)) * this%m
       this%beta = 0.d0 ! No feeding
       this%palatability = 0.5d0
 
       this%JlossPassive = cLeakage/this%r * this%m ! in units of C
   
       this%nu = 3*delta/this%r ! wall fraction
       hs = 0.0023*this%r**(0.72) !eq from ecotip report
       this%Jmax = alphaJ * this%m * (1.d0-this%nu) ! mugC/day
       this%Jresp = cR*alphaJ*this%m
       !mort = 0*0.005*(Jmax/this%m) * this%m**(-0.25);
       
     end subroutine initDiatoms
  
     subroutine calcRatesDiatoms(this, L, N, Si, gammaN, gammaSi)
       class(spectrumDiatoms), intent(inout):: this
       real(dp), intent(in):: gammaN, gammaSi
       real(dp), intent(in):: L, N, Si
       integer:: i

       !real(dp) :: dN,dN1,dN2,dN3, dSi, dSi1,dSi2,dSi3,dL ,Jlieb
       real(dp) :: dL, Jlieb
    
   
       do i = 1, this%n
          !
          ! Uptakes
          !
          this%JN(i) =   gammaN * this%AN(i)*N ! Diffusive nutrient uptake in units of N/time
          !this%JDOC(i) = gammaDOC * this%AN(i)*DOC ! Diffusive DOC uptake, units of C/time
          this%JSi(i) = gammaSi * this%AN(i)*Si! Diffusive Si uptake, units of Si/time
          this%JL(i) =   epsilonL * this%AL(i)*L  ! Photoharvesting
          !
          ! calculate Jtot Ecotip eqs(4.31) for three different cases of resourse limitation
          !
       
          !
          ! N-limited case 
          !
          dL = 1. 
          !dL=-1.
          !dN = max(0.,min(1.,(rhoCSi * (this%JL(i)*(1-bL)-Jresp(i) ) ) / &
          !( this%JN(i)*(rhoCN*bSi+rhoCSi*bN+rhoCN*rhoCsi ) )))
          !dSi = max(0.,min(1.,( rhoCN *(this%JL(i)*(1-bL)-Jresp(i) ) ) / &
          !( this%JSi(i)* (rhoCN*bSi+rhoCSi*bN+rhoCN*rhoCsi ))))

          !if ((dN .gt. 1.) .or. (dSi .gt. 1.)) then
               ! 
               ! Si-limited case
               !
               
           !    dL = max(0.,(Jresp(i)*rhoCN + this%JSi(i)*rhoCSi*rhoCN + this%JSi(i)*bSi*rhoCN + this%JSi(i)*bN*rhoCSi) &
           !    /(rhoCN*(this%JL(i)*(1 -bL))))
           !    dN = max(0.,min(1.,this%JSi(i)*rhoCSi / this%JN(i)*rhoCN))
           !    dSi = 1.

            !   if ((dL .gt. 1.) .or. (dN .gt. 1.)) then
                    !
                    ! C limited
                    !
                    !dL = ... stuff ... ! Only needed for DOC dynamics
             !       dN = max(0.,min(1.,(rhoCSi * (this%JL(i)*(1-bL)-Jresp(i) ) ) / &
             !            ( this%JN(i)*(rhoCN*bSi+rhoCSi*bN+rhoCN*rhoCsi ) )))
             !  end if
          !end if
          !
          !if ((dL .gt. 1) .or. (dN .gt. 1) .or. (dSi .gt. 1) )then
          !     write(*,*) this%JN(i)
          !     write(*,*) this%JL(i)
          !     write(*,*) this%JSi(i)
          !     write(*,*) '-----'
          ! end if  

          !
          ! Estimate the limiting growth rate:
          !
          Jlieb= min( this%JN(i)*rhoCN, this%JL(i)-this%Jresp(i), this%JSi(i)*rhoCSi )
          !    
          ! Account for possible carbon limitation
          !
          Jlieb = min( Jlieb, this%JL(i)-this%Jresp(i) - bN*Jlieb/rhoCN - bSi*Jlieb/rhoCSi)
          !
          !
          ! Liebig's Law
          !
          !Jlieb = max(0., min( this%JL(i)*(1-bL)-Jresp(i) - dN*bN*this%JN(i) - dSi*bSi*this%JSi(i), &
          !     dN*this%JN(i)*rhoCN, dSi*this%JSi(i)*rhoCSi ) )
          !
          !BEFORE THIS I need to find JCtot 
          !      
          ! Jtot saturation
          ! CHECK THIS BELOW: IT DOES NOT HANDLE NEGATIVE VALUES OF JMAX AND JTOT WELL
          this%Jtot(i)=this%Jmax(i)*Jlieb / ( Jlieb + this%Jmax(i) )
          !
          ! Fix units of JN and JSi:
          !
          this%JN(i) = this%JN(i)*rhoCN
          this%JSi(i) = this%Jsi(i)*rhoCSi
       end do
 
     end subroutine calcRatesDiatoms
   
     subroutine calcDerivativesDiatoms(this, u, dNdt, dSidt, dudt)
       class(spectrumDiatoms), intent(inout):: this
       real(dp), intent(in):: u(this%n)
       real(dp), intent(inout):: dNdt, dSidt, dudt( this%n )
       real(dp):: mortloss
       integer:: i
   
       this%mort2 = this%mort2constant*u
       do i = 1, this%n
         mortloss = u(i)*(remin2*this%mort2(i) +reminHTL* this%mortHTL(i))
         ! update rates based on Ecotip (5.7)- no JxLoss
         !
         ! Update nitrogen:
         !
         dNdt = dNdt  & 
         - ((this%Jtot(i))*u(i)/this%m(i) &
         !+ this%JNloss(i))*u(i)/this%m(i) &
         + mortloss)/rhoCN
 
         !
         ! Update DOC:
         !
         !this%dudt(idxDOC) = 0.d0 !this%dudt(idxDOC) & ! NOTE: DON'T SET dDOCdt TO ZERO!!!
              !+ (-this%JDOC(i) &
              !+ this%JCloss(i))*u(i)/this%m(i) &
              !+ mortloss
         !
         ! Update Si:
         !
         dSidt = dSidt &
              + ((-this%Jtot(i))*u(i)/this%m(i) &
              + mortloss)/rhoCSi
         !
         ! Update the diatoms:
         !
         dudt(i) = (this%Jtot(i)/this%m(i)  &
              !- mort(i) &
              - this%mortpred(i) &
              - this%mort2(i) &
              - this%mortHTL(i))*u(i)
 
      end do
    end subroutine calcDerivativesDiatoms
      
    subroutine printRatesDiatoms(this)
      class(spectrumDiatoms), intent(in):: this
 
      write(*,*) "Diatoms with ", this%n, " size classes:"
      call this%printRatesUnicellular()
 
      99 format (a10, 20f10.6)
  
      write(*,99) "jSi:", this%JSi / this%m
    end subroutine printRatesDiatoms

   end module diatoms