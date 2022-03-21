!
! Module to handle diatoms
! AMALIA,  Jtot assuming co-limitation
! with Down-regulation 
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
     !real(dp), parameter:: cm = 6d-7 ! rhoCmem 6*d.-7 ! Carbon content in cell membrane mugC(mum^-3)
     !real(dp), parameter:: cb = 1.25d-7 ! rhoCcyt 1.25*d.-7 ! Carbon content in cell cytoplasm mugC(mum^-3)  
     real(dp), parameter:: v=0.6 ! Vacuole fraction
     !
     ! Light uptake:
     !
     real(dp), parameter:: epsilonL = 0.9 ! Light uptake efficiency
     real(dp), parameter:: alphaL = 0.206
     real(dp), parameter:: rLstar = 8.25
     !real(dp), parameter :: cL=0.18 ! photosynthetic affinity constant
     !
     ! Andy's eqs
     !
     !real(dp), parameter :: y=0.05 !photosynthetic quantum yield (mugCmE^-1)
     !real(dp), parameter :: kappa=0.2 ! light investment coefficient (0.2 mum^-1)
     !real(dp), parameter :: a=2270 !cross-sectional area of a chromophore ( m^2(mol chromophore)^-1 ) 
     !real(dp), parameter :: dcr=40 !chromophore number density in cytoplasm up to~ ( mol chromophore m^-3 )
     !
     ! Costs
     !
     real(dp), parameter :: bL=0.08 ! cost of light harvesting mugC(mugC)^-1
     real(dp), parameter :: bN=0.45 ! cost of N uptake mugC(mugN)^-1
     real(dp), parameter :: bSi=0.45 ! cost of Si uptake mugC(mugSi)^-1
     real(dp), parameter :: bg=0.1 ! cost of biosynthsesis -- parameter from literature pending
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
       procedure :: getNbalanceDiatoms
       procedure :: getCbalanceDiatoms
       procedure :: getSibalanceDiatoms
     end type spectrumDiatoms

     public spectrumDiatoms, initDiatoms, calcRatesDiatoms, calcDerivativesDiatoms
     public printRatesDiatoms, getNbalanceDiatoms, getCbalanceDiatoms, getSiBalanceDiatoms

   contains
       
     subroutine initDiatoms(this, n, mMax)
       class(spectrumDiatoms):: this
       real(dp), intent(in):: mMax
       integer, intent(in):: n
       real(dp), parameter:: mMin = 3.1623d-9
       real(dp), parameter:: rho = 0.57*1d-6
       !real(dp) :: fl

       call this%initUnicellular(n, mMin, mMax)
       allocate(this%JSi(this%n))
       !
       ! Radius:
       !
       this%r = (threequarters/pi * this%m/rho/(1-v))**onethird  ! Andy's approximation

       this%AN = alphaN * this%r**(-2.) / (1.+(this%r/rNstar)**(-2.)) * this%m
       this%AL = alphaL/this%r * (1-exp(-this%r/rLstar)) * this%m * (1.d0-this%nu)
       this%AF = 0.d0
       this%JFmax = 0.d0
 
       this%JlossPassive = cLeakage/this%r * this%m ! in units of C
   
       !nu = c * this%m**(-onethird)
       this%nu = 6**twothirds*pi**onethird*delta * (this%m/rho)**(-onethird) * &
         (v**twothirds + (1.+v)**twothirds)
 
       this%Jmax = alphaJ * this%m * (1.d0-this%nu) ! mugC/day
       this%Jresp = cR*alphaJ*this%m
   
       this%beta = 0.d0 ! No feeding
       this%palatability = 0.5d0 ! Lower risk of predation

       !write(*,*) 'AN:' , this%AN
       !write(*,*) 'r:', this%r

     end subroutine initDiatoms
  
     subroutine calcRatesDiatoms(this, L, N, Si,DOC, gammaN, gammaSi,gammaDOC)
       class(spectrumDiatoms), intent(inout):: this
       real(dp), intent(in):: gammaN, gammaSi, gammaDOC
       real(dp), intent(in):: L, N, Si, DOC
       real(dp):: f, JmaxT
       !real(dp):: dL(this%n)
       real(dp) :: dL, dN, dDOC, dSi, Jnetp, Jnet,Jlim
       integer:: i

    
   
       do i = 1, this%n
          !
          ! Uptakes
          !
          this%JN(i) =  ftemp15* gammaN * this%AN(i)*N*rhoCN ! Diffusive nutrient uptake in units of C/time
          this%JDOC(i) = ftemp15*gammaDOC * this%AN(i)*DOC ! Diffusive DOC uptake, units of C/time
          this%JSi(i) = ftemp15* gammaSi * this%AN(i)*Si*rhoCSi! Diffusive Si uptake, units of C/time
          this%JL(i) =   epsilonL * this%AL(i)*L  ! Photoharvesting
          JmaxT = fTemp2*this%Jmax(i)
          !
          ! Potential net uptake
          Jnetp=this%JL(i)*(1-bL)+this%JDOC(i)*(1-bN)-this%Jresp(i)
          !
          ! Calculation of down-regulation factors for
          ! N-uptake
          dN = max(0.,min(1., 1./this%JN(i)*Jnetp/(1+bg +bN+bSi)))
          ! Si-uptake
          dSi = max(0.,min(1., 1./this%JSi(i)*Jnetp/(1+bg +bN+bSi)))
          !
           if (dN==1 .and. dSi==1) then
             if ( this%JN(i)<this%JSi(i) ) then
                 jlim=this%JN(i)
                 dSi= jlim/this%JSi(i)
             else if ( this%JSi(i)<this%JN(i) ) then
                 jlim= this%JSi(i)
                 dN= jlim/this%JN(i) 
             end if
           end if   
           
           if  ( dSi==1 .and. dN<1 ) then
              jlim=this%JSi(i)
              dN= jlim/this%jN(i)
           else if ( dN==1 .and. dSi<1 ) then
              jlim=this%JN(i)
              dSi= jlim/this%JSi(i)
           end if

           if ( dSi<1 .and. dN<1 ) then
              dL = 1 ;
           else
              dL  = min(1., 1./(this%JL(i)*(1-bL))*(jlim*(1+bg + bSi + bN ) &
              + this%Jresp(i) - (1-bN)*this%jDOC(i)) ) 

              !write(*,*) 'dL:', dL
           end if
                          
           if ( dL<0 ) then
              dL=0
              dDOC =min(1., 1./(this%JDOC(i)*(1-bN))*(jlim*(1+bg+bSi+bN)+this%Jresp(i)))
           end if
           !write(*,*) 'dL=', i ,':', dL(i)
           write(*,*) 'dN=', i ,':', dN
           write(*,*) 'dSi=', i ,':', dSi

          !
          ! Down-regulated fluxes:
          !
          this%JN(i) = dN * this%JN(i) 
          this%JDOC(i) = dDOC * this%JDOC(i)
          this%JLreal(i) = dL * this%JL(i)
          this%JSi(i) = dSi * this%JSi(i)
           !Jnetp  = this%JL(i)*(1-bL)+this%jDOC(i)*(1-bN)-this%Jresp(i)   
           !Jnet = max(0., 1./(1+bg)*(Jnetp - (bN*this%JN(i)+bSi*this%JSi(i)))) 
          !
          ! update jnetp with delta's
          !
          Jnetp  = dL*this%JL(i)*(1-bL)+dDOC*this%jDOC(i)*(1-bN)-this%Jresp(i)   
          Jnet = max(0., 1./(1+bg)*(Jnetp - (bN*dN*this%JN(i)+bSi*dSi*this%JSi(i)))) 
          f = (Jnet - this%JlossPassive(i))/(Jnet - this%JlossPassive(i) + JmaxT)
          this%Jtot(i)= f * JmaxT
          !
          !                                                        add this back to N pool    
          !write(*,*) 'N budget', i,':',(this%JN(i) - (this%JN(i)- this%JlossPassive(i) -this%Jtot(i)) &
          !- this%Jtot(i) - this%JlossPassive(i))/this%m(i)
         !write(*,*) 'Si budget', i,':',(this%JSi(i)-this%JlossPassive(i) - (this%JSi(i)-this%JlossPassive(i)-this%Jtot(i)) &
         !- this%Jtot(i))/this%m(i)
         !write(*,*) 'C budget', i,':',(this%JLreal(i) + this%JDOC(i) -(this%JLreal(i) + this%JDOC(i) &
          !   -this%JlossPassive(i) - this%Jtot(i)-fTemp2*this%Jresp(i)) &
          !   -this%JlossPassive(i)-fTemp2*this%Jresp(i) &
          !   - this%Jtot(i))/this%m(i)
          !write(*,*) 'Jnetp',i,':', Jnetp/this%m(i)
          ! write(*,*) 'Jtot',i,':', this%Jtot/this%m(i)

        end do
     end subroutine calcRatesDiatoms
   
     subroutine calcDerivativesDiatoms(this, u, dNdt, dDOCdt, dSidt, dudt)
       class(spectrumDiatoms), intent(inout):: this
       real(dp), intent(in):: u(this%n)
       real(dp), intent(inout):: dNdt, dDOCdt,dSidt, dudt( this%n )
       real(dp):: mortloss
       integer:: i
   
       this%mort2 = this%mort2constant*u
       do i = 1, this%n
         mortloss = u(i)*(remin2*this%mort2(i) +reminHTL* this%mortHTL(i))
         !
         ! Update nitrogen:
         !
         dNdt = dNdt  &
         + ((-this%JN(i) + (this%JN(i)- this%JlossPassive(i) -this%Jtot(i)) &
         +  this%JlossPassive(i))/this%m(i) &
         + this%mort2(i) &
         + reminHTL*this%mortHTL(i)) * u(i)/rhoCN
         !
         ! Update DOC:
         !
         dDOCdt = dDOCdt &
         + ((-this%JDOC(i) &
         +   this%JlossPassive(i))/this%m(i) &
         +  remin2*this%mort2(i) &
         +  reminHTL*this%mortHTL(i)) * u(i)
         ! update Si
         dSidt = dSidt &
         + ((-this%JSi(i) &
         +   this%JlossPassive(i))/this%m(i) &
         +  remin2*this%mort2(i) &
         +  reminHTL*this%mortHTL(i)) * u(i)
         !
         ! Update the diatoms:
         !
         dudt(i) = (this%Jtot(i)/this%m(i)  &
         !- mort(i) &
         - this%mortpred(i) &
         - this%mort2(i) &
         - this%mortHTL(i))*u(i)

         !write(*,*) 'Nbalance =',i, ':', (dNdt + sum( dudt &
         !+ (1-reminHTL)*this%mortHTL*u &
         !+ (1-1)*this%mort2*u & ! full N remineralization of viral mortality
         !   )/rhoCN) 

            !write(*,*) 'Nbalance =',i, ':', (dNdt + sum( dudt(i:this%n) &
            !+ (1-reminHTL)*this%mortHTL(i:this%n)*u(i:this%n) &
            !+ (1-1)*this%mort2(i:this%n)*u(i:this%n) & ! full N remineralization of viral mortality
             !  )/rhoCN) 
      end do
    end subroutine calcDerivativesDiatoms
      
    subroutine printRatesDiatoms(this)
      class(spectrumDiatoms), intent(in):: this
 
      write(*,*) "Diatoms with ", this%n, " size classes:"
      call this%printRatesUnicellular()
 
      99 format (a10, 20f10.6)
  
      write(*,99) "jSi:", this%JSi / this%m
    end subroutine printRatesDiatoms

    function getNbalanceDiatoms(this, N, dNdt, u, dudt) result(Nbalance)
      real(dp):: Nbalance
      class(spectrumDiatoms), intent(in):: this
      real(dp), intent(in):: N,dNdt, u(this%n), dudt(this%n)
  
      Nbalance = (dNdt + sum( dudt &
      + (1-reminHTL)*this%mortHTL*u &
      + (1-1)*this%mort2*u & ! full N remineralization of viral mortality
         )/rhoCN)/N 
    end function getNbalanceDiatoms 

    function getSibalanceDiatoms(this, Si, dSidt, u, dudt) result(Sibalance)
      real(dp):: Sibalance
      class(spectrumDiatoms), intent(in):: this
      real(dp), intent(in):: Si,dSidt, u(this%n), dudt(this%n)
  
      Sibalance = (dSidt + sum( dudt &
      + (1-reminHTL)*this%mortHTL*u &
      + (1-1)*this%mort2*u & ! full Si remineralization of viral mortality
      )/rhoCSi)/Si 
    end function getSibalanceDiatoms 
  
    function getCbalanceDiatoms(this, DOC, dDOCdt, u, dudt) result(Cbalance)
      real(dp):: Cbalance
      class(spectrumDiatoms), intent(in):: this
      real(dp), intent(in):: DOC, dDOCdt, u(this%n), dudt(this%n)
  
      Cbalance = (dDOCdt + sum(dudt &
      + (1-reminHTL)*this%mortHTL*u &
      + (1-remin2)*this%mort2*u &
      - this%JLreal*u/this%m & ! ??
      + fTemp2*this%Jresp*u/this%m &
      )) / DOC
    end function getCbalanceDiatoms 
    
  
    function getProdBactDiatoms(this, u) result(ProdBact)
      real(dp):: ProdBact
      class(spectrumDiatoms), intent(in):: this
      real(dp), intent(in):: u(this%n)
      integer:: i
  
      ProdBact = 0.d0
      do i = 1, this%n
        ProdBact = ProdBact + max(0.d0, this%JDOC(i) - ftemp2*this%Jresp(i))*u(i)/this%m(i)
      enddo
  
    end function getProdBactDiatoms
   end module diatoms