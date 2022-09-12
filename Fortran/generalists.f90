!
! Module to handle generalist unicellulars
!
module generalists
  use globals
  use spectrum
  implicit none

  private 

  !real(dp), parameter:: rhoCN = 5.68
  !
  ! Light uptake:
  !
  real(dp), parameter:: epsilonL = 0.8 ! Light uptake efficiency
  real(dp), parameter:: alphaL = 0.13 ! 0.206
  real(dp), parameter:: rLstar = 7.5 !8.25
  !
  ! Dissolved nutrient uptake:
  !
  real(dp), parameter:: alphaN = 0.972 !0.682 ! L/d/mugC/mum^2
  real(dp), parameter:: rNstar = 2 ! mum
  !
  ! Phagotrophy:
  !
  real(dp), parameter:: epsilonF = 0.8 ! Assimilation efficiency
  real(dp), parameter:: alphaF = 0.018 
  real(dp), parameter:: cF = 30.
  real(dp), parameter:: beta = 500.d0
  real(dp), parameter:: sigma = 1.3d0
  !
  ! Costs
  !
  real(dp), parameter :: bL=0.08 ! cost of light harvesting mugC(mugC)^-1
  real(dp), parameter :: bN=0.45 ! cost of N uptake mugC(mugN)^-1
  real(dp), parameter :: bDOC=0.45 ! cost of DOC uptake mugC(mugN)^-1
  real(dp), parameter :: bF=0.35 ! cost of food uptake mugC(mugSi)^-1
  real(dp), parameter :: bg=0.2 ! cost of biosynthsesis -- parameter from literature pending
  !
  ! Metabolism
  !
  real(dp), parameter:: cLeakage = 0.03 ! passive leakage of C and N
  real(dp), parameter:: delta = 0.05 ! Thickness of cell wall in mum
            !The constant is increased a bit to limit the lower cell size
  real(dp), parameter:: alphaJ = 1.5 ! Constant for jmax.  per day
  real(dp), parameter:: cR = 0.1
  !
  ! Biogeo:
  !
  real(dp), parameter:: remin = 0.0 ! fraction of mortality losses reminerilized to DOC
  real(dp), parameter:: remin2 = 0.5d0 ! fraction of virulysis remineralized to DOC
  real(dp), parameter:: reminF = 0.1d0
  real(dp), parameter:: reminHTL = 0.d0 ! fraction of HTL mortality remineralized to N and DOC

  type, extends(spectrumUnicellular) :: spectrumGeneralists
    real(dp), allocatable :: JFreal(:)
    
  contains
    procedure, pass :: initGeneralists
    procedure :: calcRates => calcRatesGeneralists
    procedure :: calcDerivativesGeneralists
    procedure :: printRates => printRatesGeneralists
    procedure :: getNbalanceGeneralists
    procedure :: getCbalanceGeneralists
    procedure :: getProdBact => getProdBactGeneralists
  end type spectrumGeneralists
 
  public initGeneralists, spectrumGeneralists, calcRatesGeneralists, calcDerivativesGeneralists
  public printRatesGeneralists, getNbalanceGeneralists, getCbalanceGeneralists

contains

  subroutine initGeneralists(this, n, mMax)
    class(spectrumGeneralists):: this
    real(dp), intent(in):: mMax
    integer, intent(in):: n
    integer:: i
    real(dp), parameter:: mMin = 3.1623d-9 !2.78d-8!
    real(dp), parameter:: rho = 0.4*1d6*1d-12

    call this%initUnicellular(n, mMin, mMax)
    allocate(this%JFreal(n))

    this%beta = beta
    this%sigma = sigma
    this%epsilonF = epsilonF

    this%r = (3./(4.*pi)*this%m/rho)**onethird
    
    this%nu = 3*delta/this%r
    do i = 1,this%n
      this%nu(i) = min(1.d0, this%nu(i))
    enddo

    this%AN = alphaN * this%r**(-2.) / (1.+(this%r/rNstar)**(-2.)) * this%m
    this%AL = alphaL/this%r * (1-exp(-this%r/rLstar)) * this%m * (1.d0-this%nu)
    this%AF = alphaF*this%m
    this%JFmax = cF/this%r * this%m
    
    this%JlossPassive =cLeakage/this%r * this%m ! in units of C

    this%Jmax = alphaJ * this%m * (1.d0-this%nu) ! mugC/day
    this%Jresp =cR*alphaJ*this%m

   ! this%AL = this%AL * (1.d0 - this%nu)
  end subroutine initGeneralists

  subroutine calcRatesGeneralists(this, L, N, DOC, gammaN, gammaDOC)
    class(spectrumGeneralists), intent(inout):: this
    real(dp), intent(in):: gammaN, gammaDOC
    real(dp), intent(in):: L, N, DOC
    real(dp):: f, JmaxT
    real(dp):: dL(this%n), dN(this%n), dDOC(this%n), Jnetp(this%n), Jnet(this%n),Jlim(this%n)
    integer:: i
       !this%JF(10)= 0.02

    do i = 1, this%n
       !
       ! Uptakes
       !
       this%JN(i) =  gammaN * fTemp15 * this%AN(i)*N*rhoCN ! Diffusive nutrient uptake in units of C/time
       this%JDOC(i) =gammaDOC * fTemp15 * this%AN(i)*DOC ! Diffusive DOC uptake, units of C/time
       this%JL(i) =  epsilonL * this%AL(i)*L  ! Photoharvesting
       JmaxT = fTemp2*this%Jmax(i)
       !this%JF(i)= 0.
       !
       ! Potential net uptake
       !
       Jnetp(i)=this%JL(i)*(1-bL)+this%JDOC(i)*(1-bDOC)+this%JF(i)*(1-bF)-ftemp2*this%Jresp(i) ! I think we don't need it
       !
       ! Calculation of down-regulation factors for
       ! N-uptake  
       dN(i) = min(1., 1./this%JN(i)*(Jnetp(i)-this%JF(i)*(1+bg))/(1+bg+bN))
       ! If synthesis-limited then down-regulate feeding:
       ! Photosynthesis
       if (this%JL(i) .gt. 0) then ! Needed to avoid the risk of division with zero if JL = 0
         dL(i) = min(1.,1./(this%JL(i)*(1-bL))*((this%JN(i)+this%JF(i))*(1+bg)-this%JDOC(i)*(1-bDOC)&
             -this%JF(i)*(1-bF)+ftemp2*this%Jresp(i)+bN*this%JN(i))) !+bN*this%JN(i)
       else
         dL(i) = -1.
       endif
       !*************************************************************************
       !********************************** test *********************************
       !            MAKES NO DIFFERENCE
       !if (dN(i).lt.1) then
       !dL(i)=1
       !endif
       !******************************* end of test ******************************
       !**************************************************************************
       !
       if (dL(i).lt.0) then
        dL(i)=0
        dDOC(i) =min(1.,1/(this%JDOC(i)*(1-bDOC))*((this%JN(i)+this%JF(i))*(1+bg)-this%JF(i)*(1-bF)&
                + bN*this%JN(i)+ftemp2*this%Jresp(i)))
       else 
        dDOC(i)=1.
       end if

       if (dN(i).lt.0) then ! check if N leaks out of the cell
        dN(i)=0.
        Jnet =  1./(1+bg)*(dDOC(i)*this%JDOC(i)*(1.-bDOC)+dL(i)*this%JL(i)*(1-bL) &
        + this%JF(i)*(1-bF)- fTemp2*this%Jresp(i) -bN*dN(i)*this%JN(i)) ! was with max(O.,..) 
        f = (Jnet(i) )/(Jnet(i) + JmaxT)
        this%JNlossLiebig(i) = (1-f)*this%JF(i)-f*JmaxT!-1/(1+bg)*(bN*dN(i)*this%JN(i))!this%JF(i)-!this%JF(i)-Jnet(i)!f*JmaxT
       else 
        Jnet =  1./(1+bg)*(dDOC(i)*this%JDOC(i)*(1.-bDOC)+dL(i)*this%JL(i)*(1-bL) &
        + this%JF(i)*(1-bF)- fTemp2*this%Jresp(i) -bN*dN(i)*this%JN(i)) ! was with max(O.,..) 
        f = (Jnet(i) )/(Jnet(i) + JmaxT)
        end if
       !Jnet =  1./(1+bg)*(dDOC(i)*this%JDOC(i)*(1.-bDOC)+dL(i)*this%JL(i)*(1-bL) &
       ! + this%JF(i)*(1-bF)- fTemp2*this%Jresp(i) -bN*dN(i)*this%JN(i)) ! was with max(O.,..) 
       !
       ! Saturation of net growth
       !
       !f = (Jnet(i) )/(Jnet(i) + JmaxT)
       if ((Jnet(i) + JmaxT).eq.0) then
        f=0.
       end if
        this%JCtot(i) = & 
        !
        (1-f)*(dDOC(i)*this%JDOC(i)+dL(i)*this%JL(i) + this%JF(i) )&
        - (1-f)*fTemp2*this%Jresp(i) &
        - ( (1-f)*(bDOC*dDOC(i)*this%JDOC(i)+dL(i)*this%JL(i)*bL+ this%JF(i)*bF+bN*dN(i)*this%JN(i))&
        +(1-f)*bg*Jnet(i) )
      !
       ! Apply saturation to uptake rates
       !
       this%JNreal(i)=dN(i)*(1-f)*this%JN(i)
       this%JDOCreal(i)=dDOC(i)*(1-f)*this%JDOC(i)
       this%JLreal(i)=dL(i)*(1-f)*this%JL(i)
       this%JFreal(i)=(1-f)*this%JF(i)
       this%Jtot(i)= f * JmaxT-(1-f)*this%JlossPassive(i)
      !        
      ! Actual uptakes:
      !
      this%JNtot(i) = &
            this%JNreal(i) + &
            this%JFreal(i)
      !write(*,*) f     
      !
      ! Losses:
      !
      this%JCloss_feeding(i) = (1.-epsilonF)/epsilonF*this%JFreal(i) ! Incomplete feeding (units of carbon per time)
      this%JCloss_photouptake(i) = (1.-epsilonL)/epsilonL * this%JLreal(i)
      this%Jresptot(i)= (1-f)*(fTemp2*this%Jresp(i)+bDOC*dDOC(i)*this%JDOC(i)+dL(i)*this%JL(i)*bL+ &
                        bN*dN(i)*this%JN(i))+(1-f)*bg*Jnet(i)
      !
      ! Test for conservation budget. Should be close to zero:
      !write(*,*) 'N budget', i,':',(this%JNreal(i)+this%JFreal(i)-this%JNlossLiebig(i)-(1-f)*this%JlossPassive(i) &
      !- this%Jtot(i))/this%m(i)

      write(*,*) 'C budget', i,':',(this%JCtot(i) -(1-f)*this%JlossPassive(i)&
      - this%Jtot(i))/this%m(i) !this works only if we take the negative values of jnet

      !write(*,*) 'C budget', i,':',((1-f)*Jnet(i)- (1-f)*this%JlossPassive(i)-this%Jtot(i))/this%m(i)
      this%f(i)=f
    end do
  end subroutine calcRatesGeneralists

  subroutine calcDerivativesGeneralists(this, u, dNdt, dDOCdt, dudt)
    class(spectrumGeneralists), intent(inout):: this
    real(dp), intent(in):: u(this%n)
    real(dp), intent(inout) :: dNdt, dDOCdt, dudt(this%n)
    real(dp):: mortloss
    integer:: i
    !
    ! To make mass balance check:
    !
    this%mort2 = this%mort2constant*u
    this%jPOM = 0*(1-remin2)*this%mort2 ! non-remineralized mort2 => POM

    do i = 1, this%n
      !
      ! Update nitrogen:
      !
      dNdt = dNdt  &
           + ((-this%JNreal(i) &
           +  (1-this%f(i))*this%JlossPassive(i) &
           +  this%JNlossLiebig(i) &     ! N leakage due to excess food
           +  reminF*this%JCloss_feeding(i))/this%m(i) & !reminF
           +  remin2*this%mort2(i) & 
           + reminHTL*this%mortHTL(i)) * u(i)/rhoCN
      !
      ! Update DOC:
      !
      dDOCdt = dDOCdt &
           + ((-this%JDOCreal(i) &
           + (1-this%f(i))*this%JlossPassive(i) &
           +   this%JCloss_photouptake(i) &
           +   reminF*this%JCloss_feeding(i))/this%m(i) &
           +   remin2*this%mort2(i) & 
           +  reminHTL*this%mortHTL(i)) * u(i)
      !
      ! Update the generalists:
      !
      dudt(i) = (this%Jtot(i)/this%m(i)  &
           !- mort(i) &
           - this%mortpred(i) &
           - this%mort2(i) &
           - this%mortHTL(i))*u(i)
           !write(*,*) 'u',i,':', u(i)
   end do
  
 end subroutine calcDerivativesGeneralists

subroutine printRatesGeneralists(this)
  class(spectrumGeneralists), intent(in):: this

  write(*,*) "Generalists with ", this%n, " size classes:"
  call this%printRatesUnicellular()
end subroutine printRatesGeneralists
 
  function getNbalanceGeneralists(this, N, dNdt, u, dudt) result(Nbalance)
    real(dp):: Nbalance
    class(spectrumGeneralists), intent(in):: this
    real(dp), intent(in):: N,dNdt, u(this%n), dudt(this%n)

    Nbalance = (dNdt + sum( dudt &
    + (1-reminHTL)*this%mortHTL*u &
    + (1-remin2)*this%mort2*u & ! full N remineralization of viral mortality
    + (1-reminF)*this%JCloss_feeding/this%m * u &
       )/rhoCN)/N ! full N remineralization of feeding losses
  end function getNbalanceGeneralists 

  function getCbalanceGeneralists(this, DOC, dDOCdt, u, dudt) result(Cbalance)
    real(dp):: Cbalance
    class(spectrumGeneralists), intent(in):: this
    real(dp), intent(in):: DOC, dDOCdt, u(this%n), dudt(this%n)

    Cbalance = (dDOCdt + sum(dudt &
    + (1-reminHTL)*this%mortHTL*u &
    + (1-remin2)*this%mort2*u &
    - this%JLreal*u/this%m &
    - this%JCloss_photouptake*u/this%m & !saturation effect??
    + this%Jresptot*u/this%m & !plus uptake costs
    + (1-reminF)*this%JCloss_feeding/this%m * u &
    )) / DOC
  end function getCbalanceGeneralists 
  

  function getProdBactGeneralists(this, u) result(ProdBact)
    real(dp):: ProdBact
    class(spectrumGeneralists), intent(in):: this
    real(dp), intent(in):: u(this%n)
    integer:: i

    ProdBact = 0.d0
    do i = 1, this%n
      ProdBact = ProdBact + max(0.d0, this%JDOC(i) - ftemp2*this%Jresp(i))*u(i)/this%m(i)
    enddo

  end function getProdBactGeneralists

end module generalists
