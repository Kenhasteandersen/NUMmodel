!
! Module to handle generalist unicellulars
!
module generalists
  use globals
  use spectrum
  implicit none

  real(dp), parameter:: rhoCN = 5.68
  !
  ! Light uptake:
  !
  real(dp), parameter:: epsilonL = 0.8 ! Light uptake efficiency
  real(dp), parameter:: alphaL = 0.206
  real(dp), parameter:: rLstar = 8.25
  !
  ! Dissolved nutrient uptake:
  !
  real(dp), parameter:: alphaN = 0.682 ! L/d/mugC/mum^2
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
  end type spectrumGeneralists
 
contains

  subroutine initGeneralists(this, n, mMax)
    class(spectrumGeneralists):: this
    real(dp), intent(in):: mMax
    integer, intent(in):: n
    real(dp), parameter:: mMin = 3.1623d-9
    real(dp), parameter:: rho = 0.57*1d6*1d-12

    call this%initUnicellular(n, mMin, mMax)
    allocate(this%JFreal(n))

    this%beta = beta
    this%sigma = sigma
    this%epsilonF = epsilonF

    this%r = (3./(4.*pi)*this%m/rho)**onethird
    
    this%AN = alphaN * this%r**(-2.) / (1.+(this%r/rNstar)**(-2.)) * this%m
    this%AL = alphaL/this%r * (1-exp(-this%r/rLstar)) * this%m
    this%AF = alphaF*this%m
    this%JFmax = cF/this%r * this%m
    
    this%JlossPassive = cLeakage/this%r * this%m ! in units of C

    !nu = c * this%m**(-onethird)
    this%nu = 3*delta/this%r
    this%Jmax = alphaJ * this%m * (1.d0-this%nu) ! mugC/day
    this%Jresp = cR*alphaJ*this%m
    !mort = 0*0.005*(Jmax/this%m) * this%m**(-0.25);
    !this%mort2constant = 0.0002*this%n
  end subroutine initGeneralists

  subroutine calcRatesGeneralists(this, L, N, DOC, gammaN, gammaDOC)
    class(spectrumGeneralists), intent(inout):: this
    real(dp), intent(in):: gammaN, gammaDOC
    real(dp), intent(in):: L, N, DOC
    real(dp):: f, JmaxT
    integer:: i

    do i = 1, this%n
       !
       ! Uptakes
       !
       this%JN(i) =   gammaN * fTemp15 * this%AN(i)*N*rhoCN ! Diffusive nutrient uptake in units of C/time
       this%JDOC(i) = gammaDOC * fTemp15 * this%AN(i)*DOC ! Diffusive DOC uptake, units of C/time
       this%JL(i) =   epsilonL * this%AL(i)*L  ! Photoharvesting
       ! Total nitrogen uptake:
       this%JNtot(i) = this%JN(i)+this%JF(i)-this%Jlosspassive(i) ! In units of C
       ! Total carbon uptake
       this%JCtot(i) = this%JL(i)+this%JF(i)+this%JDOC(i) & 
                        - fTemp2*this%Jresp(i)-this%JlossPassive(i)
       ! Liebig + synthesis limitation:
       this%Jtot(i) = min( this%JNtot(i), this%JCtot(i) )
       !f = rates%Jtot(ix)/(rates%Jtot(ix) + JmaxT)
       ! If synthesis-limited then down-regulate feeding:
       JmaxT = fTemp2*this%Jmax(i)
       f = this%Jtot(i)/(this%Jtot(i) + max(0.,JmaxT))
       if (this%Jtot(i) .gt. 0) then
        this%JFreal(i) = max(0.d0, min(JmaxT, this%JF(i) - (this%Jtot(i)-f*JmaxT)))
        !rates%Jtot(ix) = f * JmaxT
       else
        this%JFreal(i) = max(0.d0, this%JF(i))
       end if
      this%Jtot(i) = f * JmaxT ! Apply limitation
      
      this%JLreal(i) = this%JL(i) - max( 0.d0, &
            min((this%JCtot(i) - (this%JF(i)-this%JFreal(i))-this%Jtot(i)), this%JL(i)))

      ! Actual uptakes:
      this%JCtot(i) = &
            + this%JLreal(i)  &
            + this%JDOC(i)  &
            + this%JFreal(i)  &
            - fTemp2*this%Jresp(i)  &
            - this%JlossPassive(i)
      this%JNtot(i) = &
            this%JN(i) + &
            this%JFreal(i) - &
            this%JlossPassive(i)
      !
      ! Losses:
      !
      this%JCloss_feeding(i) = (1.-epsilonF)/epsilonF*this%JFreal(i) ! Incomplete feeding (units of carbon per time)
      this%JCloss_photouptake(i) = (1.-epsilonL)/epsilonL * this%JLreal(i)
      this%JNlossLiebig(i) = max( 0.d0, this%JNtot(i)-this%Jtot(i))  ! In units of C
      this%JClossLiebig(i) = max( 0.d0, this%JCtot(i)-this%Jtot(i)) ! C losses from Liebig, not counting losses from photoharvesting

      this%JNloss(i) = &
            this%JCloss_feeding(i) + &
            this%JNlossLiebig(i) +&
            this%JlossPassive(i) ! In units of C
      this%JCloss(i) = &
            this%JCloss_feeding(i) + &
            this%JCloss_photouptake(i) + &
            this%JClossLiebig(i) +&
            this%JlossPassive(i)
      this%JF(i) = this%JFreal(i)
      !
      ! Test for conservation budget. Should be close to zero:
      !
      !write(*,*) 'N budget', i,':',(rates%JN(ix)+JFreal(i)-JlossPassive(i) &
      !    - rates%JNlossLiebig(ix)  - rates%Jtot(ix))/this%m(i)
      !write(*,*) 'C budget', i,':',(rates%JLreal(ix) + rates%JDOC(ix)+JFreal(i) &
      !    -JlossPassive(i)-fTemp2*Jresp(i) &
      !    - rates%JClossLiebig(ix)  - rates%Jtot(ix))/this%m(i)
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
    !rates%dudt = 0*rates%dudt
    this%mort2 = this%mort2constant*u
    do i = 1, this%n
      mortloss = u(i)*(remin2*this%mort2(i) + reminHTL*this%mortHTL(i))
      !
      ! Update nitrogen:
      !
      dNdt = dNdt  &
           + ((-this%JN(i) &
           +  this%JlossPassive(i) &
           +  this%JNlossLiebig(i) &
           +  this%JCloss_feeding(i))/this%m(i) &
           + this%mort2(i) &
           + reminHTL*this%mortHTL(i)) * u(i)/rhoCN
      !
      ! Update DOC:
      !
      dDOCdt = dDOCdt &
           + ((-this%JDOC(i) &
           +   this%JlossPassive(i) &
           +   this%JClossLiebig(i) &
           +   reminF*this%JCloss_feeding(i))/this%m(i) &
           +  remin2*this%mort2(i) &
           +  reminHTL*this%mortHTL(i)) * u(i)
      !
      ! Update the generalists:
      !
      dudt(i) = (this%Jtot(i)/this%m(i)  &
           !- mort(i) &
           - this%mortpred(i) &
           - this%mort2(i) &
           - this%mortHTL(i))*u(i)
   end do
  
 end subroutine calcDerivativesGeneralists

subroutine printRatesGeneralists(this)
  class(spectrumGeneralists), intent(in):: this

  write(*,*) "Generalists with ", this%n, " size classes:"
  call this%printRatesUnicellular()
end subroutine printRatesGeneralists
 
  ! function getNbalanceGeneralists(this, u, rates) result(Nbalance)
  !   real(dp):: Nbalance
  !   class(spectrumGeneralists), intent(in):: this
  !   type(typeRates), intent(in):: rates
  !   real(dp), intent(in):: u(this%n)

  !   Nbalance = (rates%dudt(idxN) + sum(rates%dudt(1+this%ixOffset:this%ixOffset+this%n) &
  !   + (1-reminHTL)*rates%mortHTL(1+this%ixOffset:this%ixOffset+this%n)*u(1:this%n) &
  !   + (1-1)*mort2*u(1:this%n)**2 & ! full N remineralization of viral mortality
  !   + (1-1)*rates%JCloss_feeding(1+this%ixOffset:this%ixOffset+this%n)/this%m(1:this%n)&
  !      * u(1:this%n))/rhoCN)/u(idxN) ! full N remineralization of feeding losses
  ! end function getNbalanceGeneralists 

  ! function getCbalanceGeneralists(this, u, rates) result(Cbalance)
  !   real(dp):: Cbalance
  !   class(spectrumGeneralists), intent(in):: this
  !   type(typeRates), intent(in):: rates
  !   real(dp), intent(in):: u(this%n)

  !   Cbalance = (rates%dudt(idxDOC) + sum(rates%dudt(1+this%ixOffset:this%ixOffset+this%n) &
  !   + (1-reminHTL)*rates%mortHTL(1+this%ixOffset:this%ixOffset+this%n)*u(1:this%n) &
  !   + (1-remin2)*mort2*u(1:this%n)**2 &
  !   - rates%JLreal(1+this%ixOffset:this%ixOffset+this%n)*u(1:this%n)/this%m(1:this%n) &
  !   + fTemp2*Jresp(1:this%n)*u(1:this%n)/this%m(1:this%n) &
  !   + (1-reminF)*rates%JCloss_feeding(1+this%ixOffset:this%ixOffset+this%n)/this%m(1:this%n)&
  !   * u(1:this%n) ))/u(idxDOC)
  ! end function getCbalanceGeneralists 
  
end module generalists
