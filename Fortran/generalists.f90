!
! Module to handle generalist unicellulars
!
module generalists
  use globals
  use spectrum
  use debug
  implicit none

  private
  real(dp), parameter:: rhoCN = 5.68
  !
  ! Light uptake:
  !
  real(dp), parameter:: epsilonL = 0.95 ! Light uptake efficiency
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
  real(dp), parameter:: remin = 0.0 ! fraction of mortality losses reminerilized to N and DOC
  real(dp), parameter:: remin2 = 0.5d0 ! fraction of virulysis remineralized to N and DOC
  real(dp), parameter:: reminF = 0.1d0
  real(dp), parameter:: reminHTL = 0.d0 ! fraction of HTL mortality remineralized

  real(dp),  dimension(:), allocatable:: AN(:), AL(:), Jmax(:),  JlossPassive(:)
  real(dp),  dimension(:), allocatable:: nu(:), mort(:)
  real(dp),  dimension(:), allocatable:: JN(:), JL(:), Jresp(:), JFreal(:)
  real(dp):: mort2

  public initGeneralists, calcRatesGeneralists, calcDerivativesGeneralists, getProdNetGeneralists
contains

  function initGeneralists(n, ixOffset, mMax) result(this)
    type(typeSpectrum):: this
    real(dp), intent(in):: mMax
    integer, intent(in):: n, ixOffset
    real(dp), parameter:: mMin = 3.1623d-9
    real(dp):: r(n)
    real(dp), parameter:: rho = 0.57*1d6*1d-12

    this = initSpectrum(typeGeneralist, n, ixOffset, mMin, mMax)

    if ( allocated(AN) ) then
       deallocate(AN)
       deallocate(AL)
       deallocate(Jresp)
       deallocate(JlossPassive)
       deallocate(nu)
       deallocate(mort)
       
       deallocate(JN)
       deallocate(JL)
       deallocate(JFreal)
    end if
    
    allocate(AN(n))
    allocate(AL(n))
    allocate(Jresp(n))
    allocate(JlossPassive(n))
    allocate(nu(n))
    allocate(mort(n))

    allocate(JN(n))
    allocate(JL(n))
    allocate(JFreal(n))

    this%beta = beta
    this%sigma = sigma
    this%epsilonF = epsilonF

    r = (3./(4.*pi)*this%m/rho)**onethird
    
    AN = alphaN * r**(-2.) / (1.+(r/rNstar)**(-2.)) * this%m
    AL = alphaL/r * (1-exp(-r/rLstar)) * this%m
    this%AF = alphaF*this%m
    this%JFmax = cF/r * this%m
    
    JlossPassive = cLeakage/r * this%m ! in units of C

    !nu = c * this%m**(-onethird)
    nu = 3*delta/r
    Jmax = alphaJ * this%m * (1.d0-nu) ! mugC/day
    Jresp = cR*alphaJ*this%m
    mort = 0*0.005*(Jmax/this%m) * this%m**(-0.25);
    mort2 = 0.0002*n
  end function initGeneralists

  subroutine calcRatesGeneralists(this, rates, L, N, DOC, gammaN, gammaDOC)
    type(typeSpectrum), intent(in):: this
    real(dp), intent(in):: gammaN, gammaDOC
    type(typeRates), intent(inout):: rates
    real(dp), intent(in):: L, N, DOC
    real(dp):: f, JmaxT
    integer:: ix, i

    do i = 1, this%n
       ix = i+this%ixOffset
       !
       ! Uptakes
       !
       rates%JN(ix) =   gammaN * fTemp15 * AN(i)*N*rhoCN ! Diffusive nutrient uptake in units of C/time
       rates%JDOC(ix) = gammaDOC * fTemp15 * AN(i)*DOC ! Diffusive DOC uptake, units of C/time
       rates%JL(ix) =   epsilonL * AL(i)*L  ! Photoharvesting
       ! Total nitrogen uptake:
       rates%JNtot(ix) = rates%JN(ix)+rates%JF(ix)-Jlosspassive(i) ! In units of C
       ! Total carbon uptake
       rates%JCtot(ix) = rates%JL(ix)+rates%JF(ix)+rates%JDOC(ix)-fTemp2*Jresp(i)-JlossPassive(i)
       ! Liebig + synthesis limitation:
       JmaxT = fTemp2*Jmax(i)
       rates%Jtot(ix) = min( rates%JNtot(ix), rates%JCtot(ix) )
       f = rates%Jtot(ix)/(rates%Jtot(ix) + JmaxT)
       ! If synthesis-limited then down-regulate feeding:
       if (rates%Jtot(ix) .gt. 0) then
        f = rates%Jtot(ix)/(rates%Jtot(ix) + max(0.,JmaxT))
        JFreal(i) = max(0.d0, min(JmaxT, rates%JF(ix) - (rates%Jtot(ix)-f*JmaxT)))
        rates%Jtot(ix) = f * JmaxT
      else
        JFreal(i) = max(0.d0, rates%JF(ix))
      end if
      
      rates%Jtot(ix) = f * Jmax(i)
       rates%JLreal(ix) = rates%JL(ix) - max( 0.d0, &
            min((rates%JCtot(ix) - (rates%JF(ix)-JFreal(i))-rates%Jtot(ix)), rates%JL(ix)))

       ! Actual uptakes:
       rates%JCtot(ix) = &
            + rates%JLreal(ix)  &
            + rates%JDOC(ix)  &
            + JFreal(i)  &
            - fTemp2*Jresp(i)  &
            - JlossPassive(i)
       rates%JNtot(ix) = &
            rates%JN(ix) + &
            JFreal(i) - &
            JlossPassive(i)
       !
       ! Losses:
       !
       rates%JCloss_feeding(ix) = (1.-epsilonF)/epsilonF*JFreal(i) ! Incomplete feeding (units of carbon per time)
       rates%JCloss_photouptake(ix) = (1.-epsilonL)/epsilonL * rates%JLreal(ix)
       rates%JNlossLiebig(ix) = max( 0.d0, rates%JNtot(ix)-rates%Jtot(ix))  ! In units of C
       rates%JClossLiebig(ix) = max( 0.d0, rates%JCtot(ix)-rates%Jtot(ix)) ! C losses from Liebig, not counting losses from photoharvesting

       rates%JNloss(ix) = &
            rates%JCloss_feeding(ix) + &
            rates%JNlossLiebig(ix) +&
            JlossPassive(i) ! In units of C
       rates%JCloss(ix) = &
            rates%JCloss_feeding(ix) + &
            rates%JCloss_photouptake(ix) + &
            rates%JClossLiebig(ix) +&
            JlossPassive(i)
       rates%JF(ix) = JFreal(i)
    end do
  end subroutine calcRatesGeneralists

  subroutine calcDerivativesGeneralists(this, u, rates)
    type(typeSpectrum), intent(in):: this
    type(typeRates), intent(inout):: rates
    real(dp), intent(in):: u(this%n)
    real(dp):: mortloss
    integer:: i, ix

    do i = 1, this%n
      ix = i+this%ixOffset
      mortloss = u(i)*(remin2*mort2*u(i) + reminHTL*rates%mortHTL(ix))
      !
      ! Update nitrogen:
      !
      rates%dudt(idxN) = rates%dudt(idxN)  &
           + ((-rates%JN(ix) &
           +  JlossPassive(i) &
           +  rates%JNlossLiebig(ix) &
           +  reminF*rates%JCloss_feeding(ix))/this%m(i) &
           + remin2*mort2*u(i) &
           + reminHTL*rates%mortHTL(ix)) * u(i)/rhoCN
      !
      ! Update DOC:
      !
      rates%dudt(idxDOC) = rates%dudt(idxDOC) &
           + ((-rates%JDOC(ix) &
           +   JlossPassive(i) &
           +   rates%JClossLiebig(ix) &
           +   rates%JCloss_photouptake(ix) &
           +   reminF*rates%JCloss_feeding(ix))/this%m(i) &
           +  remin2*mort2*u(i) &
           +  reminHTL*rates%mortHTL(ix)) * u(i)
      !
      ! Update the generalists:
      !
      rates%dudt(ix) = (rates%Jtot(ix)/this%m(i)  &
           - mort(i) &
           - rates%mortpred(ix) &
           - mort2*u(i) &
           - rates%mortHTL(ix))*u(i)
   end do

 end subroutine calcDerivativesGeneralists

 function getProdNetGeneralists(this, u, rates) result(ProdNet)
   real(dp):: ProdNet
   type(typeSpectrum), intent(in):: this
   type(typeRates), intent(in):: rates
   real(dp), intent(in):: u(this%n)
   integer:: i

   ProdNet = 0.d0
   do i = 1, this%n
      ProdNet = ProdNet + max( 0.d0, (rates%JLreal(i+this%ixOffset)-Jresp(i))*u(i)/this%m(i) )
    end do
  end function getProdNetGeneralists
 
end module generalists
