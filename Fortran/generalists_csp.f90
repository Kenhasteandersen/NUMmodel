!
! Module to handle the protists part of the model in Serra-Pompei et al (2020).
! In comparison to generalists.f90, this model does not involve DOC uptake.
!
module generalists_csp
  use globals
  use spectrum
  implicit none

  private
  real(dp), parameter:: rhoCN = 5.68
  real(dp), parameter:: epsilonL = 0.9 ! Light uptake efficiency
  real(dp), parameter:: epsilonF = 0.8 ! Assimilation efficiency
  real(dp), parameter:: cLeakage = 0.00015 ! passive leakage of C and N
  real(dp), parameter:: c = 0.0015 ! Parameter for cell wall fraction of mass.
            !The constant is increased a bit to limit the lower cell size
  real(dp), parameter:: alphaN = 3.75E-5 !0.00004 % 0.000162 % Mathilde.  (2.5e-3 l/d/cm) * (1e6 mug/g)^(-1/3) / 1.5 (g/cm); Andersen et al 2015
  real(dp), parameter:: cN = 0.1
  real(dp), parameter:: alphaL = 0.000914 ! if using Andys shading formula for non-diatoms
  real(dp), parameter:: cL = 21 ! if using Andys shading formula for non-diatoms
  real(dp), parameter:: alphaF = 0.0024 !  Fits to TK data for protists
  real(dp), parameter:: cF = 0.1514! From paper (0.0016)
  real(dp), parameter:: cmaxN = 2.3757 ! From paper (or use 3.98?)
  real(dp), parameter:: cmaxN2 = 1.18
  real(dp), parameter:: alphaJ = 1.5 ! Constant for jmax.  per day
  real(dp), parameter:: cR = 0.1
  real(dp), parameter:: remin = 0.0 ! fraction of mortality losses reminerilized to N and DOC
  real(dp), parameter:: remin2 = 1.d0 ! fraction of virulysis remineralized to N and DOC
  real(dp), parameter:: reminHTL = 1.d0 ! fraction of HTL mortality remineralized
  real(dp), parameter:: beta = 500.d0
  real(dp), parameter:: sigma = 1.3d0
  real(dp), parameter:: pvol = 7.6E-7
  real(dp), parameter:: pvol2 = 1.22
  real(dp), parameter:: pvol3 = 2.96e7
  real(dp), parameter:: pl11 = 0.85
  real(dp), parameter:: pl12 = 0.3
  real(dp), parameter:: pl21 = 8.71
  real(dp), parameter:: pl22 = -0.2
  real(dp), parameter:: rhomup1 = 0.024
  real(dp), parameter:: rhomup2 = 1.1
  real(dp), parameter:: Qmup1 = 0.032
  real(dp), parameter:: Qmup2 = 0.76
  real(dp), parameter:: mu_infp1 = 4.7
  real(dp), parameter:: mu_infp2 = -0.26


 real(dp),  dimension(:), allocatable:: AN(:), AL(:), JNmax(:), JLmax(:), Jmax(:), volu(:), JlossPassive(:)
  real(dp),  dimension(:), allocatable:: nu(:), mort(:), mort2(:)
  real(dp),  dimension(:), allocatable:: JN(:), JL(:), Jresp(:), JFreal(:)
  real(dp),  dimension(:), allocatable:: Vol2(:), rhomu(:),Qmu(:),mu_inf(:),mu_max(:)
  !real(dp):: mort2

  public initGeneralists_csp, calcRatesGeneralists_csp, calcDerivativesGeneralists_csp
contains

  function initGeneralists_csp(n, ixOffset, mMax) result(this)
    type(typeSpectrum):: this
    real(dp), intent(in):: mMax
    integer, intent(in):: n, ixOffset
    real(dp), parameter:: mMin = 3.1623d-9

    this = initSpectrum(typeGeneralist_csp, n, ixOffset, mMin, mMax)

    if (allocated(AN)) then
       deallocate(AN)
       deallocate(AL)
       deallocate(JNmax)
       deallocate(JLmax)
       deallocate(volu)
       deallocate(Jresp)
       deallocate(JlossPassive)
       deallocate(nu)
       deallocate(mort)
       deallocate(mort2)
       deallocate(Vol2)
       deallocate(rhomu)
       deallocate(Qmu)
       deallocate(mu_inf)
       deallocate(mu_max)
       
       deallocate(JN)
       deallocate(JL)
       deallocate(JFreal)
    endif
    allocate(AN(n))
    allocate(AL(n))
    allocate(JNmax(n))
    allocate(JLmax(n))
    allocate(volu(n))
    allocate(Jresp(n))
    allocate(JlossPassive(n))
    allocate(nu(n))
    allocate(mort(n))
    allocate(mort2(n))
    allocate(Vol2(n))
    allocate(rhomu(n))
    allocate(Qmu(n))
    allocate(mu_inf(n))
    allocate(mu_max(n))

    allocate(JN(n))
    allocate(JL(n))
    allocate(JFreal(n))   
 
    this%beta = beta
    this%sigma = sigma
    this%epsilonF = epsilonF

    AN = alphaN*this%m**onethird
    AL = alphaL*this%m**twothirds * (1-exp(-cL*this%m**onethird ))  ! shading formula
    this%AF = alphaF*this%m**threequarters
    JlossPassive =0.d0! cLeakage * this%m**twothirds ! in units of C
    this%JFmax = cF*this%m**twothirds
    JNmax = cmaxN*this%m**cmaxN2
    volu=(this%m/pvol)**pvol2
    JLmax = min(pl11*volu**pl12 , pl21*volu**pl22)*this%m
    !nu = c * this%m**(-onethird)
    !Jmax = 0.d0*this%m! alphaJ * this%m * (1.d0-nu) ! mugC/day

    Vol2=pvol3*this%m**pvol2
    rhomu=rhomup1*Vol2**rhomup2
    Qmu=Qmup1*Vol2**Qmup2
    mu_inf=mu_infp1*Vol2**(mu_infp2)
    mu_max=mu_inf*rhomu/(mu_inf*Qmu + rhomu)
    mu_max=mu_max + mu_max*0.17

    Jresp = 0.2*mu_max*this%m !cR*alphaJ*this%m

    mort = 0.d0 !0*0.005*(Jmax/this%m) * this%m**(-0.25);
    mort2 = mu_max*0.03/(this%z) !0.0002*n
  end function initGeneralists_csp

  subroutine calcRatesGeneralists_csp(this, rates, L, N, gammaN)
    type(typeSpectrum), intent(in):: this
    real(dp), intent(in):: gammaN
    type(typeRates), intent(inout):: rates
    real(dp), intent(in):: L, N
    integer:: ix, i

    do i = 1, this%n
       ix = i+this%ixOffset
       !
       ! Uptakes
       !
       !rates%JN(ix) =   gammaN * AN(i)*N*rhoCN ! Diffusive nutrient uptake in units of C/time
       rates%JN(ix) =   gammaN *JNmax(i) * ftemp15*AN(i) * N * rhoCN / &
        (JNmax(i) +gammaN * ftemp15*AN(i) * N * rhoCN)  

       if (N .lt. 0.d0) then
         write(*,*) "negative N"
       end if

       rates%JDOC(ix) = 0.d0 ! Diffusive DOC uptake, units of C/time
       rates%JL(ix) =  JLmax(i) * AL(i) * L /(JLmax(i) + AL(i) * L)
       !rates%JL(ix) =   epsilonL * AL(i)*L  ! Photoharvesting
       ! Total nitrogen uptake:
       rates%JNtot(ix) = rates%JN(ix)+rates%JF(ix) ! In units of C
       ! Total carbon uptake
       rates%JCtot(ix) = rates%JL(ix)+rates%JF(ix)+rates%JDOC(ix)-Jresp(i)-JlossPassive(i)
       ! Liebig + synthesis limitation:
       rates%Jtot(ix) = min( rates%JNtot(ix), rates%JCtot(ix))
       !f = rates%Jtot(ix)/(rates%Jtot(ix) + Jmax(i))
       ! If synthesis-limited then down-regulate feeding:
       !if (rates%Jtot(ix) .gt. 0) then
      !    JFreal(i) = max(0.d0, rates%JF(ix) - (rates%Jtot(ix)-f*Jmax(i)))
       !else
      !    JFreal(i) = max(0.d0, rates%JF(ix))
       !end if
       JFreal(i) = rates%JF(ix)
       !rates%Jtot(ix) = f * Jmax(i)
       rates%JLreal(ix) = rates%JL(ix) !- max( 0.d0, &
            !min((rates%JCtot(ix) - (rates%JF(ix)-JFreal(i))-rates%Jtot(ix)), rates%JL(ix)))

       ! Actual uptakes:
       !rates%JCtot(ix) = &
        !    + rates%JLreal(ix)  &
        !    + rates%JDOC(ix)  &
        !    + JFreal(i)  &
        !    - Jresp(i)  &
        !    - JlossPassive(i)
      ! rates%JNtot(ix) = &
        !    rates%JN(ix) + &
        !    JFreal(i) - &
        !    JlossPassive(i)
       !
       ! Losses:
       !
       rates%JCloss_feeding(ix) =0.d0! (1.-epsilonF)/epsilonF*JFreal(i) ! Incomplete feeding (units of carbon per time)
       rates%JCloss_photouptake(ix) = 0.d0! (1.-epsilonL)/epsilonL * rates%JLreal(ix)
       rates%JNlossLiebig(ix) = max( 0.d0, rates%JN(ix) - rates%JL(ix) + Jresp(i))  ! In units of C
       rates%JClossLiebig(ix) =0.d0! max( 0.d0, rates%JCtot(ix)-rates%Jtot(ix)) ! C losses from Liebig, not counting losses from photoharvesting

       rates%JNloss(ix) = rates%JNlossLiebig(ix)
            !rates%JCloss_feeding(ix) + &
            !rates%JNlossLiebig(ix) +&
            !JlossPassive(i) ! In units of C
       rates%JCloss(ix) = 0.d0
            !rates%JCloss_feeding(ix) + &
            !rates%JCloss_photouptake(ix) + &
            !rates%JClossLiebig(ix) +&
            !JlossPassive(i)
       rates%JF(ix) = JFreal(i)
    end do
  end subroutine calcRatesGeneralists_csp

  subroutine calcDerivativesGeneralists_csp(this, u, rates)
    type(typeSpectrum), intent(in):: this
    type(typeRates), intent(inout):: rates
    real(dp), intent(in):: u(this%n)
    real(dp):: mortloss
    integer:: i, ix

    do i = 1, this%n
      ix = i+this%ixOffset
      mortloss = u(i)*(remin2*ftemp2*mort2(i)*u(i) +reminHTL* ftemp2*rates%mortHTL(ix))
      !
      ! Update nitrogen:
      !
!!$      rates%dudt(idxN) = rates%dudt(idxN)  &
!!$           + (-rates%JN(ix) &
!!$           + rates%JNloss(ix))*u(i)/this%m(i) &
!!$           + (remin2*mort2*u(i)*u(i) &
!!$           + remin*mortloss)/rhoCN
      rates%dudt(idxN) = rates%dudt(idxN)  &
           + ((-rates%JN(ix) &
           + rates%JNloss(ix))*u(i)/this%m(i) &
           + mortloss)/rhoCN
      !
      ! Update DOC:
      !
!!$      rates%dudt(idxDOC) = rates%dudt(idxDOC) &
!!$           + (-rates%JDOC(ix) &
!!$           + rates%JCloss(ix))*u(i)/this%m(i) &
!!$           + remin2*mort2*u(i)*u(i) &
!!$           + remin*mortloss
      rates%dudt(idxDOC) = 0.d0 !rates%dudt(idxDOC) &
           !+ (-rates%JDOC(ix) &
           !+ rates%JCloss(ix))*u(i)/this%m(i) &
           !+ mortloss
      !
      ! Update the generalists:
      !
      rates%dudt(ix) = (rates%Jtot(ix)/this%m(i)  &
           - mort(i) &
           - rates%mortpred(ix) &
           - mort2(i)*u(i) &
           - ftemp2*rates%mortHTL(ix))*u(i)
   end do
 end subroutine calcDerivativesGeneralists_csp

end module generalists_csp
