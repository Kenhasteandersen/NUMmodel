!
! Module to handle the protists part of the model in Serra-Pompei et al (2020).
! In comparison to generalists.f90, this model does not involve DOC uptake.
!
module generalists_csp
  use globals
  use spectrum
  implicit none

  private

  real(dp), parameter:: Q15corr = 1.5**(-0.5) ! Corrections because parameters are given at Tref = 15
  real(dp), parameter:: Q20corr = 2.5**(-0.5) ! wheras NUM uses Tref = 10

  !real(dp), parameter:: rhoCN = 5.68
  real(dp), parameter:: epsilonL = 0.9 ! Light uptake efficiency
  real(dp), parameter:: epsilonF = 0.8 ! Assimilation efficiency
  real(dp), parameter:: cLeakage = 0.00015 ! passive leakage of C and N
  real(dp), parameter:: c = 0.0015 ! Parameter for cell wall fraction of mass.
            !The constant is increased a bit to limit the lower cell size
  real(dp), parameter:: alphaN = 3.75E-5 *Q15corr !0.00004 % 0.000162 % Mathilde.  (2.5e-3 l/d/cm) * (1e6 mug/g)^(-1/3) / 1.5 (g/cm); Andersen et al 2015
  real(dp), parameter:: cN = 0.1
  real(dp), parameter:: alphaL = 0.000914 ! if using Andys shading formula for non-diatoms
  real(dp), parameter:: cL = 21 ! if using Andys shading formula for non-diatoms
  real(dp), parameter:: alphaF = 0.0024 *Q15corr !  Fits to TK data for protists
  real(dp), parameter:: cF = 0.1514! From paper (0.0016)
  real(dp), parameter:: cmaxN = 2.3757 *Q20corr ! From paper (or use 3.98?)
  real(dp), parameter:: cmaxN2 = 1.18 *Q20corr
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

  !real(dp):: mort2

  type, extends(spectrumUnicellular) :: spectrumGeneralists_csp
   real(dp),  dimension(:), allocatable:: JNmax(:), JLmax(:), JFreal(:)
   real(dp),  dimension(:), allocatable:: Vol2(:), rhomu(:),Qmu(:),mu_inf(:),mu_max(:), volu(:) 
  contains
    procedure, pass :: initGeneralists_csp
    procedure :: calcRates => calcRatesGeneralists_csp
    procedure :: calcDerivativesGeneralists_csp
    procedure :: printRates => printRatesGeneralists_csp
  end type spectrumGeneralists_csp

  public initGeneralists_csp, calcRatesGeneralists_csp, calcDerivativesGeneralists_csp
  public printRatesGeneralists_csp, spectrumGeneralists_csp
contains

  subroutine initGeneralists_csp(this, n, mMax)
    class(spectrumGeneralists_csp), intent(inout):: this
    real(dp), intent(in):: mMax
    integer, intent(in):: n
    real(dp), parameter:: mMin = 10**(-6.7) ! 3.1623d-9

    call this%initUnicellular(n, mMin, mMax)

    if (allocated(this%JLmax)) then
       deallocate(this%JLmax)
       deallocate(this%volu)
       deallocate(this%Vol2)
       deallocate(this%rhomu)
       deallocate(this%Qmu)
       deallocate(this%mu_inf)
       deallocate(this%mu_max)
       deallocate(this%JFreal)
    endif
    allocate(this%JNmax(n))
    allocate(this%JLmax(n))
    allocate(this%JFreal(n))
    allocate(this%volu(n))
    allocate(this%Vol2(n))
    allocate(this%rhomu(n))
    allocate(this%Qmu(n))
    allocate(this%mu_inf(n))
    allocate(this%mu_max(n))
 
    this%beta = beta
    this%sigma = sigma
    this%epsilonF = epsilonF

    this%AN = alphaN*this%m**onethird
    this%AL = alphaL*this%m**twothirds * (1-exp(-cL*this%m**onethird ))  ! shading formula
    this%AF = alphaF*this%m**threequarters
    this%JlossPassive =0.d0! cLeakage * this%m**twothirds ! in units of C
    this%JFmax = cF*this%m**twothirds * Q20corr
    this%JNmax = cmaxN*this%m**cmaxN2 * Q20corr
    this%volu=(this%m/pvol)**pvol2
    this%JLmax = min(pl11*this%volu**pl12 , pl21*this%volu**pl22)*this%m * Q20corr
    !nu = c * this%m**(-onethird)
    !Jmax = 0.d0*this%m! alphaJ * this%m * (1.d0-nu) ! mugC/day

    this%Vol2=pvol3*this%m**pvol2
    this%rhomu=rhomup1*this%Vol2**rhomup2
    this%Qmu=Qmup1*this%Vol2**Qmup2
    this%mu_inf = mu_infp1*this%Vol2**(mu_infp2)
    this%mu_max = this%mu_inf*this%rhomu/(this%mu_inf*this%Qmu + this%rhomu)
    this%mu_max = (this%mu_max + this%mu_max*0.17) * Q20corr
    
    this%Jresp = 0.2*this%mu_max*this%m *Q20corr !cR*alphaJ*this%m

    !this%mort2 = this%mu_max*0.03/(this%z) !0.0002*n
    !write(*,*) mort2
  end subroutine initGeneralists_csp

  subroutine calcRatesGeneralists_csp(this, L, N, F, gammaN)
    class(spectrumGeneralists_csp), intent(inout):: this
    real(dp), intent(in):: gammaN
    real(dp), intent(in):: L, N, F(this%n)
    integer:: i

    do i = 1, this%n
       !
       ! Uptakes
       !
       !this%JN(i) =   gammaN * AN(i)*N*rhoCN ! Diffusive nutrient uptake in units of C/time
       this%JN(i) =   gammaN * this%JNmax(i) * ftemp15*this%AN(i) * N * rhoCN / &
          (ftemp2*this%JNmax(i) +gammaN * ftemp15*this%AN(i) * N * rhoCN)  

       if (N .lt. 0.d0) then
         write(*,*) "negative N"
       end if

       this%JDOC(i) = 0.d0 ! Diffusive DOC uptake, units of C/time
       this%JL(i) =  this%JLmax(i) * this%AL(i) * L /(ftemp2*this%JLmax(i) + this%AL(i) * L)

        ! Recalc feeding because we use a different Q10:
       this%flvl(i) = ftemp15*this%AF(i)*F(i) / (ftemp15*this%AF(i)*F(i)+fTemp2*this%JFmax(i))
       this%JF(i) = this%flvl(i) * fTemp2*this%JFmax(i)
       !this%JL(i) =   epsilonL * AL(i)*L  ! Photoharvesting
       ! Total nitrogen uptake:
       this%JNtot(i) = this%JN(i)+this%JF(i) ! In units of C
       ! Total carbon uptake
       this%JCtot(i) = this%JL(i)+this%JF(i)+this%JDOC(i)-this%Jresp(i)-this%JlossPassive(i)
       ! Liebig + synthesis limitation:
       this%Jtot(i) = min( this%JNtot(i), this%JCtot(i))
       !f = this%Jtot(i)/(this%Jtot(i) + Jmax(i))
       ! If synthesis-limited then down-regulate feeding:
       !if (this%Jtot(i) .gt. 0) then
      !    JFreal(i) = max(0.d0, this%JF(i) - (this%Jtot(i)-f*Jmax(i)))
       !else
      !    JFreal(i) = max(0.d0, this%JF(i))
       !end if
       this%JFreal(i) = this%JF(i)
       !this%Jtot(i) = f * Jmax(i)
       this%JLreal(i) = this%JL(i) !- max( 0.d0, &
            !min((this%JCtot(i) - (this%JF(i)-JFreal(i))-this%Jtot(i)), this%JL(i)))

       ! Actual uptakes:
       !this%JCtot(i) = &
        !    + this%JLreal(i)  &
        !    + this%JDOC(i)  &
        !    + JFreal(i)  &
        !    - Jresp(i)  &
        !    - JlossPassive(i)
      ! this%JNtot(i) = &
        !    this%JN(i) + &
        !    JFreal(i) - &
        !    JlossPassive(i)
       !
       ! Losses:
       !
       this%JCloss_feeding(i) =0.d0! (1.-epsilonF)/epsilonF*JFreal(i) ! Incomplete feeding (units of carbon per time)
       this%JCloss_photouptake(i) = 0.d0! (1.-epsilonL)/epsilonL * this%JLreal(i)
       this%JNlossLiebig(i) = max( 0.d0, this%JN(i) - this%JL(i) + this%Jresp(i))  ! In units of C
       this%JClossLiebig(i) =0.d0! max( 0.d0, this%JCtot(i)-this%Jtot(i)) ! C losses from Liebig, not counting losses from photoharvesting

       this%JNloss(i) = this%JNlossLiebig(i)
            !this%JCloss_feeding(i) + &
            !this%JNlossLiebig(i) +&
            !JlossPassive(i) ! In units of C
       this%JCloss(i) = 0.d0
            !this%JCloss_feeding(i) + &
            !this%JCloss_photouptake(i) + &
            !this%JClossLiebig(i) +&
            !JlossPassive(i)
       this%JF(i) = this%JFreal(i)
    end do

    ! write(*,*) '----'
    ! write(*,*) 'log10(m):',log10(this%m)
    ! write(*,*) 'aN:',ftemp15*AN/this%m
    ! write(*,*) 'aL:',AL/this%m
    ! write(*,*) 'aF:',ftemp15*this%AF/this%m
    ! write(*,*) 'jN:', this%JN(3:12)/this%m
    ! write(*,*) 'jL:', this%JL(3:12)/this%m
    ! write(*,*) 'jF:', this%JF(3:12)/this%m

  end subroutine calcRatesGeneralists_csp

  subroutine calcDerivativesGeneralists_csp(this, u, dNdt, dudt)
    class(spectrumGeneralists_csp), intent(inout):: this
    real(dp), intent(in) :: u(this%n)
    real(dp), intent(inout):: dNdt, dudt(this%n)
    real(dp):: mortloss
    integer:: i

    this%mort2 = this%mu_max*0.03/(this%z) * u!0.0002*n
    do i = 1, this%n
      mortloss = u(i)*(remin2*ftemp2*this%mort2(i) +reminHTL* ftemp2*this%mortHTL(i))
      !
      ! Update nitrogen:
      !
!!$      this%dudt(idxN) = this%dudt(idxN)  &
!!$           + (-this%JN(i) &
!!$           + this%JNloss(i))*u(i)/this%m(i) &
!!$           + (remin2*mort2*u(i)*u(i) &
!!$           + remin*mortloss)/rhoCN
      dNdt = dNdt   &
           + ((-this%JN(i) &
           + this%JNloss(i))*u(i)/this%m(i) &
           + mortloss)/rhoCN
      !
      ! Update DOC:
      !
!!$      this%dudt(idxDOC) = this%dudt(idxDOC) &
!!$           + (-this%JDOC(i) &
!!$           + this%JCloss(i))*u(i)/this%m(i) &
!!$           + remin2*mort2*u(i)*u(i) &
!!$           + remin*mortloss
      !this%dudt(idxDOC) = 0.d0 !this%dudt(idxDOC) &
           !+ (-this%JDOC(i) &
           !+ this%JCloss(i))*u(i)/this%m(i) &
           !+ mortloss
      !
      ! Update the generalists:
      !
      dudt(i) = (this%Jtot(i)/this%m(i)  &
          ! - mort(i) &
           - this%mortpred(i) &
           - this%mort2(i) &
           - ftemp2*this%mortHTL(i))*u(i)
   end do
 end subroutine calcDerivativesGeneralists_csp

 subroutine printRatesGeneralists_csp(this)
   class(spectrumGeneralists_csp), intent(in):: this

   write(*,*) "Generalists_csp with ", this%n, " size classes:"
   call this%printRatesUnicellular()
 end subroutine printRatesGeneralists_csp

end module generalists_csp
