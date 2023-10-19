!
! Module to handle generalist unicellulars
!
module generalists_simple
  use globals
  use spectrum
  use read_input_module
  implicit none

  private 
  
  real(dp) :: epsilonL, remin2, reminF

  type, extends(spectrumUnicellular) :: spectrumGeneralistsSimple
    real(dp), allocatable :: JFreal(:)
    
  contains
    procedure, pass :: initGeneralistsSimple
    procedure :: calcRates => calcRatesGeneralistsSimple
    procedure :: calcDerivativesGeneralistsSimple
    procedure :: printRates => printRatesGeneralistsSimple
    procedure :: getProdBact => getProdBactGeneralistsSimple
  end type spectrumGeneralistsSimple
 
  public initGeneralistsSimple, spectrumGeneralistsSimple, calcRatesGeneralistsSimple, calcDerivativesGeneralistsSimple
  public printRatesGeneralistsSimple

contains

 subroutine initGeneralistsSimple(this, n,errorio,errorstr)
    use iso_c_binding, only: c_char
    class(spectrumGeneralistsSimple):: this
    integer, intent(in):: n
    logical(1), intent(out):: errorio 
    character(c_char), dimension(*), intent(out) :: errorstr
    integer:: i
    real(dp), parameter:: rho = 0.4*1d6*1d-12
    real(dp) :: mMinGeneralist, mMaxGeneralist
    real(dp) :: alphaN, rNstar  
    real(dp) :: alphaL, rLstar
    real(dp) :: alphaF, cF  
    real(dp) :: cLeakage, delta, alphaJ, cR  
    
    ! no errors to begin with
    errorio=.false.

    print*, 'Loading parameter for generalist simple from ', inputfile, ':'
    call read_input(inputfile,'generalists_simple','mMinGeneralist',mMinGeneralist,errorio,errorstr)
    call read_input(inputfile,'generalists_simple','mMaxGeneralist',mMaxGeneralist,errorio,errorstr)
    call this%initUnicellular(n, mMinGeneralist, mMaxGeneralist)
    call read_input(inputfile,'generalists_simple','alphaN',alphaN,errorio,errorstr)
    call read_input(inputfile,'generalists_simple','rNstar',rNstar,errorio,errorstr)
    call read_input(inputfile,'generalists_simple','alphaL',alphaL,errorio,errorstr)
    call read_input(inputfile,'generalists_simple','alphaF',alphaF,errorio,errorstr)
    call read_input(inputfile,'generalists_simple','cF',cF,errorio,errorstr)
    call read_input(inputfile,'generalists_simple','cLeakage',cLeakage,errorio,errorstr)
    call read_input(inputfile,'generalists_simple','alphaJ',alphaJ,errorio,errorstr)
    call read_input(inputfile,'generalists_simple','cR',cR,errorio,errorstr)
    call read_input(inputfile,'generalists_simple','rLstar',rLstar,errorio,errorstr)
    call read_input(inputfile,'generalists_simple','delta',delta,errorio,errorstr)
    call read_input(inputfile,'generalists_simple','epsilonL',epsilonL,errorio,errorstr)
    call read_input(inputfile,'generalists_simple','remin2',remin2,errorio,errorstr)
    call read_input(inputfile,'generalists_simple','reminF',reminF,errorio,errorstr)
    call read_input(inputfile,'generalists_simple','beta',this%beta,errorio,errorstr)
    call read_input(inputfile,'generalists_simple','sigma',this%sigma,errorio,errorstr)
    call read_input(inputfile,'generalists_simple','epsilonF',this%epsilonF,errorio,errorstr)
    
    
    allocate(this%JFreal(n))

    this%r = (3./(4.*pi)*this%m/rho)**onethird
    
    this%nu = 3*delta/this%r
    do i = 1,this%n
      this%nu(i) = min(1.d0, this%nu(i))
    enddo

    this%AN = alphaN * this%r**(-2.) / (1.+(this%r/rNstar)**(-2.)) * this%m
    this%AL = alphaL/this%r * (1-exp(-this%r/rLstar)) * this%m * (1.d0-this%nu)
    this%AF = alphaF*this%m
    this%JFmax = cF/this%r * this%m
    
    this%JlossPassive = cLeakage/this%r * this%m ! in units of C

    this%Jmax = alphaJ * this%m * (1.d0-this%nu) ! mugC/day
    this%Jresp = cR*alphaJ*this%m

  end subroutine initGeneralistsSimple

  subroutine calcRatesGeneralistsSimple(this, L, N, DOC, gammaN, gammaDOC)
    class(spectrumGeneralistsSimple), intent(inout):: this
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
       !f = this%Jtot(ix)/(this%Jtot(ix) + JmaxT)
       ! If synthesis-limited then down-regulate feeding:
       JmaxT = fTemp2*this%Jmax(i)

       if (this%Jtot(i) .gt. 0) then
        f = this%Jtot(i)/(this%Jtot(i) + max(0.,JmaxT))
        this%JFreal(i) = max(0.d0, min(JmaxT, this%JF(i) - (this%Jtot(i)-f*JmaxT)))
        this%Jtot(i) = f * JmaxT
       else
        !f = this%Jtot(i) / max(,JmaxT)
        this%JFreal(i) = max(0.d0, this%JF(i))
       end if
      
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
      this%JCloss_feeding(i) = (1.-this%epsilonF)/this%epsilonF*this%JFreal(i) ! Incomplete feeding (units of carbon per time)
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
      this%Jresptot(i) = fTemp2*this%Jresp(i)
      this%JDOCreal(i) = this%JDOC(i) - this%JClossLiebig(i)
    end do
    !
    ! Test for conservation budget. Should be close to zero:
    !
    !write(*,*) 'N budget:',(-this%Jtot+this%JN+this%JFreal & ! Gains
    !  -this%JNlossLiebig-this%JlossPassive)/this%m           ! Losses
    !write(*,*) 'C budget:',(-this%Jtot+this%JLreal+this%JDOC+this%JFreal & ! Gains
    !  -fTemp2*this%Jresp - this%JClossLiebig - this%JlossPassive)/this%m   ! Losses
end subroutine calcRatesGeneralistsSimple

  subroutine calcDerivativesGeneralistsSimple(this, u, dNdt, dDOCdt, dudt)
    class(spectrumGeneralistsSimple), intent(inout):: this
    real(dp), intent(in):: u(this%n)
    real(dp), intent(inout) :: dNdt, dDOCdt, dudt(this%n)
    !real(dp):: mortloss
    integer:: i

    this%mort2 = this%mort2constant*u ! "quadratic" mortality
    this%jPOM = (1-remin2)*this%mort2 ! non-remineralized mort2 => POM

    do i = 1, this%n
      !
      ! Update nitrogen:
      !
      dNdt = dNdt  &
           + ((-this%JN(i) &
           +  this%JlossPassive(i) &
           +  this%JNlossLiebig(i) &
           +  this%JCloss_feeding(i))/this%m(i) & ! All feeding losses are reminineralized
           +  remin2*this%mort2(i) & 
           ) * u(i)/rhoCN
      !
      ! Update DOC:
      !
      dDOCdt = dDOCdt &
           + ((-this%JDOC(i) &
           +   this%JlossPassive(i) &
           +   this%JClossLiebig(i) &
           +   this%JCloss_photouptake(i) &
           +   reminF*this%JCloss_feeding(i))/this%m(i) &
           +   remin2*this%mort2(i) & 
           ) * u(i)
      !
      ! Update the generalists:
      !
      dudt(i) = (this%Jtot(i)/this%m(i)  &
           - this%mortpred(i) &
           - this%mort2(i) &
           - this%mortHTL(i))*u(i)
   end do
  
 end subroutine calcDerivativesGeneralistsSimple

subroutine printRatesGeneralistsSimple(this)
  class(spectrumGeneralistsSimple), intent(in):: this

  write(*,*) "Generalists Simple with ", this%n, " size classes:"
  call this%printRatesUnicellular()
end subroutine printRatesGeneralistsSimple

  function getProdBactGeneralistsSimple(this, u) result(ProdBact)
    real(dp):: ProdBact
    class(spectrumGeneralistsSimple), intent(in):: this
    real(dp), intent(in):: u(this%n)
    integer:: i

    ProdBact = 0.d0
    do i = 1, this%n
      ProdBact = ProdBact + max(0.d0, this%JDOC(i) - ftemp2*this%Jresp(i))*u(i)/this%m(i)
    enddo

  end function getProdBactGeneralistsSimple

end module generalists_simple
