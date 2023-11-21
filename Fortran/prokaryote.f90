!
! Module to handle prokaryote. Based on generalist_simple
!
module prokaryote
  use globals
  use spectrum
  use read_input_module
  implicit none

  private 
  
  real(dp) :: epsilonL, remin2, reminF

  type, extends(spectrumUnicellular) :: spectrumProkaryote
    real(dp), allocatable :: JFreal(:)
    
  contains
    procedure, pass :: initProkaryote
    procedure :: calcRates => calcRatesProkaryote
    procedure :: calcDerivativesProkaryote
    procedure :: printRates => printRatesProkaryote
    procedure :: getProdBact => getProdBactProkaryote
  end type spectrumProkaryote
 
  public initProkaryote, spectrumProkaryote, calcRatesProkaryote, calcDerivativesProkaryote
  public printRatesProkaryote
  public initProkaryoteX


contains

 subroutine initProkaryote(this, n,errorio,errorstr)
    use iso_c_binding, only: c_char
    class(spectrumProkaryote):: this
    integer, intent(in):: n
    logical(1), intent(out):: errorio 
    character(c_char), dimension(*), intent(out) :: errorstr
    integer:: i
    real(dp), parameter:: rho = 0.4*1d6*1d-12
    real(dp) :: ii
    real(dp) :: mMinProkaryote, mMaxProkaryote
    real(dp) :: alphaN, rNstar  
    real(dp) :: alphaL, rLstar
    real(dp) :: alphaF, cF  
    real(dp) :: cLeakage, delta, alphaJ, cR  
    
    ! no errors to begin with
    errorio=.false.

    print*, 'Loading parameter for prokaryote simple from ', inputfile, ':'
    call read_input(inputfile,'prokaryote','mMinProkaryote',mMinProkaryote,errorio,errorstr)
    call read_input(inputfile,'prokaryote','mMaxProkaryote',mMaxProkaryote,errorio,errorstr)
    call read_input(inputfile,'prokaryote','alphaN',alphaN,errorio,errorstr)
    call read_input(inputfile,'prokaryote','rNstar',rNstar,errorio,errorstr)
    call read_input(inputfile,'prokaryote','alphaL',alphaL,errorio,errorstr)
    call read_input(inputfile,'prokaryote','alphaF',alphaF,errorio,errorstr)
    call read_input(inputfile,'prokaryote','cF',cF,errorio,errorstr)
    call read_input(inputfile,'prokaryote','cLeakage',cLeakage,errorio,errorstr)
    call read_input(inputfile,'prokaryote','alphaJ',alphaJ,errorio,errorstr)
    call read_input(inputfile,'prokaryote','cR',cR,errorio,errorstr)
    call read_input(inputfile,'prokaryote','rLstar',rLstar,errorio,errorstr)
    call read_input(inputfile,'prokaryote','delta',delta,errorio,errorstr)
    call read_input(inputfile,'prokaryote','epsilonL',epsilonL,errorio,errorstr)
    call read_input(inputfile,'prokaryote','remin2',remin2,errorio,errorstr)
    call read_input(inputfile,'prokaryote','reminF',reminF,errorio,errorstr)
    call read_input(inputfile,'prokaryote','beta',this%beta,errorio,errorstr)
    call read_input(inputfile,'prokaryote','sigma',this%sigma,errorio,errorstr)
    call read_input(inputfile,'prokaryote','epsilonF',this%epsilonF,errorio,errorstr)
    
    
    
    call this%initUnicellular(n, mMinProkaryote, mMaxProkaryote)
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


  end subroutine initProkaryote
  
  subroutine initProkaryoteX(this, n,k,errorio,errorstr)
    use iso_c_binding, only: c_char
    class(spectrumProkaryote):: this
    integer, intent(in):: n,k
    logical(1), intent(out):: errorio 
    character(c_char), dimension(*), intent(out) :: errorstr
    integer:: i
    real(dp), parameter:: rho = 0.4*1d6*1d-12
    real(dp) :: ii
    real(dp) :: mMinProkaryote, mMaxProkaryote
    real(dp) :: alphaN, rNstar  
    real(dp) :: alphaL, rLstar
    real(dp) :: alphaF, cF  
    real(dp) :: cLeakage, delta, alphaJ, cR, palatability  
    
    ! no errors to begin with
    errorio=.false.
    inputfile='../input/input_generalists_simpleX.h'

    print*, 'Loading parameter for prokaryite from ', inputfile, ':'
    call read_inputX(inputfile,'prokaryote','mMinProkaryote',mMinProkaryote,k,errorio,errorstr)
    print*, '    *mMinProkaryote=',mMinProkaryote
    call read_inputX(inputfile,'prokaryote','mMaxProkaryote',mMaxProkaryote,k,errorio,errorstr)
    print*, '    *mMaxProkaryote=',mMaxProkaryote
    call read_inputX(inputfile,'prokaryote','alphaN',alphaN,k,errorio,errorstr)
    print*, '    *alphaN=',alphaN
    call read_inputX(inputfile,'prokaryote','rNstar',rNstar,k,errorio,errorstr)
    print*, '    *rNstar=',rNstar
    call read_inputX(inputfile,'prokaryote','alphaL',alphaL,k,errorio,errorstr)
    print*, '    *alphaL=',alphaL
    call read_inputX(inputfile,'prokaryote','alphaF',alphaF,k,errorio,errorstr)
    print*, '    *alphaF=',alphaF
    call read_inputX(inputfile,'prokaryote','cF',cF,k,errorio,errorstr)
    print*, '    *cF=',cF
    call read_inputX(inputfile,'prokaryote','cLeakage',cLeakage,k,errorio,errorstr)
    print*, '    *cLeakage=',cLeakage
    call read_inputX(inputfile,'prokaryote','alphaJ',alphaJ,k,errorio,errorstr)
    print*, '    *alphaJ=',alphaJ
    call read_inputX(inputfile,'prokaryote','cR',cR,k,errorio,errorstr)
    print*, '    *cR=',cR
    call read_inputX(inputfile,'prokaryote','rLstar',rLstar,k,errorio,errorstr)
    print*, '    *rLstar=',rLstar
    call read_inputX(inputfile,'prokaryote','delta',delta,k,errorio,errorstr)
    print*, '    *delta=',delta
    call read_inputX(inputfile,'prokaryote','epsilonL',epsilonL,k,errorio,errorstr)
    print*, '    *epsilonL=',epsilonL
    call read_inputX(inputfile,'prokaryote','remin2',remin2,k,errorio,errorstr)
    print*, '    *remin2=',remin2
    call read_inputX(inputfile,'prokaryote','reminF',reminF,k,errorio,errorstr)
    print*, '    *reminF=',reminF
    call read_inputX(inputfile,'prokaryote','beta',this%beta,k,errorio,errorstr)
    print*, '    *beta=',this%beta
    call read_inputX(inputfile,'prokaryote','sigma',this%sigma,k,errorio,errorstr)
    print*, '    *sigma=',this%sigma
    call read_inputX(inputfile,'prokaryote','epsilonF',this%epsilonF,k,errorio,errorstr)
    print*, '    *epsilonF=',this%epsilonF
    call read_inputX(inputfile,'prokaryote','palatability',palatability,k,errorio,errorstr)
    
    
    
    call this%initUnicellular(n, mMinProkaryote, mMaxProkaryote)
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
    this%palatability = palatability


  end subroutine initProkaryoteX

  subroutine calcRatesProkaryote(this, L, N, DOC, gammaN, gammaDOC)
    class(spectrumProkaryote), intent(inout):: this
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
end subroutine calcRatesProkaryote

  subroutine calcDerivativesProkaryote(this, u, dNdt, dDOCdt, dudt)
    class(spectrumProkaryote), intent(inout):: this
    real(dp), intent(in):: u(this%n)
    real(dp), intent(inout) :: dNdt, dDOCdt, dudt(this%n)
    !real(dp):: mortloss
    integer:: i
    
    print*, 'uType is:',this%uType
    

    !this%mort2 = 0*u !this%mort2constant*u ! "quadratic" mortality
    this%mort2 = this%mort2constant*this%uType ! "quadratic" mortality
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
      ! Update the Prokaryote:
      !
      dudt(i) = (this%Jtot(i)/this%m(i)  &
           - this%mortpred(i) &
           - this%mort2(i) &
           - this%mortHTL(i))*u(i)
   end do
  
 end subroutine calcDerivativesProkaryote

subroutine printRatesProkaryote(this)
  class(spectrumProkaryote), intent(in):: this

  write(*,*) "Prokaryote with ", this%n, " size classes:"
  call this%printRatesUnicellular()
end subroutine printRatesProkaryote

  function getProdBactProkaryote(this, u) result(ProdBact)
    real(dp):: ProdBact
    class(spectrumProkaryote), intent(in):: this
    real(dp), intent(in):: u(this%n)
    integer:: i

    ProdBact = 0.d0
    do i = 1, this%n
      ProdBact = ProdBact + max(0.d0, this%JDOC(i) - ftemp2*this%Jresp(i))*u(i)/this%m(i)
    enddo

  end function getProdBactProkaryote

end module prokaryote
