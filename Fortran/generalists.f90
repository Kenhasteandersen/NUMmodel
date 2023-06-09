!
! Module to handle generalist unicellulars
!
module generalists
  use globals
  use spectrum
  use read_input_module
  use read_input_module2
  implicit none

  private 

  type, extends(spectrumUnicellular) :: spectrumGeneralists
    real(dp), allocatable :: JFreal(:)
    !
    ! Regulation factors:
    !
    real(dp), allocatable :: dL(:), dN(:), dDOC(:)
    
  contains
    procedure, pass :: initGeneralists
    procedure :: calcRates => calcRatesGeneralists
    procedure :: calcDerivativesGeneralists
    procedure :: printRates => printRatesGeneralists
    procedure :: getProdBact => getProdBactGeneralists 
  end type spectrumGeneralists
 
  public initGeneralists, spectrumGeneralists, calcRatesGeneralists, calcDerivativesGeneralists
  public printRatesGeneralists

contains
  
  subroutine initGeneralists(this, n,errorio,errorstr)
    use iso_c_binding, only: c_char
    class(spectrumGeneralists):: this
    integer, intent(in):: n
    logical(1), intent(out):: errorio 
    character(c_char), dimension(*), intent(out) :: errorstr
    integer:: i
    real(dp), parameter:: mMin = 3.1623d-9 !2.78d-8!
    real(dp), parameter:: rho = 0.4*1d6*1d-12
    real(dp) :: mMinGeneralist, mMaxGeneralist
    real(dp) :: alphaL, rLstar !Light uptake
    real(dp) :: alphaN, rNstar !osmotrophic uptake
    real(dp) :: epsilonF, alphaF, cF, beta, sigma !Phagotrophy
    real(dp) :: cLeakage, delta, alphaJ, cR !Metabolism
    real(dp) ::remin2, reminF 
    
    ! no errors to begin with
    errorio=.false.
    
    print*, 'Loading parameter for generalist from ', inputfile, ':'
    call read_input2(inputfile,'generalists','mMinGeneralist',mMinGeneralist,errorio,errorstr)
    call read_input2(inputfile,'generalists','mMaxGeneralist',mMaxGeneralist,errorio,errorstr)
    call read_input2(inputfile,'generalists','alphaL',alphaL,errorio,errorstr)
    call read_input2(inputfile,'generalists','rLstar',rLstar,errorio,errorstr)
    call read_input2(inputfile,'generalists','alphaN',alphaN,errorio,errorstr)
    call read_input2(inputfile,'generalists','rNstar',rNstar,errorio,errorstr)
    call read_input2(inputfile,'generalists','alphaF',alphaF,errorio,errorstr)
    call read_input2(inputfile,'generalists','cF',cF,errorio,errorstr)
    call read_input2(inputfile,'generalists','cLeakage',cLeakage,errorio,errorstr)
    call read_input2(inputfile,'generalists','delta',delta,errorio,errorstr)
    call read_input2(inputfile,'generalists','alphaJ',alphaJ,errorio,errorstr)
    call read_input2(inputfile,'generalists','cR',cR,errorio,errorstr)
   
    call read_input2(inputfile,'generalists','epsilonL',this%epsilonL,errorio,errorstr)
    call read_input2(inputfile,'generalists','bL',this%bL,errorio,errorstr)
    call read_input2(inputfile,'generalists','bN',this%bN,errorio,errorstr)
    call read_input2(inputfile,'generalists','bDOC',this%bDOC,errorio,errorstr)
    call read_input2(inputfile,'generalists','epsilonF',this%epsilonF,errorio,errorstr)
    call read_input2(inputfile,'generalists','beta',this%beta,errorio,errorstr)
    call read_input2(inputfile,'generalists','sigma',this%sigma,errorio,errorstr)
    call read_input2(inputfile,'generalists','bF',this%bF,errorio,errorstr)
    call read_input2(inputfile,'generalists','bg',this%bg,errorio,errorstr)
    call read_input2(inputfile,'generalists_simple','remin2',this%remin2,errorio,errorstr)
    call read_input2(inputfile,'generalists_simple','reminF',this%reminF,errorio,errorstr)
    
    !call read_input(inputfile,'generalists')
    call this%initUnicellular(n, mMinGeneralist, mMaxGeneralist)
    allocate(this%JFreal(n))

    this%beta = beta
    this%sigma = sigma
    !this%epsilonF = epsilonF

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
    this%Jresp =cR*alphaJ*this%m ! decrease basal resp. according to Ken

    allocate(this%dL(n))
    allocate(this%dN(n))
    allocate(this%dDOC(n))

  end subroutine initGeneralists

  subroutine calcRatesGeneralists(this, L, N, DOC, gammaN, gammaDOC)
    class(spectrumGeneralists), intent(inout):: this
    real(dp), intent(in):: gammaN, gammaDOC
    real(dp), intent(in):: L, N, DOC
    real(dp):: f, JmaxT
    real(dp):: Jnetp(this%n), Jnet(this%n),Jlim(this%n)
    integer:: i

    do i = 1, this%n
       !
       ! Uptakes
       !
       this%JN(i) = gammaN * fTemp15 * this%AN(i)*N*rhoCN ! Diffusive nutrient uptake in units of C/time
       this%JDOC(i) =gammaDOC * fTemp15 * this%AN(i)*DOC ! Diffusive DOC uptake, units of C/time
       this%JL(i) =  epsilonL * this%AL(i)*L  ! Photoharvesting
       JmaxT = fTemp2*this%Jmax(i)
       !write(*,*) this%JL(i)/this%m(i), L, this%AL(i)
       !this%JF(i)= 0.
       !
       ! Potential net uptake
       !
       Jnetp(i)=this%JL(i)*(1-this%bL)+this%JDOC(i)*(1-this%bDOC)+this%JF(i)*(1-this%bF)-ftemp2*this%Jresp(i) ! I think we don't need it
       !
       ! Calculation of down-regulation factors for
       ! N-uptake  
       this%dN(i) = min(1., 1./this%JN(i)*(Jnetp(i)-this%JF(i)*(1+this%bg))/(1+this%bg+this%bN))
       ! If synthesis-limited then down-regulate feeding:
       ! Photosynthesis
       if (this%JL(i) .gt. 0) then ! Needed to avoid the risk of division with zero if JL = 0
        this%dL(i) = min(1.,1./(this%JL(i)*(1-this%bL))*((this%JN(i)+this%JF(i))*(1+this%bg)-this%JDOC(i)*(1-this%bDOC)&
             -this%JF(i)*(1-this%bF)+ftemp2*this%Jresp(i)+this%bN*this%JN(i))) !+this%bN*this%JN(i)
       else
        this%dL(i) = -1.
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
       if (this%dL(i).lt.0) then
        this%dL(i)=0
        this%dDOC(i) =min(1.,1/(this%JDOC(i)*(1-this%bDOC))*((this%JN(i)+this%JF(i))*(1+this%bg)-this%JF(i)*(1-this%bF)&
                + this%bN*this%JN(i)+ftemp2*this%Jresp(i)))
       else 
        this%dDOC(i)=1.
       end if
       !write(*,*) dN(i)
       if (this%dN(i).lt.0) then ! check if N leaks out of the cell
        this%dN(i)=0.
        Jnet =  1./(1+this%bg)*(this%dDOC(i)*this%JDOC(i)*(1.-this%bDOC)+this%dL(i)*this%JL(i)*(1-this%bL) &
        + this%JF(i)*(1-this%bF)- fTemp2*this%Jresp(i) -this%bN*this%dN(i)*this%JN(i)) ! was with max(O.,..) 
        f = (Jnet(i) )/(Jnet(i) + JmaxT)
        !if (f .gt. 0) then

        this%JNlossLiebig(i) = (1-f)*this%JF(i)-f*JmaxT!-1/(1+this%bg)*(this%bN*dN(i)*this%JN(i))!this%JF(i)-!this%JF(i)-Jnet(i)!f*JmaxT
        !endif
       else 
        Jnet =  1./(1+this%bg)*(this%dDOC(i)*this%JDOC(i)*(1.-this%bDOC)+this%dL(i)*this%JL(i)*(1-this%bL) &
        + this%JF(i)*(1-this%bF)- fTemp2*this%Jresp(i) -this%bN*this%dN(i)*this%JN(i)) ! was with max(O.,..) 
        f = (Jnet(i) )/(Jnet(i) + JmaxT)
        this%JNlossLiebig(i) = 0
        end if
       !Jnet =  1./(1+this%bg)*(dDOC(i)*this%JDOC(i)*(1.-this%bDOC)+dL(i)*this%JL(i)*(1-this%bL) &
       ! + this%JF(i)*(1-this%bF)- fTemp2*this%Jresp(i) -this%bN*dN(i)*this%JN(i)) ! was with max(O.,..) 
       !
       ! Saturation of net growth
       !
       !f = (Jnet(i) )/(Jnet(i) + JmaxT)
       if ((Jnet(i) + JmaxT).eq.0) then
        f=0.
       end if
        this%JCtot(i) = & 
        (1-f)*(this%dDOC(i)*this%JDOC(i)+this%dL(i)*this%JL(i) + this%JF(i) )&
        - (1-f)*fTemp2*this%Jresp(i) &
        - ( (1-f)&
        *(this%bDOC*this%dDOC(i)*this%JDOC(i)+this%dL(i)*this%JL(i)*this%bL+ this%JF(i)*this%bF+this%bN*this%dN(i)*this%JN(i))&
        +(1-f)*this%bg*Jnet(i) )
      !
       ! Apply saturation to uptake rates
       !
       this%JNreal(i)=this%dN(i)*(1-f)*this%JN(i)
       this%JDOCreal(i)=this%dDOC(i)*(1-f)*this%JDOC(i)
       this%JLreal(i)=this%dL(i)*(1-f)*this%JL(i)
       this%JFreal(i)=(1-f)*this%JF(i)
       this%Jtot(i)= f * JmaxT-(1-f)*( this%JlossPassive(i))
      !        
      ! Actual uptakes:
      !
      this%JNtot(i) = &
            this%JNreal(i) + &
            this%JFreal(i)
      !write(*,*) f,1-f,Jnet(i)   
      !write(*,*) dN(i),dL(i),this%JNlossLiebig(i)
      !
      ! Losses:
      !
      this%JCloss_feeding(i) = (1.-epsilonF)/epsilonF*this%JFreal(i) ! Incomplete feeding (units of carbon per time)
      this%JCloss_photouptake(i) = (1.-epsilonL)/epsilonL * this%JLreal(i)
      this%Jresptot(i)= (1-f)*fTemp2*this%Jresp(i) + &
            (1-f)*(this%bDOC*this%dDOC(i)*this%JDOC(i) + &
                   this%bL*this%dL(i)*this%JL(i) + &
                   this%bN*this%dN(i)*this%JN(i) + &
                   this%bF*this%JF(i) + &
                   this%bg*Jnet(i))

      !write(*,*) this%Jresptot(i)                 
      ! 
      !-------------------------------------------------------
      ! Test for conservation budget. Should be close to zero:
      !-------------------------------------------------------
      !write(*,*) 'N budget', i,':',(this%JNreal(i)+this%JFreal(i)-this%JNlossLiebig(i)-(1-f)*this%JlossPassive(i) &
      !- this%Jtot(i))/this%m(i)

      !write(*,*) 'C budget', i,':',(this%JCtot(i) -(1-f)*this%JlossPassive(i)&
      !- this%Jtot(i))/this%m(i) !this works only if we take the negative values of jnet
      !this%JF(i) = this%JFreal(i)

      this%f(i)=f
      !this%JF(i) = this%JFreal(i)
    end do
    this%jN = this%jNreal  ! Needed to get the the actual uptakes out with "getRates"
    this%jDOC = this%jDOCreal
    this%JF = this%JFreal
    !
    ! Test for conservation budget. Should be close to zero:
    !
    !write(*,*) 'N budget:',(-this%Jtot+this%JN+this%JFreal & ! Gains
    !  -this%JNlossLiebig-this%JlossPassive)/this%m           ! Losses
    !write(*,*) 'C budget:',(-this%Jtot+this%JLreal+this%JDOC+this%JFreal & ! Gains
    !  -fTemp2*this%Jresp - this%JClossLiebig - this%JlossPassive)/this%m   ! Losses
end subroutine calcRatesGeneralists

  subroutine calcDerivativesGeneralists(this, u, dNdt, dDOCdt, dudt)
    class(spectrumGeneralists), intent(inout):: this
    real(dp), intent(in):: u(this%n)
    real(dp), intent(inout) :: dNdt, dDOCdt, dudt(this%n)
    !real(dp):: mortloss
    integer:: i

    this%mort2 = this%mort2constant*u ! "quadratic" mortality
    this%jPOM = (1-remin2)*this%mort2  &! non-remineralized mort2 => POM
      + (1-reminF)*this%JCloss_feeding/this%m

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
           ) * u(i)/rhoCN
      !
      ! Update DOC:
      !
      dDOCdt = dDOCdt &
           + ((-this%JDOCreal(i) &
           + (1-this%f(i))*this%JlossPassive(i) &
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
  
 end subroutine calcDerivativesGeneralists

subroutine printRatesGeneralists(this)
  class(spectrumGeneralists), intent(in):: this

  write(*,*) "Generalists with ", this%n, " size classes:"
  call this%printRatesUnicellular()

  99 format (a10, 20f10.6)
  
  write(*,99) "jFreal:", this%JFreal / this%m
  write(*,99) "jResptot:", this%Jresptot / this%m
  write(*,99) "deltaL:", this%dL
  write(*,99) "deltaN:", this%dN
  write(*,99) "deltaDOC:", this%dDOC
end subroutine printRatesGeneralists

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
