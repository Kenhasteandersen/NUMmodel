!
! Module to handle generalist unicellulars
!
module generalists
  use globals
  use spectrum
  use read_input_module
  implicit none

  private 
  
  real(dp) :: epsilonL,bL,bN,bDOC,bF,bg,remin2,reminF

  type, extends(spectrumUnicellular) :: spectrumGeneralists
    real(dp), allocatable :: JFreal(:)
    !
    ! Regulation factors:
    !
    real(dp), allocatable :: dL(:), dN(:), dDOC(:), jNet(:)
    
  contains
    procedure, pass :: initGeneralists
    procedure :: calcRates => calcRatesGeneralists
    procedure :: calcDerivativesGeneralists
    procedure :: printRates => printRatesGeneralists
    procedure :: getProdNet => getProdNetGeneralists
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
    !real(dp), parameter:: rho != 0.4*1d6*1d-12
    real(dp) :: mMinGeneralist, mMaxGeneralist, rho
    real(dp) :: alphaL, rLstar !Light uptake
    real(dp) :: alphaN, rNstar !osmotrophic uptake
    real(dp) :: alphaF, cF !Phagotrophy
    real(dp) :: cLeakage, delta, alphaJ, cR !Metabolism

    ! no errors to begin with
    errorio=.false.
    
    print*, 'Loading parameter for generalist from ', inputfile, ':'
    call read_input(inputfile,'generalists','mMinGeneralist',mMinGeneralist,errorio,errorstr)
    call read_input(inputfile,'generalists','mMaxGeneralist',mMaxGeneralist,errorio,errorstr)
    call this%initUnicellular(n, mMinGeneralist, mMaxGeneralist)
    call read_input(inputfile,'generalists','alphaL',alphaL,errorio,errorstr)
    call read_input(inputfile,'generalists','rLstar',rLstar,errorio,errorstr)
    call read_input(inputfile,'generalists','alphaN',alphaN,errorio,errorstr)
    call read_input(inputfile,'generalists','rNstar',rNstar,errorio,errorstr)
    call read_input(inputfile,'generalists','alphaF',alphaF,errorio,errorstr)
    call read_input(inputfile,'generalists','cF',cF,errorio,errorstr)
    call read_input(inputfile,'generalists','cLeakage',cLeakage,errorio,errorstr)
    call read_input(inputfile,'generalists','delta',delta,errorio,errorstr)
    call read_input(inputfile,'generalists','alphaJ',alphaJ,errorio,errorstr)
    call read_input(inputfile,'generalists','cR',cR,errorio,errorstr)
    call read_input(inputfile,'generalists','epsilonL',epsilonL,errorio,errorstr)
    call read_input(inputfile,'generalists','bL',bL,errorio,errorstr)
    call read_input(inputfile,'generalists','bN',bN,errorio,errorstr)
    call read_input(inputfile,'generalists','bDOC',bDOC,errorio,errorstr)
    call read_input(inputfile,'generalists','bF',bF,errorio,errorstr)
    call read_input(inputfile,'generalists','bg',bg,errorio,errorstr)
    call read_input(inputfile,'generalists','remin2',remin2,errorio,errorstr)
    call read_input(inputfile,'generalists','reminF',reminF,errorio,errorstr)
    call read_input(inputfile,'generalists','rho',rho,errorio,errorstr)
    call read_input(inputfile,'generalists','epsilonF',this%epsilonF,errorio,errorstr)
    call read_input(inputfile,'generalists','beta',this%beta,errorio,errorstr)
    call read_input(inputfile,'generalists','sigma',this%sigma,errorio,errorstr)

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
    this%Jresp =cR*alphaJ*this%m ! decrease basal resp. according to Ken

    allocate(this%dL(n))
    allocate(this%dN(n))
    allocate(this%dDOC(n))
    allocate(this%Jnet(n))

  end subroutine initGeneralists

  subroutine calcRatesGeneralists(this, L, N, DOC, gammaN, gammaDOC)
    class(spectrumGeneralists), intent(inout):: this
    real(dp), intent(in):: gammaN, gammaDOC
    real(dp), intent(in):: L, N, DOC
    real(dp):: f, JmaxT, tmp
    real(dp):: Jnetp(this%n)
    integer:: i

    do i = 1, this%n
       !
       ! Encounters:
       !
       this%JN(i)   = gammaN * fTemp15 * this%AN(i)*N*rhoCN ! Diffusive nutrient uptake in units of C/time
       this%JDOC(i) = gammaDOC * fTemp15 * this%AN(i)*DOC ! Diffusive DOC uptake, units of C/time
       this%JL(i)   = epsilonL * this%AL(i)*L  ! Photoharvesting
       JmaxT = fTemp2*this%Jmax(i)
       !
       ! Potential net uptake
       !
       Jnetp(i) = this%JL(i)*(1-bL) + this%JDOC(i)*(1-bDOC) + this%JF(i)*(1-bF) - ftemp2*this%Jresp(i)
       !
       ! Calculation of down-regulation factors for N-uptake and the net uptake:
       !
       if (Jnetp(i) .lt. 0) then
         this%Jnet(i) = Jnetp(i)  ! Severe carbon limitations => negative growth
         this%dN(i) = 0.d0
       else
         if (this%JN(i) .eq. 0) then
           this%dN(i) = 1.d0
         else
           this%dN(i) = max( 0.d0, min( 1.d0, (Jnetp(i) - this%JF(i)*(bg+1))/(this%JN(i)*(1+bg+bN)) ) )
         endif
         this%Jnet(i) = min( (Jnetp(i)-bg*(this%dN(i)*this%JN(i)))/(1+bg) , & ! Carbon limitation
                            this%JF(i) + this%dN(i)*this%JN(i))              ! N limitation
       endif 
       !
       ! Synthesis limitation:
       !
       if ( this%Jnet(i) .gt. this%JlossPassive(i) ) then ! Apply FR only if net growth is positive
         this%Jnet(i) = JmaxT * this%Jnet(i) / ( this%Jnet(i) + JmaxT )
       endif
       this%Jtot(i) = this%Jnet(i) - this%JlossPassive(i)

       ! Take up N only to the degree that is is not supplied by feeding (ie priotize feeding):
       this%JNreal(i) = max( 0.d0, this%Jnet(i) - this%jF(i) )
       
       !
       ! Regulate carbon uptakes for growth + respiration towards lowered jNet.
       !

       ! First carbon from F assuming no uptakes from DOC and L:
       this%JFreal(i) = min( this%JF(i), &
           (this%Jnet(i) + bg*max(0.d0, this%Jnet(i)) + ftemp2*this%Jresp(i) + bN*this%Jnreal(i) )/(1-bF)) 
       ! Then prioritize DOC:
       tmp = this%Jnet(i) + bg*max(0.d0, this%Jnet(i)) + bN*this%JNreal(i) + ftemp2*this%Jresp(i) - this%JFreal(i)*(1-bF)
       this%jDOCreal(i) = min( this%JDOC(i), tmp/(1-bDOC) )
       ! And finally light:
       this%JLreal(i) = min( this%JL(i), (tmp - this%jDOCreal(i)*(1-bDOC))/(1-bL) )
      
       ! Exude surplus N:
       this%JNlossLiebig(i) = max( 0.d0, this%Jnreal(i) + this%Jfreal(i) - this%JlossPassive(i) - this%Jtot(i) )
       this%JClossLiebig(i) = 0.d0 ! There are never surplus C uptakes
      !        
      ! Actual uptakes:
      !
      this%JNtot(i) = this%JNreal(i) + this%JFreal(i)
      !
      ! Losses:
      !
      this%JCloss_feeding(i)     = (1.-this%epsilonF)/this%epsilonF * this%JFreal(i) ! Incomplete feeding (units of carbon per time)
      this%JCloss_photouptake(i) = (1.-epsilonL)/epsilonL * this%JLreal(i)
      this%Jresptot(i)= &
            fTemp2*this%Jresp(i) + &
            bDOC*this%JDOCreal(i) + &
            bL*this%JLreal(i) + &
            bN*this%JNreal(i) + &
            bF*this%JFreal(i) + &
            max(0.d0,bg*this%Jnet(i))
 
      this%f(i) = f   
    end do

    ! Needed to get the right rates with getRates:
    this%jN = this%jNreal  
    this%jDOC = this%jDOCreal
    this%JF = this%JFreal
    this%JNloss = this%JNlossLiebig
    !
    ! Test for conservation budget. Should be close to zero:
    !
    !write(*,*) 'N budget:',(-this%Jtot+this%JN+this%JFreal & ! Gains
    !  -this%JNlossLiebig-this%JlossPassive)/this%m           ! Losses
    !write(*,*) 'C budget:',(this%JLreal+this%JDOCreal+this%JFreal & ! Gains
    !  -this%Jtot-this%Jresptot - this%JClossLiebig - this%JlossPassive)/this%m   ! Losses
end subroutine calcRatesGeneralists

  subroutine calcDerivativesGeneralists(this, u, dNdt, dDOCdt, dudt)
    class(spectrumGeneralists), intent(inout):: this
    real(dp), intent(in):: u(this%n)
    real(dp), intent(inout) :: dNdt, dDOCdt, dudt(this%n)
    !real(dp):: mortloss
    integer:: i

    this%mort2 = this%mort2constant*u ! "quadratic" mortality
    this%jPOM = (1-remin2)*this%mort2  &! non-remineralized mort2 => POM
              + (1-reminF)*this%JCloss_feeding/this%m ! Feeding losses

    do i = 1, this%n
      !
      ! Update nitrogen:
      !
      dNdt = dNdt  &
           + ((-this%JNreal(i) &
           +  this%JlossPassive(i) &
           +  this%JNlossLiebig(i) &     ! N leakage due to excess food
           +  reminF*this%JCloss_feeding(i))/this%m(i) & ! Remineralized feeding losses
           +  remin2*this%mort2(i) & ! Remineralized viral lysis
           ) * u(i)/rhoCN
      !
      ! Update DOC:
      !
      dDOCdt = dDOCdt &
           + ((-this%JDOCreal(i) &
           +   this%JlossPassive(i) &
           +   this%JCloss_photouptake(i) &
           +   reminF*this%JCloss_feeding(i))/this%m(i) & ! Remineralized feeding losses
           +   remin2*this%mort2(i) & ! Remineralized viral lysis
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
  !write(*,99) "deltaL:", this%JLreal/this%JL
  write(*,99) "deltaN:", this%dN
  !write(*,99) "deltaDOC:", this%JDOCreal/this%JDOC
end subroutine printRatesGeneralists

!
  ! Returns the net primary production calculated as the total amount of carbon fixed
  ! by photsynthesis minus the respiration due to basal respiration,
  ! photosynthesis, nutrint uptake, and growth. Units: mugC/day/m3
  ! (See Andersen and Visser (2023) table 5)
  !
function getProdNetGeneralists(this, u) result(ProdNet)
  real(dp):: ProdNet
  class(spectrumGeneralists), intent(in):: this
  real(dp), intent(in):: u(this%n)
  integer:: i
  real(dp):: resp, tmp, tmp2

  ProdNet = 0.d0
  do i = 1, this%n
    if ( (this%JLreal(i) + this%JDOCreal(i)) .ne. 0.d0 ) then
      tmp = this%JLreal(i) / (this%JLreal(i) + this%JDOCreal(i))
    else
      tmp = 0.d0
    endif

    if ( (this%JLreal(i) + this%JDOCreal(i)+ this%JFreal(i)) .ne. 0.d0 ) then
      tmp2 = this%JLreal(i) / (this%JLreal(i) + this%JDOCreal(i) + this%JFreal(i))
    else
      tmp2 = 0.d0
    endif

    resp = &
      fTemp2*this%Jresp(i) + & ! Basal metabolism
      bL*this%JLreal(i) + &    ! Light uptake metabolism
      bN*this%JNreal(i) * tmp + &  ! The fraction of N uptake that is not associated to DOC uptake  
      bg*this%Jnet(i) * tmp2 ! The fraction of growth not associated with DOC or feeding
    ProdNet = ProdNet + max( 0.d0, (this%JLreal(i) - resp) * u(i)/this%m(i) )

  end do
end function getProdNetGeneralists

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
