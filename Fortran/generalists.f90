!
! Module to handle generalist unicellulars
!
module generalists
  use globals
  use spectrum
  use read_input_module
  implicit none

  private 
  
  real(dp) :: bL,bN,bDOC,bF,bg,remin2,reminF

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
   ! procedure :: getProdNet => getProdNetGeneralists
   ! procedure :: getProdBact => getProdBactGeneralists 
  end type spectrumGeneralists
 
  public initGeneralists, spectrumGeneralists, calcRatesGeneralists, calcDerivativesGeneralists
  public printRatesGeneralists
  public ratesGeneralistsCore, derivativesGeneralistsCore, getParametersGeneralists

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
    real(dp) :: alphaL, rLstar, mUpperAlphaL !Light uptake
    real(dp) :: alphaN, rNstar !osmotrophic uptake
    real(dp) :: alphaF, cF !Phagotrophy
    real(dp) :: cLeakage, delta, alphaJ, cR !Metabolism

    ! no errors to begin with
    errorio=.false.
    
    call read_input(inputfile,'generalists','mMinGeneralist',mMinGeneralist,errorio,errorstr)
    call read_input(inputfile,'generalists','mMaxGeneralist',mMaxGeneralist,errorio,errorstr)
    call this%initUnicellular(n, mMinGeneralist, mMaxGeneralist)
    call read_input(inputfile,'generalists','alphaL',alphaL,errorio,errorstr)
    call read_input(inputfile,'generalists','mUpperAlphaL',mUpperAlphaL,errorio,errorstr)
    call read_input(inputfile,'generalists','rLstar',rLstar,errorio,errorstr)
    call read_input(inputfile,'generalists','alphaN',alphaN,errorio,errorstr)
    call read_input(inputfile,'generalists','rNstar',rNstar,errorio,errorstr)
    call read_input(inputfile,'generalists','alphaF',alphaF,errorio,errorstr)
    call read_input(inputfile,'generalists','cF',cF,errorio,errorstr)
    call read_input(inputfile,'generalists','cLeakage',cLeakage,errorio,errorstr)
    call read_input(inputfile,'generalists','delta',delta,errorio,errorstr)
    call read_input(inputfile,'generalists','alphaJ',alphaJ,errorio,errorstr)
    call read_input(inputfile,'generalists','cR',cR,errorio,errorstr)
    call read_input(inputfile,'generalists','epsilonL',this%epsilonL,errorio,errorstr)
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
    do i = 1, n
      if (this%m(i) .gt. mUpperAlphaL) then
        this%AL(i) = 0.d0
      end if
    end do
    this%AF = alphaF*this%m
    this%JFmax = cF/this%r * this%m
    
    this%JlossPassive = cLeakage/this%r * this%m ! in units of C

    this%Jmax = alphaJ * this%m * (1.d0-this%nu) ! mugC/day
    this%Jresp =cR*alphaJ*this%m 

    allocate(this%dL(n))
    allocate(this%dN(n))
    allocate(this%dDOC(n))
    allocate(this%Jnet(n))

  end subroutine initGeneralists

  subroutine calcRatesGeneralists(this, L, N, DOC, gammaN, gammaDOC)
    class(spectrumGeneralists), intent(inout):: this
    real(dp), intent(in):: gammaN, gammaDOC
    real(dp), intent(in):: L, N, DOC

    call ratesGeneralistsCore(L, N, DOC, gammaN, gammaDOC, fTemp2, fTemp15, rhoCN, &
         bL, bN, bDOC, bF, bg, this%epsilonL, this%epsilonF, &
         this%AN, this%AL, this%Jmax, this%Jresp, this%JlossPassive, this%JF, &
         this%JN, this%JDOC, this%JL, this%Jnet, this%dN, this%f, this%Jtot, &
         this%JNreal, this%JFreal, this%JDOCreal, this%JLreal, &
         this%JNlossLiebig, this%JClossLiebig, this%JNtot, &
         this%JCloss_feeding, this%JCloss_photouptake, this%Jresptot)

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

  ! -----------------------------------------------
  ! Rates of one size class of generalists.
  ! The physiology is here, separated from the object, so that the same code
  ! is used by the object-oriented library (called with arrays over size
  ! classes) and by the flat GPU kernel in NUMmodel_offload (called with scalars).
  ! -----------------------------------------------
  elemental subroutine ratesGeneralistsCore(L, N, DOC, gammaN, gammaDOC, fTemp2, fTemp15, rhoCN, &
       bL, bN, bDOC, bF, bg, epsilonL, epsilonF, &
       AN, AL, Jmax, Jresp, JlossPassive, JF, &
       JN, JDOC, JL, Jnet, dN, f, Jtot, JNreal, JFreal, JDOCreal, JLreal, &
       JNlossLiebig, JClossLiebig, JNtot, JCloss_feeding, JCloss_photouptake, Jresptot)
    !$omp declare target
    real(dp), intent(in):: L, N, DOC, gammaN, gammaDOC, fTemp2, fTemp15, rhoCN
    real(dp), intent(in):: bL, bN, bDOC, bF, bg, epsilonL, epsilonF
    real(dp), intent(in):: AN, AL, Jmax, Jresp, JlossPassive
    real(dp), intent(in):: JF ! Available food (from calcFeeding)
    real(dp), intent(out):: JN, JDOC, JL, Jnet, dN, f, Jtot, JNreal, JFreal, JDOCreal, JLreal
    real(dp), intent(out):: JNlossLiebig, JClossLiebig, JNtot, JCloss_feeding, JCloss_photouptake, Jresptot
    real(dp):: JmaxT, Jnetp, tmp
    !
    ! Encounters:
    !
    JN   = gammaN * fTemp15 * AN*N*rhoCN ! Diffusive nutrient uptake in units of C/time
    JDOC = gammaDOC * fTemp15 * AN*DOC ! Diffusive DOC uptake, units of C/time
    JL   = epsilonL * AL*L  ! Photoharvesting
    JmaxT = fTemp2*Jmax
    !
    ! Potential net uptake
    !
    Jnetp = JL*(1-bL) + JDOC*(1-bDOC) + JF*(1-bF) - ftemp2*Jresp
    !
    ! Calculation of down-regulation factors for N-uptake and the net uptake:
    !
    if (Jnetp .lt. 0) then
      Jnet = Jnetp  ! Severe carbon limitations => negative growth
      dN = 0.d0
    else
      if (JN .eq. 0) then
        dN = 1.d0
      else
        dN = max( 0.d0, min( 1.d0, (Jnetp - JF*(bg+1))/(JN*(1+bg+bN)) ) )
      endif
      Jnet = min( (Jnetp-bN*(dN*JN))/(1+bg) , & ! Carbon limitation
                  JF + dN*JN)                   ! N limitation
    endif 
    !
    ! Synthesis limitation:
    !
    f = 0
    if ( Jnet .gt. JlossPassive ) then ! Apply FR only if net growth is positive
      f = Jnet / ( Jnet + JmaxT )
      Jnet = JmaxT * f
    endif
    Jtot = Jnet - JlossPassive

    ! Take up N only to the degree that is is not supplied by feeding:
    JNreal = max( 0.d0, Jnet - JF )
    !
    ! Regulate carbon uptakes for growth + respiration towards lowered jNet.
    !

    ! First carbon from F assuming no uptakes from DOC and L:
    JFreal = min( JF, &
        (Jnet + bg*max(0.d0, Jnet) + ftemp2*Jresp + bN*JNreal )/(1-bF)) 
    ! Then divide evenly btw DOC and L:
    tmp = ( (1 - bDOC)*JDOC + (1 - bL)*JL  )
    if (tmp .eq. 0.0d0) then
      JDOCreal = 0.0d0
      JLreal = 0.0d0
    else
      tmp = ( Jnet + bg*max(0.d0, Jnet) + bN*JNreal + ftemp2*Jresp - &
             JFreal*(1 - bF) ) / tmp
      JDOCreal = tmp * JDOC
      JLreal = tmp * JL
    endif       
    ! Exude surplus N:
    JNlossLiebig = max( 0.d0, JNreal + JFreal - JlossPassive - Jtot )
    JClossLiebig = 0.d0 ! There are never surplus C uptakes
    !        
    ! Actual uptakes:
    !
    JNtot = JNreal + JFreal
    !
    ! Losses:
    !
    JCloss_feeding     = (1.-epsilonF)/epsilonF * JFreal ! Incomplete feeding (units of carbon per time)
    JCloss_photouptake = (1.-epsilonL)/epsilonL * JLreal
    Jresptot = &
         fTemp2*Jresp + &
         bDOC*JDOCreal + &
         bL*JLreal + &
         bN*JNreal + &
         bF*JFreal + &
         max(0.d0,bg*Jnet)
  end subroutine ratesGeneralistsCore

  subroutine calcDerivativesGeneralists(this, u, dNdt, dDOCdt, dudt)
    class(spectrumGeneralists), intent(inout):: this
    real(dp), intent(in):: u(this%n)
    real(dp), intent(inout) :: dNdt, dDOCdt, dudt(this%n)
    real(dp):: dNcontrib(this%n), dDOCcontrib(this%n)
    integer:: i

    call derivativesGeneralistsCore(u, this%m, this%mort2constant, remin2, reminF, rhoCN, &
         this%JNreal, this%JlossPassive, this%JNlossLiebig, this%JCloss_feeding, &
         this%JDOCreal, this%JCloss_photouptake, this%Jtot, this%mortpred, this%mortHTL, &
         this%mort2, this%jPOM, dNcontrib, dDOCcontrib, dudt)
    do i = 1, this%n
      dNdt = dNdt + dNcontrib(i)
      dDOCdt = dDOCdt + dDOCcontrib(i)
    end do
  end subroutine calcDerivativesGeneralists

  ! -----------------------------------------------
  ! Derivatives of one size class of generalists and its contributions
  ! to the derivatives of N and DOC (shared with the flat GPU kernel).
  ! -----------------------------------------------
  elemental subroutine derivativesGeneralistsCore(u, m, mort2constant, remin2, reminF, rhoCN, &
       JNreal, JlossPassive, JNlossLiebig, JCloss_feeding, JDOCreal, JCloss_photouptake, &
       Jtot, mortpred, mortHTL, mort2, jPOM, dNcontrib, dDOCcontrib, dudt)
    !$omp declare target
    real(dp), intent(in):: u, m, mort2constant, remin2, reminF, rhoCN
    real(dp), intent(in):: JNreal, JlossPassive, JNlossLiebig, JCloss_feeding
    real(dp), intent(in):: JDOCreal, JCloss_photouptake, Jtot, mortpred, mortHTL
    real(dp), intent(out):: mort2, jPOM, dNcontrib, dDOCcontrib, dudt

    mort2 = mort2constant*u ! "quadratic" mortality
    jPOM = (1-remin2)*mort2  &! non-remineralized mort2 => POM
         + (1-reminF)*JCloss_feeding/m ! Feeding losses
    !
    ! Nitrogen:
    !
    dNcontrib = ((-JNreal &
         +  JlossPassive &
         +  JNlossLiebig &     ! N leakage due to excess food
         +  reminF*JCloss_feeding)/m & ! Remineralized feeding losses
         +  remin2*mort2 & ! Remineralized viral lysis
         ) * u/rhoCN
    !
    ! DOC:
    !
    dDOCcontrib = ((-JDOCreal &
         +   JlossPassive &
         +   JCloss_photouptake &
         +   reminF*JCloss_feeding)/m & ! Remineralized feeding losses
         +   remin2*mort2 & ! Remineralized viral lysis
         ) * u
    !
    ! The generalists:
    !
    dudt = (Jtot/m  &
         - mortpred &
         - mort2 &
         - mortHTL)*u
  end subroutine derivativesGeneralistsCore

  ! -----------------------------------------------
  ! Access to the module parameters (used by NUMmodel_offload)
  ! -----------------------------------------------
  subroutine getParametersGeneralists(bL_, bN_, bDOC_, bF_, bg_, remin2_, reminF_)
    real(dp), intent(out):: bL_, bN_, bDOC_, bF_, bg_, remin2_, reminF_
    bL_ = bL; bN_ = bN; bDOC_ = bDOC; bF_ = bF; bg_ = bg
    remin2_ = remin2; reminF_ = reminF
  end subroutine getParametersGeneralists

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
! function getProdNetGeneralists(this, u) result(ProdNet)
!   real(dp):: ProdNet
!   class(spectrumGeneralists), intent(in):: this
!   real(dp), intent(in):: u(this%n)
!   integer:: i
!   real(dp):: resp, tmp, tmp2

!   ProdNet = 0.d0
!   do i = 1, this%n
!     if ( (this%JLreal(i) + this%JDOCreal(i)) .ne. 0.d0 ) then
!       tmp = this%JLreal(i) / (this%JLreal(i) + this%JDOCreal(i))
!     else
!       tmp = 0.d0
!     endif

!     if ( (this%JLreal(i) + this%JDOCreal(i)+ this%JFreal(i)) .ne. 0.d0 ) then
!       tmp2 = this%JLreal(i) / (this%JLreal(i) + this%JDOCreal(i) + this%JFreal(i))
!     else
!       tmp2 = 0.d0
!     endif

!     resp = &
!       fTemp2*this%Jresp(i) + & ! Basal metabolism
!       bL*this%JLreal(i) + &    ! Light uptake metabolism
!       bN*this%JNreal(i) * tmp + &  ! The fraction of N uptake that is not associated to DOC uptake  
!       bg*this%Jnet(i) * tmp2 ! The fraction of growth not associated with DOC or feeding
!     ProdNet = ProdNet + max( 0.d0, (this%JLreal(i) - resp) * u(i)/this%m(i) )

!   end do
! end function getProdNetGeneralists

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
