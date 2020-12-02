
module NUMmodel
  implicit none
  private
  integer, parameter :: dp=kind(0.d0) ! double precision

  type typeParameters
     integer:: n
     real(dp), dimension(:), allocatable:: m(:), AN(:), AL(:), AF(:), Jmax(:), JFmax(:), Jresp(:)
     real(dp), dimension(:,:), allocatable:: theta(:,:)
     real(dp), dimension(:), allocatable:: mort(:), mortHTL(:)
     real(dp):: rhoCN, epsilonL, epsilonF, mort2, remin, remin2, cLeakage
  end type typeParameters

  type typeRates
     real(dp), dimension(:), allocatable:: JN, JDOC, JL, JF, F, JFreal
     real(dp), dimension(:), allocatable:: JNtot, JLreal, JCtot, Jtot;
     real(dp), dimension(:), allocatable:: JCloss_feeding, JCloss_photouptake, JNlossLiebig, JClossLiebig
     real(dp), dimension(:), allocatable:: JNloss, JCloss
     real(dp), dimension(:), allocatable:: mortpred
  end type typeRates

  type(typeParameters):: p
  type(typeRates):: rates
  integer, parameter:: unitDebug = 10
  real(dp), dimension(:), allocatable:: B, ANmT, JmaxT, JFmaxT, JrespT
 
  public setParameters, calcRates
contains

  subroutine openDebug
    open(unit=unitDebug, FILE='debug.out', status='replace')
    write(unitDebug,*) '------------------------------------'
  end subroutine openDebug
  
  subroutine setParameters(n, m, rhoCN, epsilonL, epsilonF, AN, AL, AF, Jmax, JFmax, Jresp, &
         theta, mort, mort2, mortHTL, remin, remin2, cLeakage)
    integer, intent(in):: n
    real(dp), intent(in), dimension(:):: m(:)
    real(dp), intent(in):: rhoCN, epsilonL, epsilonF, mort2, remin, remin2, cLeakage
    real(dp), intent(in), dimension(:):: AN(:), AL(:), AF(:), Jmax(:), JFmax(:), Jresp(:), mort(:), mortHTL(:)
    real(dp), intent(in), dimension(:,:):: theta(:,:)
    
    allocate(p%m(n))
    allocate(p%AN(n))
    allocate(p%AL(n))
    allocate(p%AF(n))
    allocate(p%Jmax(n))
    allocate(p%JFmax(n))
    allocate(p%Jresp(n))
    allocate(p%theta(n,n))
    allocate(p%mort(n))
    allocate(p%mortHTL(n))

    p%n = n
    p%m = m
    p%rhoCN = rhoCN
    p%epsilonL = epsilonL
    p%epsilonF = epsilonF
    p%AN = AN
    p%AL = AL    
    p%AF = AF
    p%Jmax = Jmax
    p%JFmax = JFmax
    p%Jresp = Jresp
    p%theta = theta
    p%mort = mort
    p%mort2 = mort2
    p%mortHTL = mortHTL
    p%remin = remin
    p%remin2 = remin2
    p%cLeakage = cLeakage
    
    !call openDebug
    !write(unitDebug,*) p%cLeakage    
    !close(unitDebug)
    !
    ! Init private temps:
    !
    allocate(B(n))
    allocate(ANmT(n))
    allocate(JmaxT(n))
    allocate(JFmaxT(n))
    allocate(JrespT(n))
    !
    !  Init rates:
    !
    allocate(rates%JN(n))
    allocate(rates%JDOC(n))
    allocate(rates%JL(n))
    allocate(rates%JF(n))
    allocate(rates%JFreal(n))
    allocate(rates%F(n))
    allocate(rates%JNtot(n))
    allocate(rates%JLreal(n))
    allocate(rates%JCtot(n))
    allocate(rates%Jtot(n))
    allocate(rates%JCloss_feeding(n))
    allocate(rates%JCloss_photouptake(n))
    allocate(rates%JNlossLiebig(n))
    allocate(rates%JClossLiebig(n))
    allocate(rates%JNloss(n))
    allocate(rates%JCloss(n))
    allocate(rates%mortpred(n))

  end subroutine setParameters
    
  function calcRates(T, L, u, dudt, gammaN, gammaDOC) result(dudt)
    real(dp), intent(in):: u(:)
    real(dp):: dudt(size(u))

    dudt = 0.01*u
  end function calcRates

  void calcRates(const double& T, const double& L, const double* u, 
               double* dudt, const double gammaN, const double gammaDOC) {
  int i;
  double f;
  static double Tsaved;
  static double f15, f20;
  /*
   * Make sure values of B and DOC are positive:
   */
  for (i=0; i<p.n; i++)
    B[i] = max(0, u[idxB+i]); 
  double DOC = max(0, u[idxDOC]);
  double N = max(0, u[idxN]);

end module NUMmodel


