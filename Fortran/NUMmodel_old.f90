
module NUMmodel
  implicit none
  private
  integer, parameter :: dp=kind(0.d0) ! double precision
  integer, parameter :: idxN = 1
  integer, parameter :: idxDOC = 2
  integer, parameter :: idxB = 3
  
  type typeParameters
     integer:: n
     real(dp), dimension(:), allocatable:: m(:), AN(:), AL(:), AF(:), Jmax(:), JFmax(:), Jresp(:), JlossPassive(:)
     real(dp), dimension(:,:), allocatable:: theta(:,:)
     real(dp), dimension(:), allocatable:: mort(:), mortHTL(:)
     real(dp):: rhoCN, epsilonL, epsilonF, mort2, remin, remin2, cLeakage
  end type typeParameters


  type(typeParameters):: p
  type(typeRates):: rates

  real(dp), dimension(:), allocatable:: B, ANmT, JmaxT, JFmaxT, JrespT
 
  public setParameters, calcRates
contains

  subroutine openDebug
    open(unit=unitDebug, FILE='debug.out', status='replace')
    write(unitDebug,*) '------------------------------------'
  end subroutine openDebug
  
  subroutine setParameters(n, m, rhoCN, epsilonL, epsilonF, AN, AL, AF, Jmax, JFmax, Jresp, JlossPassive, &
         theta, mort, mort2, mortHTL, remin, remin2, cLeakage)
    integer, intent(in):: n
    real(dp), intent(in), dimension(:):: m(:)
    real(dp), intent(in):: rhoCN, epsilonL, epsilonF, mort2, remin, remin2, cLeakage
    real(dp), intent(in), dimension(:):: AN(:), AL(:), AF(:), Jmax(:), JFmax(:), Jresp(:), JlossPassive(:), mort(:), mortHTL(:)
    real(dp), intent(in), dimension(:,:):: theta(:,:)
    
    allocate(p%m(n))
    allocate(p%AN(n))
    allocate(p%AL(n))
    allocate(p%AF(n))
    allocate(p%Jmax(n))
    allocate(p%JFmax(n))
    allocate(p%Jresp(n))
    allocate(p%JlossPassive(n))
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
    p%JlossPassive = JlossPassive
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

  subroutine printRates
    call openDebug
    !write(unitDebug, *) 'N: ', N
    !write(unitDebug, *) 'DOC: ', DOC
    write(unitDebug, *) 'B: ', B
    write(unitDebug, *) 'JN/m: ', rates%Jn/p%m
    write(unitDebug, *) 'JL/m: ', rates%JL/p%m
    write(unitDebug, *) 'F/m: ', rates%F/p%m
    write(unitDebug, *) 'JF/m: ', rates%JF/p%m
    write(unitDebug, *) 'JFreal/m: ', rates%JFreal/p%m
    write(unitDebug, *) 'Jtot/m: ', rates%Jtot/p%m
    write(unitDebug, *) 'mortpred: ', rates%mortpred
    
    close(unitDebug)
  end subroutine printRates

  subroutine printU(u)
    real(dp), intent(in):: u(:)
    call openDebug
    write(unitDebug, *) 'u: ', u
    close(unitDebug)
  end subroutine printU

  function fTemp(Q10, T) result(f)
    real(dp), intent(in), value:: Q10, T
    real(dp):: f

    f = Q10**(T/10.-1.)
  end function fTemp
  
  function calcRates(T, L, u, gammaN, gammaDOC) result(dudt)
    real(dp), intent(in):: T, L, gammaN, gammaDOC, u(:)
    real(dp):: dudt(size(u))
    real(dp):: N, DOC
    real(dp):: f, mortloss
    real(dp), save:: Tsaved, f15, f20
    integer:: i, j
    !
    ! Make sure values are positive
    !
    N = max(0., u(idxN))
    DOC = max(0., u(idxDOC))
    do i = 1, p%n
       B(i) = max(0., u(idxB+i-1))
    end do
    !
    ! Temperature corrections:
    !
    if (T .ne. Tsaved) then
       f15 = fTemp(1.5d0, T);
       f20 = fTemp(2.0d0, T);
       Tsaved = T
       do i=1, p%n
          ANmT(i) = p%AN(i)*f15;
          JmaxT(i) = p%Jmax(i)*f20;
          JrespT(i) = p%Jresp(i)*f20;
          JFmaxT(i) = p%JFmax(i)*f20;
       end do
    end if
    !
    ! Uptakes
    !
    do i=1, p%n
       !
       ! Uptakes:
       ! 
       
       ! Nutrients:
       rates%JN(i) = gammaN * ANmT(i)*N*p%rhoCN    
       ! DOC:
       rates%JDOC(i) = gammaDOC * ANmT(i)*DOC
       ! Light:
       rates%JL(i) = p%epsilonL * p%AL(i)*L
       ! Feeding:
       rates%F(i) = 0.d0
       do j = 1, p%n
          rates%F(i) = rates%F(i) + p%theta(i,j)*B(j)
       enddo
       rates%JF(i) = p%epsilonF * JFmaxT(i) * p%AF(i)*rates%F(i) &
            / (p%AF(i)*rates%F(i) + JFmaxT(i))
    
       ! Combined N and C fluxes:
       rates%JNtot(i) = rates%JN(i) + rates%JF(i) - p%JlossPassive(i)
       rates%JCtot(i) = rates%JL(i) + rates%JF(i) + rates%JDOC(i)  &
            - JrespT(i) - p%JlossPassive(i)
    
       ! Synthesis:
       rates%Jtot(i) = min( rates%JNtot(i), rates%JCtot(i) ) 
       f = rates%Jtot(i)/(rates%Jtot(i) + JmaxT(i)) ! Feeding level
    
       !If synthesis-limited then down-regulate feeding:
       if (rates%Jtot(i) .gt. 0) then
          rates%JFreal(i) = max(0.d0, rates%JF(i) - (rates%Jtot(i)-f*JmaxT(i)))
       else
          rates%JFreal(i) = max(0.d0, rates%JF(i))
       end if
       rates%Jtot(i) = f*JmaxT(i)
       rates%JLreal(i) = rates%JL(i) & 
            - max(0.d0, min((rates%JCtot(i) - (rates%JF(i)-rates%JFreal(i))-rates%Jtot(i)), rates%JL(i)))
    
       ! Actual uptakes:
       rates%JCtot(i) = rates%JLreal(i) + rates%JDOC(i) &
            + rates%JFreal(i) - p%Jresp(i) - p%JlossPassive(i) 
       rates%JNtot(i) = rates%JN(i) + rates%JFreal(i) - p%JlossPassive(i)
    
       ! Losses:
       rates%JCloss_feeding(i) = (1-p%epsilonF)/p%epsilonF*rates%JFreal(i)
       rates%JCloss_photouptake(i) = (1-p%epsilonL)/p%epsilonL*rates%JLreal(i)
       rates%JNlossLiebig(i) = max(0.d0, rates%JNtot(i)-rates%Jtot(i))
       rates%JClossLiebig(i) = max(0.d0, rates%JCtot(i)-rates%Jtot(i)) 
    
       rates%JNloss(i) = rates%JCloss_feeding(i) &
            +rates%JNlossLiebig(i) &
            +p%JlossPassive(i) 
       rates%JCloss(i) = rates%JCloss_feeding(i) &
            + rates%JCloss_photouptake(i) &
            + rates%JClossLiebig(i) &
            + p%JlossPassive(i)
    enddo
    !
    ! Mortality:
    !
    do i=1, p%n
       rates%mortpred(i) = 0
       do j=1, p%n
          if (rates%F(j) .ne. 0) then
             rates%mortpred(i) = rates%mortpred(i) + &
                  p%theta(j,i) * rates%JFreal(j)*B(j)/(p%epsilonF*p%m(j)*rates%F(j))
          end if
       end do
    end do
    !
    ! Reaction terms:
    !
    dudt(idxN) = 0
    dudt(idxDOC) = 0
    do i = 1, p%n
       mortloss = B(i)*((1-p%remin2)*p%mort2*B(i) + p%mortHTL(i))
       dudt(idxN) = dudt(idxN) + &
            ((-rates%JN(i) &
            + rates%JNloss(i))*B(i)/p%m(i) &
            + p%remin2*p%mort2*B(i)*B(i) &
            + p%remin*mortloss)/p%rhoCN
       dudt(idxDOC) = dudt(idxDOC) + &
            (-rates%JDOC(i) &
            + rates%JCloss(i))*B(i)/p%m(i) &
            + p%remin2*p%mort2*B(i)*B(i) &
            + p%remin*mortloss
       
       dudt(idxB+i-1) = (rates%Jtot(i)/p%m(i)  &
            - (p%mort(i) &
            + rates%mortpred(i) &
            + p%mort2*B(i) &
            + p%mortHTL(i)))*B(i)
    end do
       
    !call printRates
  end function calcRates

  
end module NUMmodel


