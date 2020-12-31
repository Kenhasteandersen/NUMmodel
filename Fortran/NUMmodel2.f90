
module NUMmodel
  implicit none
  private
  integer, parameter :: dp=kind(0.d0) ! double precision
  integer, parameter :: idxN = 1
  integer, parameter :: idxDOC = 2
  integer, parameter :: idxB = 3
  integer, parameter :: unitDebug = 10
  
  integer:: n, nGroups
 
  type typeRates
     real(dp), dimension(:), allocatable:: JN, JDOC, JL, JF, F, JFreal
     real(dp), dimension(:), allocatable:: JNtot, JLreal, JCtot, Jtot;
     real(dp), dimension(:), allocatable:: JCloss_feeding, JCloss_photouptake, JNlossLiebig, JClossLiebig
     real(dp), dimension(:), allocatable:: JNloss, JCloss
     real(dp), dimension(:), allocatable:: mortpred
  end type typeRates

  real(dp), dimension(:), allocatable:: B, ANmT, JmaxT, JFmaxT, JrespT
  
  public
  real(dp), dimension(:), allocatable:: u(:), m(:), AF(:),  JFmax(:), Jresp(:), epsilonF(:)
  real(dp), dimension(:,:), allocatable:: theta(:,:)
  real(dp), dimension(:), allocatable:: mortHTL(:)
  type(typeRates):: rates

  parametersInit, setParameters, calcRates, calcGrid
contains

  subroutine openDebug
    open(unit=unitDebug, FILE='debug.out', status='replace')
    write(unitDebug,*) '------------------------------------'
  end subroutine openDebug

  subroutine calcGrid(mMin, mMax, n, ixStart)
    real(dp), intent(in):: mMin, mMax
    integer, intent(in):: n, ix
    real(dp):: deltax
    integer:: i, ix

    deltax = (log(mMax)-log(mMin))/(n-1)x(2)-x(1)
    do i=1,n
       ix = ixStart+i-1
       x(i) = log(mMin) + (i-1)*deltax
       m(ix) = exp(x(i))
       mLower(ix) = exp(x(i)-0.5*deltax)
       mDelta(ix) =  exp(x(i)+0.5*deltax)-mLower(ix)
       z(ix) = mLower(ix)/(mLower(ix)+mDelta(ix))
    end do
  end subroutine calcGrid
  
  subroutine parametersInit(n)
    integer, intent(in):: n

    nGroups = 0
    allocate(m(n))
    allocate(AF(n))
    allocate(JFmax(n))
    allocate(Jresp(n))
    allocate(epsilonF(n))
  end subroutine parametersInit

  subroutine parametersAddGroup(typeGroup, n, mMax)
    integer, intent(in):: typeGroup, n
    real(dp), intent(in):: mMax

    p%nGroups = p%nGroups + 1
    p%typeGroups(p%nGroups) = typeGroup
    !
    ! Add the groups
    !
    select case (typeGroup)
    case (typeGeneralist)
       parametersAddGeneralists(n)
    end select
  end subroutine parametersAddGroup

  subroutine calcDerivatives
    !
    ! Calc uptakes of food
    !

    !
    ! Calc other uptakes of unicellular groups
    !
    calcRatesGeneralists()
    !
    ! Calc predation mortality
    !

    !
    ! Assemble derivatives:
    !
    rates%dudt(idxN) = 0
    rates%dudt(idxDOC) = 0
    select case (typeGroup)
    case (typeGeneralist)
       calcDerivativesGeneralists()
    end select
  end subroutine calcDerivatives
  
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
  
  
end module NUMmodel


