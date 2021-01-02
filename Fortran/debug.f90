module debug
  use globals
  implicit none
  
  private
  integer, parameter:: unitDebug = 6

  
  public openDebug, closeDebug, printU, printRates
contains

subroutine openDebug
    open(unit=unitDebug, FILE='debug.out', status='replace')
    write(unitDebug,*) '------------------------------------'
  end subroutine openDebug

  subroutine closeDebug
    close(unit=unitDebug)
  end subroutine closeDebug

  subroutine printRates(rates)
    type(typeRates), intent(in):: rates
    
    !call openDebug
    write(unitDebug, *) 'm: ', m(idxB:nGrid)
!1 format(10E1.4)
    !write(unitDebug, *) 'N: ', N
    !write(unitDebug, *) 'DOC: ', DOC
    !write(unitDebug, *) 'B: ', B
    write(unitDebug, *) 'JN/m: ', rates%JN(idxB:nGrid)/m(idxB:nGrid)
    write(unitDebug, *) 'JL/m: ', rates%JL(idxB:nGrid)/m(idxB:nGrid)
    write(unitDebug, *) 'JDOC/m: ', rates%JDOC(idxB:nGrid)/m(idxB:nGrid)
    
    write(unitDebug, *) 'F/m: ', rates%F(idxB:nGrid)/m(idxB:nGrid)
    write(unitDebug, *) 'JF/m: ', rates%JF(idxB:nGrid)/m(idxB:nGrid)
    write(unitDebug, *) 'f: ', rates%flvl(idxB:nGrid)

    write(unitDebug, *) 'JLreal/m: ', rates%JLreal(idxB:nGrid)/m(idxB:nGrid)
    write(unitDebug, *) 'JNtot/m: ', rates%JNtot(idxB:nGrid)/m(idxB:nGrid)
    write(unitDebug, *) 'JCtot/m: ', rates%JCtot(idxB:nGrid)/m(idxB:nGrid)
    
   ! write(unitDebug, *) 'JFreal/m: ', rates%JFreal/p%m
    write(unitDebug, *) 'Jtot/m: ', rates%Jtot(idxB:nGrid)/m(idxB:nGrid)
    write(unitDebug, *) 'mortpred: ', rates%mortpred(idxB:nGrid)
    write(unitDebug, *) 'JNloss/m: ', rates%JNloss(idxB:nGrid)/m(idxB:nGrid)
    write(unitDebug, *) 'JCloss/m: ', rates%JCloss(idxB:nGrid)/m(idxB:nGrid)

    write(unitDebug, *) 'dudt', rates%dudt
    
    !close(unitDebug)
  end subroutine printRates

  subroutine printU(u)
    real(dp), intent(in):: u(:)
    !call openDebug
    write(unitDebug, *) 'u: ', u
    write(unitDebug, *) '---'
    !close(unitDebug)
  end subroutine printU

end module debug
