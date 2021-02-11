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

  subroutine printRates(m, rates)
    type(typeRates), intent(in):: rates
    real(dp):: m(:)
    
    !call openDebug
    write(unitDebug, *) 'm: ', m
!1 format(10E1.4)
    !write(unitDebug, *) 'N: ', N
    !write(unitDebug, *) 'DOC: ', DOC
    !write(unitDebug, *) 'B: ', B
    write(unitDebug, *) 'JN/m(idxB:nGrid): ', rates%JN(idxB:nGrid)/m(idxB:nGrid)
    write(unitDebug, *) 'JL/m(idxB:nGrid): ', rates%JL(idxB:nGrid)/m(idxB:nGrid)
    write(unitDebug, *) 'JDOC/m(idxB:nGrid): ', rates%JDOC(idxB:nGrid)/m(idxB:nGrid)
    
    write(unitDebug, *) 'F/m(idxB:nGrid): ', rates%F(idxB:nGrid)/m(idxB:nGrid)
    write(unitDebug, *) 'JF/m(idxB:nGrid): ', rates%JF(idxB:nGrid)/m(idxB:nGrid)
    write(unitDebug, *) 'f: ', rates%flvl(idxB:nGrid)

    write(unitDebug, *) 'JLreal/m(idxB:nGrid): ', rates%JLreal(idxB:nGrid)/m(idxB:nGrid)
    write(unitDebug, *) 'JNtot/m(idxB:nGrid): ', rates%JNtot(idxB:nGrid)/m(idxB:nGrid)
    write(unitDebug, *) 'JCtot/m(idxB:nGrid): ', rates%JCtot(idxB:nGrid)/m(idxB:nGrid)
    
   ! write(unitDebug, *) 'JFreal/m(idxB:nGrid): ', rates%JFreal/p%m
    write(unitDebug, *) 'Jtot/m(idxB:nGrid): ', rates%Jtot(idxB:nGrid)/m(idxB:nGrid)
    write(unitDebug, *) 'mortpred: ', rates%mortpred(idxB:nGrid)
    write(unitDebug, *) 'JNloss/m(idxB:nGrid): ', rates%JNloss(idxB:nGrid)/m(idxB:nGrid)
    write(unitDebug, *) 'JCloss/m(idxB:nGrid): ', rates%JCloss(idxB:nGrid)/m(idxB:nGrid)

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
