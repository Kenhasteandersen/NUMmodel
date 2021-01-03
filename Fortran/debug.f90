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
    write(unitDebug, *) 'JN/m: ', rates%JN/m
    write(unitDebug, *) 'JL/m: ', rates%JL/m
    write(unitDebug, *) 'JDOC/m: ', rates%JDOC/m
    
    write(unitDebug, *) 'F/m: ', rates%F/m
    write(unitDebug, *) 'JF/m: ', rates%JF/m
    write(unitDebug, *) 'f: ', rates%flvl

    write(unitDebug, *) 'JLreal/m: ', rates%JLreal/m
    write(unitDebug, *) 'JNtot/m: ', rates%JNtot/m
    write(unitDebug, *) 'JCtot/m: ', rates%JCtot/m
    
   ! write(unitDebug, *) 'JFreal/m: ', rates%JFreal/p%m
    write(unitDebug, *) 'Jtot/m: ', rates%Jtot/m
    write(unitDebug, *) 'mortpred: ', rates%mortpred
    write(unitDebug, *) 'JNloss/m: ', rates%JNloss/m
    write(unitDebug, *) 'JCloss/m: ', rates%JCloss/m

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
