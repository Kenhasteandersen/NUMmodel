        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar  4 16:51:55 2021
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE F_CALCDERIVATIVES__genmod
          INTERFACE 
            SUBROUTINE F_CALCDERIVATIVES(NN,U,L,DT,DUDT)
              INTEGER(KIND=4), INTENT(IN) :: NN
              REAL(KIND=8), INTENT(IN) :: U(NN)
              REAL(KIND=8), INTENT(IN) :: L
              REAL(KIND=8), INTENT(IN) :: DT
              REAL(KIND=8), INTENT(OUT) :: DUDT(NN)
            END SUBROUTINE F_CALCDERIVATIVES
          END INTERFACE 
        END MODULE F_CALCDERIVATIVES__genmod
