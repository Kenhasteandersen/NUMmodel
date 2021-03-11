        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar  4 16:51:55 2021
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE F_GETRATES__genmod
          INTERFACE 
            SUBROUTINE F_GETRATES(JN,JDOC,JL,JF,JTOT,JMAX,JR,           &
     &JLOSSPASSIVE,JNLOSS,JLREAL,MORTPRED,MORTHTL,MORT2,MORT)
              USE GLOBALS
              REAL(KIND=8), INTENT(OUT) :: JN(NGRID)
              REAL(KIND=8), INTENT(OUT) :: JDOC(NGRID)
              REAL(KIND=8), INTENT(OUT) :: JL(NGRID)
              REAL(KIND=8), INTENT(OUT) :: JF(NGRID)
              REAL(KIND=8), INTENT(OUT) :: JTOT(NGRID)
              REAL(KIND=8), INTENT(OUT) :: JMAX(NGRID)
              REAL(KIND=8), INTENT(OUT) :: JR(NGRID)
              REAL(KIND=8), INTENT(OUT) :: JLOSSPASSIVE(NGRID)
              REAL(KIND=8), INTENT(OUT) :: JNLOSS(NGRID)
              REAL(KIND=8), INTENT(OUT) :: JLREAL(NGRID)
              REAL(KIND=8), INTENT(OUT) :: MORTPRED(NGRID)
              REAL(KIND=8), INTENT(OUT) :: MORTHTL(NGRID)
              REAL(KIND=8), INTENT(OUT) :: MORT2(NGRID)
              REAL(KIND=8), INTENT(OUT) :: MORT(NGRID)
            END SUBROUTINE F_GETRATES
          END INTERFACE 
        END MODULE F_GETRATES__genmod
