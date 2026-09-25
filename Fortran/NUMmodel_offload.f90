!
! Prototype: OpenMP target offload (GPU) of the Euler integration over many
! grid cells for the "generalists only" setup.
!
! The object-oriented library (NUMmodel) remains the source of truth:
!  - Parameters are packed from the group objects into flat arrays at each call
!    (setupOffloadGeneralists), so setup, setHTL etc. work as usual.
!  - The physiology of a size class is in elemental procedures shared with the
!    OO code (ratesGeneralistsCore, derivativesGeneralistsCore, feedingCore).
!  - Only the assembly of the community (food availability, predation mortality,
!    Euler correction of nutrient uptakes, HTL losses) is repeated here in flat form.
!    It mirrors calcDerivatives in NUMmodel.f90 and must be kept in sync with it.
!
! Without an offload compiler/device the target region runs on the host.
!
module NUMmodel_offload
  use globals
  use spectrum
  use generalists
  use NUMmodel
  implicit none

  private
  public setupOffloadGeneralists, simulateEulerCellsGeneralists

  ! Max number of state variables per cell (size of per-thread work arrays on the device):
  integer, parameter:: nMaxGrid = 64

  ! Flat copies of the parameters (size classes of the generalists group):
  integer:: nSize = 0
  real(dp), allocatable:: pm(:), pAN(:), pAL(:), pAF(:), pJFmax(:), pJmax(:), pJresp(:), &
       pJlossPassive(:), ppHTL(:), ptheta(:,:)
  ! Scalars: epsilonL, epsilonF, mort2constant, bL, bN, bDOC, bF, bg, remin2, reminF, rhoCN, fracHTL_to_N
  real(dp):: ps(12)

contains

  ! -----------------------------------------------
  ! Pack parameters from the objects into flat arrays.
  ! Requires a generalists-only setup (setupGeneralistsOnly).
  ! -----------------------------------------------
  subroutine setupOffloadGeneralists()
    integer:: n

    if (nGroups .ne. 1 .or. nNutrients .ne. 2 .or. idxPOM .ne. 0) &
       stop 'NUMmodel_offload: only implemented for setupGeneralistsOnly'
    if (nGrid .gt. nMaxGrid) stop 'NUMmodel_offload: nGrid > nMaxGrid'

    select type (spec => group(1)%spec)
    type is (spectrumGeneralists)
       n = spec%n
       if (nSize .ne. n) then
          if (allocated(pm)) deallocate(pm, pAN, pAL, pAF, pJFmax, pJmax, pJresp, &
               pJlossPassive, ppHTL, ptheta)
          allocate(pm(n), pAN(n), pAL(n), pAF(n), pJFmax(n), pJmax(n), pJresp(n), &
               pJlossPassive(n), ppHTL(n), ptheta(n,n))
          nSize = n
       end if
       pm = spec%m
       pAN = spec%AN
       pAL = spec%AL
       pAF = spec%AF
       pJFmax = spec%JFmax
       pJmax = spec%Jmax
       pJresp = spec%Jresp
       pJlossPassive = spec%JlossPassive
       ps(1) = spec%epsilonL
       ps(2) = spec%epsilonF
       ps(3) = spec%mort2constant
    class default
       stop 'NUMmodel_offload: only implemented for setupGeneralistsOnly'
    end select

    call getParametersGeneralists(ps(4), ps(5), ps(6), ps(7), ps(8), ps(9), ps(10))
    ps(11) = rhoCN
    ps(12) = fracHTL_to_N
    ppHTL = pHTL(idxB:nGrid)
    ptheta = theta(idxB:nGrid, idxB:nGrid)
  end subroutine setupOffloadGeneralists

  ! -----------------------------------------------
  ! Integrate all cells over tEnd with Euler steps dt.
  ! Same interface as simulateEulerCells in NUMmodel.
  ! -----------------------------------------------
  subroutine simulateEulerCellsGeneralists(nCells, nGridIn, u, L, T, tEnd, dt)
    integer, intent(in):: nCells, nGridIn
    real(dp), intent(inout):: u(nCells, nGridIn)
    real(dp), intent(in):: L(nCells), T(nCells), tEnd, dt

    call setupOffloadGeneralists()
    call eulerCellsKernel(nCells, nGridIn, nSize, u, L, T, tEnd, dt, bQuadraticHTL, &
         pm, pAN, pAL, pAF, pJFmax, pJmax, pJresp, pJlossPassive, ppHTL, ptheta, ps)
  end subroutine simulateEulerCellsGeneralists

  subroutine eulerCellsKernel(nCells, nG, n, u, L, T, tEnd, dt, bQuadratic, &
       m, AN, AL, AF, JFmax, Jmax, Jresp, JlossPassive, pHTLb, thetab, ps)
    integer, intent(in):: nCells, nG, n
    real(dp), intent(inout):: u(nCells, nG)
    real(dp), intent(in):: L(nCells), T(nCells), tEnd, dt
    logical, intent(in):: bQuadratic
    real(dp), intent(in):: m(n), AN(n), AL(n), AF(n), JFmax(n), Jmax(n), Jresp(n), &
         JlossPassive(n), pHTLb(n), thetab(n,n), ps(12)
    real(dp):: ucell(nMaxGrid)
    integer:: k

    !$omp target teams distribute parallel do private(ucell) &
    !$omp&  map(tofrom: u) map(to: L, T, m, AN, AL, AF, JFmax, Jmax, Jresp, JlossPassive, pHTLb, thetab, ps)
    do k = 1, nCells
       ucell(1:nG) = u(k,1:nG)
       call eulerCellGeneralists(ucell, n, L(k), T(k), tEnd, dt, bQuadratic, &
            m, AN, AL, AF, JFmax, Jmax, Jresp, JlossPassive, pHTLb, thetab, ps)
       u(k,1:nG) = ucell(1:nG)
    end do
    !$omp end target teams distribute parallel do
  end subroutine eulerCellsKernel

  ! -----------------------------------------------
  ! Euler integration of one cell. State: u(1)=N, u(2)=DOC, u(3:n+2) = generalists.
  ! Mirrors calcDerivatives/simulateEuler in NUMmodel.f90 for this setup.
  ! -----------------------------------------------
  subroutine eulerCellGeneralists(u, n, L, T, tEnd, dt, bQuadratic, &
       m, AN, AL, AF, JFmax, Jmax, Jresp, JlossPassive, pHTLb, thetab, ps)
    !$omp declare target
    integer, intent(in):: n
    real(dp), intent(inout):: u(nMaxGrid)
    real(dp), intent(in):: L, T, tEnd, dt
    logical, intent(in):: bQuadratic
    real(dp), intent(in):: m(n), AN(n), AL(n), AF(n), JFmax(n), Jmax(n), Jresp(n), &
         JlossPassive(n), pHTLb(n), thetab(n,n), ps(12)
    integer, parameter:: iN = 1, iDOC = 2, iB = 3
    real(dp):: epsilonL, epsilonF, mort2constant, bL, bN, bDOC, bF, bg, remin2, reminF, rhoCN_, fracHTL
    real(dp):: fT2, fT15, gammaN, gammaDOC
    real(dp):: dudt(nMaxGrid), upos(nMaxGrid)
    real(dp), dimension(nMaxGrid):: F, flvl, JF, JN, JDOC, JL, Jnet, dN, Jtot, JNreal, JFreal, fSynth, &
         JDOCreal, JLreal, JNlossLiebig, JClossLiebig, JNtot, JCloss_feeding, JCloss_photouptake, &
         Jresptot, mortpred, mortHTL, mort2, jPOM, dNc, dDOCc
    integer:: istep, i, j, nG

    epsilonL = ps(1); epsilonF = ps(2); mort2constant = ps(3)
    bL = ps(4); bN = ps(5); bDOC = ps(6); bF = ps(7); bg = ps(8)
    remin2 = ps(9); reminF = ps(10); rhoCN_ = ps(11); fracHTL = ps(12)
    nG = n + 2

    fT2 = fTemp(2.d0, T)
    fT15 = fTemp(1.5d0, T)

    do istep = 1, floor(tEnd/dt)
       do i = 1, nG
          upos(i) = max(0.d0, u(i))
       end do
       !
       ! Available food and feeding:
       !
       do i = 1, n
          F(i) = 0.d0
          do j = 1, n
             F(i) = F(i) + thetab(i,j)*upos(iB-1+j)
          end do
       end do
       call feedingCore(epsilonF, AF(1:n), JFmax(1:n), fT2, F(1:n), flvl(1:n), JF(1:n))
       !
       ! HTL mortality:
       !
       if (bQuadratic) then
          mortHTL(1:n) = pHTLb(1:n) * upos(iB:nG)
       else
          mortHTL(1:n) = pHTLb(1:n)
       end if
       !
       ! Predictor step:
       !
       gammaN = 1.d0
       gammaDOC = 1.d0
       call unicellulars()
       !
       ! Correction if nutrients would become negative (as in calcDerivatives):
       !
       if ((u(iN) + dudt(iN)*dt) .lt. 0) then
          gammaN = max(0.d0, min(1.d0, -u(iN)/(dudt(iN)*dt)))
       end if
       if ((u(iDOC) + dudt(iDOC)*dt) .lt. 0) then
          gammaDOC = max(0.d0, min(1.d0, u(iDOC)/(dudt(iDOC)*dt)))
       end if
       if ((gammaN .lt. 1.d0) .or. (gammaDOC .lt. 1.d0)) then
          call unicellulars()
       end if
       !
       ! Some HTL mortality ends up as nutrients:
       !
       dudt(iN) = dudt(iN) + fracHTL * sum( upos(iB:nG) * mortHTL(1:n) )/rhoCN_
       !
       ! Euler step:
       !
       do i = 1, nG
          u(i) = u(i) + dudt(i)*dt
       end do
    end do

  contains

    subroutine unicellulars()
      call ratesGeneralistsCore(L, upos(iN), upos(iDOC), gammaN, gammaDOC, fT2, fT15, rhoCN_, &
           bL, bN, bDOC, bF, bg, epsilonL, epsilonF, &
           AN(1:n), AL(1:n), Jmax(1:n), Jresp(1:n), JlossPassive(1:n), JF(1:n), &
           JN(1:n), JDOC(1:n), JL(1:n), Jnet(1:n), dN(1:n), fSynth(1:n), Jtot(1:n), &
           JNreal(1:n), JFreal(1:n), JDOCreal(1:n), JLreal(1:n), &
           JNlossLiebig(1:n), JClossLiebig(1:n), JNtot(1:n), &
           JCloss_feeding(1:n), JCloss_photouptake(1:n), Jresptot(1:n))
      JF(1:n) = JFreal(1:n) ! as in calcRatesGeneralists
      !
      ! Predation mortality:
      !
      do i = 1, n
         mortpred(i) = 0.d0
         do j = 1, n
            if (F(j) .gt. 0.d0) then
               mortpred(i) = mortpred(i) &
                    + thetab(j,i) * JF(j)*upos(iB-1+j) / (epsilonF*m(j)*F(j))
            end if
         end do
      end do
      !
      ! Derivatives:
      !
      dudt(iN) = 0.d0
      dudt(iDOC) = 0.d0
      call derivativesGeneralistsCore(upos(iB:nG), m(1:n), mort2constant, remin2, reminF, rhoCN_, &
           JNreal(1:n), JlossPassive(1:n), JNlossLiebig(1:n), JCloss_feeding(1:n), &
           JDOCreal(1:n), JCloss_photouptake(1:n), Jtot(1:n), mortpred(1:n), mortHTL(1:n), &
           mort2(1:n), jPOM(1:n), dNc(1:n), dDOCc(1:n), dudt(iB:nG))
      do i = 1, n
         dudt(iN) = dudt(iN) + dNc(i)
         dudt(iDOC) = dudt(iDOC) + dDOCc(i)
      end do
    end subroutine unicellulars

  end subroutine eulerCellGeneralists

end module NUMmodel_offload
