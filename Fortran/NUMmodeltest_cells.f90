!
! Test of simulations over many grid cells:
!   1) reference: serial loop over cells calling simulateEuler (as matlab does now)
!   2) simulateEulerCells: OpenMP threads over cells, object-oriented code
!   3) (generalists only) simulateEulerCellsGeneralists: flat kernel for GPU offload
! Run from a directory one level below the repository root (input file is ../input/input.yaml).
!
program NUMmodeltest_cells
  use globals
  use NUMmodel
  use NUMmodel_offload
  !$ use omp_lib
  implicit none

  integer, parameter:: nCells = 4000, nRep = 20
  real(dp), parameter:: tEnd = 0.5d0, dt = 0.1d0
  character(len=200):: errorstr
  logical(1):: errorio
  integer:: nThreads

  errorio = .false.
  nThreads = 1
  !$ nThreads = omp_get_max_threads()
  write(*,'(a,i0,a,i0)') 'Threads: ', nThreads, '   cells: ', nCells

  write(*,*) '--- Generalists only (10 size classes) ---'
  call setupGeneralistsOnly(10, errorio, errorstr)
  call runTests(.true.)

  write(*,*) '--- Generalists only, quadratic HTL mortality ---'
  call setHTL(0.01d0, 0.1d0, logical(.true.,1), logical(.false.,1), logical(.false.,1))
  call runTests(.true.)

  write(*,*) '--- Full NUM model (generalists, copepods, POM) ---'
  call setupNUMmodel(10, 6, 1, (/0.1d0/), (/1.d0, 10.d0/), errorio, errorstr)
  call runTests(.false.)

contains

  subroutine runTests(bFlat)
    logical, intent(in):: bFlat
    real(dp), allocatable:: u0(:,:), uRef(:,:), uCells(:,:), uFlat(:,:), L(:), T(:)
    real(dp):: ucell(nGrid), t0, tRef, tCells, tFlat
    integer:: k, iRep

    allocate(u0(nCells,nGrid), uRef(nCells,nGrid), uCells(nCells,nGrid), uFlat(nCells,nGrid))
    allocate(L(nCells), T(nCells))
    !
    ! Cells with a range of light, temperature and nutrients:
    !
    do k = 1, nCells
       L(k) = 300.d0 * modulo(0.618034d0*k, 1.d0)
       T(k) = -2.d0 + 32.d0 * modulo(0.414214d0*k, 1.d0)
       u0(k,:) = 0.5d0 + modulo(0.732051d0*k, 1.d0)
       u0(k,idxN) = 1.d0 + 150.d0 * modulo(0.236068d0*k, 1.d0)
       u0(k,idxDOC) = 0.01d0 + 10.d0 * modulo(0.3166d0*k, 1.d0)
       if (nNutrients .gt. 2) u0(k,idxSi) = 10.d0
    end do
    !
    ! Reference: serial loop (as in matlab now):
    !
    uRef = u0
    t0 = wtime()
    do iRep = 1, nRep
       do k = 1, nCells
          ucell = uRef(k,:)
          call simulateEuler(ucell, L(k), T(k), tEnd, dt)
          uRef(k,:) = ucell
       end do
    end do
    tRef = wtime() - t0
    !
    ! Threaded, object oriented:
    !
    uCells = u0
    t0 = wtime()
    do iRep = 1, nRep
       call simulateEulerCells(nCells, uCells, L, T, tEnd, dt)
    end do
    tCells = wtime() - t0
    write(*,'(a,f8.3,a)')       '  serial reference      : ', tRef, ' s'
    write(*,'(a,f8.3,a,es10.2)') '  simulateEulerCells    : ', tCells, ' s   max rel diff: ', maxRelDiff(uCells, uRef)
    !
    ! Flat kernel:
    !
    if (bFlat) then
       uFlat = u0
       t0 = wtime()
       call setupOffloadGeneralists()
       do iRep = 1, nRep
          call simulateEulerCellsGeneralists(nCells, nGrid, uFlat, L, T, tEnd, dt)
       end do
       tFlat = wtime() - t0
       write(*,'(a,f8.3,a,es10.2)') '  flat kernel (offload) : ', tFlat, ' s   max rel diff: ', maxRelDiff(uFlat, uRef)
    end if
  end subroutine runTests

  function wtime() result(t)
    real(dp):: t
    integer(8):: count, rate
    call system_clock(count, rate)
    t = real(count, dp)/real(rate, dp)
  end function wtime

  function maxRelDiff(a, b) result(d)
    real(dp), intent(in):: a(:,:), b(:,:)
    real(dp):: d
    d = maxval( abs(a-b) / (abs(b) + 1d-10) )
  end function maxRelDiff

end program NUMmodeltest_cells
