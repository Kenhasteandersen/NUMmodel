!
! Benchmark of the transport matrix step w = Aimp*(Aexp*u), written as a threaded
! CSR sparse-matrix times dense-matrix product, to compare against the matlab
! version. This is a measurement tool and not part of the library. Build it with
!
!   gfortran -O3 -fopenmp -o spmmbench Fortran/spmmbench.f90
!
! and run it through testTransportSpeed.m, which writes the files it reads.
!
! -O3 matters here: at -O2 the inner loop is not vectorized and the tracer-major
! kernel is twice as slow.
!
! Three layouts of the state are timed:
!   cell-major    as matlab holds it, u(nb,nGrid)
!   tracer-major  transposed, u(nGrid,nb), so that the tracers of one grid cell
!                 lie next to each other and the inner loop becomes a
!                 vectorizable axpy. This is the layout matlab cannot use.
!   tracer single the same, in single precision
!
program spmmbench
  use omp_lib
  implicit none
  integer, parameter:: dp = kind(0.d0), sp = kind(0.e0)

  integer:: nb, nG, nnzE, nnzI, nt, i, r, nmax
  integer, allocatable:: ptrE(:), colE(:), ptrI(:), colI(:)
  real(dp), allocatable:: valE(:), valI(:), u(:,:), w(:,:), v(:,:), ref(:,:)
  real(dp), allocatable:: uT(:,:), wT(:,:), vT(:,:)
  real(sp), allocatable:: valEs(:), valIs(:), uTs(:,:), wTs(:,:), vTs(:,:)
  real(dp):: t0, tbest
  character(len=512):: dir

  call get_command_argument(1, dir)
  if (len_trim(dir) == 0) then
     write(*,*) 'usage: spmmbench <data directory written by testTransportSpeed.m>'
     stop 1
  end if
  dir = trim(dir)//'/'

  call readints(trim(dir)//'dims.txt', nb, nG)
  call readint (trim(dir)//'Aexp_nnz.txt', nnzE)
  call readint (trim(dir)//'Aimp_nnz.txt', nnzI)
  nmax = omp_get_max_threads()
  write(*,'(a,i0,a,i0,a,f8.0,a)') 'nb = ', nb, ', nGrid = ', nG, &
       ', state = ', nb*real(nG,dp)*8/2.d0**20, ' MB'
  write(*,'(a,i0,a,i0,a,i0)') 'nnz(Aexp) = ', nnzE, ', nnz(Aimp) = ', nnzI, &
       ', threads available = ', nmax

  allocate(ptrE(nb+1), colE(nnzE), valE(nnzE), ptrI(nb+1), colI(nnzI), valI(nnzI))
  call readi4(trim(dir)//'Aexp_rowptr.bin', ptrE, nb+1)
  call readi4(trim(dir)//'Aexp_col.bin',    colE, nnzE)
  call readr8(trim(dir)//'Aexp_val.bin',    valE, nnzE)
  call readi4(trim(dir)//'Aimp_rowptr.bin', ptrI, nb+1)
  call readi4(trim(dir)//'Aimp_col.bin',    colI, nnzI)
  call readr8(trim(dir)//'Aimp_val.bin',    valI, nnzI)
  allocate(u(nb,nG), w(nb,nG), v(nb,nG), ref(nb,nG))
  call readr8(trim(dir)//'u.bin',   u,   nb*nG)
  call readr8(trim(dir)//'ref.bin', ref, nb*nG)
  allocate(uT(nG,nb), wT(nG,nb), vT(nG,nb))
  uT = transpose(u)

  write(*,*) ''
  write(*,'(a)') 'seconds per transport step, best of 3'
  write(*,'(a)') ' threads   cell-major  tracer-major  tracer single'
  do i = 1, 10
     nt = threadCount(i, nmax)
     if (nt < 0) exit
     call omp_set_num_threads(nt)
     write(*,'(i6,a)', advance='no') nt, '   '

     tbest = huge(1.d0)
     do r = 1, 3
        t0 = wclock()
        call spmmCell(nb, nG, ptrE, colE, valE, u, v)
        call spmmCell(nb, nG, ptrI, colI, valI, v, w)
        tbest = min(tbest, wclock()-t0)
     end do
     write(*,'(f10.3,a)', advance='no') tbest, '   '

     tbest = huge(1.d0)
     do r = 1, 3
        t0 = wclock()
        call spmmTracer(nb, nG, ptrE, colE, valE, uT, vT)
        call spmmTracer(nb, nG, ptrI, colI, valI, vT, wT)
        tbest = min(tbest, wclock()-t0)
     end do
     write(*,'(f11.3,a)', advance='no') tbest, '   '

     if (.not. allocated(valEs)) then
        allocate(valEs(nnzE), valIs(nnzI), uTs(nG,nb), wTs(nG,nb), vTs(nG,nb))
        valEs = real(valE,sp); valIs = real(valI,sp); uTs = real(uT,sp)
     end if
     tbest = huge(1.d0)
     do r = 1, 3
        t0 = wclock()
        call spmmTracerSingle(nb, nG, ptrE, colE, valEs, uTs, vTs)
        call spmmTracerSingle(nb, nG, ptrI, colI, valIs, vTs, wTs)
        tbest = min(tbest, wclock()-t0)
     end do
     write(*,'(f11.3)') tbest
  end do

  write(*,*) ''
  write(*,'(a,es10.2)') 'max relative difference from matlab, double = ', &
       maxval(abs(transpose(wT)-ref))/maxval(abs(ref))
  write(*,'(a,es10.2)') 'max relative difference from matlab, single = ', &
       maxval(abs(real(transpose(wTs),dp)-ref))/maxval(abs(ref))

contains

  !
  ! 1, 2, 4, 8 ... up to the number of threads available, and then that number
  !
  function threadCount(i, nmax) result(nt)
    integer, intent(in):: i, nmax
    integer:: nt

    nt = 2**(i-1)
    if (nt > nmax) then
       if (2**(i-2) < nmax) then
          nt = nmax ! the last step, e.g. 24 after 16
       else
          nt = -1   ! done
       end if
    end if
  end function threadCount

  !
  ! The state as matlab holds it: the tracers of one grid cell are nb apart, so
  ! each nonzero touches nGrid separate cache lines.
  !
  subroutine spmmCell(nb, nG, ptr, col, val, x, y)
    integer, intent(in):: nb, nG, ptr(:), col(:)
    real(dp), intent(in):: val(:), x(nb,nG)
    real(dp), intent(out):: y(nb,nG)
    real(dp):: acc(nG), a
    integer:: i, k, j, c

    !$omp parallel do schedule(static) private(i,k,j,c,a,acc)
    do i = 1, nb
       acc = 0.d0
       do k = ptr(i)+1, ptr(i+1)
          j = col(k)
          a = val(k)
          do c = 1, nG
             acc(c) = acc(c) + a*x(j,c)
          end do
       end do
       do c = 1, nG
          y(i,c) = acc(c)
       end do
    end do
    !$omp end parallel do
  end subroutine spmmCell

  !
  ! The state transposed: x(:,j) is contiguous, so the inner loop is an axpy over
  ! nGrid elements and vectorizes.
  !
  subroutine spmmTracer(nb, nG, ptr, col, val, x, y)
    integer, intent(in):: nb, nG, ptr(:), col(:)
    real(dp), intent(in):: val(:), x(nG,nb)
    real(dp), intent(out):: y(nG,nb)
    real(dp):: acc(nG), a
    integer:: i, k, j

    !$omp parallel do schedule(static) private(i,k,j,a,acc)
    do i = 1, nb
       acc = 0.d0
       do k = ptr(i)+1, ptr(i+1)
          j = col(k)
          a = val(k)
          acc = acc + a*x(:,j)
       end do
       y(:,i) = acc
    end do
    !$omp end parallel do
  end subroutine spmmTracer

  subroutine spmmTracerSingle(nb, nG, ptr, col, val, x, y)
    integer, intent(in):: nb, nG, ptr(:), col(:)
    real(sp), intent(in):: val(:), x(nG,nb)
    real(sp), intent(out):: y(nG,nb)
    real(sp):: acc(nG), a
    integer:: i, k, j

    !$omp parallel do schedule(static) private(i,k,j,a,acc)
    do i = 1, nb
       acc = 0.e0
       do k = ptr(i)+1, ptr(i+1)
          j = col(k)
          a = val(k)
          acc = acc + a*x(:,j)
       end do
       y(:,i) = acc
    end do
    !$omp end parallel do
  end subroutine spmmTracerSingle

  function wclock() result(tt)
    integer(8):: c, r
    real(dp):: tt

    call system_clock(c, r)
    tt = real(c,dp)/real(r,dp)
  end function wclock

  subroutine readints(f, a, b)
    character(len=*), intent(in):: f
    integer, intent(out):: a, b

    open(10, file=f, status='old')
    read(10,*) a, b
    close(10)
  end subroutine readints

  subroutine readint(f, a)
    character(len=*), intent(in):: f
    integer, intent(out):: a

    open(10, file=f, status='old')
    read(10,*) a
    close(10)
  end subroutine readint

  subroutine readi4(f, a, n)
    character(len=*), intent(in):: f
    integer, intent(in):: n
    integer, intent(out):: a(n)

    open(11, file=f, form='unformatted', access='stream', status='old')
    read(11) a
    close(11)
  end subroutine readi4

  subroutine readr8(f, a, n)
    character(len=*), intent(in):: f
    integer, intent(in):: n
    real(dp), intent(out):: a(n)

    open(11, file=f, form='unformatted', access='stream', status='old')
    read(11) a
    close(11)
  end subroutine readr8

end program spmmbench
