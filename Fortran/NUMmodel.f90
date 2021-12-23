
!
! Module to handle the NUM model framework
!
module NUMmodel
  use globals
  use spectrum
  use generalists
  !use generalists_csp
  !use diatoms
  !use diatoms_simple
  !use copepods
  !use debug
  implicit none

  !
  ! Variables that contain the size spectrum groups
  !
  integer:: nGroups ! Number of groups
  integer:: iCurrentGroup ! The current group to be added
  integer:: nNutrients ! Number of nutrient state variables
  integer:: idxB ! First index into non-nutrient groups (=nNutrients+1)
  type(spectrumContainer), allocatable :: group(:) ! Structure pointing to each group
  integer, dimension(:), allocatable :: ixStart, ixEnd ! Indices into u for each group

  real(dp), dimension(:,:), allocatable:: theta ! Interaction matrix
  real(dp), dimension(:), allocatable:: dudt
  real(dp), dimension(:), allocatable:: upositive ! State variable constrained to be positive
  real(dp), dimension(:), allocatable:: F ! Available food
  !
  ! Variables for HTL mortalities:
  !
  !real(dp), dimension(:), allocatable:: pHTL ! Selectivity function for HTL mortality
  !real(dp) :: mortHTL ! Level of HTL mortality (see below for override)
  !real(dp):: gammaHTL ! Parameter for quadratic HTL mortality
  !logical:: bQuadraticHTL ! Boolean flag to signify whether mortality is standard or "quadratic"
  ! Defaults to true; can be overridden but parameters must then be set with a
  ! call to parametersFinalize()

contains

  ! ======================================
  !  Various model setups
  ! ======================================

  ! -----------------------------------------------
  ! A basic setup with only generalists
  ! -----------------------------------------------
  subroutine setupGeneralistsOnly(n)
    integer, intent(in):: n
    call parametersInit(1, n, 2) ! 1 group, n size classes (excl nutrients and DOC)
    call parametersAddGroup(typeGeneralist, n, 1.d0) ! generalists with n size classes
    call parametersFinalize(0.2d0, .false.) ! Use standard "linear" mortality
  end subroutine setupGeneralistsOnly

  ! -----------------------------------------------
  ! A basic setup with only generalists -- (Serra-Pompei et al 2020 version)
  ! -----------------------------------------------
!   subroutine setupGeneralistsOnly_csp()
!     call parametersInit(1, 10, 2) ! 1 group, 10 size classes (excl nutrients and DOC)
!     call parametersAddGroup(typeGeneralist_csp, 10, 10.d0**(-1.3d0)) ! generalists with 10 size classes
!     call parametersFinalize(0.003d0, .true.) ! Serra-Pompei (2020))
!   end subroutine setupGeneralistsOnly_csp

!   ! -----------------------------------------------
!   ! A basic setup with only diatoms:
!   ! -----------------------------------------------
!   subroutine setupDiatomsOnly(n)
!    integer, intent(in):: n
!    call parametersInit(1, n, 3) ! 1 group, n size classes (excl nutrients)
!    call parametersAddGroup(typeDiatom, n, 1.d0) ! diatoms with n size classes
!    call parametersFinalize(0.2d0, .false.)
!  end subroutine setupDiatomsOnly

!  ! -----------------------------------------------
!   ! A basic setup with only simple diatoms:
!   ! -----------------------------------------------
!  subroutine setupDiatoms_simpleOnly(n)
!    integer, intent(in):: n
!    call parametersInit(1, n, 3) ! 1 group, n size classes (excl nutrients)
!    call parametersAddGroup(typeDiatom_simple, n, 1.d0) ! diatoms with n size classes
!    call parametersFinalize(0.2d0, .false.)
!  end subroutine setupDiatoms_simpleOnly
 
!   ! -----------------------------------------------
!   ! Generalists and diatoms:
!   ! -----------------------------------------------
!    subroutine setupGeneralistsDiatoms(n)
!       integer, intent(in):: n
!       call parametersInit(2, 2*n, 3)
!       call parametersAddGroup(typeGeneralist, n, 1.d0) ! generalists with n size classes
!       call parametersAddGroup(typeDiatom, n, 1.d0) ! diatoms with n size classes
!       call parametersFinalize(.2d0, .false.)
!    end subroutine setupGeneralistsDiatoms
 
!    subroutine setupGeneralistsDiatoms_simple(n)
!       integer, intent(in):: n
!       call parametersInit(2, 2*n, 3)
!       call parametersAddGroup(typeGeneralist, n, 1.d0) ! generalists with n size classes
!       call parametersAddGroup(typeDiatom_simple, n, 1.d0) ! diatoms with n size classes
!       call parametersFinalize(0.2d0, .false.)
!    end subroutine setupGeneralistsDiatoms_simple
 
!   ! -----------------------------------------------
!   ! A basic setup with generalists and 1 copepod
!   ! -----------------------------------------------
!   subroutine setupGeneralistsCopepod()
!     call parametersInit(2, 20, 2)
!     call parametersAddGroup(typeGeneralist, 10, 0.1d0)
!     call parametersAddGroup(typeCopepod, 10, .1d0) ! add copepod with adult mass .1 mugC
!     call parametersFinalize(0.003d0, .true.) ! Use quadratic mortality
!   end subroutine setupGeneralistsCopepod

!   ! -----------------------------------------------
!   ! A generic setup with generalists and a number of copepod species
!   ! -----------------------------------------------
!   subroutine setupGeneric(mAdult)
!     real(dp), intent(in):: mAdult(:)
!     integer, parameter:: n = 10 ! number of size classes in each group
!     integer:: iCopepod

!     call parametersInit(size(mAdult)+1, n*(size(mAdult)+1), 2)
!     call parametersAddGroup(typeGeneralist, n, 0.1d0)
!     if ( size(mAdult) .eq. 0) then
!        bQuadraticHTL = .false. ! Use standard "linear" mortality
!     else
!        do iCopepod = 1, size(mAdult)
!           call parametersAddGroup(typeCopepod, n, mAdult(iCopepod)) ! add copepod
!        end do
!     end if
!     call parametersFinalize(0.003d0, .true.)
!   end subroutine setupGeneric

!   ! -----------------------------------------------
!   ! A generic setup with generalists and a number of copepod species
!   ! -----------------------------------------------
!   subroutine setupGeneric_csp(mAdult)
!     real(dp), intent(in):: mAdult(:)
!     integer, parameter:: n = 10 ! number of size classes in each group
!     integer:: iCopepod

!     call parametersInit(size(mAdult)+1, n*(size(mAdult)+1), 2)
!     call parametersAddGroup(typeGeneralist_csp, n, 0.1d0)
!     do iCopepod = 1, size(mAdult)
!        call parametersAddGroup(typeCopepod, n, mAdult(iCopepod)) ! add copepod
!     end do
!     call parametersFinalize(0.003d0, .true.)
!   end subroutine setupGeneric_csp



  ! ======================================
  !  Model initialization stuff:
  ! ======================================

  ! -----------------------------------------------
  ! Initialize parameters
  ! In:
  !    nnGroups: number of size spectrum groups
  !    nnGrid: total length of the grid (excl nnNutrients points for N, DOC, etc.)
  ! -----------------------------------------------
  subroutine parametersInit(nnGroups, nnGrid, nnNutrients)
    integer, intent(in):: nnGrid, nnGroups, nnNutrients
    !
    ! Set groups:
    !
    nGroups = nnGroups
    iCurrentGroup = 0
    nNutrients = nnNutrients
    nGrid = nnGrid+nnNutrients
    idxB = nNutrients + 1
    !
    ! Allocate variables:
    !
    if (allocated(upositive)) then
       deallocate(group)
       deallocate(ixStart)
       deallocate(ixEnd)
       deallocate(upositive)
       deallocate(dudt)
       deallocate(F)
       deallocate(theta)
    end if

    allocate( group(nGroups) )
    allocate(ixStart(nGroups))
    allocate(ixEnd(nGroups))
    allocate(upositive(nGrid))
    allocate(dudt(nGrid))
    allocate(F(nGrid))
    allocate(theta(nGrid,nGrid))     ! Interaction matrix:
 
  end subroutine parametersInit
  ! -----------------------------------------------
  !  Add a size spectrum group
  !  In:
  !    typeGroup: the group type (see definitions in Globals.f90
  !    n: number of grid points
  !    mMax: the maximum size (mid-point in grid cell)
  ! -----------------------------------------------
  subroutine parametersAddGroup(typeGroup, n, mMax)
    integer, intent(in):: typeGroup, n
    real(dp), intent(in):: mMax

    type(spectrumGeneralists) :: generalists
    !
    ! Find the group number and grid location:
    !
    iCurrentGroup = iCurrentGroup + 1

    if (iCurrentGroup.eq.1) then
      ixStart(iCurrentGroup) = idxB
    else
      ixStart(iCurrentGroup) = ixEnd(iCurrentGroup-1)+1
    end if
    ixEnd(iCurrentGroup) = ixStart(iCurrentGroup)+n-1
    !
    ! Add the group
    !
    select case (typeGroup)
    case (typeGeneralist)
      call initGeneralists(generalists, n, mMax)
      allocate( group( iCurrentGroup )%spec, source=generalists )
        !group(iCurrentGroup) = spectrumContainer( )
      !call group(iCurrentGroup)%init( n, ixStart-1, mMax)
   !  case (typeGeneralist_csp)
   !    group(iCurrentGroup) = initGeneralists_csp(n, group(iCurrentGroup)%ixStart-1, mMax)
   !  case (typeDiatom)
   !    group(iCurrentGroup) = initDiatoms(n, group(iCurrentGroup)%ixStart-1, mMax)
   !    !palatability(group(iCurrentGroup)%ixStart:group(iCurrentGroup)%ixEnd) = 0.5d0 ! Lower palatability for diatoms
   !  case (typeDiatom_simple)
   !    group(iCurrentGroup) = initDiatoms_simple(n, group(iCurrentGroup)%ixStart-1, mMax)
   !    !palatability(group(iCurrentGroup)%ixStart:group(iCurrentGroup)%ixEnd) = 0.5d0 ! Lower palatability for diatoms
   !  case(typeCopepod)
   !     group(iCurrentGroup) = initCopepod(n, group(iCurrentGroup)%ixStart-1, mMax)
    end select
  end subroutine parametersAddGroup
  ! -----------------------------------------------
  !  Finalize the setting of parameters. Must be called when
  !  all groups have been added.
  ! -----------------------------------------------
  subroutine parametersFinalize(mortHTL, bQuadraticHTL)
    real(dp), intent(in):: mortHTL
    logical, intent(in):: bQuadraticHTL
    integer:: i,j, iGroup, jGroup
    real(dp),parameter :: betaHTL = 500.d0
    real(dp):: mHTL
    !
    ! Calc theta:
    !
    do iGroup = 1, nGroups
      do i = 1, group(iGroup)%spec%n !group(iGroup)%ixStart, group(iGroup)%ixEnd
         do jGroup = 1, nGroups
            do j = 1, group(jGroup)%spec%n!group(jGroup)%ixStart, group(jGroup)%ixEnd
               theta(i+ixStart(jGroup)-1, j+ixStart(jGroup)-1) = &
                  group(jGroup)%spec%palatability * &
                  calcPhi(group(iGroup)%spec%m(i)/group(jGroup)%spec%m(j), &
                     group(iGroup)%spec%beta, group(iGroup)%spec%sigma, &
                     group(iGroup)%spec%z(i))
            end do
         end do
      end do
   end do
   !
   ! Set HTL mortality
   !

   ! Find the largest mass
   mHTL = 0.d0
   do iGroup = 1, nGroups
      mHTL = max( mHTL, group(iGroup)%spec%m( group(iGroup)%spec%n ) )
   end do
   ! Calc the mass where HTL mortality is 50%
   mHTL = mHTL/betaHTL**1.5

   call setHTL(mHTL, mortHTL, bQuadraticHTL)

  contains
    !
    ! Calculate the interaction coefficient between two size groups.
    ! In:
    !   z : The predator:prey body mass ratio between the two groups
    !   beta: preferred predator:prey body mass ratio
    !   sigma: width of selection
    !   Delta: ratio between upper and lower body mass in size groups
    !
    function calcPhi(z, beta,sigma, Delta) result(res)
      real(dp), intent(in):: z,beta,sigma,Delta
      real(dp):: res, s

      if (beta .eq. 0.d0) then
         res = 0.d0 ! If beta = 0 it is interpreted as if the group is not feeding
      else
         s = 2*sigma*sigma
         res = max(0.d0, &
         (Sqrt(Delta)*(((exp(-Log((beta*Delta)/z)**2/s) - 2/exp(Log(z/beta)**2/s) + &
         exp(-Log((Delta*z)/beta)**2/s))*s)/2. - &
         (Sqrt(Pi)*Sqrt(s)*(Erf((-Log(beta*Delta) + Log(z))/Sqrt(s))*Log((beta*Delta)/z) + &
         2*Erf(Log(z/beta)/Sqrt(s))*Log(z/beta) + &
         Erf((Log(beta) - Log(Delta*z))/Sqrt(s))*Log((Delta*z)/beta)))/2.))/ &
         ((-1 + Delta)*Log(Delta)) )
      end if
    end function calcPhi

  end subroutine parametersFinalize

  subroutine setHTL(mHTL, mortalityHTL, boolQuadraticHTL)
    real(dp), intent(in):: mHTL ! The size where HTL is 50% of max
    real(dp), intent(in):: mortalityHTL ! The level of HTL mortality
    logical, intent(in):: boolQuadraticHTL ! Whether to use "quadratic" mortality
    real(dp), parameter:: betaHTL = 500.
    integer:: iGroup
 
    !     
    ! Calc htl mortality
    !
    if (.not. boolQuadraticHTL) then
      !
      ! Standard HTL mortality. In this case "pHTL" represent the selectivity of HTL mortality
      !
         do iGroup = 1, nGroups
         group(iGroup)%spec%mortHTL = &
               mortalityHTL * (1 / (1+(group(iGroup)%spec%m/mHTL)**(-2)) )
      end do
    else
      !
      ! Linear HTL mortality (commonly referred to as "quadratic")
      ! In this case "pHTL" represents p_HTL*m^-1/4 from eq (16) in Serra-Pompei (2020)
      !
     !  gammaHTL = 0.2 ! ibid
     !  mMax = maxval(m(idxB:nGrid))
     !  do i=idxB, nGrid
     !     if (m(i) .lt. (mMax/betaHTL)) then
     !        pHTL(i) = exp( -(log(m(i)*betaHTL/mMax))**2/4.)
     !     else
     !        pHTL(i) = 1.
     !     end if
     !  end do
     !  pHTL(idxB:nGrid) = pHTL(idxB:nGrid) * m(idxB:nGrid)**(-0.25)
    end if
    end subroutine setHTL
  
  ! ======================================
  !  Calculate rates and derivatives:
  ! ======================================

 
  ! -----------------------------------------------
  !  Calculate the derivatives for all groups:
  !  In:
  !    u: the vector of state variables (nutrients and biomasses)
  !    L: light level
  !    T: temperature
  !    dt: time step for predictor-corrector
  ! 
  !  Uses a simple predictor-corrector scheme.
  !  If one of the nutrients would become negative after an Euler 
  !  time step dt, then the uptake of said nutrient is reduced
  !  by a factor gamma to avoid the nutrient becoming negative.
  ! -----------------------------------------------
  subroutine calcDerivatives(u, L, T, dt)
    real(dp), intent(in):: L, T, dt, u(:)
    integer:: i, j, iGroup
    real(dp):: gammaN, gammaDOC, gammaSi

    !
    ! Use only the positive part of biomasses for calculation of derivatives:
    !
    do i = 1, nGrid
       upositive(i) = max( 0.d0, u(i) )
    end do
    !
    ! Update temperature corrections (in global.f90):
    !
    call updateTemperature(T)
    !
    ! Calc available food:
    !
    do i = idxB, nGrid
       F(i) = 0.d0
       do j = idxB, nGrid
          F(i) = F(i) + theta(i,j)*upositive(j)
       end do
    end do
    ! Calculate feeding for each group:
    do iGroup = 1, nGroups
      call calcFeeding(group(iGroup)%spec, F( ixStart(iGroup):ixEnd(iGroup) ))
    end do 
    !
    ! Calc HTL mortality:
    !
    !if (bQuadraticHTL) then
    !   do i = idxB, nGrid
    !      rates%mortHTL(i) = calcHTL(upositive, i)*upositive(i)
    !   end do
    !else
    !   rates%mortHTL(idxB:nGrid) = mortHTL*pHTL(idxB:nGrid)
    !end if
    !
    ! Calc derivatives of unicellular groups (predictor step)
    !
    gammaN = 1.d0
    gammaDOC = 1.d0
    gammaSi = 1.d0

    call calcDerivativesUnicellulars(upositive, L, gammaN, gammaDOC, gammaSi)
    !
    ! Make a correction if nutrient fields will become less than zero:
    !
    if ((u(idxN) + dudt(idxN)*dt) .lt. 0) then
       gammaN = max(0.d0, min(1.d0, -u(idxN)/(dudt(idxN)*dt)))
    end if
    if ((u(idxDOC) + dudt(idxDOC)*dt) .lt. 0) then
       gammaDOC = max(0.d0, min(1.d0, -u(idxDOC)/(dudt(idxDOC)*dt)))
    end if
    if ((u(idxSi) + dudt(idxSi)*dt) .lt. 0) then
      gammaSi = max(0.d0, min(1.d0, -u(idxSi)/(dudt(idxSi)*dt)))
    end if
    if ((gammaN .lt. 1.d0) .or. (gammaDOC .lt. 1.d0) .or. (gammaSi .lt. 1.d0)) then
       call calcDerivativesUnicellulars(upositive, L, gammaN, gammaDOC, gammaSi)
    end if
    !
    ! Calc derivatives of multicellular groups:
    !
   !  do iGroup = 1, nGroups
   !    if (group(iGroup)%type .eq. typeCopepod) then
   !       call calcDerivativesCopepod(group(iGroup), &
   !          upositive(group(iGroup)%ixStart:group(iGroup)%ixEnd), &
   !          rates)
   !    end if
   !  end do

    contains

     !
  ! Calculate derivatives for unicellular groups
  ! In:
  !   gammaN and gammaDOC are reduction factors [0...1] of uptakes of N and DOC,
  !   used for correction of Euler integration. If no correction is used, just set to 1.0
  !   This correction procedure is needed for correct Euler integration.
  subroutine calcDerivativesUnicellulars(upositive, L, gammaN, gammaDOC, gammaSi)
   real(dp), intent(in):: upositive(:), L, gammaN, gammaDOC, gammaSi
   integer:: i,j, iGroup, jGroup, ixj, ixi
   !
   ! Calc uptakes of all unicellular groups:
   !
   do iGroup = 1, nGroups
      select type (spectrum => group(iGroup)%spec)
      type is (spectrumGeneralists)
         call calcRatesGeneralists(spectrum, &
                     L, upositive(idxN), upositive(idxDOC), gammaN, gammaDOC)
        ! call calcRatesGeneralists(group(iGroup), &
        !      rates, L, upositive(idxN), upositive(idxDOC), gammaN, gammaDOC)
      ! case(typeGeneralist_csp)
      !    call calcRatesGeneralists_csp(group(iGroup), &
      !         rates, L, upositive(idxN),  gammaN)
      !      case(typeDiatom)
      !         call calcRatesDiatoms(group(iGroup), &
      !         rates, L, upositive(idxN), upositive(idxDOC) , upositive(idxSi), &
      !         gammaN, gammaDOC, gammaSi)
      !      case(typeDiatom_simple)
      !         call calcRatesDiatoms_simple(group(iGroup), &
      !         rates, L, upositive(idxN), upositive(idxDOC) , upositive(idxSi), &
      !         gammaN, gammaDOC, gammaSi) 
      end select
   end do
   !
   ! Calc predation mortality
   !
   do iGroup = 1, nGroups
      group(iGroup)%spec%mortpred = 0.d0
      do i = ixStart(iGroup), ixEnd(iGroup)
         ixi = i-ixStart(iGroup)+1
         do jGroup = 1, nGroups
            do j = ixStart(jGroup), ixEnd(jGroup)
               ixj = j-ixStart(jGroup)+1
               if (F(j) .gt. 0.d0) then
                  group(iGroup)%spec%mortpred(ixi) = group(iGroup)%spec%mortpred(ixi) &
                     + theta(j,i) * group(jGroup)%spec%JF(ixj)*upositive(j) &
                     / (group(jGroup)%spec%epsilonF*group(jGroup)%spec%m(ixj)*F(j))
               end if
            end do
         end do
      end do
   end do 
   !
   ! Assemble derivatives:
   !
   dudt(1:(idxB-1)) = 0.d0 ! Set derivatives of nutrients to zero
   
   do iGroup = 1, nGroups
      select type (spectrum => group(iGroup)%spec)
      type is (spectrumGeneralists)
         call calcDerivativesGeneralists(spectrum, &
              upositive(ixStart(iGroup):ixEnd(iGroup)), &
              dudt(idxN), dudt(idxDOC), dudt(ixStart(iGroup):ixEnd(iGroup)))
      ! case (typeGeneralist_csp)
      !    call calcDerivativesGeneralists_csp(group(iGroup),&
      !         upositive(group(iGroup)%ixStart:group(iGroup)%ixEnd), &
      !         rates)
      !    case (typeDiatom)
      !         call calcDerivativesDiatoms(group(iGroup),&
      !         upositive(group(iGroup)%ixStart:group(iGroup)%ixEnd), &
      !         rates)
      !   case (typeDiatom_simple)
      !         call calcDerivativesDiatoms_simple(group(iGroup),&
      !         upositive(group(iGroup)%ixStart:group(iGroup)%ixEnd), &
      !         rates)
                 
      end select
   end do
 end subroutine calcDerivativesUnicellulars

   
  end subroutine calcDerivatives


  ! ======================================
  !  Simulate models:
  ! ======================================

  ! -----------------------------------------------
  ! Simulate a chemostat with Euler integration
  ! Ndeep is a vector with the concentrations of
  ! nutrients in the deep layer.
  ! -----------------------------------------------
  subroutine simulateChemostatEuler(u, L, T, Ndeep, diff, tEnd, dt)
    real(dp), intent(inout):: u(:) ! Initial conditions and result after integration
    real(dp), intent(in):: L      ! Light level
    real(dp), intent(in):: T ! Temperature
    real(dp), intent(in):: Ndeep(nNutrients) ! Nutrients in the deep layer
    real(dp), intent(in):: diff      ! Diffusivity
    real(dp), intent(in):: tEnd ! Time to simulate
    real(dp), intent(in):: dt    ! time step
    integer:: i, iEnd

    iEnd = floor(tEnd/dt)
   
    do i=1, iEnd
       call calcDerivatives(u, L, T, dt)
       dudt(idxN) = dudt(idxN) + diff*(Ndeep(idxN)-u(idxN))
       dudt(idxDOC) = dudt(idxDOC) + diff*(Ndeep(idxDOC) - u(idxDOC))
       if (idxB .gt. idxSi) then
         dudt(idxSi) = dudt(idxSi) + diff*(Ndeep(idxSi) - u(idxSi))
       end if  
       !
       ! Note: should not be done for copepods:
       !
       dudt(idxB:nGrid) = dudt(idxB:nGrid) + diff*(0.d0 - u(idxB:nGrid))
       u = u + dudt*dt
    end do
  end subroutine simulateChemostatEuler

  ! -----------------------------------------------
  ! Simulate with Euler integration
  ! -----------------------------------------------
  subroutine simulateEuler(u, L, T, tEnd, dt)
    real(dp), intent(inout):: u(:) ! Initial conditions and result after integration
    real(dp), intent(in):: L      ! Light level
    real(dp), intent(in):: T ! Temperature
    real(dp), intent(in):: tEnd ! Time to simulate
    real(dp), intent(in):: dt    ! time step
    integer:: i, iEnd

    iEnd = floor(tEnd/dt)

    do i=1, iEnd
       call calcDerivatives(u, L, T, dt)
       u = u + dudt*dt
    end do
  end subroutine simulateEuler


  !=========================================
  ! Diagnostic functions
  !=========================================

   subroutine printRates()
      integer :: iGroup
      99 format (a10, 20f10.6)

      do iGroup = 1, nGroups
         call group(iGroup)%spec%printRates()
      end do
      write(*,99) "dudt", dudt

   end subroutine printRates

  
  function calcN(u) result(N)
   real(dp), intent(in):: u(:)
   integer:: i
   real(dp):: N

   N = 0
   N = u(idxN)
   do i = 1, nGrid
      N = N + u(nNutrients+i)/5.68
   end do
 end function calcN
 
 
 subroutine getMass(m, mDelta)
   real(dp), intent(inout):: m(nGrid), mDelta(nGrid)
   integer:: i

   do i = 1,nGroups
      m(ixStart(i):ixEnd(i)) = group(i)%spec%m
      mDelta(ixStart(i):ixEnd(i)) = group(i)%spec%mDelta
   end do
   end subroutine getMass
  
  ! ---------------------------------------------------
  ! Get the ecosystem functions as calculated from the last call
  ! to calcDerivatives
  ! ---------------------------------------------------
!   subroutine getFunctions(ProdGross, ProdNet,ProdHTL,eHTL,Bpico,Bnano,Bmicro)
!     real(dp), intent(out):: ProdGross, ProdNet,ProdHTL,eHTL,Bpico,Bnano,Bmicro
!     real(dp) :: conversion
!     real(dp) :: ESD(nGrid)
!     integer:: i

!     ProdGross = 0.d0
!     ProdNet = 0.d0
!     ProdHTL = 0.d0
!     Bpico = 0.d0
!     Bnano = 0.d0
!     Bmicro = 0.d0
    
!     conversion = 365.*1d-6*1000. ! Convert to gC/yr/m3
!     do i = 1, nGroups
!        select type (spectrum => group(i)%spec)
!           type is (spectrumGeneralists)
!             ProdGross = ProdGross + conversion * &
!                sum(  spectrum%JL * upositive(idxB:nGrid) / spectrum%m )
          
!           ProdNet = ProdNet + conversion * &
!                getProdNetGeneralists(spectrum,  upositive(group(i)%spec%ixStart:group(i)%spec%ixEnd))
!        end select
!     end do

!     ESD = 10000. * 1.5 * (m*1d-6)**onethird
!     conversion = 1d-6*1000 ! Convert to gC/m3
!     do i = idxB, nGrid
!        if (ESD(i) .le. 2.) then
!           Bpico = Bpico + conversion*upositive(i)
!        endif
       
!        if ((ESD(i).gt.2.) .and. (ESD(i) .le. 20.)) then
!           Bnano = Bnano + conversion*upositive(i)
!        endif

!        if (ESD(i) .gt. 20.) then
!           Bmicro = Bmicro + conversion*upositive(i)
!        endif

!        ProdHTL = ProdHTL + 365*conversion*mortHTL(i)*upositive(i)
!     end do

!     eHTL = eHTL / ProdNet
!     if (eHTL .gt. 1) then
!        eHTL = -1.
!     end if
!   end subroutine getFunctions

  ! ---------------------------------------------------
  ! Returns mass conservation calculated from last call to calcDerivatives
  ! ---------------------------------------------------
!   subroutine getBalance(Nbalance,Cbalance)
!    real(dp), intent(out):: Nbalance, Cbalance
!    integer:: i
   
!    i = 1 ! Do it only for the first group, whihc we assume are generalists
!       if (group(i)%type .eq. typeGeneralist) then
!          Nbalance = getNbalanceGeneralists(group(i),  upositive(ixStart(i):ixEnd(i)))  
!          Cbalance = getCbalanceGeneralists(group(i),  upositive(ixStart(i):ixEnd(i)))
    
!       end if   
 
! end subroutine getBalance
!   ! ---------------------------------------------------
!   ! Returns the rates calculated from last call to calcDerivatives
!   ! ---------------------------------------------------
  subroutine getRates(jN, jDOC, jL, jSi, jF, jFreal,&
    jTot, jMax, jFmax, jR, jLossPassive, &
    jNloss,jLreal, &
    mortpred, mortHTL, mort2, mort)
    use globals
    real(dp), intent(out):: jN(nGrid-nNutrients), jDOC(nGrid-nNutrients), jL(nGrid-nNutrients)
    real(dp), intent(out):: jSi(nGrid-nNutrients)
    real(dp), intent(out):: jF(nGrid-nNutrients), jFreal(nGrid-nNutrients)
    real(dp), intent(out):: jTot(nGrid-nNutrients), jMax(nGrid-nNutrients), jFmax(nGrid-nNutrients)
    real(dp), intent(out):: jR(nGrid-nNutrients)
    real(dp), intent(out):: jLossPassive(nGrid-nNutrients), jNloss(nGrid-nNutrients), jLreal(nGrid-nNutrients)
    real(dp), intent(out):: mortpred(nGrid-nNutrients), mortHTL(nGrid-nNutrients)
    real(dp), intent(out):: mort2(nGrid-nNutrients), mort(nGrid-nNutrients)
   integer :: iGroup, i1, i2

   do iGroup = 1, nGroups
      i1 = ixStart(iGroup)-nNutrients
      i2 = ixEnd(iGroup)-nNutrients
      ! Extract common fields:
      jF( i1:i2 ) = group(iGroup)%spec%flvl * group(iGroup)%spec%JFmax / group(iGroup)%spec%m
      jFreal( i1:i2 ) = group(iGroup)%spec%JF / group(iGroup)%spec%m
      jFmax( i1:i2 ) = group(iGroup)%spec%JFmax / group(iGroup)%spec%m
      Jtot( i1:i2 ) = group(iGroup)%spec%Jtot / group(iGroup)%spec%m
      mortpred( i1:i2 ) = group(iGroup)%spec%mortpred
      mortHTL( i1:i2 ) = group(iGroup)%spec%mortHTL
      mort2( i1:i2 ) = group(iGroup)%spec%mort2
      jNloss( i1:i2 ) = group(iGroup)%spec%JNloss / group(iGroup)%spec%m
      jR( i1:i2 ) = group(iGroup)%spec%Jresp / group(iGroup)%spec%m

      select type (spectrum => group(iGroup)%spec)
      class is (spectrumUnicellular)
        jN( i1:i2 ) = spectrum%JN / spectrum%m
        jDOC( i1:i2 ) = spectrum%JDOC / spectrum%m
        jL( i1:i2 ) = spectrum%JL / spectrum%m
        jMax( i1:i2 ) = spectrum%Jmax / spectrum%m
        jLossPassive( i1:i2 ) = spectrum%JlossPassive / spectrum%m
        jLreal( i1:i2 ) = spectrum%JLreal / spectrum%m

      end select
   end do
  end subroutine getRates

end module NUMmodel
