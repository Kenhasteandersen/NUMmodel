
!
! Module to handle the NUM model framework
!
module NUMmodel
  use globals
  use spectrum
  use generalists
  use generalists_csp
  use diatoms
  use diatoms_simple
  use copepods
  use debug
  implicit none

  !
  ! Variables that contain the size spectrum groups
  !
  integer:: nGroups ! Number of groups
  integer:: iCurrentGroup ! The current group to be added
  integer:: nNutrients ! Number of nutrient state variables
  integer:: idxB ! First index into non-nutrient groups (=nNutrients+1)
  type(typeSpectrum), dimension(:), allocatable:: group ! Structure for each group

  real(dp), dimension(:,:), allocatable:: theta ! Interaction matrix
  real(dp), dimension(:), allocatable:: palatability ! Palatability of each size group
  real(dp), dimension(:), allocatable:: upositive ! State variable constrained to be positive
  !
  ! Variables for HTL mortalities:
  !
  real(dp), dimension(:), allocatable:: pHTL ! Selectivity function for HTL mortality
  real(dp) :: mortHTL ! Level of HTL mortality (see below for override)
  real(dp):: gammaHTL ! Parameter for quadratic HTL mortality
  logical:: bQuadraticHTL ! Boolean flag to signify whether mortality is standard or "quadratic"
  ! Defaults to true; can be overridden but parameters must then be set with a
  ! call to parametersFinalize()

  type(typeRates):: rates

  real(dp), dimension(:), allocatable:: m, z, beta, sigma, AF, JFmax, epsilonF ! Feeding parameters

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
  subroutine setupGeneralistsOnly_csp()
    call parametersInit(1, 10, 2) ! 1 group, 10 size classes (excl nutrients and DOC)
    call parametersAddGroup(typeGeneralist_csp, 10, 0.1d0) ! generalists with 10 size classes
    call parametersFinalize(0.003d0, .true.) ! Serra-Pompei (2020))
  end subroutine setupGeneralistsOnly_csp

  ! -----------------------------------------------
  ! A basic setup with only diatoms:
  ! -----------------------------------------------
  subroutine setupDiatomsOnly(n)
   integer, intent(in):: n
   call parametersInit(1, n, 3) ! 1 group, n size classes (excl nutrients)
   call parametersAddGroup(typeDiatom, n, 1.d0) ! diatoms with n size classes
   call parametersFinalize(0.2d0, .false.)
 end subroutine setupDiatomsOnly

 ! -----------------------------------------------
  ! A basic setup with only simple diatoms:
  ! -----------------------------------------------
 subroutine setupDiatoms_simpleOnly(n)
   integer, intent(in):: n
   call parametersInit(1, n, 3) ! 1 group, n size classes (excl nutrients)
   call parametersAddGroup(typeDiatom_simple, n, 1.d0) ! diatoms with n size classes
   call parametersFinalize(0.2d0, .false.)
 end subroutine setupDiatoms_simpleOnly
 
  ! -----------------------------------------------
  ! Generalists and diatoms:
  ! -----------------------------------------------
   subroutine setupGeneralistsDiatoms(n)
      integer, intent(in):: n
      call parametersInit(2, 2*n, 3)
      call parametersAddGroup(typeGeneralist, n, 1.d0) ! generalists with n size classes
      call parametersAddGroup(typeDiatom, n, 1.d0) ! diatoms with n size classes
      call parametersFinalize(.2d0, .false.)
   end subroutine setupGeneralistsDiatoms
 
   subroutine setupGeneralistsDiatoms_simple(n)
      integer, intent(in):: n
      call parametersInit(2, 2*n, 3)
      call parametersAddGroup(typeGeneralist, n, 1.d0) ! generalists with n size classes
      call parametersAddGroup(typeDiatom_simple, n, 1.d0) ! diatoms with n size classes
      call parametersFinalize(0.2d0, .false.)
   end subroutine setupGeneralistsDiatoms_simple
 
  ! -----------------------------------------------
  ! A basic setup with generalists and 1 copepod
  ! -----------------------------------------------
  subroutine setupGeneralistsCopepod()
    call parametersInit(2, 20, 2)
    call parametersAddGroup(typeGeneralist, 10, 0.1d0)
    call parametersAddGroup(typeCopepod, 10, .1d0) ! add copepod with adult mass .1 mugC
    call parametersFinalize(0.003d0, .true.) ! Use quadratic mortality
  end subroutine setupGeneralistsCopepod

  ! -----------------------------------------------
  ! A generic setup with generalists and a number of copepod species
  ! -----------------------------------------------
  subroutine setupGeneric(mAdult)
    real(dp), intent(in):: mAdult(:)
    integer, parameter:: n = 10 ! number of size classes in each group
    integer:: iCopepod

    call parametersInit(size(mAdult)+1, n*(size(mAdult)+1), 2)
    call parametersAddGroup(typeGeneralist, n, 0.1d0)
    if ( size(mAdult) .eq. 0) then
       bQuadraticHTL = .false. ! Use standard "linear" mortality
    else
       do iCopepod = 1, size(mAdult)
          call parametersAddGroup(typeCopepod, n, mAdult(iCopepod)) ! add copepod
       end do
    end if
    call parametersFinalize(0.003d0, .true.)
  end subroutine setupGeneric

  ! -----------------------------------------------
  ! A generic setup with generalists and a number of copepod species
  ! -----------------------------------------------
  subroutine setupGeneric_csp(mAdult)
    real(dp), intent(in):: mAdult(:)
    integer, parameter:: n = 10 ! number of size classes in each group
    integer:: iCopepod

    call parametersInit(size(mAdult)+1, n*(size(mAdult)+1), 2)
    call parametersAddGroup(typeGeneralist_csp, n, 0.1d0)
    do iCopepod = 1, size(mAdult)
       call parametersAddGroup(typeCopepod, n, mAdult(iCopepod)) ! add copepod
    end do
    call parametersFinalize(0.003d0, .true.)
  end subroutine setupGeneric_csp



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
    if (allocated(m)) then
       deallocate(m)
       deallocate(palatability)
       deallocate(z)
       deallocate(beta)
       deallocate(sigma)
       deallocate(AF)
       deallocate(JFmax)
       deallocate(epsilonF)

       deallocate(group)

       deallocate(upositive)

       ! Interaction matrix:
       deallocate(theta)
       !
       ! Deallocate rates:
       !
       deallocate(rates%dudt)

       deallocate(rates%F)
       deallocate(rates%flvl)
       deallocate(rates%JF)
       deallocate(rates%JEnc)

       deallocate(rates%JN)
       deallocate(rates%JL)
       deallocate(rates%JDOC)
       deallocate(rates%JSi)
       deallocate(rates%JNtot)
       deallocate(rates%JCtot)
       deallocate(rates%Jtot)
       deallocate(rates%JLreal)
       deallocate(rates%JCloss_feeding)
       deallocate(rates%JCloss_photouptake)
       deallocate(rates%JNlossLiebig)
       deallocate(rates%JClossLiebig)
       deallocate(rates%JNloss)
       deallocate(rates%JCloss)
       deallocate(rates%JSiloss)

       deallocate(rates%mortpred)
       deallocate(rates%mortHTL)
       deallocate(pHTL)

       deallocate(rates%g)
       deallocate(rates%mortStarve)
       deallocate(rates%mort)

    end if

    allocate(m(nGrid))
    m = 0.d0
    allocate(palatability(nGrid))
    palatability = 1.d0 ! Default
    allocate(z(nGrid))
    allocate(beta(nGrid))
    allocate(sigma(nGrid))
    allocate(AF(nGrid))
    allocate(JFmax(nGrid))
    allocate(epsilonF(nGrid))

    allocate(group(nGroups))

    allocate(upositive(nGrid))

    ! Interaction matrix:
    allocate(theta(nGrid,nGrid))
    !
    ! Allocate rates:
    !
    allocate(rates%dudt(nGrid))

    allocate(rates%F(nGrid))
    allocate(rates%flvl(nGrid))
    allocate(rates%JF(nGrid))
    allocate(rates%JEnc(nGrid))

    allocate(rates%JN(nGrid))
    allocate(rates%JL(nGrid))
    allocate(rates%JDOC(nGrid))
    allocate(rates%JSi(nGrid))
    allocate(rates%JNtot(nGrid))
    allocate(rates%JCtot(nGrid))
    allocate(rates%Jtot(nGrid))
    allocate(rates%JLreal(nGrid))
    allocate(rates%JCloss_feeding(nGrid))
    allocate(rates%JCloss_photouptake(nGrid))
    allocate(rates%JNlossLiebig(nGrid))
    allocate(rates%JClossLiebig(nGrid))
    allocate(rates%JNloss(nGrid))
    allocate(rates%JCloss(nGrid))
    allocate(rates%JSiloss(nGrid))

    allocate(rates%mortpred(nGrid))
    allocate(rates%mortHTL(nGrid))
    allocate(pHTL(nGrid))

    allocate(rates%g(nGrid))
    allocate(rates%mortStarve(nGrid))
    allocate(rates%mort(nGrid))

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

    !
    ! Find the group number and grid location:
    !
    iCurrentGroup = iCurrentGroup + 1
    if (iCurrentGroup.eq.1) then
       group(iCurrentGroup)%ixStart = idxB
    else
       group(iCurrentGroup)%ixStart = group(iCurrentGroup-1)%ixEnd+1
    end if
    group(iCurrentGroup)%ixEnd = group(iCurrentGroup)%ixStart+n-1

    if (group(iCurrentGroup)%ixEnd .gt. nGrid) then
       write(6,*) 'Attempting to add more grid points than allocated', &
           group(iCurrentGroup)%ixEnd, nGrid
       stop 1
    end if
    !
    ! Add the group
    !
    select case (typeGroup)
    case (typeGeneralist)
      group(iCurrentGroup) = initGeneralists(n, group(iCurrentGroup)%ixStart-1, mMax)
    case (typeGeneralist_csp)
      group(iCurrentGroup) = initGeneralists_csp(n, group(iCurrentGroup)%ixStart-1, mMax)
    case (typeDiatom)
      group(iCurrentGroup) = initDiatoms(n, group(iCurrentGroup)%ixStart-1, mMax)
      palatability(group(iCurrentGroup)%ixStart:group(iCurrentGroup)%ixEnd) = 0.5d0 ! Lower palatability for diatoms
    case (typeDiatom_simple)
      group(iCurrentGroup) = initDiatoms_simple(n, group(iCurrentGroup)%ixStart-1, mMax)
      palatability(group(iCurrentGroup)%ixStart:group(iCurrentGroup)%ixEnd) = 0.5d0 ! Lower palatability for diatoms
    case(typeCopepod)
       group(iCurrentGroup) = initCopepod(n, group(iCurrentGroup)%ixStart-1, mMax)
    end select
    group(iCurrentGroup)%type = typeGroup
    !
    ! Import grid to globals:
    !
    m(group(iCurrentGroup)%ixStart : group(iCurrentGroup)%ixEnd) = group(iCurrentGroup)%m
    z(group(iCurrentGroup)%ixStart : group(iCurrentGroup)%ixEnd) = group(iCurrentGroup)%z
    !
    ! Import feeding parameters:
    !
    beta(group(iCurrentGroup)%ixStart : group(iCurrentGroup)%ixEnd) = group(iCurrentGroup)%beta
    sigma(group(iCurrentGroup)%ixStart : group(iCurrentGroup)%ixEnd) = group(iCurrentGroup)%sigma
    AF(group(iCurrentGroup)%ixStart : group(iCurrentGroup)%ixEnd) = group(iCurrentGroup)%AF
    JFmax(group(iCurrentGroup)%ixStart : group(iCurrentGroup)%ixEnd) = group(iCurrentGroup)%JFmax
    epsilonF(group(iCurrentGroup)%ixStart : group(iCurrentGroup)%ixEnd) = group(iCurrentGroup)%epsilonF

  end subroutine parametersAddGroup
  ! -----------------------------------------------
  !  Finalize the setting of parameters. Must be called when
  !  all groups have been added.
  ! -----------------------------------------------
  subroutine parametersFinalize(mortHTL, bQuadraticHTL)
    real(dp), intent(in):: mortHTL
    logical, intent(in):: bQuadraticHTL
    integer:: i,j
    real(dp):: betaHTL, mHTL, mMax
    !
    ! Calc theta:
    !
    do i = idxB, nGrid
       do j = idxB, nGrid
          !theta(i,j) = exp( -(log(m(i)/m(j)/beta(i)))**2/(2*sigma(i)**2))
          theta(i,j) = palatability(j) * calcPhi(m(i)/m(j), beta(i), sigma(i), z(i))
       end do
    end do
    !
    ! Set HTL mortality
    !
    call parametersHTL(mortHTL, bQuadraticHTL)

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

  !
  ! Sets the HTL mortality. 
  ! In:
  !  mortHTL - the level of mortality
  !  bQuadraticHTL - whether to use a "quadratic" mortality
  !
  subroutine parametersHTL(mortalityHTL, boolQuadraticHTL)
   real(dp), intent(in):: mortalityHTL
   logical, intent(in):: boolQuadraticHTL
   real(dp):: betaHTL, mMax, mHTL
   integer:: i

   mortHTL = mortalityHTL
   bQuadraticHTL = boolQuadraticHTL
   !     
   ! Calc htl mortality
   !
   betaHTL = 500.
   if (.not. bQuadraticHTL) then
     !
     ! Standard HTL mortality. In this case "pHTL" represent the selectivity of HTL mortality
     !
     mHTL = m(nGrid)/betaHTL**1.5  ! Bins affected by HTL mortality
     pHTL(idxB:nGrid) = (1/(1+(m(idxB:nGrid)/mHTL)**(-2)))
   else
     !
     ! Linear HTL mortality (commonly referred to as "quadratic")
     ! In this case "pHTL" represents p_HTL*m^-1/4 from eq (16) in Serra-Pompei (2020)
     !
     gammaHTL = 0.2 ! ibid
     mMax = maxval(m(idxB:nGrid))
     do i=idxB, nGrid
        if (m(i) .lt. (mMax/betaHTL)) then
           pHTL(i) = exp( -(log(m(i)*betaHTL/mMax))**2/4.)
        else
           pHTL(i) = 1.
        end if
     end do
     pHTL(idxB:nGrid) = pHTL(idxB:nGrid) * m(idxB:nGrid)**(-0.25)
   end if
   end subroutine parametersHTL

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
  !  time step
  !  with length dt, then the uptake of said nutrient is reduced
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
    ! Calc uptakes of food
    !
    do i = idxB, nGrid
       rates%F(i) = 0.d0
       do j = idxB, nGrid
          rates%F(i) = rates%F(i) + theta(i,j)*upositive(j)
       end do
    end do

    rates%flvl(idxB:nGrid) = AF(idxB:nGrid)*rates%F(idxB:nGrid) / (AF(idxB:nGrid)*rates%F(idxB:nGrid) + fTemp2*JFmax(idxB:nGrid))
    rates%JF(idxB:nGrid) = rates%flvl(idxB:nGrid) * fTemp2*JFmax(idxB:nGrid)
    !
    ! Calc HTL mortality:
    !
    if (bQuadraticHTL) then
       do i = idxB, nGrid
          rates%mortHTL(i) = calcHTL(upositive, i)*upositive(i)
       end do
    else
       rates%mortHTL(idxB:nGrid) = mortHTL*pHTL(idxB:nGrid)
    end if
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
    if ((u(idxN) + rates%dudt(idxN)*dt) .lt. 0) then
       gammaN = max(0.d0, min(1.d0, -u(idxN)/(rates%dudt(idxN)*dt)))
    end if
    if ((u(idxDOC) + rates%dudt(idxDOC)*dt) .lt. 0) then
       gammaDOC = max(0.d0, min(1.d0, -u(idxDOC)/(rates%dudt(idxDOC)*dt)))
    end if
    if ((u(idxSi) + rates%dudt(idxSi)*dt) .lt. 0) then
      gammaSi = max(0.d0, min(1.d0, -u(idxSi)/(rates%dudt(idxSi)*dt)))
    end if
    if ((gammaN .lt. 1.d0) .or. (gammaDOC .lt. 1.d0) .or. (gammaSi .lt. 1.d0)) then
       call calcDerivativesUnicellulars(upositive, L, gammaN, gammaDOC, gammaSi)
    end if
    !
    ! Calc derivatives of multicellular groups:
    !
    do iGroup = 1, nGroups
      if (group(iGroup)%type .eq. typeCopepod) then
         call calcDerivativesCopepod(group(iGroup), &
            upositive(group(iGroup)%ixStart:group(iGroup)%ixEnd), &
            rates)
      end if
    end do

    contains

     !
  ! Calculate derivatives for unicellular groups
  ! In:
  !   gammaN and gammaDOC are reduction factors [0...1] of uptakes of N and DOC,
  !   used for correction of Euler integration. If no correction is used, just set to 1.0
  !   This correction procedure is needed for correct Euler integration.
  subroutine calcDerivativesUnicellulars(upositive, L, gammaN, gammaDOC, gammaSi)
   real(dp), intent(in):: upositive(:), L, gammaN, gammaDOC, gammaSi
   integer:: i,j, iGroup
   !
   ! Calc uptakes of all unicellular groups:
   !
   do iGroup = 1, nGroups
      select case (group(iGroup)%type)
      case (typeGeneralist)
         call calcRatesGeneralists(group(iGroup), &
              rates, L, upositive(idxN), upositive(idxDOC), gammaN, gammaDOC)
      case(typeGeneralist_csp)
         call calcRatesGeneralists_csp(group(iGroup), &
              rates, L, upositive(idxN),  gammaN)
           case(typeDiatom)
              call calcRatesDiatoms(group(iGroup), &
              rates, L, upositive(idxN), upositive(idxDOC) , upositive(idxSi), &
              gammaN, gammaDOC, gammaSi)
           case(typeDiatom_simple)
              call calcRatesDiatoms_simple(group(iGroup), &
              rates, L, upositive(idxN), upositive(idxDOC) , upositive(idxSi), &
              gammaN, gammaDOC, gammaSi) 
      end select
   end do
   !
   ! Calc predation mortality
   !
   do i=idxB, nGrid
      rates%mortpred(i) = 0.d0
      do j=idxB, nGrid
         if (rates%F(j) .ne. 0.d0) then
            rates%mortpred(i) = rates%mortpred(i)  &
                 + theta(j,i) * rates%JF(j)*upositive(j)/(epsilonF(j)*m(j)*rates%F(j))
         end if
      end do
   end do
   !
   ! Assemble derivatives:
   !
   rates%dudt(1:(idxB-1)) = 0.d0 ! Set derivatives of nutrients to zero
   
   do iGroup = 1, nGroups
      select case (group(iGroup)%type)
      case (typeGeneralist)
         call calcDerivativesGeneralists(group(iGroup),&
              upositive(group(iGroup)%ixStart:group(iGroup)%ixEnd), &
              rates)
      case (typeGeneralist_csp)
         call calcDerivativesGeneralists_csp(group(iGroup),&
              upositive(group(iGroup)%ixStart:group(iGroup)%ixEnd), &
              rates)
         case (typeDiatom)
              call calcDerivativesDiatoms(group(iGroup),&
              upositive(group(iGroup)%ixStart:group(iGroup)%ixEnd), &
              rates)
        case (typeDiatom_simple)
              call calcDerivativesDiatoms_simple(group(iGroup),&
              upositive(group(iGroup)%ixStart:group(iGroup)%ixEnd), &
              rates)
                 
      end select
   end do
 end subroutine calcDerivativesUnicellulars

   
  end subroutine calcDerivatives
  !
  ! Returns the htl mortality divided by h for size bin i
  !
  function calcHTL(u, i) result(mHTL)
    real(dp) :: mHTL, B
    real(dp), intent(in):: u(:)
    integer, intent(in):: i
    integer:: j

    !
    ! First calc the biomass within a size range:
    !
    B = 0.d0
    do j = 1, nGrid
       if (m(j)>m(i)/3.16 .and. m(j)<m(i)*3.16) then
          B = B + u(j)
       end if
    end do

    mHTL = pHTL(i)*mortHTL/z(i)*u(i)**gammaHTL*B**(1.-gammaHTL)
  end function calcHTL

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
       rates%dudt(idxN) = rates%dudt(idxN) + diff*(Ndeep(idxN)-u(idxN))
       rates%dudt(idxDOC) = rates%dudt(idxDOC) + diff*(Ndeep(idxDOC) - u(idxDOC))
       if (idxB .gt. idxSi) then
         rates%dudt(idxSi) = rates%dudt(idxSi) + diff*(Ndeep(idxSi) - u(idxSi))
       end if  
       !
       ! Note: should not be done for copepods:
       !
       rates%dudt(idxB:nGrid) = rates%dudt(idxB:nGrid) + diff*(0.d0 - u(idxB:nGrid))
       u = u + rates%dudt*dt
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
       u = u + rates%dudt*dt
    end do
  end subroutine simulateEuler


  !=========================================
  ! Diagnostic functions
  !=========================================

  
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
 
 
 subroutine getMass(m_, mDelta)
   real(dp), intent(inout):: m_(nGrid), mDelta(nGrid)
   integer:: i

   m_ = m;
   do i = 1,nGroups
      mDelta(group(i)%ixStart:group(i)%ixEnd) = group(i)%mDelta
   end do
   end subroutine getMass
  
  ! ---------------------------------------------------
  ! Get the ecosystem functions as calculated from the last call
  ! to calcDerivatives
  ! ---------------------------------------------------
  subroutine getFunctions(ProdGross, ProdNet,ProdHTL,eHTL,Bpico,Bnano,Bmicro)
    real(dp), intent(out):: ProdGross, ProdNet,ProdHTL,eHTL,Bpico,Bnano,Bmicro
    real(dp) :: conversion
    real(dp) :: ESD(nGrid)
    integer:: i

    ProdGross = 0.d0
    ProdNet = 0.d0
    ProdHTL = 0.d0
    Bpico = 0.d0
    Bnano = 0.d0
    Bmicro = 0.d0
    
    conversion = 365.*1d-6*1000. ! Convert to gC/yr/m3
    do i = 1, nGroups
       if (group(i)%type .eq. typeGeneralist) then
          ProdGross = ProdGross + conversion * &
               sum(  rates%JL(idxB:nGrid) * upositive(idxB:nGrid) / m(idxB:nGrid) )
          
          ProdNet = ProdNet + conversion * &
               getProdNetGeneralists(group(i),  upositive(group(i)%ixStart:group(i)%ixEnd), rates)
       end if
    end do

    ESD = 10000. * 1.5 * (m*1d-6)**onethird
    conversion = 1d-6*1000 ! Convert to gC/m3
    do i = idxB, nGrid
       if (ESD(i) .le. 2.) then
          Bpico = Bpico + conversion*upositive(i)
       endif
       
       if ((ESD(i).gt.2.) .and. (ESD(i) .le. 20.)) then
          Bnano = Bnano + conversion*upositive(i)
       endif

       if (ESD(i) .gt. 20.) then
          Bmicro = Bmicro + conversion*upositive(i)
       endif

       ProdHTL = ProdHTL + 365*conversion*rates%mortHTL(i)*upositive(i)
    end do

    eHTL = eHTL / ProdNet
    if (eHTL .gt. 1) then
       eHTL = -1.
    end if
  end subroutine getFunctions

  ! ---------------------------------------------------
  ! Returns the rates calculated from last call to calcDerivatives
  ! ---------------------------------------------------
  subroutine getRates(jN, jDOC, jL, jSi, jF, jFreal,&
    jTot, jMax, jFmaxx, jR, jLossPassive, &
    jNloss,jLreal, &
    mortpred, mortHTL, mort2, mort)
    use globals
    real(dp), intent(out):: jN(nGrid-nNutrients), jDOC(nGrid-nNutrients), jL(nGrid-nNutrients)
    real(dp), intent(out):: jSi(nGrid-nNutrients)
    real(dp), intent(out):: jF(nGrid-nNutrients), jFreal(nGrid-nNutrients)
    real(dp), intent(out):: jTot(nGrid-nNutrients), jMax(nGrid-nNutrients), jFmaxx(nGrid-nNutrients)
    real(dp), intent(out):: jR(nGrid-nNutrients)
    real(dp), intent(out):: jLossPassive(nGrid-nNutrients), jNloss(nGrid-nNutrients), jLreal(nGrid-nNutrients)
    real(dp), intent(out):: mortpred(nGrid-nNutrients), mortHTL(nGrid-nNutrients)
    real(dp), intent(out):: mort2(nGrid-nNutrients), mort(nGrid-nNutrients)

    jN = rates%JN(idxB:nGrid) / m(idxB:nGrid)
    jDOC = rates%JDOC(idxB:nGrid) / m(idxB:nGrid)
    jL = rates%JL(idxB:nGrid) / m(idxB:nGrid)
    jSi = rates%JSi(idxB:nGrid) / m(idxB:nGrid)
    jF = rates%flvl(idxB:nGrid) * JFmax(idxB:nGrid)/ m(idxB:nGrid)
    jFreal = rates%JF(idxB:nGrid) / m(idxB:nGrid)
    jTot = rates%Jtot(idxB:nGrid) / m(idxB:nGrid)
    jMax = 1.5 + 0*m(idxB:nGrid) ! NOTE: HARDCODED. Should be taken from generalists
    jFmaxx = JFmax(idxB:nGrid) / m(idxB:nGrid)
    jR = 1.5*0.1 + 0*m(idxB:nGrid)  ! NOTE: HARDCODED. Should be taken from generalists
    jLossPassive = 0* m(idxB:nGrid)  ! NOTE: HARDCODED. Should be taken from generalists
    jNloss = rates%JNloss(idxB:nGrid) / m(idxB:nGrid)
    jLreal = rates%JLreal(idxB:nGrid) / m(idxB:nGrid)
    mortpred = rates%mortpred(idxB:nGrid)
    mortHTL = rates%mortHTL(idxB:nGrid)
    mort2 = 0.0002*(nGrid-nNutrients)*upositive(idxB:nGrid)  ! NOTE: HARDCODED. Should be taken from generalists
    mort = 0*m(idxB:nGrid)  ! NOTE: HARDCODED. Should be taken from generalists
  end subroutine getRates

end module NUMmodel
