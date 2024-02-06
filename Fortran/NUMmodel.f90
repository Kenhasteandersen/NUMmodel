!
! The NUM model framework
! 
! This module consists of three parts:
!  1) Setting up a new simulation. This handled by all the subroutines starting with
!     setupXXX and parametersXXX.
!  2) The core routines to make a derivative and simulate with Euler integration
!  3) The remaining subroutines are getting information out about parameters or the simulation
!
module NUMmodel
  use iso_c_binding, only: c_char, c_null_char
  use globals
  use spectrum
  use generalists
  use generalists_simple
  use diatoms_simple
  use copepods
  use diatoms
  use POM
  use read_input_module
  implicit none

  ! Indices into the state-variable vector:
  integer, parameter :: idxN = 1
  integer, parameter :: idxDOC = 2
  integer, parameter :: idxSi = 3

  ! Types of spectra:
  integer, parameter :: typeGeneralistSimple = 1
  integer, parameter :: typeGeneralist = 5  
  integer, parameter :: typeDiatom = 3
  integer, parameter :: typeDiatom_simple = 4
  integer, parameter :: typeCopepodActive = 10
  integer, parameter :: typeCopepodPassive = 11
  integer, parameter :: typePOM = 100
  !
  ! Variables that contain the size spectrum groups
  !
  integer:: nGroups ! Number of groups
  integer:: iCurrentGroup ! The current group to be added (used in the parametersXX subroutines)
  integer:: nNutrients ! Number of nutrient state variables
  integer:: idxB ! First index into non-nutrient groups (=nNutrients+1)
  integer:: nGrid ! Total number of state variables in u incl. variables for nutrients
  type(spectrumContainer), allocatable :: group(:) ! Structure pointing to each group
  integer, dimension(:), allocatable :: ixStart, ixEnd ! Indices into u for each group

  real(dp), dimension(:,:), allocatable:: theta ! Interaction matrix
  real(dp), dimension(:), allocatable:: upositive ! State variable constrained to be positive
  real(dp), dimension(:), allocatable:: F ! Available food
  !
  ! Variables for HTL mortalities:
  !
  real(dp), dimension(:), allocatable:: pHTL ! Selectivity function for HTL mortality
  logical:: bQuadraticHTL ! Boolean flag to signify whether mortality is standard or "quadratic"
                          ! Defaults to true; can be overridden but parameters must then be set with a
                          ! call to parametersFinalize()
  !
  ! Variables for POM:
  !
  integer :: idxPOM = 0     ! Index to the POM group:
  integer, dimension(:), allocatable:: thetaPOM ! Which POM size class does
                                       ! each size class in u deliver POM to

contains

  ! ======================================
  !  Various model setups
  ! ======================================

  ! -----------------------------------------------
  ! A basic setup with only simple generalists
  ! -----------------------------------------------
  subroutine setupGeneralistsSimpleOnly(n,errorio,errorstr)
    integer, intent(in):: n
    logical(1), intent(out):: errorio ! Whether to losses to the deep
    character(c_char), dimension(*) :: errorstr
    
    call parametersInit(1, n, 2,errorio,errorstr) ! 1 group, n size classes (excl nutrients and DOC)
      IF ( errorio ) RETURN 
    call parametersAddGroup(typeGeneralistSimple, n, 0.d0,errorio,errorstr) ! generalists with n size classes
    call parametersFinalize(0.1d0, .false., .false.) ! Use standard "linear" mortality
  end subroutine setupGeneralistsSimpleOnly
  
  ! -----------------------------------------------
  ! A basic setup with only generalists
  ! -----------------------------------------------
  subroutine setupGeneralistsOnly(n,errorio,errorstr)
    integer, intent(in):: n
    logical(1), intent(out):: errorio ! Whether to losses to the deep
    character(c_char), dimension(*) :: errorstr
    call parametersInit(1, n, 2,errorio,errorstr) ! 1 group, n size classes (excl nutrients and DOC)
      IF ( errorio ) RETURN 
    call parametersAddGroup(typeGeneralist, n, 0.0d0,errorio,errorstr) ! generalists with n size classes
    call parametersFinalize(0.1d0, .false., .false.) ! Use standard "linear" mortality
  end subroutine setupGeneralistsOnly

  ! -----------------------------------------------
  ! A basic setup with generalists simple and POM
  ! -----------------------------------------------
  subroutine setupGeneralistsSimplePOM(n, nPOM,errorio,errorstr)
    integer, intent(in):: n, nPOM
    logical(1), intent(out):: errorio ! Whether to losses to the deep
    character(c_char), dimension(*) :: errorstr
    call parametersInit(2, n+nPOM, 2,errorio,errorstr) ! 2 groups, n+nPOM size classes (excl nutrients and DOC)
      IF ( errorio ) RETURN 
    call parametersAddGroup(typeGeneralistSimple, n, 0.d0,errorio,errorstr) ! generalists with n size classes
      IF ( errorio ) RETURN 
    call parametersAddGroup(typePOM, nPOM, 1.0d0,errorio,errorstr) ! POM with nPOM size classes and max size 1 ugC
    call parametersFinalize(0.1d0, .false., .false.) ! Use standard "linear" mortality
  end subroutine setupGeneralistsSimplePOM
  
  
    ! -----------------------------------------------
  ! A basic setup with generalists and POM
  ! -----------------------------------------------
  subroutine setupGeneralistsPOM(n, nPOM,errorio,errorstr)
    integer, intent(in):: n, nPOM
    logical(1), intent(out):: errorio ! Whether to losses to the deep
    character(c_char), dimension(*) :: errorstr
    call parametersInit(2, n+nPOM, 2,errorio,errorstr) ! 2 groups, n+nPOM size classes (excl nutrients and DOC)
      IF ( errorio ) RETURN 
    call parametersAddGroup(typeGeneralist, n, 0.d0,errorio,errorstr) ! generalists with n size classes
      IF ( errorio ) RETURN 
    call parametersAddGroup(typePOM, nPOM, 1.0d0,errorio,errorstr) ! POM with nPOM size classes and max size 1 ugC
    call parametersFinalize(0.1d0, .false., .false.) ! Use standard "linear" mortality
  end subroutine setupGeneralistsPOM

  ! -----------------------------------------------
  ! A basic setup with only diatoms:
  ! -----------------------------------------------
  subroutine setupDiatomsOnly(n,errorio,errorstr)
    integer, intent(in):: n
    logical(1), intent(out):: errorio ! Whether to losses to the deep
    character(c_char), dimension(*) :: errorstr
    call parametersInit(1, n, 3,errorio,errorstr) ! 1 group, n size classes (excl nutrients)
      IF ( errorio ) RETURN 
    call parametersAddGroup(typeDiatom, n, 1.d0,errorio,errorstr) ! diatoms with n size classes
    call parametersFinalize(0.1d0, .false., .false.)
  end subroutine setupDiatomsOnly

  ! -----------------------------------------------
  ! A basic setup with only simple diatoms:
  ! -----------------------------------------------
  subroutine setupDiatoms_simpleOnly(n,errorio,errorstr)
   integer, intent(in):: n
   logical(1), intent(out):: errorio ! Whether to losses to the deep
   character(c_char), dimension(*) :: errorstr
   call parametersInit(1, n, 3,errorio,errorstr) ! 1 group, n size classes (excl nutrients)
     IF ( errorio ) RETURN 
   call parametersAddGroup(typeDiatom_simple, n, 1.d0,errorio,errorstr) ! diatoms with n size classes
   call parametersFinalize(0.1d0, .false., .false.)
  end subroutine setupDiatoms_simpleOnly
 
  ! -----------------------------------------------
  ! Generalists and diatoms:
  ! -----------------------------------------------
   subroutine setupGeneralistsDiatoms(n,errorio,errorstr)
      integer, intent(in):: n
      logical(1), intent(out):: errorio ! Whether to losses to the deep
      character(c_char), dimension(*) :: errorstr
      call parametersInit(2, 2*n, 3,errorio,errorstr)
        IF ( errorio ) RETURN 
      call parametersAddGroup(typeGeneralist, n, 0.0d0,errorio,errorstr) ! generalists with n size classes
        IF ( errorio ) RETURN 
      call parametersAddGroup(typeDiatom, n, 0.d0,errorio,errorstr) ! diatoms with n size classes
      call parametersFinalize(.1d0, .false., .false.)
   end subroutine setupGeneralistsDiatoms
 
   subroutine setupGeneralistsDiatoms_simple(n,errorio,errorstr)
      integer, intent(in):: n
      logical(1), intent(out):: errorio ! Whether to losses to the deep
      character(c_char), dimension(*) :: errorstr
      call parametersInit(2, 2*n, 3,errorio,errorstr)
        IF ( errorio ) RETURN 
      call parametersAddGroup(typeGeneralistSimple, n, 0.d0,errorio,errorstr) ! generalists with n size classes
        IF ( errorio ) RETURN 
      call parametersAddGroup(typeDiatom_simple, n, 1.d0,errorio,errorstr) ! diatoms with n size classes
      call parametersFinalize(0.1d0, .false., .false.)
   end subroutine setupGeneralistsDiatoms_simple
 
  ! -----------------------------------------------
  ! A basic setup with generalists and 1 copepod
  ! -----------------------------------------------
  subroutine setupGeneralistsSimpleCopepod(errorio,errorstr)
    logical(1), intent(out):: errorio ! Whether to losses to the deep
    character(c_char), dimension(*) :: errorstr
    call parametersInit(2, 20, 2,errorio,errorstr)
      IF ( errorio ) RETURN 
    call parametersAddGroup(typeGeneralistSimple, 10, 0.0d0,errorio,errorstr)
      IF ( errorio ) RETURN 
    call parametersAddGroup(typeCopepodActive, 10, .1d0,errorio,errorstr) ! add copepod with adult mass .1 mugC
    call parametersFinalize(0.003d0, .true., .true.) ! Use quadratic mortality
  end subroutine setupGeneralistsSimpleCopepod

  ! -----------------------------------------------
  ! A generic setup with generalists and a number of copepod species
  ! -----------------------------------------------
  subroutine setupGeneric(mAdult,errorio,errorstr)
    real(dp), intent(in):: mAdult(:)
    logical(1), intent(out):: errorio ! Whether to losses to the deep
    character(c_char), dimension(*) :: errorstr
    integer, parameter:: n = 10 ! number of size classes in each group
    integer:: iCopepod

    call parametersInit(size(mAdult)+1, n*(size(mAdult)+1), 2,errorio,errorstr)
      IF ( errorio ) RETURN 
    call parametersAddGroup(typeGeneralistSimple, n, 0.0d0,errorio,errorstr)
      IF ( errorio ) RETURN 
    if ( size(mAdult) .eq. 0) then
       call parametersFinalize(0.1d0, .true., .true.)
    else
       do iCopepod = 1, size(mAdult)
          call parametersAddGroup(typeCopepodActive, n, mAdult(iCopepod),errorio,errorstr) ! add copepod
            IF ( errorio ) RETURN 
       end do
       call parametersFinalize(0.001d0, .true., .true.)
    end if
  end subroutine setupGeneric
  ! -----------------------------------------------
  ! Full NUM model setup with generalists, copepods, and POM
  ! -----------------------------------------------
  subroutine setupNUMmodel(n, nCopepod, nPOM, mAdultPassive, mAdultActive,errorio,errorstr)
   integer, intent(in):: n, nCopepod, nPOM ! number of size classes in each group
   real(dp), intent(in):: mAdultPassive(:), mAdultActive(:)
   logical(1), intent(out):: errorio ! Whether to losses to the deep
   character(c_char), dimension(*) :: errorstr
   integer:: iCopepod
 
   call parametersInit(size(mAdultActive)+size(mAdultPassive)+3, 2*n + nPOM &
           + nCopepod*(size(mAdultPassive)+size(mAdultActive)), 3,errorio,errorstr)
     IF ( errorio ) RETURN 
   call parametersAddGroup(typeGeneralist, n, 0.0d0,errorio,errorstr)
     IF ( errorio ) RETURN 
   call parametersAddGroup(typeDiatom, n, 1.0d0,errorio,errorstr)
   IF ( errorio ) RETURN 

   do iCopepod = 1, size(mAdultPassive)
      call parametersAddGroup(typeCopepodPassive, nCopepod, mAdultPassive(iCopepod),errorio,errorstr) ! add copepod
        IF ( errorio ) RETURN 
   end do
   
   do iCopepod = 1, size(mAdultActive)
      call parametersAddGroup(typeCopepodActive, nCopepod, mAdultActive(iCopepod),errorio,errorstr) ! add copepod
        IF ( errorio ) RETURN 
   end do
   
   call parametersAddGroup(typePOM, nPOM, maxval(group(nGroups-1)%spec%mPOM),errorio,errorstr) ! POM with nPOM size classes and max size 1 ugC
   call parametersFinalize(0.07d0, .true., .false.)

  end subroutine setupNUMmodel

  ! -----------------------------------------------
  ! Full NUM model setup with generalistsSImple, copepods, and POM
  ! -----------------------------------------------
    subroutine setupNUMmodelSimple(n, nCopepod, nPOM, mAdult,errorio,errorstr)
   integer, intent(in):: n, nCopepod, nPOM ! number of size classes in each group
   real(dp), intent(in):: mAdult(:)
   logical(1), intent(out):: errorio ! Whether to losses to the deep
   character(c_char), dimension(*) :: errorstr
   integer:: iCopepod
 
   call parametersInit(size(mAdult)+3, 2*n + nPOM + nCopepod*size(mAdult), 3,errorio,errorstr)
     IF ( errorio ) RETURN 
   call parametersAddGroup(typeGeneralistSimple, n, 0.0d0,errorio,errorstr)
     IF ( errorio ) RETURN 
   call parametersAddGroup(typeDiatom_simple, n, 1.0d0,errorio,errorstr)
     IF ( errorio ) RETURN 

   do iCopepod = 1, size(mAdult)
      call parametersAddGroup(typeCopepodActive, nCopepod, mAdult(iCopepod),errorio,errorstr) ! add copepod
        IF ( errorio ) RETURN 
   end do
   call parametersAddGroup(typePOM, nPOM, maxval(group(nGroups-1)%spec%mPOM),errorio,errorstr) ! POM with nPOM size classes and max size 1 ugC
   call parametersFinalize(0.001d0, .true., .true.)

  end subroutine setupNUMmodelSimple
  
  ! -------------------------------------------------------
  ! A generic setup with generalists, diatoms and copepods
  ! -------------------------------------------------------
  subroutine setupGenDiatCope(n, nCopepod, nPOM, mAdult,errorio,errorstr)
   integer, intent(in):: n, nCopepod, nPOM ! number of size classes in each group
   real(dp), intent(in):: mAdult(:)
   logical(1), intent(out):: errorio ! Whether to losses to the deep
   character(c_char), dimension(*) :: errorstr
   integer:: iCopepod
 
   call parametersInit(size(mAdult)+3, 2*n + nPOM + nCopepod*size(mAdult), 3,errorio,errorstr)
     IF ( errorio ) RETURN 
   call parametersAddGroup(typeDiatom, n, 0.0d0,errorio,errorstr)
     IF ( errorio ) RETURN 
   call parametersAddGroup(typeGeneralist, n, 0.0d0,errorio,errorstr)
     IF ( errorio ) RETURN 

   do iCopepod = 1, size(mAdult)
      call parametersAddGroup(typeCopepodActive, nCopepod, mAdult(iCopepod),errorio,errorstr) ! add copepod
        IF ( errorio ) RETURN 
   end do
   call parametersAddGroup(typePOM, nPOM, maxval(group(nGroups-1)%spec%mPOM),errorio,errorstr) ! POM with nPOM size classes and max size 1 ugC
   call parametersFinalize(0.007d0, .true., .false.)

  end subroutine setupGenDiatCope

  ! ======================================
  !  Model initialization stuff:
  ! ======================================

  ! -----------------------------------------------
  ! Initialize parameters
  ! In:
  !    nnGroups: number of size spectrum groups
  !    nnGrid: total length of the grid (excl nnNutrients points for N, DOC, etc.)
  ! -----------------------------------------------
  subroutine parametersInit(nnGroups, nnGrid, nnNutrients,errorio,errorstr)
    integer, intent(in):: nnGrid, nnGroups, nnNutrients
    logical(1), intent(out):: errorio 
    character(c_char), dimension(*), intent(out) :: errorstr
    !
    ! read general input parameters:
    !
    errorio=.false.
    call read_input(inputfile,'general','rhoCN',rhoCN,errorio,errorstr)
    call read_input(inputfile,'general','fracHTL_to_N',fracHTL_to_N,errorio,errorstr)
    call read_input(inputfile,'general','fracHTL_to_POM',fracHTL_to_POM,errorio,errorstr)
    !call read_input(inputfile,'general')
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
       deallocate(F)
       deallocate(theta)
       deallocate(pHTL)
       if (allocated(thetaPOM)) then
         deallocate(thetaPOM)
         idxPOM = 0
       endif
    end if

    allocate( group(nGroups) )
    allocate(ixStart(nGroups))
    allocate(ixEnd(nGroups))
    allocate(upositive(nGrid))
    allocate(F(nGrid))
    allocate(pHTL(nGrid))
    allocate(theta(nGrid,nGrid))     ! Interaction matrix:
  end subroutine parametersInit
  
 

  ! -----------------------------------------------
  !  Add a size spectrum group
  !  In:
  !    typeGroup: the group type (see definitions in Globals.f90
  !    n: number of grid points
  !    mMax: the maximum size (upper size of a grid cell)
  ! -----------------------------------------------
  subroutine parametersAddGroup(typeGroup, n, mMax,errorio,errorstr)
    integer, intent(in):: typeGroup, n
    real(dp), intent(in):: mMax
    logical(1), intent(out) :: errorio
    character(c_char), dimension(*) :: errorstr

    type(spectrumGeneralists) :: specGeneralists
    type(spectrumGeneralistsSimple) :: specGeneralistsSimple
    type(spectrumDiatoms_simple):: specDiatoms_simple
    type(spectrumDiatoms):: specDiatoms
    type(spectrumCopepod):: specCopepod
    type(spectrumPOM):: specPOM
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
    case (typeGeneralistSimple)
      call initGeneralistsSimple(specGeneralistsSimple, n, errorio, errorstr)
      allocate( group( iCurrentGroup )%spec, source=specGeneralistsSimple )
    case (typeGeneralist)
      call initGeneralists(specGeneralists, n, errorio, errorstr)
      allocate( group( iCurrentGroup )%spec, source=specGeneralists )
    case (typeDiatom_simple)
      call initDiatoms_simple(specDiatoms_simple, n, mMax, errorio, errorstr)
      allocate( group( iCurrentGroup )%spec, source=specDiatoms_simple )
    case (typeDiatom)
      call initDiatoms(specDiatoms, n, errorio, errorstr)
      allocate( group( iCurrentGroup )%spec, source=specDiatoms )
   case(typeCopepodPassive)
      call initCopepod(specCopepod, passive, n, mMax, errorio, errorstr)
      allocate (group( iCurrentGroup )%spec, source=specCopepod)
   case(typeCopepodActive)
      call initCopepod(specCopepod, active, n, mMax, errorio, errorstr)
      allocate (group( iCurrentGroup )%spec, source=specCopepod)
   case(typePOM)
      call initPOM(specPOM, n, mMax, errorio, errorstr)
      allocate (group( iCurrentGroup )%spec, source=specPOM)
      idxPOM = iCurrentGroup 
   end select

  end subroutine parametersAddGroup
  ! -----------------------------------------------
  !  Finalize the setting of parameters. Must be called when
  !  all groups have been added.
  ! -----------------------------------------------
  subroutine parametersFinalize(mortHTL, boolQuadraticHTL, boolDecliningHTL)
    real(dp), intent(in):: mortHTL
    logical, intent(in):: boolQuadraticHTL, boolDecliningHTL
    integer:: i,j, iGroup, jGroup
    real(dp),parameter :: betaHTL = 500.d0
    real(dp):: mHTL
    !
    ! Calc theta:
    !
    theta = 0.d0
    do iGroup = 1, nGroups
      do i = 1, group(iGroup)%spec%n !group(iGroup)%ixStart, group(iGroup)%ixEnd
         do jGroup = 1, nGroups
            do j = 1, group(jGroup)%spec%n!group(jGroup)%ixStart, group(jGroup)%ixEnd
               theta(i+ixStart(iGroup)-1, j+ixStart(jGroup)-1) = &
                  group(jGroup)%spec%palatability * &
                  calcPhi(group(iGroup)%spec%m(i)/group(jGroup)%spec%m(j), &
                     group(iGroup)%spec%beta, group(iGroup)%spec%sigma, &
                     group(iGroup)%spec%z(i))
               !
               ! Passive copepods have a lower preference for feeding on diatoms:
               !
               select type (spec => group(iGroup)%spec) ! Predator group
                  type is (spectrumCopepod)
                    select type (specPrey => group(jGroup)%spec) ! Prey group
                      type is (spectrumDiatoms)
                          theta(i+ixStart(iGroup)-1, j+ixStart(jGroup)-1) = &
                            spec%DiatomsPreference * theta(i+ixStart(iGroup)-1, j+ixStart(jGroup)-1)
                      type is (spectrumDiatoms_simple)
                           theta(i+ixStart(iGroup)-1, j+ixStart(jGroup)-1) = &
                            spec%DiatomsPreference * theta(i+ixStart(iGroup)-1, j+ixStart(jGroup)-1)
                     end select
               end select
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

   call setHTL(mHTL, mortHTL, boolQuadraticHTL, boolDecliningHTL)
   !
   ! If there is a POM group then calculate the interactions with other groups
   !
   if ( idxPOM .ne. 0) then
      allocate( thetaPOM(nGrid) )
      thetaPOM = 0
      do iGroup = 1, nGroups
         if (iGroup .ne. idxPOM) then
           do i = 1, group(iGroup)%spec%n
              ! Find the size class in POM corresponding to mPOM:
              j = 1
              do while ( (group(iGroup)%spec%mPOM(i) .gt. &
                      (group(idxPOM)%spec%mLower(j)+group(idxPOM)%spec%mDelta(j))) &
                 .and. (j .lt. group(idxPOM)%spec%n))
                 j = j + 1
              end do
              thetaPOM( ixStart(iGroup)+i-1 ) = j
           end do
         end if
      end do
   end if

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
         res = 0.d0 ! beta = 0 is interpreted as if the group is not feeding
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
  ! Set the "HTL" mortality experienced by the largest groups.
  !
  ! IN:
  !  mHTL : The mass where the HTL mortality begins to act
  !  mortalityHTL : the level of the morality (see below)
  !  boolQuadraticHTL : whether to use a constant mortality (false) or a mortality
  !                     that is proportional to the biomass density (true)
  !  boolDecliningHTL : whether the mortality declines with size as mass^-1/4 (true)
  !                     or is constant (false)
  !
  ! If the mortality is constant (boolQuadraticHTL=false) then the mortality is at a level 
  ! of mortalityHTL at mHTL. The mortality may decline from there is boolDecliningHTL = true.
  !
  ! If the mortality is "quadratic" (boolQuadratic=true) then the mortality is
  ! at a level mortalityHTL at mHTL if the biomass density is 1 (ugC/l) / ugC.
  !
  ! The decline in mortality (boolDecliningHTL=true) is set such that the mortality
  ! is mortalityHTL at a mass mRef = .1 ugC.
  !
  subroutine setHTL(mHTL, mortalityHTL, boolQuadraticHTL, boolDecliningHTL)
    real(dp), intent(in):: mHTL ! The size where HTL is 50% of max
    real(dp), intent(in):: mortalityHTL ! The level of HTL mortality (at a reference size of 1 ugC
                                        ! B/z = 1/l )
    logical, intent(in):: boolQuadraticHTL ! Whether to use "quadratic" mortality
    logical, intent(in):: boolDecliningHTL ! Whether the mortality declines with size
    real(dp), parameter:: mRef = .1d0 ! Reference mass (in ugC)
    !real(dp), parameter:: betaHTL = 500.
    integer:: iGroup
 
    !     
    ! Calc htl mortality
    !

    pHTL = 0.d0
    ! Find the selectivity:
    do iGroup = 1, nGroups
      if (iGroup .ne. idxPOM) then ! POM is unaffected by HTL mortality
         pHTL( ixStart(iGroup):ixEnd(iGroup) ) = &
             (1 / (1+(group(iGroup)%spec%m/mHTL)**(-2))) ! The size selectivity switch around mHTL
         if (boolDecliningHTL) then
            pHTL( ixStart(iGroup):ixEnd(iGroup) ) = pHTL( ixStart(iGroup):ixEnd(iGroup) ) &
                * (group(iGroup)%spec%m/mHTL)**(-0.25)
         end if
      end if
    end do

    if (.not. boolQuadraticHTL) then
      !
      ! Standard HTL mortality that is constant over time:
      !
      do iGroup = 1, nGroups
         group(iGroup)%spec%mortHTL = mortalityHTL * pHTL( ixStart(iGroup):ixEnd(iGroup) )
      end do
    else
      !
      ! Linear HTL mortality (commonly referred to as "quadratic")
      ! The selectivity is now normalized by the width of the size classes
      !
      do iGroup = 1, nGroups
        pHTL( ixStart(iGroup):ixEnd(iGroup) ) = mortalityHTL &
           * pHTL( ixStart(iGroup):ixEnd(iGroup) ) &
           / log(1/group(iGroup)%spec%z)
      end do
    end if

    bQuadraticHTL = boolQuadraticHTL ! Set the global type of HTL mortality
  end subroutine setHTL
  !
  ! Routine for the user to override the calculation of HTL mortality done by setHTL:
  !
  subroutine setMortHTL(mortHTL)
   real(dp), intent(in):: mortHTL(nGrid-idxB+1)
   integer:: iGroup

   do iGroup = 1, nGroups
      group(iGroup)%spec%mortHTL = mortHTL( (ixStart(iGroup)-idxB+1):(ixEnd(iGroup)-idxB+1) )
   end do
  end subroutine setMortHTL

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
  !    dudt: vector to hold the derivative (which is returned)
  ! 
  !  Uses a simple predictor-corrector scheme.
  !  If one of the nutrients would become negative after an Euler 
  !  time step dt, then the uptake of said nutrient is reduced
  !  by a factor gamma to avoid the nutrient becoming negative.
  ! -----------------------------------------------
  subroutine calcDerivatives(u, L, T, dt, dudt)
    real(dp), intent(in):: L, T, dt, u(nGrid)
    real(dp), intent(inout) :: dudt(nGrid)
    integer:: i, j, iGroup, ix
    real(dp):: gammaN, gammaDOC, gammaSi

    dudt = 0.d0
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
    if (bQuadraticHTL) then
       do iGroup = 1, nGroups
         group(iGroup)%spec%mortHTL = pHTL(ixStart(iGroup):ixEnd(iGroup)) &
            * upositive(ixStart(iGroup):ixEnd(iGroup))
       end do
     end if
    !
    ! Calc derivatives of unicellular groups (predictor step)
    !
    gammaN = 1.d0
    gammaDOC = 1.d0
    gammaSi = 1.d0

    call calcDerivativesUnicellulars()
    !
    ! Make a correction if nutrient fields will become less than zero:
    !
    if ((u(idxN) + dudt(idxN)*dt) .lt. 0) then
       gammaN = max(0.d0, min(1.d0, -u(idxN)/(dudt(idxN)*dt)))
    end if
    if ((u(idxDOC) + dudt(idxDOC)*dt) .lt. 0) then
       gammaDOC = max(0.d0, min(1.d0, u(idxDOC)/(dudt(idxDOC)*dt)))
    end if
    if (nNutrients .gt. 2) then
      if ((u(idxSi) + dudt(idxSi)*dt) .lt. 0) then
        gammaSi = max(0.d0, min(1.d0, -u(idxSi)/(dudt(idxSi)*dt)))
      end if
    end if
    if ((gammaN .lt. 1.d0) .or. (gammaDOC .lt. 1.d0) .or. (gammaSi .lt. 1.d0)) then
       call calcDerivativesUnicellulars()
    end if
    !
    ! Calc derivatives of multicellular groups:
    !
    do iGroup = 1, nGroups
      select type (spec => group(iGroup)%spec)
      type is (spectrumCopepod)
         call calcDerivativesCopepod(spec, &
            upositive(ixStart(iGroup):ixEnd(iGroup)), &
            dudt(idxN), &
            dudt(ixStart(iGroup):ixEnd(iGroup)))
      end select
    end do
    !
    ! Transfer POM to POM group:
    !
    if (idxPOM .ne. 0) then
      do iGroup = 1, nGroups ! run over all groups
         if (iGroup .ne. idxPOM) then
            do i = 1, group(iGroup)%spec%n ! run over all size classes
               ix = ixStart(iGroup)+i-1
               if (thetaPOM(ix) .ne. 0) then
                  j = ixStart(idxPOM)+thetaPOM(ix)-1 ! find the size class that it delivers POM to
                  dudt(j) = dudt(j) + group(iGroup)%spec%jPOM(i)*upositive(ix) 
               end if
               ! Throw a fraction of HTL production into the largest POM group:
               dudt(ixEnd(idxPOM)) = dudt(ixEnd(idxPOM)) + &
                 (1-fracHTL_to_N) * upositive(ix) * group(iGroup)%spec%mortHTL(i)
            end do
         end if
      end do
      !
      ! Update POM
      ! 
      select type (spec => group(idxPOM)%spec)
      type is (spectrumPOM)
         call calcDerivativesPOM(spec, &
           upositive(ixStart(idxPOM):ixEnd(idxPOM)), &
           dudt(idxN), dudt(idxDOC), dudt(ixStart(idxPOM):ixEnd(idxPOM)))
      end select
    end if
    !
    ! Some HTL mortality ends up as nutrients:
    !
    do iGroup = 1, nGroups
      if (iGroup .ne. idxPOM) then
        dudt(idxN) = dudt(idxN) + &
            fracHTL_to_N * sum( upositive(ixStart(iGroup):ixEnd(iGroup)) * group(iGroup)%spec%mortHTL )/rhoCN
      end if
    end do
    !
    ! Check: Should be close to zero
    !
    !Nbalance = dudt(idxN) + sum(dudt(idxB:nGrid))/rhoCN
    !do iGroup = 1, nGroups
    !  Nbalance = Nbalance  &
    !    + (1.d0-fracHTL_to_N) * sum( u(ixStart(iGroup):ixEnd(iGroup)) * group(iGroup)%spec%mortHTL )/rhoCN &
    !    + sum(u(ixStart(iGroup):ixEnd(iGroup)) * group(iGroup)%spec%jPOM) / rhoCN
    !end do
    !write(*,*) 'N balance:', Nbalance

    contains

  !
  ! Calculate derivatives for unicellular groups
  ! In:
  !   gammaN and gammaDOC are reduction factors [0...1] of uptakes of N and DOC,
  !   used for correction of Euler integration. If no correction is used, just set to 1.0
  !   This correction procedure is needed for correct Euler integration.
  subroutine calcDerivativesUnicellulars()
   integer :: jGroup, ixj, ixi
   !
   ! Calc uptakes of all unicellular groups:
   !
   do iGroup = 1, nGroups
      select type (spec => group(iGroup)%spec)
      type is (spectrumGeneralistsSimple)
         call calcRatesGeneralistsSimple(spec, &
                     L, upositive(idxN), upositive(idxDOC), gammaN, gammaDOC)
      type is (spectrumGeneralists)
         call calcRatesGeneralists(spec, &
                     L, upositive(idxN), upositive(idxDOC), gammaN, gammaDOC)
      type is (spectrumDiatoms_simple)
         call calcRatesDiatoms_simple(spec, &
                     L, upositive(idxN), upositive(idxSi), gammaN, gammaSi)
      type is (spectrumDiatoms)
         call calcRatesDiatoms(spec, &
                     L, upositive(idxN),  upositive(idxDOC), upositive(idxSi),gammaN, gammaDOC, gammaSi)
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
      select type (spec => group(iGroup)%spec)
      type is (spectrumGeneralistsSimple)
         call calcDerivativesGeneralistsSimple(spec, &
              upositive(ixStart(iGroup):ixEnd(iGroup)), &
              dudt(idxN), dudt(idxDOC), dudt(ixStart(iGroup):ixEnd(iGroup)))
      type is (spectrumGeneralists)
         call calcDerivativesGeneralists(spec, &
              upositive(ixStart(iGroup):ixEnd(iGroup)), &
              dudt(idxN), dudt(idxDOC), dudt(ixStart(iGroup):ixEnd(iGroup)))
      type is (spectrumDiatoms_simple)
         call calcDerivativesDiatoms_simple(spec, &
              upositive(ixStart(iGroup):ixEnd(iGroup)), &
              dudt(idxN), dudt(idxSi), dudt(ixStart(iGroup):ixEnd(iGroup)))
      type is (spectrumDiatoms)
         call calcDerivativesDiatoms(spec, &
              upositive(ixStart(iGroup):ixEnd(iGroup)), &
              dudt(idxN), dudt(idxDOC), dudt(idxSi), dudt(ixStart(iGroup):ixEnd(iGroup)))
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
  subroutine simulateChemostatEuler(u, L, T, Ndeep, diff, tEnd, dt, bLosses)
    real(dp), intent(inout):: u(:) ! Initial conditions and result after integration
    real(dp), intent(in):: L      ! Light level
    real(dp), intent(in):: T ! Temperature
    real(dp), intent(in):: Ndeep(nNutrients) ! Nutrients in the deep layer
    real(dp), intent(in):: diff      ! Diffusivity
    real(dp), intent(in):: tEnd ! Time to simulate
    real(dp), intent(in):: dt    ! time step
    logical(1), intent(in):: bLosses ! Whether to losses to the deep
    real(dp) :: dudt(nGrid)
    integer:: i, iEnd, iGroup

    iEnd = floor(tEnd/dt)
   
    do i=1, iEnd
       call calcDerivatives(u, L, T, dt, dudt)
       dudt(idxN) = dudt(idxN) + diff*(Ndeep(idxN)-u(idxN))
       dudt(idxDOC) = dudt(idxDOC) + diff*(Ndeep(idxDOC) - u(idxDOC))
       if (idxB .gt. idxSi) then
         dudt(idxSi) = dudt(idxSi) + diff*(Ndeep(idxSi) - u(idxSi))
       end if  
       !
       ! Mixing:
       !
       if (bLosses) then
         do iGroup = 1, nGroups
            if ( (group(iGroup)%spec%type .ne. typeCopepodActive) .and. (group(iGroup)%spec%type .ne. typeCopepodPassive) ) then
               dudt( ixStart(iGroup):ixEnd(iGroup) ) = dudt( ixStart(iGroup):ixEnd(iGroup) ) + diff*(0.d0 - u(idxB:nGrid))
            end if
         end do
       end if
       !
       ! Sinking:
       !
       do iGroup = 1, nGroups
         dudt( ixStart(iGroup):ixEnd(iGroup) ) = dudt( ixStart(iGroup):ixEnd(iGroup) ) - &
            group(iGroup)%spec%velocity*u( ixStart(iGroup):ixEnd(iGroup) )
       end do

       u = u + dudt*dt ! Euler update
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
    real(dp) :: dudt(nGrid)
    integer:: i

    do i = 1, floor(tEnd/dt)
       call calcDerivatives(u, L, T, dt, dudt)
       u = u + dudt*dt
    end do
  end subroutine simulateEuler

  ! -----------------------------------------------
  ! Simulate with Euler integration and also return functions
  ! -----------------------------------------------
  subroutine simulateEulerFunctions(u, L, T, tEnd, dt, &
            ProdGross, ProdNet,ProdHTL,prodBact,eHTL,Bpico,Bnano,Bmicro)
    real(dp), intent(inout):: u(:) ! Initial conditions and result after integration
    real(dp), intent(in):: L      ! Light level
    real(dp), intent(in):: T ! Temperature
    real(dp), intent(in):: tEnd ! Time to simulate
    real(dp), intent(in):: dt    ! time step
    real(dp), intent(out):: ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro

    call simulateEuler(u, L, T, tEnd, dt)
    call getFunctions(u, ProdGross, ProdNet,ProdHTL,prodBact,eHTL,Bpico,Bnano,Bmicro)
 end subroutine simulateEulerFunctions

  !=========================================
  ! Diagnostic functions
  !=========================================

   subroutine printRates()
      integer :: iGroup

      do iGroup = 1, nGroups
         call group(iGroup)%spec%printRates()
      end do
   end subroutine printRates

  function calcN(u) result(N)
   real(dp), intent(in):: u(:)
   integer:: i
   real(dp):: N

   N = 0
   N = u(idxN)
   do i = 1, nGrid
      N = N + u(nNutrients+i)/rhoCN
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

   subroutine getSinking(velocity)
      real(dp), intent(inout):: velocity(nGrid)
      integer:: iGroup

      velocity(1:(idxB-1)) = 0.d0 ! No sinking of nutrient groups
      do iGroup = 1,nGroups
         velocity( ixStart(iGroup):ixEnd(iGroup) ) = group(iGroup)%spec%velocity
      end do
     
   end subroutine getSinking
   !
   ! Set the sinking velocity of all groups. Typically used to change the velocity of POM
   !
   subroutine setSinking(velocity)
      real(dp), intent(in) :: velocity(nGrid)
      integer:: iGroup
  
      do iGroup = 1, nGroups
         group(iGroup)%spec%velocity = velocity( (ixStart(iGroup)):(ixEnd(iGroup)) )
      end do
    end subroutine setSinking
  
  ! ---------------------------------------------------
  ! Get the ecosystem functions as calculated from the last call
  ! to calcDerivatives
  ! ---------------------------------------------------
  subroutine getFunctions(u, ProdGross, ProdNet,ProdHTL,prodBact,eHTL,Bpico,Bnano,Bmicro)
    real(dp), intent(in):: u(nGrid)
    real(dp), intent(out):: ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro
    real(dp) :: conversion
    real(dp) :: ESD(nGrid)
    real(dp):: m(nGrid), mDelta(nGrid)
    integer:: i

    ProdGross = 0.d0
    ProdNet = 0.d0
    ProdHTL = 0.d0
    ProdBact = 0.d0
    Bpico = 0.d0
    Bnano = 0.d0
    Bmicro = 0.d0
    !
    ! Get primary production only from unicellular spectra:
    !
    !conversion = 365.*1d-6*1000. ! Convert to gC/yr/m3
    conversion = 1d-3*1000. ! Convert to mgC/d/m3
    do i = 1, nGroups
       select type (spec => group(i)%spec)
          class is (spectrumUnicellular)
            ProdGross = ProdGross + conversion * &
               sum(  spec%JLreal * u( ixStart(i):ixEnd(i) ) / spec%m )
          
            ProdNet = ProdNet + conversion * &
               spec%getProdNet(u( ixStart(i):ixEnd(i) ))

            ProdBact = ProdBact + conversion * &
               spec%getProdBact(u( ixStart(i):ixEnd(i) ))
       end select
    end do
    !
    ! Make a rough estimate of pico-nano-micro plankton biomasses:
    !
    call getMass(m, mDelta )
    ESD = 10000. * 1.5 * (m*1d-6)**onethird
    !conversion = 1d-6*1000 ! Convert to gC/m3
    conversion = 1d-3*1000 ! Convert to mgC/m3

    do i = idxB, nGrid
       if (ESD(i) .le. 2.) then
          Bpico = Bpico + conversion*u(i)
       endif
       
       if ((ESD(i).gt.2.) .and. (ESD(i) .le. 20.)) then
          Bnano = Bnano + conversion*u(i)
       endif

       if (ESD(i) .gt. 20.) then
          Bmicro = Bmicro + conversion*u(i)
       endif
    end do
    !
    ! HTL production:
    !  
    do i = 1, nGroups
       !ProdHTL = ProdHTL + 365*conversion* &
        ProdHTL = ProdHTL + conversion* &
          sum(group(i)%spec%mortHTL*u( ixStart(i):ixEnd(i) ))
    end do

    eHTL = ProdHTL / ProdNet
 
  end subroutine getFunctions

  ! ---------------------------------------------------
  ! Returns mass conservation calculated from last call to calcDerivatives.
  ! The numbers returned are the changes in carbon, N, and Si relative to the 
  ! concentrations of N, DOC, and Si, in units of 1/day
  ! 
  ! ISSUE:
  ! Does not deal with remin2 losses
  !
  ! ---------------------------------------------------
subroutine getBalance(u, dudt, Cbalance,Nbalance,SiBalance)
   real(dp), intent(in):: u(nGrid), dudt(nGrid)
   real(dp), intent(out):: Cbalance, Nbalance, Sibalance
   real(dp) :: Clost, Nlost, SiLost
   integer :: iGroup

   call getLost(u, Clost, Nlost, SiLost) ! Find what is lost from the system
   !
   ! The balance is the change in the nutrient + the change in biomass + whatever is lost:
   !
   Cbalance = dudt(idxDOC) + sum( dudt(idxB:nGrid) )       + Clost
   Nbalance = dudt(idxN)   + sum( dudt(idxB:nGrid) )/rhoCN + Nlost 
   ! Only diatom groups for silicate
   Sibalance = dudt(idxSi) + SiLost 

   do iGroup = 1, nGroups
      select type ( spec => group(iGroup)%spec )
         type is (spectrumDiatoms)
            Sibalance = Sibalance + sum( dudt( ixStart(iGroup):ixEnd(iGroup) ) )/3.4d0 ! rhoC:Si hardcoded here
            
         type is (spectrumDiatoms_simple)
            Sibalance = Sibalance + sum( dudt( ixStart(iGroup):ixEnd(iGroup) ) )/3.4d0 ! rhoC:Si hardcoded here
      end select
   end do
   !
   ! Normalize by the concentrations:
   !
   Cbalance = Cbalance / u(idxDOC)
   Nbalance = Nbalance / u(idxN)
   if (nNutrients .gt. 2) then
     Sibalance = Sibalance / u(idxSi)
   else
     Sibalance = 0.d0
   endif

end subroutine getBalance
!
! Return what is lost (or gained) of carbon, nutrients and silicate from the system:
!
subroutine getLost(u, Clost, Nlost, SiLost)
   real(dp), intent(in):: u(nGrid)
   real(dp), intent(out):: Clost, Nlost, SiLost
   integer:: iGroup
   real(dp) :: HTLloss, POMloss
   real(dp),parameter :: rhoCSi = 3.4d0 ! Hard coded here

   Clost = 0.d0
   Nlost = 0.d0
   SiLost = 0.d0

   do iGroup = 1, nGroups
      !
      ! Get carbon losses from each group:
      !
      Clost = Clost + group(iGroup)%spec%getClost( &
         u(ixStart(iGroup):ixEnd(iGroup)))
      !
      ! Deal with higher trophic level losses and POM:
      !
      HTLloss = sum( group(iGroup)%spec%mortHTL * u(ixStart(iGroup):ixEnd(iGroup) ))
      POMloss = sum( group(iGroup)%spec%jPOM * u(ixStart(iGroup):ixEnd(iGroup)) )

      if (idxPOM .eq. 0) then
         !
         ! If there is no POM, then POM is lost from the system:
         !
         Nlost = Nlost + (1-fracHTL_to_N) * HTLloss/rhoCN ! The part of POM which is not respired
         Clost = Clost +  HTLloss  ! Both respiration losses and POM losses
      
         Nlost = Nlost + POMloss / rhoCN
         Clost = Clost + POMloss   
      else
         !
         ! If there is POM then just account for the respiration
         !
         Clost = Clost +  fracHTL_to_N * HTLloss  ! Respiration losses from HTL
      end if
      !
      ! Deal with silicate balance. Clunky and does not deal with remin2-losses
      !
      select type ( spec => group(iGroup)%spec )
         type is (spectrumDiatoms)
            ! Consumed silicate is considered lost:
            SiLost = SiLost + ( &
               sum( spec%mortpred*u(ixStart(iGroup):ixEnd(iGroup)) ) & ! The silicate from consumed diatoms is lost
               + HTLloss & ! All HTL silicate is lost. rhoCSi is hard-coded here
               + POMloss) / rhoCSi ! POM from diatoms is lost:
            
         type is (spectrumDiatoms_simple)
           ! NOT IMPLEMENTED
            SiLost = SiLost + ( &
               sum( spec%mortpred*u(ixStart(iGroup):ixEnd(iGroup)) ) & ! The silicate from consumed diatoms is lost
               + HTLloss & ! All HTL silicate is lost. rhoCSi is hard-coded here
               + POMloss) / rhoCSi ! POM from diatoms is lost:
      
      end select
   end do
end subroutine getLost

  ! ---------------------------------------------------
  ! Returns the rates calculated from last call to calcDerivatives
  ! ---------------------------------------------------
  subroutine getRates(jN, jDOC, jL, jSi, jF, jFreal, ff,&
    jTot, jMax, jFmax, jR, jResptot, jLossPassive, &
    jNloss,jLreal, jPOM, &
    mortpred, mortHTL, mort2, mort)
    use globals
    real(dp), intent(out):: jN(nGrid-nNutrients), jDOC(nGrid-nNutrients), jL(nGrid-nNutrients)
    real(dp), intent(out):: jSi(nGrid-nNutrients)
    real(dp), intent(out):: jF(nGrid-nNutrients), jFreal(nGrid-nNutrients), ff(nGrid-nNutrients)
    real(dp), intent(out):: jTot(nGrid-nNutrients), jMax(nGrid-nNutrients), jFmax(nGrid-nNutrients)
    real(dp), intent(out):: jR(nGrid-nNutrients), jResptot(nGrid-nNutrients)
    real(dp), intent(out):: jLossPassive(nGrid-nNutrients), jNloss(nGrid-nNutrients), jLreal(nGrid-nNutrients)
    real(dp), intent(out):: jPOM(nGrid-nNutrients)
    real(dp), intent(out):: mortpred(nGrid-nNutrients), mortHTL(nGrid-nNutrients)
    real(dp), intent(out):: mort2(nGrid-nNutrients), mort(nGrid-nNutrients)
   integer :: iGroup, i1, i2

   do iGroup = 1, nGroups
      i1 = ixStart(iGroup)-nNutrients
      i2 = ixEnd(iGroup)-nNutrients
      ! Extract common fields:
      jF( i1:i2 ) = group(iGroup)%spec%flvl * group(iGroup)%spec%JFmax / group(iGroup)%spec%m
      jFreal( i1:i2 ) = group(iGroup)%spec%JF / group(iGroup)%spec%m
      jFmax( i1:i2 ) = fTemp2 * group(iGroup)%spec%JFmax / group(iGroup)%spec%m
      Jtot( i1:i2 ) = group(iGroup)%spec%Jtot / group(iGroup)%spec%m
      mortpred( i1:i2 ) = group(iGroup)%spec%mortpred
      mortHTL( i1:i2 ) = group(iGroup)%spec%mortHTL
      mort2( i1:i2 ) = group(iGroup)%spec%mort2
      jNloss( i1:i2 ) = group(iGroup)%spec%JNloss / group(iGroup)%spec%m
      jR( i1:i2 ) = fTemp2 * group(iGroup)%spec%Jresp / group(iGroup)%spec%m
      jResptot( i1:i2 ) = group(iGroup)%spec%JResptot / group(iGroup)%spec%m
      jPOM( i1:i2 ) = group(iGroup)%spec%jPOM

      select type (spectrum => group(iGroup)%spec)
      class is (spectrumUnicellular)
        jN( i1:i2 ) = fTemp15 * spectrum%JN / spectrum%m
        jDOC( i1:i2 ) = fTemp15 * spectrum%JDOC / spectrum%m
        jL( i1:i2 ) = spectrum%JL / spectrum%m
        jMax( i1:i2 ) = fTemp2 * spectrum%Jmax / spectrum%m
        jLossPassive( i1:i2 ) = spectrum%JlossPassive / spectrum%m
        jLreal( i1:i2 ) = spectrum%JLreal / spectrum%m
        ff( i1:i2 ) = group(iGroup)%spec%f
      class is (spectrumMulticellular)
        jN( i1:i2 ) = 0.d0
        jDOC( i1:i2 ) = 0.d0
        jL( i1:i2 ) = 0.d0
        jMax( i1:i2 ) = 0.d0
        jLossPassive( i1:i2 ) = 0.d0
        jLreal( i1:i2 ) = 0.d0
        ff( i1:i2 ) = 0.d0
      end select

      select type (spectrum => group(iGroup)%spec)
      class is (spectrumDiatoms_simple)
        jSi( i1:i2 ) = spectrum%JSi / spectrum%m
        jDOC( i1:i2 ) = 0.d0 ! Diatoms simple don't take up DOC
      class is (spectrumDiatoms)
        jSi( i1:i2 ) = spectrum%JSi / spectrum%m
      end select

      mort = 0.d0 ! For odd reasons this gives a segfault when called from R

   end do
  end subroutine getRates

end module NUMmodel
