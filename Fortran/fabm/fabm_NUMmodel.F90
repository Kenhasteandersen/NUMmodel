#include "fabm_driver.h"
!
! FABM coupling for the NUMmodel size-structured plankton model.
! Based on setupNUMmodel: generalists + diatoms + passive/active copepods + POM.
!
! Units used by this module:
!   Dissolved nitrogen    : µgN/L
!   Dissolved silicate    : µgSi/L
!   Dissolved org. carbon : µgC/L
!   Plankton biomass      : µgC/L
!   Rates from NUMmodel   : 1/day  --> converted to 1/s for FABM
!   Sinking velocities    : m/day  --> converted to m/s for FABM (negative = downward)
!   Light                 : W/m²  (downwelling PAR from FABM standard variables)
!   Temperature           : °C
!
! IMPORTANT LIMITATION:
!   NUMmodel uses module-level (global) workspace arrays (F, upositive, etc.).
!   Only a single FABM instance of this model is supported per executable.
!   Thread-safety (OpenMP over spatial points) is not guaranteed.
!
module fabm_num_model
  use fabm_types
  use iso_c_binding, only: c_char
  use NUMmodel, only: &
      setupNUMmodel, calcDerivatives, getSinking, getFunctions, &
      nGrid, nNutrients, nGroups, ixStart, ixEnd, &
      idxN, idxDOC, idxSi, idxB, idxPOM, group
  use globals, only: dp

  implicit none

  private
  public :: type_num_model

  real(rk), parameter :: secs_per_day = 86400._rk
  ! Nominal timestep (days) passed to calcDerivatives for its nutrient-limiter
  ! predictor-corrector.  One hour is a reasonable ocean model timestep.
  real(dp), parameter :: dt_nominal = 1._dp / 24._dp

  type, extends(type_base_model) :: type_num_model

    ! ---------------------------------------------------------------
    ! Configuration parameters (read from fabm.yaml)
    ! ---------------------------------------------------------------
    integer :: n_size    ! size classes per generalist/diatom group
    integer :: n_copepod ! size classes per copepod group
    integer :: n_pom     ! POM size classes
    integer :: n_passive ! number of passive copepod groups
    integer :: n_active  ! number of active copepod groups

    ! ---------------------------------------------------------------
    ! State variable IDs
    ! ---------------------------------------------------------------
    type(type_state_variable_id) :: id_N    ! dissolved nitrogen  [µgN/L]
    type(type_state_variable_id) :: id_DOC  ! dissolved org. C    [µgC/L]
    type(type_state_variable_id) :: id_Si   ! dissolved silicate  [µgSi/L]

    ! One entry per biomass size class (flat index 1..nGrid-nNutrients)
    type(type_state_variable_id), allocatable :: id_B(:)

    ! ---------------------------------------------------------------
    ! Environmental dependencies
    ! ---------------------------------------------------------------
    type(type_dependency_id) :: id_T    ! temperature [°C]
    type(type_dependency_id) :: id_PAR  ! downwelling PAR [W/m²]

    ! ---------------------------------------------------------------
    ! Diagnostic variable IDs
    ! ---------------------------------------------------------------
    type(type_diagnostic_variable_id) :: id_ProdGross
    type(type_diagnostic_variable_id) :: id_ProdNet
    type(type_diagnostic_variable_id) :: id_ProdHTL
    type(type_diagnostic_variable_id) :: id_Bpico
    type(type_diagnostic_variable_id) :: id_Bnano
    type(type_diagnostic_variable_id) :: id_Bmicro

  contains
    procedure :: initialize
    procedure :: do
    procedure :: get_vertical_movement
  end type type_num_model

contains

  ! ================================================================
  subroutine initialize(self, configunit)
  ! ================================================================
    class(type_num_model), intent(inout), target :: self
    integer,               intent(in)            :: configunit

    integer  :: i, iGroup, ig, n_biomass
    real(rk) :: mAdult_tmp
    real(dp), allocatable :: mAdultPassive(:), mAdultActive(:)
    real(dp), allocatable :: velocity(:)
    logical(1)                             :: errorio
    character(kind=c_char), dimension(256) :: errorstr
    character(len=32) :: varname, longname
    character(len=8)  :: prefix

    ! ------------------------------------------------------------------
    ! Read parameters from fabm.yaml
    ! ------------------------------------------------------------------
    call self%get_parameter(self%n_size,    'n_size',    '-', &
        'size classes per generalist/diatom group',    default=10)
    call self%get_parameter(self%n_copepod, 'n_copepod', '-', &
        'size classes per copepod group',              default=10)
    call self%get_parameter(self%n_pom,     'n_pom',     '-', &
        'POM size classes',                            default=1)
    call self%get_parameter(self%n_passive, 'n_passive', '-', &
        'number of passive-feeding copepod groups',    default=2)
    call self%get_parameter(self%n_active,  'n_active',  '-', &
        'number of active-feeding copepod groups',     default=2)

    allocate(mAdultPassive(self%n_passive))
    allocate(mAdultActive(self%n_active))

    ! Passive copepod adult masses (default: 10, 1000 µgC)
    do i = 1, self%n_passive
      write(varname, '(a,i0)') 'mAdultPassive', i
      call self%get_parameter(mAdult_tmp, trim(varname), 'ug C', &
          'adult mass of passive copepod group '//trim(varname), &
          default=10._rk**(real(-1 + 2*i, rk)))
      mAdultPassive(i) = real(mAdult_tmp, dp)
    end do

    ! Active copepod adult masses (default: 10, 1000 µgC)
    do i = 1, self%n_active
      write(varname, '(a,i0)') 'mAdultActive', i
      call self%get_parameter(mAdult_tmp, trim(varname), 'ug C', &
          'adult mass of active copepod group '//trim(varname), &
          default=10._rk**(real(-1 + 2*i, rk)))
      mAdultActive(i) = real(mAdult_tmp, dp)
    end do

    ! ------------------------------------------------------------------
    ! Initialise the NUMmodel
    ! ------------------------------------------------------------------
    call setupNUMmodel(self%n_size, self%n_copepod, self%n_pom, &
                       mAdultPassive, mAdultActive, errorio, errorstr)
    deallocate(mAdultPassive)
    deallocate(mAdultActive)

    if (errorio) then
      call self%fatal_error('fabm_NUMmodel', 'setupNUMmodel returned an error')
      return
    end if

    ! ------------------------------------------------------------------
    ! Register nutrient state variables
    ! ------------------------------------------------------------------
    self%id_N = self%register_state_variable( &
        'N', 'ug N L-1', 'dissolved inorganic nitrogen', &
        minimum=0._rk, initial_value=14._rk)

    self%id_DOC = self%register_state_variable( &
        'DOC', 'ug C L-1', 'dissolved organic carbon', &
        minimum=0._rk, initial_value=0._rk)

    self%id_Si = self%register_state_variable( &
        'Si', 'ug Si L-1', 'dissolved silicate', &
        minimum=0._rk, initial_value=150._rk)

    ! ------------------------------------------------------------------
    ! Obtain size-class sinking velocities (m/day) to register with
    ! each biomass state variable as a constant vertical velocity.
    ! ------------------------------------------------------------------
    n_biomass = nGrid - nNutrients
    allocate(velocity(nGrid))
    call getSinking(velocity)

    ! ------------------------------------------------------------------
    ! Register biomass state variables (one per size class in every group)
    ! ------------------------------------------------------------------
    allocate(self%id_B(n_biomass))
    ig = 0
    do iGroup = 1, nGroups
      ! Choose a short prefix that identifies the group type
      select case (group(iGroup)%spec%type)
      case (5)     ; prefix = 'Gen'
      case (1)     ; prefix = 'GenS'
      case (3)     ; prefix = 'Dia'
      case (4)     ; prefix = 'DiaS'
      case (10)    ; prefix = 'ACop'
      case (11)    ; prefix = 'PCop'
      case (100)   ; prefix = 'POM'
      case default ; prefix = 'B'
      end select

      do i = 1, group(iGroup)%spec%n
        ig = ig + 1
        write(varname,  '(a,i0)') trim(prefix), i
        write(longname, '(a,a,i0)') trim(prefix), ' biomass size class ', i
        ! velocity(idxB+ig-1) is in m/day, positive = sinking downward.
        ! FABM vertical_movement is in m/s, positive = upward; hence the negation.
        self%id_B(ig) = self%register_state_variable( &
            trim(varname), 'ug C L-1', trim(longname), &
            minimum=0._rk, initial_value=1.e-4_rk, &
            vertical_movement=-real(velocity(idxB + ig - 1), rk) / secs_per_day)
      end do
    end do

    deallocate(velocity)

    ! ------------------------------------------------------------------
    ! Register environmental dependencies
    ! ------------------------------------------------------------------
    call self%register_dependency(self%id_T, &
        standard_variables%temperature)
    call self%register_dependency(self%id_PAR, &
        standard_variables%downwelling_photosynthetic_radiative_flux)

    ! ------------------------------------------------------------------
    ! Register diagnostic variables
    ! ------------------------------------------------------------------
    self%id_ProdGross = self%register_diagnostic_variable( &
        'ProdGross', 'mg C d-1 m-3', 'gross primary production')
    self%id_ProdNet   = self%register_diagnostic_variable( &
        'ProdNet',   'mg C d-1 m-3', 'net primary production')
    self%id_ProdHTL   = self%register_diagnostic_variable( &
        'ProdHTL',   'mg C d-1 m-3', &
        'production removed by higher trophic levels')
    self%id_Bpico     = self%register_diagnostic_variable( &
        'Bpico', 'mg C m-3', 'pico-plankton biomass (ESD < 2 µm)')
    self%id_Bnano     = self%register_diagnostic_variable( &
        'Bnano', 'mg C m-3', 'nano-plankton biomass (2-20 µm ESD)')
    self%id_Bmicro    = self%register_diagnostic_variable( &
        'Bmicro', 'mg C m-3', 'micro-plankton biomass (ESD > 20 µm)')

  end subroutine initialize

  ! ================================================================
  ! do: compute local source/sink rates for all state variables.
  ! Called once per spatial point per time step by FABM.
  ! All rates returned to FABM are in units of [state variable unit] / s.
  ! ================================================================
  subroutine do(self)
    class(type_num_model), intent(inout) :: self

    real(rk) :: T_rk, PAR_rk, B_tmp
    real(dp) :: u(nGrid), dudt(nGrid)
    real(dp) :: ProdGross, ProdNet, ProdHTL, ProdBact, eHTL
    real(dp) :: Bpico, Bnano, Bmicro, mHTL
    integer  :: i

    _LOOP_BEGIN_

      ! Retrieve forcing
      _GET_(self%id_T,   T_rk)
      _GET_(self%id_PAR, PAR_rk)

      ! Retrieve nutrient state variables
      _GET_(self%id_N,   B_tmp)
      u(idxN)   = real(B_tmp, dp)
      _GET_(self%id_DOC, B_tmp)
      u(idxDOC) = real(B_tmp, dp)
      _GET_(self%id_Si,  B_tmp)
      u(idxSi)  = real(B_tmp, dp)

      ! Retrieve biomass state variables
      do i = 1, nGrid - nNutrients
        _GET_(self%id_B(i), B_tmp)
        u(idxB + i - 1) = real(B_tmp, dp)
      end do

      ! Compute rates of change (dudt in 1/day)
      call calcDerivatives(u, real(PAR_rk, dp), real(T_rk, dp), &
                           dt_nominal, dudt)

      ! Provide nutrient rates to FABM (convert 1/day -> 1/s)
      _ADD_SOURCE_(self%id_N,   real(dudt(idxN),   rk) / secs_per_day)
      _ADD_SOURCE_(self%id_DOC, real(dudt(idxDOC), rk) / secs_per_day)
      _ADD_SOURCE_(self%id_Si,  real(dudt(idxSi),  rk) / secs_per_day)

      ! Provide biomass rates to FABM (convert 1/day -> 1/s)
      do i = 1, nGrid - nNutrients
        _ADD_SOURCE_(self%id_B(i), &
            real(dudt(idxB + i - 1), rk) / secs_per_day)
      end do

      ! Compute and set ecosystem-function diagnostics
      call getFunctions(u, ProdGross, ProdNet, ProdHTL, ProdBact, &
                        eHTL, Bpico, Bnano, Bmicro, mHTL)
      _SET_DIAGNOSTIC_(self%id_ProdGross, real(ProdGross, rk))
      _SET_DIAGNOSTIC_(self%id_ProdNet,   real(ProdNet,   rk))
      _SET_DIAGNOSTIC_(self%id_ProdHTL,   real(ProdHTL,   rk))
      _SET_DIAGNOSTIC_(self%id_Bpico,     real(Bpico,     rk))
      _SET_DIAGNOSTIC_(self%id_Bnano,     real(Bnano,     rk))
      _SET_DIAGNOSTIC_(self%id_Bmicro,    real(Bmicro,    rk))

    _LOOP_END_

  end subroutine do

  ! ================================================================
  ! get_vertical_movement: sinking velocities are registered as
  ! constants in initialize (vertical_movement parameter), so this
  ! procedure is a no-op placeholder.
  ! ================================================================
  subroutine get_vertical_movement(self)
    class(type_num_model), intent(inout) :: self
  end subroutine get_vertical_movement

end module fabm_num_model
