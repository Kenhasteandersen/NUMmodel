#include "fabm_driver.h"
!
! FABM coupling for the NUMmodel size-structured plankton model.
! Based on setupNUMmodel: generalists + diatoms + passive/active copepods + POM.
!
! Units:
!   Dissolved nitrogen / silicate / carbon : µgN/L, µgSi/L, µgC/L
!   Plankton biomass                        : µgC/L
!   Rates from NUMmodel (1/day)  --> /s for FABM  (divide by 86400)
!   Sinking velocities  (m/day)  --> m/s for FABM (divide by 86400, negate)
!   Light : W/m²   Temperature : °C
!
! LIMITATION: NUMmodel uses global workspace arrays; only one FABM instance
! per executable is supported.  Not OpenMP-safe over spatial points.
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
  real(dp), parameter :: dt_nominal   = 1._dp / 24._dp  ! 1-hour nominal dt (days)

  type, extends(type_base_model) :: type_num_model
    ! Configuration
    integer :: n_size, n_copepod, n_pom, n_passive, n_active

    ! Nutrient state variable IDs
    type(type_state_variable_id) :: id_N, id_DOC, id_Si

    ! Biomass state variable IDs (one per size class)
    type(type_state_variable_id), allocatable :: id_B(:)

    ! Environmental dependencies
    type(type_dependency_id) :: id_T, id_PAR

    ! Ecosystem function diagnostics
    type(type_diagnostic_variable_id) :: id_ProdGross, id_ProdNet, id_ProdHTL
    type(type_diagnostic_variable_id) :: id_Bpico, id_Bnano, id_Bmicro

  contains
    procedure :: initialize
    procedure :: do
    procedure :: get_vertical_movement
  end type type_num_model

contains

  subroutine initialize(self, configunit)
    class(type_num_model), intent(inout), target :: self
    integer,               intent(in)            :: configunit

    integer  :: i, iGroup, ig, n_biomass
    real(rk) :: mAdult_tmp
    real(dp), allocatable :: mAdultPassive(:), mAdultActive(:), velocity(:)
    logical(1)                             :: errorio
    character(kind=c_char), dimension(256) :: errorstr
    character(len=32) :: varname, longname
    character(len=8)  :: prefix

    call self%get_parameter(self%n_size,    'n_size',    '-', &
        'size classes per generalist/diatom group', default=10)
    call self%get_parameter(self%n_copepod, 'n_copepod', '-', &
        'size classes per copepod group',           default=10)
    call self%get_parameter(self%n_pom,     'n_pom',     '-', &
        'POM size classes',                         default=1)
    call self%get_parameter(self%n_passive, 'n_passive', '-', &
        'number of passive-feeding copepod groups', default=2)
    call self%get_parameter(self%n_active,  'n_active',  '-', &
        'number of active-feeding copepod groups',  default=2)

    allocate(mAdultPassive(self%n_passive))
    allocate(mAdultActive(self%n_active))

    do i = 1, self%n_passive
      write(varname, '(a,i0)') 'mAdultPassive', i
      call self%get_parameter(mAdult_tmp, trim(varname), 'ug C', &
          'adult mass of passive copepod group '//trim(varname), &
          default=10._rk**(real(-1 + 2*i, rk)))
      mAdultPassive(i) = real(mAdult_tmp, dp)
    end do

    do i = 1, self%n_active
      write(varname, '(a,i0)') 'mAdultActive', i
      call self%get_parameter(mAdult_tmp, trim(varname), 'ug C', &
          'adult mass of active copepod group '//trim(varname), &
          default=10._rk**(real(-1 + 2*i, rk)))
      mAdultActive(i) = real(mAdult_tmp, dp)
    end do

    call setupNUMmodel(self%n_size, self%n_copepod, self%n_pom, &
                       mAdultPassive, mAdultActive, errorio, errorstr)
    deallocate(mAdultPassive)
    deallocate(mAdultActive)

    if (errorio) then
      call self%fatal_error('fabm_NUMmodel', 'setupNUMmodel returned an error')
      return
    end if

    ! Nutrients
    call self%register_state_variable(self%id_N, &
        'N', 'ug N L-1', 'dissolved inorganic nitrogen', &
        minimum=0._rk, initial_value=14._rk)
    call self%register_state_variable(self%id_DOC, &
        'DOC', 'ug C L-1', 'dissolved organic carbon', &
        minimum=0._rk, initial_value=0._rk)
    call self%register_state_variable(self%id_Si, &
        'Si', 'ug Si L-1', 'dissolved silicate', &
        minimum=0._rk, initial_value=150._rk)

    ! Sinking velocities (m/day -> m/s, positive=up so negate)
    n_biomass = nGrid - nNutrients
    allocate(velocity(nGrid))
    call getSinking(velocity)

    ! Biomass state variables
    allocate(self%id_B(n_biomass))
    ig = 0
    do iGroup = 1, nGroups
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
        call self%register_state_variable(self%id_B(ig), &
            trim(varname), 'ug C L-1', trim(longname), &
            minimum=0._rk, initial_value=1.e-4_rk, &
            vertical_movement=-real(velocity(idxB + ig - 1), rk) / secs_per_day)
      end do
    end do
    deallocate(velocity)

    ! Environmental dependencies
    call self%register_dependency(self%id_T,   standard_variables%temperature)
    call self%register_dependency(self%id_PAR, &
        standard_variables%downwelling_photosynthetic_radiative_flux)

    ! Diagnostics
    call self%register_diagnostic_variable(self%id_ProdGross, &
        'ProdGross', 'mg C d-1 m-3', 'gross primary production')
    call self%register_diagnostic_variable(self%id_ProdNet, &
        'ProdNet', 'mg C d-1 m-3', 'net primary production')
    call self%register_diagnostic_variable(self%id_ProdHTL, &
        'ProdHTL', 'mg C d-1 m-3', 'production removed by higher trophic levels')
    call self%register_diagnostic_variable(self%id_Bpico, &
        'Bpico', 'mg C m-3', 'pico-plankton biomass (ESD < 2 um)')
    call self%register_diagnostic_variable(self%id_Bnano, &
        'Bnano', 'mg C m-3', 'nano-plankton biomass (2-20 um ESD)')
    call self%register_diagnostic_variable(self%id_Bmicro, &
        'Bmicro', 'mg C m-3', 'micro-plankton biomass (ESD > 20 um)')

  end subroutine initialize

  ! Compute source/sink rates for all state variables.
  ! self is intent(in) to match the base-class interface.
  ! _ADD_SOURCE_ macros must not span line continuations.
  subroutine do(self, _ARGUMENTS_DO_)
    class(type_num_model), intent(in) :: self
    _DECLARE_ARGUMENTS_DO_

    real(rk) :: T_rk, PAR_rk, B_tmp, rate_rk
    real(dp) :: u(nGrid), dudt(nGrid)
    real(dp) :: ProdGross, ProdNet, ProdHTL, ProdBact, eHTL
    real(dp) :: Bpico, Bnano, Bmicro, mHTL
    integer  :: i

    _LOOP_BEGIN_

      _GET_(self%id_T,   T_rk)
      _GET_(self%id_PAR, PAR_rk)

      _GET_(self%id_N,   B_tmp)
      u(idxN)   = real(B_tmp, dp)
      _GET_(self%id_DOC, B_tmp)
      u(idxDOC) = real(B_tmp, dp)
      _GET_(self%id_Si,  B_tmp)
      u(idxSi)  = real(B_tmp, dp)

      do i = 1, nGrid - nNutrients
        _GET_(self%id_B(i), B_tmp)
        u(idxB + i - 1) = real(B_tmp, dp)
      end do

      call calcDerivatives(u, real(PAR_rk, dp), real(T_rk, dp), dt_nominal, dudt)

      rate_rk = real(dudt(idxN), rk) / secs_per_day
      _ADD_SOURCE_(self%id_N, rate_rk)
      rate_rk = real(dudt(idxDOC), rk) / secs_per_day
      _ADD_SOURCE_(self%id_DOC, rate_rk)
      rate_rk = real(dudt(idxSi), rk) / secs_per_day
      _ADD_SOURCE_(self%id_Si, rate_rk)

      do i = 1, nGrid - nNutrients
        rate_rk = real(dudt(idxB + i - 1), rk) / secs_per_day
        _ADD_SOURCE_(self%id_B(i), rate_rk)
      end do

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

  ! Sinking velocities are registered as constants in initialize;
  ! this no-op satisfies the FABM interface.
  subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
    class(type_num_model), intent(in) :: self
    _DECLARE_ARGUMENTS_DO_
  end subroutine get_vertical_movement

end module fabm_num_model
