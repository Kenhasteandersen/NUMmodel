#include "fabm_driver.h"
!
! FABM model factory for NUMmodel.
!
! Registers the NUMmodel coupling under the identifier "num/num_model".
! In fabm.yaml, reference the model as:
!
!   instances:
!     plankton:
!       model: num/num_model
!       parameters:
!         ...
!
module num_model_library
  use fabm_types,   only: type_base_model_factory, type_base_model
  use fabm_num_model

  implicit none

  private
  public :: type_factory

  type, extends(type_base_model_factory) :: type_factory
  contains
    procedure :: create
  end type type_factory

  type(type_factory), target, public :: num_model_factory

contains

  subroutine create(self, name, model)
    class(type_factory),     intent(in)  :: self
    character(*),            intent(in)  :: name
    class(type_base_model),  pointer     :: model

    select case (name)
    case ('num_model')
      allocate(type_num_model :: model)
    case default
      call self%type_base_model_factory%create(name, model)
    end select
  end subroutine create

end module num_model_library
