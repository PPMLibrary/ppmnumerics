!minclude ppm_header(ppm_module_interfaces)

#define __REAL 3 
#define __COMPLEX 4 
#define __integer 5 
#define __LONGINT 6 
#define __LOGICAL 7 
#define __CHAR 8 

module ppm_module_numerics_interfaces
!!! Declares all data types 
!!! The derived types are declared as abstract. They contain all the
!!! different fields as well as the interfaces to the type-bound
!!! procedures.
!!! Each type is then extended in another module, which also
!!! contains the implementations of the type-bound procedures.
!!! This allows for type-bound procedures to take arguments of any
!!! other derived type (because they only need to use the module
!!! with the abstract types and there is no risk of circular
!!! dependencies between modules).

!----------------------------------------------------------------------
!  Modules
!----------------------------------------------------------------------
use ppm_module_core


implicit none

!----------------------------------------------------------------------
! Global parameters
!----------------------------------------------------------------------

include 'ppm_numerics.h'

!----------------------------------------------------------------------
! Global variables 
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Module variables 
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Type declaration
!----------------------------------------------------------------------
  

type,abstract :: ppm_t_integrator_
    !!! Data structure for time integrators
    !!! The time integrator implements a specific integration scheme.

    integer                                         :: ID = 0
    !!! global identifier 
    character(len=ppm_char)                         :: name
    !!! string description
    integer                     :: scheme_order
    !!! order of integration scheme
    integer                     :: scheme_memsize
    !!! memory requirements for intermediate buffers
    integer                     :: scheme_nstages
    !!! number of stages for multistage time integration schemes
    integer                     :: scheme_kickoff
    !!! suitable kick off scheme
    
    procedure(ppm_p_rhsfunc_s),       pointer, nopass :: rhsfunc_s        => null()
    procedure(ppm_p_rhsfunc_d),       pointer, nopass :: rhsfunc_d        => null()

    class(ppm_v_main_abstr),          pointer         :: variables        => null()
    class(ppm_v_discr_kind),          pointer         :: discretizations  => null()
    class(ppm_v_var_discr_pair),      pointer         :: rhs_variables    => null()
    class(ppm_v_main_abstr),          pointer         :: changes          => null()
    class(ppm_v_field),               pointer         :: buffers          => null() 
    contains
    procedure(integrator_create_s_),     deferred :: create_s
    procedure(integrator_create_d_),     deferred :: create_d
    generic :: create => create_s, create_d
    procedure(integrator_destroy_),      deferred :: destroy
    procedure(integrator_step_s_),       deferred :: step_s
    procedure(integrator_step_d_),       deferred :: step_d
    generic :: step => step_s, step_d
end type ppm_t_integrator_
!minclude ppm_create_collection(integrator,integrator_,generate="abstract")
!minclude ppm_create_collection(integrator,integrator_,generate="abstract",vec=true,def_ptr=false)

type,abstract :: ppm_t_ode_
    !!! Data structure for ODEs
    !!! An ODE represents the time integration scheme and provides the
    !!! functionality for performing a time step using one of the provided
    !!! schemes.

    integer                                         :: ID = 0
    !!! global identifier 
    
    integer                     :: state
    !!! state of time ODE
    !!!
    !!! One of:
    !!! * ode_state_init 
    !!! * ode_state_kickoff
    !!! * ode_state_running
    !!! * ode_state_finished
    !!! 
    !!! WARNING: This is internal.
    
    class(ppm_t_integrator_), pointer       :: integrator => null()
    class(ppm_t_integrator_), pointer       :: kickoff => null()

    contains
    procedure(ode_create_s_),      deferred :: create_s
    procedure(ode_create_d_),      deferred :: create_d
    generic :: create => create_s, create_d
    procedure(ode_destroy_),       deferred :: destroy
    procedure(ode_step_),          deferred :: step
    procedure(ode_map_push_),      deferred :: map_push
    procedure(ode_map_pop_),       deferred :: map_pop
end type ppm_t_ode_

!----------------------------------------------------------------------
!  INTERFACES
!----------------------------------------------------------------------
abstract interface


  real(ppm_kind_single) function ppm_p_rhsfunc_s(vars_and_discr,time,changes)
  import ppm_v_main_abstr
  import ppm_kind_single
  import ppm_v_var_discr_pair
  class(ppm_v_var_discr_pair), pointer      :: vars_and_discr
  real(ppm_kind_single)                     :: time
  class(ppm_v_main_abstr),    pointer       :: changes
  end function ppm_p_rhsfunc_s
  
  real(ppm_kind_double) function ppm_p_rhsfunc_d(vars_and_discr,time,changes)
  import ppm_v_main_abstr
  import ppm_kind_double
  import ppm_v_var_discr_pair
  class(ppm_v_var_discr_pair), pointer      :: vars_and_discr
  real(ppm_kind_double)                     :: time
  class(ppm_v_main_abstr),    pointer       :: changes
  end function ppm_p_rhsfunc_d
end interface


interface

#include "ode/integrator_interfaces.f"

#include "ode/ode_interfaces.f"

end interface


end module ppm_module_numerics_interfaces
