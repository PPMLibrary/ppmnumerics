!minclude ppm_header(ppm_module_interfaces)

      MODULE ppm_module_numerics_interfaces
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
      USE ppm_module_core
      USE ppm_module_interfaces, ONLY : ppm_v_main_abstr,ppm_v_var_discr_pair,  &
      &   ppm_v_discr_kind,ppm_t_options,ppm_t_part_prop_s_,ppm_t_part_prop_d_, &
      &   ppm_t_field_discr_pair,ppm_t_var_discr_pair
      USE ppm_module_field_typedef, ONLY : ppm_v_field

      USE ppm_module_numerics_data
      IMPLICIT NONE

      !----------------------------------------------------------------------
      ! Global parameters
      !----------------------------------------------------------------------
      !----------------------------------------------------------------------
      ! Global variables
      !----------------------------------------------------------------------
      !----------------------------------------------------------------------
      ! Module variables
      !----------------------------------------------------------------------
      !----------------------------------------------------------------------
      ! Type declaration
      !----------------------------------------------------------------------
      TYPE,ABSTRACT :: ppm_t_integrator_
          !!! Data structure for time integrators
          !!! The time integrator implements a specific integration scheme.
          INTEGER                                      :: ID = 0
          !!! global identifier
          CHARACTER(LEN=ppm_char)                      :: name
          !!! string description
          INTEGER                                      :: scheme_order
          !!! order of integration scheme
          INTEGER                                      :: scheme_memsize
          !!! memory requirements for intermediate buffers
          INTEGER                                      :: scheme_nstages
          !!! number of stages for multistage time integration schemes
          INTEGER                                      :: scheme_kickoff
          !!! suitable kick off scheme

          PROCEDURE(ppm_p_rhsfunc_s),  POINTER, NOPASS :: rhsfunc_s       => NULL()
          PROCEDURE(ppm_p_rhsfunc_d),  POINTER, NOPASS :: rhsfunc_d       => NULL()

          CLASS(ppm_v_main_abstr),     POINTER         :: variables       => NULL()
          CLASS(ppm_v_discr_kind),     POINTER         :: discretizations => NULL()
          CLASS(ppm_v_var_discr_pair), POINTER         :: rhs_variables   => NULL()
          CLASS(ppm_v_main_abstr),     POINTER         :: changes         => NULL()
          CLASS(ppm_v_field),          POINTER         :: buffers         => NULL()

      CONTAINS
          PROCEDURE(integrator_create_s_), DEFERRED :: create_s
          PROCEDURE(integrator_create_d_), DEFERRED :: create_d
          GENERIC :: create => create_s, create_d

          PROCEDURE(integrator_destroy_), DEFERRED :: destroy

          PROCEDURE(integrator_step_s_), DEFERRED :: step_s
          PROCEDURE(integrator_step_d_), DEFERRED :: step_d
          GENERIC :: step => step_s, step_d

      END TYPE ppm_t_integrator_
      !minclude ppm_create_collection(integrator,integrator_,generate="abstract")
      !minclude ppm_create_collection(integrator,integrator_,generate="abstract",vec=true,def_ptr=false)

      TYPE,ABSTRACT :: ppm_t_ode_
          !!! Data structure for ODEs
          !!! An ODE represents the time integration scheme and provides the
          !!! functionality for performing a time step using one of the provided
          !!! schemes.
          INTEGER                           :: ID = 0
          !!! global identifier
          INTEGER                           :: state
          !!! state of time ODE
          !!!
          !!! One of:
          !!! * ode_state_init
          !!! * ode_state_kickoff
          !!! * ode_state_running
          !!! * ode_state_finished
          !!!
          !!! WARNING: This is internal.
          CLASS(ppm_t_integrator_), POINTER :: integrator => NULL()
          CLASS(ppm_t_integrator_), POINTER :: kickoff    => NULL()

      CONTAINS
          PROCEDURE(ode_create_s_), DEFERRED :: create_s
          PROCEDURE(ode_create_d_), DEFERRED :: create_d
          GENERIC :: create => create_s, create_d

          PROCEDURE(ode_destroy_),  DEFERRED :: destroy

          PROCEDURE(ode_step_),     DEFERRED :: step

          PROCEDURE(ode_map_push_), DEFERRED :: map_push
          PROCEDURE(ode_map_pop_),  DEFERRED :: map_pop

      END TYPE ppm_t_ode_

      !----------------------------------------------------------------------
      !  INTERFACES
      !----------------------------------------------------------------------
      ABSTRACT INTERFACE
        REAL(ppm_kind_single) FUNCTION ppm_p_rhsfunc_s(vars_and_discr,time,changes)
          IMPORT :: ppm_v_var_discr_pair
          IMPORT :: ppm_kind_single
          IMPORT :: ppm_v_main_abstr
          IMPLICIT NONE
          CLASS(ppm_v_var_discr_pair), POINTER :: vars_and_discr
          REAL(ppm_kind_single)                :: time
          CLASS(ppm_v_main_abstr),     POINTER :: changes
        END FUNCTION ppm_p_rhsfunc_s

        REAL(ppm_kind_double) FUNCTION ppm_p_rhsfunc_d(vars_and_discr,time,changes)
          IMPORT :: ppm_v_var_discr_pair
          IMPORT :: ppm_kind_double
          IMPORT :: ppm_v_main_abstr
          IMPLICIT NONE
          CLASS(ppm_v_var_discr_pair), POINTER :: vars_and_discr
          REAL(ppm_kind_double)                :: time
          CLASS(ppm_v_main_abstr),     POINTER :: changes
        END FUNCTION ppm_p_rhsfunc_d
      END INTERFACE

      INTERFACE
#include "integrator/integrator_interfaces.f"

#include "ode/ode_interfaces.f"
      END INTERFACE

      END MODULE ppm_module_numerics_interfaces
