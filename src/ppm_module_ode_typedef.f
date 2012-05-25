module ppm_module_ode_typedef

!----------------------------------------------------------------------
!  Modules
!----------------------------------------------------------------------
use ppm_module_core
use ppm_module_numerics_interfaces

implicit none
!-----------------------------------------------------------------------
! ODE states
!-----------------------------------------------------------------------
integer, parameter :: ode_state_init      = 0
integer, parameter :: ode_state_kickoff   = 1 ! currently not used
integer, parameter :: ode_state_running   = 2
integer, parameter :: ode_state_finished  = 3

!----------------------------------------------------------------------
!  Types
!----------------------------------------------------------------------
type,extends(ppm_t_ode_) :: ppm_t_ode
  contains
  procedure :: create => ode_create
  procedure :: destroy => ode_destroy
  procedure :: step => ode_step
end type ppm_t_ode
        

!----------------------------------------------------------------------
!  Type-bound procedures
!----------------------------------------------------------------------
contains

subroutine ode_create(this,scheme,fields,rhsfunc,rhs_fields_discr,info,kickoff_scheme)
  use ppm_module_integrator_typedef
  implicit none

  class(ppm_t_ode)                       :: this
  integer,          intent(in   )        :: scheme
  class(ppm_v_main_abstr), pointer       :: fields
  procedure(ppm_p_rhsfunc)               :: rhsfunc
  class(ppm_v_field_discr_pair), pointer :: rhs_fields_discr
  integer,          intent(  out)        :: info
  integer,optional, intent(in   )        :: kickoff_scheme
  !----------------------------------------------------------------------
  !  Variables
  !----------------------------------------------------------------------
  integer                          :: kickoff

  start_subroutine("ode_add_integrator")

  select case(scheme)
  case(ppm_param_ode_scheme_eulerf)
    ! allocate changes array
    allocate(ppm_t_eulerf::this%integrator,STAT=info)
    call this%integrator%create(fields,rhsfunc,rhs_fields_discr,info)
  case default
    ppm_fail("Integrator not implemented")
  end select
  
  ! use default or user kickoff scheme
  if(PRESENT(kickoff_scheme)) then
    kickoff = kickoff_scheme
  else
    kickoff = this%integrator%scheme_kickoff
  end if
  if (kickoff.eq.scheme) then
    nullify(this%kickoff)
    this%state = ode_state_init
  else
    select case(kickoff)
    case(ppm_param_ode_scheme_eulerf)
      ! allocate changes array
      allocate(ppm_t_eulerf::this%kickoff,STAT=info)
      call this%kickoff%create(fields,rhsfunc,rhs_fields_discr,info)
     or_fail("Creating eulerf failed")
    case default
      ppm_fail("Integrator not implemented")
    end select
    this%state = ode_state_kickoff
  end if
  

  end_subroutine()
end subroutine ode_create

subroutine ode_destroy(this,info)
  IMPLICIT NONE
  class(ppm_t_ode)      :: this
  integer,               intent(  out)   :: info
  start_subroutine("ode_destroy")
  
  end_subroutine()
end subroutine ode_destroy

subroutine ode_step(this,t,dt,info)
  implicit none
  class(ppm_t_ode)                                   :: this
  real(ppm_kind_double),            intent(inout)    :: t
  real(ppm_kind_double),            intent(in   )    :: dt
  integer,                          intent(  out)    :: info
  start_subroutine("ode_step")



  if (this%state.EQ.ode_state_init) then
    call this%integrator%step(t,dt,info)
    this%state = ode_state_running
  else if (this%state.EQ.ode_state_kickoff) then
    call this%kickoff%step(t,dt,info)
    this%state = ode_state_running
  else if(this%state.EQ.ode_state_running) then
    call this%integrator%step(t,dt,info)
  end if
  
  end_subroutine()

end subroutine ode_step

end module ppm_module_ode_typedef
