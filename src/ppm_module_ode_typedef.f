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
  procedure :: map_push => ode_map_push
  procedure :: map_pop => ode_map_pop
end type ppm_t_ode
        

!----------------------------------------------------------------------
!  Type-bound procedures
!----------------------------------------------------------------------
contains

subroutine ode_create(this,scheme,fields,rhsfunc,rhs_fields_discr,info,options,kickoff_scheme)
  use ppm_module_integrator_typedef
  implicit none

  class(ppm_t_ode)                       :: this
  integer,          intent(in   )        :: scheme
  class(ppm_v_main_abstr)                :: fields
  procedure(ppm_p_rhsfunc)               :: rhsfunc
  class(ppm_v_field_discr_pair)          :: rhs_fields_discr
  integer,          intent(  out)        :: info
  class(ppm_t_options),target,optional,intent(in   ) :: options
  integer,optional, intent(in   )        :: kickoff_scheme
  !----------------------------------------------------------------------
  !  Variables
  !----------------------------------------------------------------------
  integer                          :: kickoff

  start_subroutine("ode_create")

  select case(scheme)
  case(ppm_param_ode_scheme_eulerf)
    ! allocate changes array
    allocate(ppm_t_eulerf::this%integrator,STAT=info)
    call this%integrator%create(fields,rhsfunc,rhs_fields_discr,info,options)
    or_fail("Creating eulerf failed")
  case(ppm_param_ode_scheme_sts)
    ! allocate changes array
    allocate(ppm_t_sts::this%integrator,STAT=info)
    call this%integrator%create(fields,rhsfunc,rhs_fields_discr,info,options)
    or_fail("Creating STS failed")
  case(ppm_param_ode_scheme_tvdrk2)
    ! allocate changes array
    allocate(ppm_t_tvdrk2::this%integrator,STAT=info)
    call this%integrator%create(fields,rhsfunc,rhs_fields_discr,info,options)
    or_fail("Creating tvd RK2 failed")
  case(ppm_param_ode_scheme_midrk2)
    ! allocate changes array
    allocate(ppm_t_midrk2::this%integrator,STAT=info)
    call this%integrator%create(fields,rhsfunc,rhs_fields_discr,info,options)
    or_fail("Creating mid RK2 failed")
  case(ppm_param_ode_scheme_rk4)
    ! allocate changes array
    allocate(ppm_t_rk4::this%integrator,STAT=info)
    call this%integrator%create(fields,rhsfunc,rhs_fields_discr,info,options)
    or_fail("Creating RK4 failed")
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
      call this%kickoff%create(fields,rhsfunc,rhs_fields_discr,info,options)
     or_fail("Creating eulerf failed")
    case(ppm_param_ode_scheme_sts)
      ! allocate changes array
      allocate(ppm_t_sts::this%kickoff,STAT=info)
      call this%kickoff%create(fields,rhsfunc,rhs_fields_discr,info,options)
     or_fail("Creating STS failed")
    case(ppm_param_ode_scheme_tvdrk2)
      ! allocate changes array
      allocate(ppm_t_tvdrk2::this%kickoff,STAT=info)
      call this%kickoff%create(fields,rhsfunc,rhs_fields_discr,info,options)
      or_fail("Creating tvd RK2 failed")
    case(ppm_param_ode_scheme_midrk2)
      ! allocate changes array
      allocate(ppm_t_midrk2::this%kickoff,STAT=info)
      call this%kickoff%create(fields,rhsfunc,rhs_fields_discr,info,options)
      or_fail("Creating mid RK2 failed")
    case(ppm_param_ode_scheme_rk4)
      ! allocate changes array
      allocate(ppm_t_rk4::this%kickoff,STAT=info)
      call this%kickoff%create(fields,rhsfunc,rhs_fields_discr,info,options)
      or_fail("Creating RK4 failed")
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

subroutine ode_step(this,t,dt,istage,info)
  implicit none
  class(ppm_t_ode)                                   :: this
  real(ppm_kind_double),            intent(inout)    :: t
  real(ppm_kind_double),            intent(in   )    :: dt
  integer,                          intent(in   )    :: istage
  integer,                          intent(  out)    :: info
  start_subroutine("ode_step")



  if (this%state.EQ.ode_state_init) then
    call this%integrator%step(t,dt,istage,info)
    this%state = ode_state_running
  else if (this%state.EQ.ode_state_kickoff) then
    call this%kickoff%step(t,dt,istage,info)
    this%state = ode_state_running
  else if(this%state.EQ.ode_state_running) then
    call this%integrator%step(t,dt,istage,info)
  end if
  
  end_subroutine()

end subroutine ode_step

subroutine ode_map_push(this,info)
  IMPLICIT NONE
  class(ppm_t_ode)      :: this
  integer,               intent(  out)   :: info
  class(ppm_t_field_), pointer           :: buffer => null()
  class(ppm_t_discr_info_), pointer      :: di => null()
  class(ppm_t_particles_d), pointer      :: pset => null()
  start_subroutine("ode_map_push")
  
  if (this%state.EQ.ode_state_kickoff) then
    buffer => this%kickoff%buffers%begin()
    do while (associated(buffer))
      di => buffer%discr_info%begin()
      select type(disc => di%discr_ptr)
        class is (ppm_t_particles_d)
        pset => disc
        call pset%map_ghost_push(info,buffer)
      end select

      buffer => this%kickoff%buffers%next()
    end do
  else
    buffer => this%integrator%buffers%begin()
    do while (associated(buffer))
      di => buffer%discr_info%begin()
      select type(disc => di%discr_ptr)
        class is (ppm_t_particles_d)
        pset => disc
        call pset%map_ghost_push(info,buffer)
      end select

      buffer => this%integrator%buffers%next()
    end do
  end if

  end_subroutine()
end subroutine ode_map_push

subroutine ode_map_pop(this,info)
  IMPLICIT NONE
  class(ppm_t_ode)      :: this
  integer,               intent(  out)   :: info
  class(ppm_t_field_), pointer           :: buffer => null()
  class(ppm_t_discr_info_), pointer      :: di => null()
  class(ppm_t_particles_d), pointer      :: pset => null()
  start_subroutine("ode_map_pop")
    
  if (this%state.EQ.ode_state_kickoff) then
    ! traverse in reverse order!
    buffer => this%kickoff%buffers%last()
    do while (associated(buffer))
      di => buffer%discr_info%begin()
      select type(disc => di%discr_ptr)
        class is (ppm_t_particles_d)
        pset => disc
        call pset%map_ghost_pop(info,buffer)
      end select

      buffer => this%kickoff%buffers%prev()
    end do
  else
    ! traverse in reverse order!
    buffer => this%integrator%buffers%last()
    do while (associated(buffer))
      di => buffer%discr_info%begin()
      select type(disc => di%discr_ptr)
        class is (ppm_t_particles_d)
        pset => disc
        call pset%map_ghost_pop(info,buffer)
      end select

      buffer => this%integrator%buffers%prev()
    end do
  end if

  end_subroutine()
end subroutine ode_map_pop

end module ppm_module_ode_typedef
