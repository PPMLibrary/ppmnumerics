
subroutine ode_create_(this,scheme,fields,rhsfunc,rhs_fields_discr,info,options,kickoff_scheme)
  import ppm_v_main_abstr, ppm_v_field_discr_pair, ppm_t_options
  import ppm_t_ode_, ppm_p_rhsfunc
  class(ppm_t_ode_)                      :: this
  integer,          intent(in   )        :: scheme
  class(ppm_v_main_abstr), pointer       :: fields
  procedure(ppm_p_rhsfunc)               :: rhsfunc
  class(ppm_v_field_discr_pair), pointer :: rhs_fields_discr
  integer,          intent(  out)        :: info
  class(ppm_t_options),target,optional,intent(in   ) :: options
  integer,optional, intent(in   )        :: kickoff_scheme
end subroutine

subroutine ode_destroy_(this,info)
  import ppm_t_ode_
  class(ppm_t_ode_)                      :: this
  integer,          intent(  out)        :: info
end subroutine

subroutine ode_step_(this,t,dt,istage,info)
  import ppm_t_ode_, ppm_kind_double
  class(ppm_t_ode_)                      :: this
  real(ppm_kind_double), intent(inout)   :: t
  real(ppm_kind_double), intent(in   )   :: dt
  integer,               intent(in   )   :: istage
  integer,               intent(  out)   :: info
end subroutine

subroutine ode_map_push_(this,info)
  import ppm_t_ode_
  class(ppm_t_ode_)                      :: this
  integer,               intent(  out)   :: info
end subroutine

subroutine ode_map_pop_(this,info)
  import ppm_t_ode_
  class(ppm_t_ode_)                      :: this
  integer,               intent(  out)   :: info
end subroutine
