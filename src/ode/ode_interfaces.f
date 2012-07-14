
subroutine ode_create_s_(this,scheme,variables,rhsfunc,rhs_variables,info,options,kickoff_scheme)
  import ppm_v_main_abstr, ppm_v_var_discr_pair, ppm_t_options
  import ppm_t_ode_, ppm_p_rhsfunc_s
  class(ppm_t_ode_)                      :: this
  integer,          intent(in   )        :: scheme
  class(ppm_v_main_abstr)                :: variables
  procedure(ppm_p_rhsfunc_s)             :: rhsfunc
  class(ppm_v_var_discr_pair)            :: rhs_variables
  integer,          intent(  out)        :: info
  class(ppm_t_options),target,optional,intent(in   ) :: options
  integer,optional, intent(in   )        :: kickoff_scheme
end subroutine

subroutine ode_create_d_(this,scheme,variables,rhsfunc,rhs_variables,info,options,kickoff_scheme)
  import ppm_v_main_abstr, ppm_v_var_discr_pair, ppm_t_options
  import ppm_t_ode_, ppm_p_rhsfunc_d
  class(ppm_t_ode_)                      :: this
  integer,          intent(in   )        :: scheme
  class(ppm_v_main_abstr)                :: variables
  procedure(ppm_p_rhsfunc_d)             :: rhsfunc
  class(ppm_v_var_discr_pair)            :: rhs_variables
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
  integer,          intent(  out)        :: info
end subroutine

subroutine ode_map_pop_(this,info)
  import ppm_t_ode_
  class(ppm_t_ode_)                      :: this
  integer,          intent(  out)        :: info
end subroutine

