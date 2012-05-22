
subroutine ode_create_(this,scheme,fields,rhsfunc,rhs_fields,info,kickoff_scheme)
  import ppm_t_ode_,ppm_v_main_abstr,ppm_p_rhsfunc
  class(ppm_t_ode_)                      :: this
  integer,          intent(in   )        :: scheme
  class(ppm_v_main_abstr), pointer       :: fields
  procedure(ppm_p_rhsfunc)               :: rhsfunc
  class(ppm_v_main_abstr), pointer       :: rhs_fields
  integer,          intent(  out)        :: info
  integer,optional, intent(in   )        :: kickoff_scheme
end subroutine

subroutine ode_destroy_(this,info)
  import ppm_t_ode_
  class(ppm_t_ode_)                      :: this
  integer,          intent(  out)        :: info
end subroutine

subroutine ode_step_(this,t,dt,info)
  import ppm_t_ode_, ppm_kind_double
  class(ppm_t_ode_)                      :: this
  real(ppm_kind_double), intent(inout)   :: t
  real(ppm_kind_double), intent(in   )   :: dt
  integer,               intent(  out)   :: info
end subroutine

