

subroutine integrator_create_s_(this,variables,rhsfunc,rhs_variables,info,options)
  import ppm_v_main_abstr, ppm_v_var_discr_pair, ppm_t_options
  import ppm_t_integrator_, ppm_p_rhsfunc_s
  class(ppm_t_integrator_)                    :: this
  class(ppm_v_main_abstr)                     :: variables
  procedure(ppm_p_rhsfunc_s)                  :: rhsfunc
  class(ppm_v_var_discr_pair)               :: rhs_variables
  integer,                      intent(  out) :: info
  class(ppm_t_options),target,optional,intent(in   ) :: options
end subroutine

subroutine integrator_create_d_(this,variables,rhsfunc,rhs_variables,info,options)
  import ppm_v_main_abstr, ppm_v_var_discr_pair, ppm_t_options
  import ppm_t_integrator_, ppm_p_rhsfunc_d
  class(ppm_t_integrator_)                    :: this
  class(ppm_v_main_abstr)                     :: variables
  procedure(ppm_p_rhsfunc_d)                  :: rhsfunc
  class(ppm_v_var_discr_pair)               :: rhs_variables
  integer,                      intent(  out) :: info
  class(ppm_t_options),target,optional,intent(in   ) :: options
end subroutine

subroutine integrator_destroy_(this,info)
  import ppm_t_integrator_
  class(ppm_t_integrator_)           :: this
  integer,             intent(  out) :: info
end subroutine

subroutine integrator_step_s_(this,t,dt,istage,info)
  import ppm_t_integrator_, ppm_kind_single
  class(ppm_t_integrator_)                 :: this
  real(ppm_kind_single)            ,intent(inout)    :: t
  real(ppm_kind_single)            ,intent(in   )    :: dt
  integer,                          intent(in   )    :: istage
  integer,                          intent(  out)    :: info
end subroutine

subroutine integrator_step_d_(this,t,dt,istage,info)
  import ppm_t_integrator_, ppm_kind_double
  class(ppm_t_integrator_)                 :: this
  real(ppm_kind_double)            ,intent(inout)    :: t
  real(ppm_kind_double)            ,intent(in   )    :: dt
  integer,                          intent(in   )    :: istage
  integer,                          intent(  out)    :: info
end subroutine
