

subroutine integrator_create_(this,fields,rhsfunc,rhs_fields,info)
  import ppm_t_integrator_, ppm_v_main_abstr, ppm_p_rhsfunc
  class(ppm_t_integrator_)               :: this
  class(ppm_v_main_abstr), pointer       :: fields
  procedure(ppm_p_rhsfunc)               :: rhsfunc
  class(ppm_v_main_abstr), pointer       :: rhs_fields
  integer,                 intent(  out) :: info
end subroutine

subroutine integrator_destroy_(this,info)
  import ppm_t_integrator_
  class(ppm_t_integrator_)           :: this
  integer,             intent(  out) :: info
end subroutine

subroutine integrator_step_(this,t,dt,info)
  import ppm_t_integrator_, ppm_kind_double
  class(ppm_t_integrator_)                 :: this
  real(ppm_kind_double)            ,intent(inout)    :: t
  real(ppm_kind_double)            ,intent(in   )    :: dt
  integer,                          intent(  out)    :: info
end subroutine
