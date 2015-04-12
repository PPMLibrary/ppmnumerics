      SUBROUTINE integrator_create_s_(this,variables,rhsfunc,rhs_variables,info,options)
        IMPORT :: ppm_v_main_abstr, ppm_v_var_discr_pair, ppm_t_options
        IMPORT :: ppm_t_integrator_, ppm_p_rhsfunc_s
        IMPLICIT NONE
        CLASS(ppm_t_integrator_)                              :: this
        CLASS(ppm_v_main_abstr)                               :: variables
        PROCEDURE(ppm_p_rhsfunc_s)                            :: rhsfunc
        CLASS(ppm_v_var_discr_pair)                           :: rhs_variables
        INTEGER,                                INTENT(  OUT) :: info
        CLASS(ppm_t_options), OPTIONAL, TARGET, INTENT(IN   ) :: options
      END SUBROUTINE

      SUBROUTINE integrator_create_d_(this,variables,rhsfunc,rhs_variables,info,options)
        IMPORT :: ppm_v_main_abstr, ppm_v_var_discr_pair, ppm_t_options
        IMPORT :: ppm_t_integrator_, ppm_p_rhsfunc_d
        IMPLICIT NONE
        CLASS(ppm_t_integrator_)                              :: this
        CLASS(ppm_v_main_abstr)                               :: variables
        PROCEDURE(ppm_p_rhsfunc_d)                            :: rhsfunc
        CLASS(ppm_v_var_discr_pair)                           :: rhs_variables
        INTEGER,                                INTENT(  OUT) :: info
        CLASS(ppm_t_options), OPTIONAL, TARGET, INTENT(IN   ) :: options
      END SUBROUTINE

      SUBROUTINE integrator_destroy_(this,info)
        IMPORT :: ppm_t_integrator_
        IMPLICIT NONE
        CLASS(ppm_t_integrator_) :: this
        INTEGER,   INTENT(  OUT) :: info
      END SUBROUTINE

      SUBROUTINE integrator_step_s_(this,t,dt,istage,info)
        IMPORT :: ppm_t_integrator_, ppm_kind_single
        IMPLICIT NONE
        CLASS(ppm_t_integrator_)             :: this
        REAL(ppm_kind_single), INTENT(INOUT) :: t
        REAL(ppm_kind_single), INTENT(IN   ) :: dt
        INTEGER,               INTENT(IN   ) :: istage
        INTEGER,               INTENT(  OUT) :: info
      END SUBROUTINE

      SUBROUTINE integrator_step_d_(this,t,dt,istage,info)
        IMPORT :: ppm_t_integrator_, ppm_kind_double
        IMPLICIT NONE
        CLASS(ppm_t_integrator_)             :: this
        REAL(ppm_kind_double), INTENT(INOUT) :: t
        REAL(ppm_kind_double), INTENT(IN   ) :: dt
        INTEGER,               INTENT(IN   ) :: istage
        INTEGER,               INTENT(  OUT) :: info
      END SUBROUTINE
