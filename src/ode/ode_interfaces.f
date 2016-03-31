
      SUBROUTINE ode_create_s_(this,scheme,variables,rhsfunc,rhs_variables,info,options,kickoff_scheme)
        IMPORT :: ppm_v_main_abstr, ppm_v_var_discr_pair, ppm_t_options
        IMPORT :: ppm_t_ode_, ppm_p_rhsfunc_s
        IMPLICIT NONE
        CLASS(ppm_t_ode_)                                     :: this
        INTEGER,                                INTENT(IN   ) :: scheme
        CLASS(ppm_v_main_abstr)                               :: variables
        PROCEDURE(ppm_p_rhsfunc_s)                            :: rhsfunc
        CLASS(ppm_v_var_discr_pair)                           :: rhs_variables
        INTEGER,                                INTENT(  OUT) :: info
        CLASS(ppm_t_options), OPTIONAL, TARGET, INTENT(IN   ) :: options
        INTEGER,optional, INTENT(IN   )                       :: kickoff_scheme
      END SUBROUTINE

      SUBROUTINE ode_create_d_(this,scheme,variables,rhsfunc,rhs_variables,info,options,kickoff_scheme)
        IMPORT :: ppm_v_main_abstr, ppm_v_var_discr_pair, ppm_t_options
        IMPORT :: ppm_t_ode_, ppm_p_rhsfunc_d
        IMPLICIT NONE
        CLASS(ppm_t_ode_)                                     :: this
        INTEGER,                                INTENT(IN   ) :: scheme
        CLASS(ppm_v_main_abstr)                               :: variables
        PROCEDURE(ppm_p_rhsfunc_d)                            :: rhsfunc
        CLASS(ppm_v_var_discr_pair)                           :: rhs_variables
        INTEGER,                                INTENT(  OUT) :: info
        CLASS(ppm_t_options), OPTIONAL, TARGET, INTENT(IN   ) :: options
        INTEGER,optional, INTENT(IN   )        :: kickoff_scheme
      END SUBROUTINE

      SUBROUTINE ode_destroy_(this,info)
        IMPORT :: ppm_t_ode_
        IMPLICIT NONE
        CLASS(ppm_t_ode_)      :: this
        INTEGER, INTENT(  OUT) :: info
      END SUBROUTINE

      SUBROUTINE ode_step_(this,t,dt,istage,info)
        IMPORT :: ppm_t_ode_, ppm_kind_double
        IMPLICIT NONE
        CLASS(ppm_t_ode_)                    :: this
        REAL(ppm_kind_double), INTENT(INOUT) :: t
        REAL(ppm_kind_double), INTENT(IN   ) :: dt
        INTEGER,               INTENT(IN   ) :: istage
        INTEGER,               INTENT(  OUT) :: info
      END SUBROUTINE

      SUBROUTINE ode_map_push_(this,info)
        IMPORT :: ppm_t_ode_
        IMPLICIT NONE
        CLASS(ppm_t_ode_)      :: this
        INTEGER, INTENT(  OUT) :: info
      END SUBROUTINE

      SUBROUTINE ode_map_pop_(this,info)
        IMPORT :: ppm_t_ode_
        IMPLICIT NONE
        CLASS(ppm_t_ode_)      :: this
        INTEGER, INTENT(  OUT) :: info
      END SUBROUTINE

