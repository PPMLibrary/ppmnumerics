      TYPE,EXTENDS(ppm_t_integrator_) :: ppm_t_integrator
      CONTAINS
        PROCEDURE :: create_s => integrator_create_s
        PROCEDURE :: create_d => integrator_create_d
        PROCEDURE :: destroy  => integrator_destroy
        PROCEDURE :: step_s   => integrator_step_s
        PROCEDURE :: step_d   => integrator_step_d
      END TYPE ppm_t_integrator

      TYPE,EXTENDS(ppm_t_integrator) :: ppm_t_sts
        REAL(ppm_kind_single), DIMENSION(20) :: stsnu_s
        REAL(ppm_kind_double), DIMENSION(20) :: stsnu_d
        INTEGER                              :: nsts = 1
      CONTAINS
        PROCEDURE :: create_s => sts_create_s
        PROCEDURE :: create_d => sts_create_d
        PROCEDURE :: step_s   => sts_step_s
        PROCEDURE :: step_d   => sts_step_d
      END TYPE ppm_t_sts

      TYPE,EXTENDS(ppm_t_options) :: ppm_t_sts_options_s
        INTEGER               :: nsts = 1
        REAL(ppm_kind_single) :: nu   = 0.0_ppm_kind_single
      END TYPE ppm_t_sts_options_s

      TYPE,EXTENDS(ppm_t_options) :: ppm_t_sts_options_d
        INTEGER               :: nsts = 1
        REAL(ppm_kind_double) :: nu   = 0.0_ppm_kind_double
      END TYPE ppm_t_sts_options_d

      TYPE,EXTENDS(ppm_t_integrator) :: ppm_t_eulerf
      CONTAINS
        PROCEDURE :: create_s => eulerf_create_s
        PROCEDURE :: create_d => eulerf_create_d
        PROCEDURE :: step_s   => eulerf_step_s
        PROCEDURE :: step_d   => eulerf_step_d
      END TYPE ppm_t_eulerf

     TYPE,EXTENDS(ppm_t_integrator) :: ppm_t_tvdrk2
      CONTAINS
        PROCEDURE :: create_s => tvdrk2_create_s
        PROCEDURE :: create_d => tvdrk2_create_d
        PROCEDURE :: step_s   => tvdrk2_step_s
        PROCEDURE :: step_d   => tvdrk2_step_d
      END TYPE ppm_t_tvdrk2

      TYPE,EXTENDS(ppm_t_integrator) :: ppm_t_midrk2
      CONTAINS
        PROCEDURE :: create_s => midrk2_create_s
        PROCEDURE :: create_d => midrk2_create_d
        PROCEDURE :: step_s   => midrk2_step_s
        PROCEDURE :: step_d   => midrk2_step_d
      END TYPE ppm_t_midrk2

      TYPE,EXTENDS(ppm_t_integrator) :: ppm_t_rk4
      CONTAINS
        PROCEDURE :: create_s => rk4_create_s
        PROCEDURE :: create_d => rk4_create_d
        PROCEDURE :: step_s   => rk4_step_s
        PROCEDURE :: step_d   => rk4_step_d
      END TYPE ppm_t_rk4