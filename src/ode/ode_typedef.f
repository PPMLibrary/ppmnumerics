      !----------------------------------------------------------------------
      !  Types
      !----------------------------------------------------------------------
      TYPE,EXTENDS(ppm_t_ode_) :: ppm_t_ode
      CONTAINS
        PROCEDURE :: create_s => ode_create_s
        PROCEDURE :: create_d => ode_create_d

        PROCEDURE :: destroy  => ode_destroy

        PROCEDURE :: step     => ode_step

        PROCEDURE :: map_push => ode_map_push
        PROCEDURE :: map_pop  => ode_map_pop

      END TYPE ppm_t_ode