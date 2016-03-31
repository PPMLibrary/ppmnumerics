
      template <ppm_p_rhsfunc_t:[ppm_p_rhsfunc_s,ppm_p_rhsfunc_d]> nointerface suffixes [s,d]
      SUBROUTINE ode_create(this,scheme,variables,rhsfunc,rhs_variables,info,options,kickoff_scheme)
        !-----------------------------------------------------------------------
        !  Modules
        !-----------------------------------------------------------------------
        USE ppm_module_integrator_typedef
        IMPLICIT NONE

        !-----------------------------------------------------------------------
        !  Arguments
        !-----------------------------------------------------------------------
        CLASS(ppm_t_ode)                                      :: this

        INTEGER,                                INTENT(IN   ) :: scheme

        CLASS(ppm_v_main_abstr)                               :: variables

        PROCEDURE(ppm_p_rhsfunc_t)                            :: rhsfunc

        CLASS(ppm_v_var_discr_pair)                           :: rhs_variables

        INTEGER,                                INTENT(  OUT) :: info

        CLASS(ppm_t_options), OPTIONAL, TARGET, INTENT(IN   ) :: options

        INTEGER,              OPTIONAL,         INTENT(IN   ) :: kickoff_scheme

        !-----------------------------------------------------------------------
        ! Local Variables
        !-----------------------------------------------------------------------
        REAL(ppm_kind_double) :: t0

        INTEGER :: kickoff

        CHARACTER(LEN=*), PARAMETER :: caller = 'ode_create'

        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        CALL substart(caller,t0,info)

        SELECT CASE(scheme)
        CASE (ppm_param_ode_scheme_eulerf)
           ! allocate changes array
           ALLOCATE(ppm_t_eulerf::this%integrator,STAT=info)
           or_fail_alloc("this%integrator")

           CALL this%integrator%create(variables,rhsfunc,rhs_variables,info,options)
           or_fail("Creating eulerf failed")

        CASE (ppm_param_ode_scheme_sts)
           ! allocate changes array
           ALLOCATE(ppm_t_sts::this%integrator,STAT=info)
           or_fail_alloc("this%integrator")

           CALL this%integrator%create(variables,rhsfunc,rhs_variables,info,options)
           or_fail("Creating STS failed")

        CASE (ppm_param_ode_scheme_tvdrk2)
           ! allocate changes array
           ALLOCATE(ppm_t_tvdrk2::this%integrator,STAT=info)
           or_fail_alloc("this%integrator")

           CALL this%integrator%create(variables,rhsfunc,rhs_variables,info,options)
           or_fail("Creating tvd RK2 failed")

        CASE (ppm_param_ode_scheme_midrk2)
           ! allocate changes array
           ALLOCATE(ppm_t_midrk2::this%integrator,STAT=info)
           or_fail_alloc("this%integrator")

           CALL this%integrator%create(variables,rhsfunc,rhs_variables,info,options)
           or_fail("Creating mid RK2 failed")

        CASE (ppm_param_ode_scheme_rk4)
           ! allocate changes array
           ALLOCATE(ppm_t_rk4::this%integrator,STAT=info)
           or_fail_alloc("this%integrator")

           CALL this%integrator%create(variables,rhsfunc,rhs_variables,info,options)
           or_fail("Creating RK4 failed")

        CASE DEFAULT
           fail("Integrator not implemented",ppm_error=ppm_error_fatal)

        END SELECT

        ! use default or user kickoff scheme
        kickoff = MERGE(kickoff_scheme,this%integrator%scheme_kickoff,PRESENT(kickoff_scheme))

        IF (kickoff.EQ.scheme) THEN
           NULLIFY(this%kickoff)
           this%state = ode_state_init
        ELSE
           SELECT CASE (kickoff)
           CASE (ppm_param_ode_scheme_eulerf)
              ! allocate changes array
              ALLOCATE(ppm_t_eulerf::this%kickoff,STAT=info)
              or_fail_alloc("this%kickoff")

              CALL this%kickoff%create(variables,rhsfunc,rhs_variables,info,options)
              or_fail("Creating eulerf failed")

           CASE (ppm_param_ode_scheme_sts)
              ! allocate changes array
              ALLOCATE(ppm_t_sts::this%kickoff,STAT=info)
              or_fail_alloc("this%kickoff")

              CALL this%kickoff%create(variables,rhsfunc,rhs_variables,info,options)
              or_fail("Creating STS failed")

           CASE (ppm_param_ode_scheme_tvdrk2)
              ! allocate changes array
              ALLOCATE(ppm_t_tvdrk2::this%kickoff,STAT=info)
              or_fail_alloc("this%kickoff")

              CALL this%kickoff%create(variables,rhsfunc,rhs_variables,info,options)
              or_fail("Creating tvd RK2 failed")

           CASE (ppm_param_ode_scheme_midrk2)
              ! allocate changes array
              ALLOCATE(ppm_t_midrk2::this%kickoff,STAT=info)
              or_fail_alloc("this%kickoff")

              CALL this%kickoff%create(variables,rhsfunc,rhs_variables,info,options)
              or_fail("Creating mid RK2 failed")

           CASE (ppm_param_ode_scheme_rk4)
              ! allocate changes array
              ALLOCATE(ppm_t_rk4::this%kickoff,STAT=info)
              or_fail_alloc("this%kickoff")

              CALL this%kickoff%create(variables,rhsfunc,rhs_variables,info,options)
              or_fail("Creating RK4 failed")

           CASE DEFAULT
              fail("Integrator not implemented",ppm_error=ppm_error_fatal)

           END SELECT

           this%state = ode_state_kickoff
        ENDIF !(kickoff.EQ.scheme)

      9999 CONTINUE
        !-----------------------------------------------------------------------
        ! substop
        !-----------------------------------------------------------------------
        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE ode_create

      SUBROUTINE ode_destroy(this,info)
        IMPLICIT NONE

        !-----------------------------------------------------------------------
        !  Arguments
        !-----------------------------------------------------------------------
        CLASS(ppm_t_ode)       :: this

        INTEGER, INTENT(  OUT) :: info
        !!! return status

        !-----------------------------------------------------------------------
        !  Local Variables
        !-----------------------------------------------------------------------
        REAL(ppm_kind_double) :: t0

        CHARACTER(LEN=*), PARAMETER :: caller="ode_destroy"

        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        CALL substart(caller,t0,info)

        IF (ASSOCIATED(this%kickoff)) THEN
           CALL this%kickoff%destroy(info)
           or_fail("Failed to destroy kickoff!")

           DEALLOCATE(this%kickoff,STAT=info)
           or_fail_dealloc("this%kickoff")
        ENDIF

        CALL this%integrator%destroy(info)
        or_fail("Failed to destroy the integrator!")

        DEALLOCATE(this%integrator,STAT=info)
        or_fail_dealloc("this%integrator")

      9999 CONTINUE
        !-----------------------------------------------------------------------
        ! substop
        !-----------------------------------------------------------------------
        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE ode_destroy

      SUBROUTINE ode_step(this,t,dt,istage,info)
      !!! integrates one stage of a mode

        IMPLICIT NONE

        !-----------------------------------------------------------------------
        !  Arguments
        !-----------------------------------------------------------------------
        CLASS(ppm_t_ode)                     :: this

        REAL(ppm_kind_double), INTENT(INOUT) :: t
        REAL(ppm_kind_double), INTENT(IN   ) :: dt

        INTEGER,               INTENT(IN   ) :: istage
        !!! which stage were in
        INTEGER,               INTENT(  OUT) :: info
        !!! return status

        !-----------------------------------------------------------------------
        !  Local Variables
        !-----------------------------------------------------------------------
        REAL(ppm_kind_double) :: t0

        CHARACTER(LEN=*), PARAMETER :: caller="ode_step"

        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        CALL substart(caller,t0,info)

        SELECT CASE (this%state)
        CASE (ode_state_init)
          CALL this%integrator%step(t,dt,istage,info)
          or_fail("integrator%step")

          this%state = ode_state_running

        CASE (ode_state_kickoff)
          CALL this%kickoff%step(t,dt,istage,info)
          or_fail("kickoff%step")

          this%state = ode_state_running

        CASE (ode_state_running)
          CALL this%integrator%step(t,dt,istage,info)
          or_fail("integrator%step")

        END SELECT

      9999 CONTINUE
        !-----------------------------------------------------------------------
        ! substop
        !-----------------------------------------------------------------------
        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE ode_step

      SUBROUTINE ode_map_push(this,info)
      !!! pushes whats needed of the buffer
        IMPLICIT NONE

        !-----------------------------------------------------------------------
        !  Arguments
        !-----------------------------------------------------------------------
        CLASS(ppm_t_ode)       :: this

        INTEGER, INTENT(  OUT) :: info
        !!! Return status

        !-----------------------------------------------------------------------
        ! Local Variables
        !-----------------------------------------------------------------------
        CLASS(ppm_t_field_),      POINTER :: buffer
        CLASS(ppm_t_discr_info_), POINTER :: di
        CLASS(ppm_t_particles_d), POINTER :: pset

        REAL(ppm_kind_double) :: t0

        CHARACTER(LEN=*), PARAMETER :: caller = 'ode_map_push'

        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        CALL substart(caller,t0,info)

        IF (this%state.EQ.ode_state_kickoff) THEN
           buffer => this%kickoff%buffers%begin()
           DO WHILE (ASSOCIATED(buffer))

              di => buffer%discr_info%begin()
              IF (ASSOCIATED(di)) THEN
                 SELECT TYPE(disc => di%discr_ptr)
                 CLASS IS (ppm_t_particles_d)
                    pset => disc
                    CALL pset%map_push(info,buffer)
                    or_fail("pset%map_push")

                 CLASS DEFAULT
                    fail("This type is not suported right now!",ppm_error=ppm_error_fatal)

                 END SELECT
              ENDIF

              buffer => this%kickoff%buffers%next()
           ENDDO !(ASSOCIATED(buffer))
        ELSE
           buffer => this%integrator%buffers%begin()
           DO WHILE (ASSOCIATED(buffer))

              di => buffer%discr_info%begin()
              IF (ASSOCIATED(di)) THEN
                 SELECT TYPE(disc => di%discr_ptr)
                 CLASS IS (ppm_t_particles_d)
                    pset => disc
                    CALL pset%map_push(info,buffer)
                    or_fail("pset%map_push")

                 CLASS DEFAULT
                    fail("This type is not suported right now!",ppm_error=ppm_error_fatal)

                 END SELECT
              ENDIF

              buffer => this%integrator%buffers%next()
           ENDDO !(ASSOCIATED(buffer))
        ENDIF !(this%state.EQ.ode_state_kickoff)

      9999 CONTINUE
        !-----------------------------------------------------------------------
        ! substop
        !-----------------------------------------------------------------------
        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE ode_map_push


      SUBROUTINE ode_map_pop(this,info)
      !!! pops whats needed of the buffer

        IMPLICIT NONE

        !-----------------------------------------------------------------------
        !  Arguments
        !-----------------------------------------------------------------------
        CLASS(ppm_t_ode)       :: this

        INTEGER, INTENT(  OUT) :: info
        !!! Return status

        !-----------------------------------------------------------------------
        ! Local Variables
        !-----------------------------------------------------------------------
        CLASS(ppm_t_field_),      POINTER :: buffer
        CLASS(ppm_t_discr_info_), POINTER :: di
        CLASS(ppm_t_particles_d), POINTER :: pset

        REAL(ppm_kind_double) :: t0

        CHARACTER(LEN=*), PARAMETER :: caller = 'ode_map_pop'

        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        CALL substart(caller,t0,info)

        IF (this%state.EQ.ode_state_kickoff) THEN
           ! traverse in reverse order!
           buffer => this%kickoff%buffers%last()
           DO WHILE (ASSOCIATED(buffer))

              di => buffer%discr_info%begin()
              IF (ASSOCIATED(di)) THEN
                 SELECT TYPE(disc => di%discr_ptr)
                 CLASS IS (ppm_t_particles_d)
                    pset => disc
                    CALL pset%map_pop(info,buffer)
                    or_fail("pset%map_pop")

                 CLASS DEFAULT
                    fail("This type is not suported right now!",ppm_error=ppm_error_fatal)

                 END SELECT
              ENDIF

              buffer => this%kickoff%buffers%prev()
           ENDDO !(ASSOCIATED(buffer))
        ELSE
           ! traverse in reverse order!
           buffer => this%integrator%buffers%last()
           DO WHILE (ASSOCIATED(buffer))

              di => buffer%discr_info%begin()
              IF (ASSOCIATED(di)) THEN
                 SELECT TYPE(disc => di%discr_ptr)
                 CLASS IS (ppm_t_particles_d)
                    pset => disc
                    CALL pset%map_pop(info,buffer)
                    or_fail("pset%map_pop")

                 CLASS DEFAULT
                    fail("This type is not suported right now!",ppm_error=ppm_error_fatal)

                 END SELECT
              ENDIF

              buffer => this%integrator%buffers%prev()
           ENDDO
        ENDIF !(this%state.EQ.ode_state_kickoff)

      9999 CONTINUE
        !-----------------------------------------------------------------------
        ! substop
        !-----------------------------------------------------------------------
        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE ode_map_pop
