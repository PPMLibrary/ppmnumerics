        !-----------------------------------------------------------------------
        ! ppm_ode_ischeme
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_ode_ischeme, tlda, iopt, info)
        or_fail_alloc('PPM_ODE_ISCHEME',ppm_error=ppm_error_fatal)

        !-----------------------------------------------------------------------
        ! ppm_ode_adaptive
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_ode_adaptive, tlda, iopt, info)
        or_fail_alloc('PPM_ODE_ADAPTIVE',ppm_error=ppm_error_fatal)

        !-----------------------------------------------------------------------
        ! ppm_ode_stages
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_ode_stages, tlda, iopt, info)
        or_fail_alloc('PPM_ODE_STAGES',ppm_error=ppm_error_fatal)

        !-----------------------------------------------------------------------
        ! ppm_ode_state
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_ode_state, tlda, iopt, info)
        or_fail_alloc('PPM_ODE_STATE',ppm_error=ppm_error_fatal)

        !-----------------------------------------------------------------------
        ! ppm_ode_sent
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_ode_sent, tlda, iopt, info)
        or_fail_alloc('PPM_ODE_SENT',ppm_error=ppm_error_fatal)

        !-----------------------------------------------------------------------
        ! ppm_ode_bfrsize
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_ode_bfrsize, tlda, iopt, info)
        or_fail_alloc('PPM_ODE_BFRSIZE',ppm_error=ppm_error_fatal)

        !-----------------------------------------------------------------------
        ! ppm_ode_sent
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_ode_sent, tlda, iopt, info)
        or_fail_alloc('PPM_ODE_SENT',ppm_error=ppm_error_fatal)

        !-----------------------------------------------------------------------
        ! ppm_internal_mid
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_internal_mid, tlda, iopt, info)
        or_fail_alloc('PPM_INTERNAL_MID',ppm_error=ppm_error_fatal)

        !-----------------------------------------------------------------------
        ! ppm_user_mid
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_user_mid, tlda, iopt, info)
        or_fail_alloc('PPM_USER_MID',ppm_error=ppm_error_fatal)

        !-----------------------------------------------------------------------
        ! ppm_ode_kscheme
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_ode_kscheme, tlda, iopt, info)
        or_fail_alloc('PPM_ODE_KSCHEME',ppm_error=ppm_error_fatal)
