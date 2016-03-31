      template <ppm_p_rhsfunc_t:[ppm_p_rhsfunc_s,ppm_p_rhsfunc_d],rhsfunc_t:[rhsfunc_s,rhsfunc_d]> nointerface suffixes [s,d]
      SUBROUTINE integrator_create(this,variables,rhsfunc,rhs_variables,info,options)
        IMPLICIT NONE

        !----------------------------------------------------------------------
        !  Arguments
        !----------------------------------------------------------------------
        CLASS(ppm_t_integrator)                               :: this

        CLASS(ppm_v_main_abstr)                               :: variables
        !!! This vector holds all entities to be updated by the integrator the
        !!! elements may be fields or particle discretizations.
        !!!
        !!! fields must either have one discretization or wrapped together with the
        !!! intended discretization in a ppm_t_pair object

        PROCEDURE(ppm_p_rhsfunc_t)                            :: rhsfunc
        !!! The right hand side function to be executed by the eulerf::step function

        CLASS(ppm_v_var_discr_pair)                           :: rhs_variables
        !!! The fields or particle properties (discr_data) to be passed to the right
        !!! hand side function.
        !!!
        !!! The elements of the container are (var,discr) pairs where var can be
        !!! either a field or a discr_data and discr is the corresponding
        !!! discretization.

        INTEGER,                                INTENT(  OUT) :: info

        CLASS(ppm_t_options), OPTIONAL, TARGET, INTENT(IN   ) :: options

        !----------------------------------------------------------------------
        !  Variables
        !----------------------------------------------------------------------
        CLASS(ppm_t_main_abstr),       POINTER :: el
        CLASS(ppm_t_main_abstr),       POINTER :: temp
        CLASS(ppm_t_field_),           POINTER :: cfield
!         CLASS(ppm_t_field_),           POINTER :: pfield
        CLASS(ppm_t_part_prop_s_),     POINTER :: cprop_s
        CLASS(ppm_t_part_prop_d_),     POINTER :: cprop_d
        CLASS(ppm_t_field_),           POINTER :: buf
        CLASS(ppm_t_field_),           POINTER :: el_f
        CLASS(ppm_t_discr_data),       POINTER :: prop
        CLASS(ppm_t_discr_info_),      POINTER :: di
        CLASS(ppm_t_discr_kind),       POINTER :: d
        CLASS(ppm_t_field_discr_pair), POINTER :: el_p
        CLASS(ppm_t_var_discr_pair),   POINTER :: el_vp

        INTEGER :: ifield

        CHARACTER(LEN=16) :: bname
        CHARACTER(LEN=16) :: cname

        LOGICAL :: mkbuf = .TRUE.

        !-----------------------------------------------------------------------
        !  start_subroutine
        !-----------------------------------------------------------------------
        start_subroutine("integrator_create")

        mkbuf = this%scheme_memsize.NE.0

        ALLOCATE(this%variables,STAT=info)
        or_fail_alloc("this%fields")

        IF (mkbuf) THEN
           ALLOCATE(this%buffers,STAT=info)
           or_fail_alloc("this%buffers")
        ELSE
           NULLIFY(this%buffers)
        ENDIF

        ALLOCATE(this%changes,STAT=info)
        or_fail_alloc("this%changes")

        ALLOCATE(this%discretizations,STAT=info)
        or_fail_alloc("this%discretizations")

        el => variables%begin()
        ifield = 0
        DO WHILE (ASSOCIATED(el))
          ifield = ifield + 1
          write(cname,'(A,I0)') 'ode_change_',ifield
          IF (mkbuf) THEN
             ALLOCATE(ppm_t_field::buf,STAT=info)
             or_fail_alloc("buf")

             write(bname,'(A,I0)') 'ode_buffer_',ifield
          ENDIF

          SELECT TYPE (el)
          CLASS is (ppm_t_field_)
             el_f => el
             di => el_f%discr_info%begin()

             ALLOCATE(ppm_t_field::cfield,STAT=info)
             or_fail_alloc("cfield")

             CALL cfield%create(el_f%lda,info,name=cname)
             or_fail("cfield%create")

             CALL cfield%discretize_on(di%discr_ptr,info)
             or_fail("Discretizing change failed")

             temp => cfield

             CALL this%changes%push(temp,info)
             or_fail("Pushing change failed")

             CALL this%discretizations%push(di%discr_ptr,info)
             or_fail("Pushing change discretization failed")

             temp => el
             CALL this%variables%push(temp,info)
             or_fail("Pushing field failed")

             IF (mkbuf) THEN
                CALL buf%create(el_f%lda*this%scheme_memsize,info,name=bname)
                or_fail("buf%create")

                CALL buf%discretize_on(di%discr_ptr,info,with_ghosts=.FALSE.)
                or_fail("buf%discretize_on")
             ENDIF

          CLASS is (ppm_t_discr_data)
             prop => el

             SELECT TYPE (parts => prop%discr)
             CLASS is (ppm_t_particles_s)
                NULLIFY(cprop_s)
                CALL parts%create_prop(info,name=cname,lda=prop%lda,part_prop=cprop_s)
                or_fail("creating change property")

                temp => cprop_s

                CALL this%changes%push(temp,info)

             CLASS is (ppm_t_particles_d)
                NULLIFY(cprop_d)
                CALL parts%create_prop(info,name=cname,lda=prop%lda,part_prop=cprop_d)
                or_fail("creating change property")

                temp => cprop_d

                CALL this%changes%push(temp,info)

             CLASS DEFAULT
                fail("Only particle properties allowed",ppm_err_argument)

             END SELECT
             or_fail("Pushing change failed")

             CALL this%discretizations%push(prop%discr,info)
             or_fail("Pushing change discretization failed")

             temp => el

             CALL this%variables%push(temp,info)
             or_fail("Pushing field failed")

             !IF (mkbuf) THEN
             !  CALL buf%create(el_f%lda*this%scheme_memsize,info,name=bname)
             !  CALL buf%discretize_on(di%discr_ptr,info,with_ghosts=.FALSE.)
             !ENDIF

          CLASS is (ppm_t_discr_kind)
             d => el
             SELECT TYPE (parts => d)
             CLASS is (ppm_t_particles_s)
                NULLIFY(cprop_s)
                CALL parts%create_prop(info,name=cname,lda=ppm_dim,part_prop=cprop_s)
                or_fail("creating change property")

                temp => cprop_s
                CALL this%changes%push(temp,info)

             CLASS is (ppm_t_particles_d)
                NULLIFY(cprop_d)
                CALL parts%create_prop(info,name=cname,lda=ppm_dim,part_prop=cprop_d)
                or_fail("creating change property")

                temp => cprop_d
                CALL this%changes%push(temp,info)

             CLASS DEFAULT
                fail("Only particles allowed",ppm_err_argument)

             END SELECT
             or_fail("Pushing change failed")

             CALL this%discretizations%push(d,info)
             or_fail("Pushing change discretization failed")

             temp => d
             CALL this%variables%push(temp,info)
             or_fail("Pushing positions failed")

             IF (mkbuf) THEN
                CALL buf%create(ppm_dim*this%scheme_memsize,info,name=bname)
                CALL buf%discretize_on(d,info,with_ghosts=.FALSE.)
             ENDIF

          CLASS is (ppm_t_field_discr_pair)
             el_p => el
             ALLOCATE(ppm_t_field::cfield,STAT=info)

             CALL cfield%create(el_p%field%lda,info,name=cname)
             or_fail("cfield%create")

             CALL cfield%discretize_on(el_p%discr,info)
             or_fail("Discretizing change failed")

             temp => cfield
             CALL this%changes%push(temp,info)
             or_fail("Pushing change failed")

             CALL this%discretizations%push(el_p%discr,info)
             or_fail("Pushing change discretization failed")

             temp => el_p%field
             CALL this%variables%push(temp,info)
             or_fail("Pushing field failed")

             IF (mkbuf) THEN
                CALL buf%create(el_p%field%lda*this%scheme_memsize,info,name=bname)
                CALL buf%discretize_on(el_p%discr,info,with_ghosts=.FALSE.)
             ENDIF

          CLASS DEFAULT
             fail("variables should only contain fields, props, discrs and pairs",ppm_err_argument)

          END SELECT

          IF (mkbuf) THEN
             CALL this%buffers%push(buf,info)
          ENDIF

          el => variables%next()
        ENDDO !WHILE (ASSOCIATED(el))

        ! copy the rhs_variables vector
        ALLOCATE(this%rhs_variables,STAT=info)
        or_fail_alloc("this%rhs_variables")

        el_vp => rhs_variables%begin()
        DO WHILE (ASSOCIATED(el_vp))
           CALL this%rhs_variables%push(el_vp,info)
           el_vp => rhs_variables%next()
        ENDDO

        ! set the right hand side function
        this%rhsfunc_t => rhsfunc

        end_subroutine()
      END SUBROUTINE integrator_create


      SUBROUTINE integrator_destroy(this,info)
        IMPLICIT NONE

        !----------------------------------------------------------------------
        ! Arguments
        !----------------------------------------------------------------------
        CLASS(ppm_t_integrator)               :: this
        INTEGER,                INTENT(  OUT) :: info

        !----------------------------------------------------------------------
        ! Local Variables
        !----------------------------------------------------------------------
        CLASS(ppm_t_main_abstr), POINTER :: change
        CLASS(ppm_t_main_abstr), POINTER :: buffer

        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        start_subroutine("eulerf_destroy")


        !TODO
        !TOCHECK
        !This part seems to be corrected !!!
        change => this%changes%begin()
        DO WHILE (ASSOCIATED(change))
           !CALL change%destroy(info)  !FIXME
           ! or_fail("Destroying change failed")
           change => this%changes%next()
        ENDDO

        DEALLOCATE(this%changes,STAT=info)
        or_fail_dealloc("this%changes")

        !TODO
        !TOCHECK
        !This part seems to be corrected !!!
        IF (ASSOCIATED(this%buffers)) THEN
           buffer => this%buffers%begin()
           DO WHILE (ASSOCIATED(buffer))
              !CALL buffer%destroy(info) !FIXME
              ! or_fail("Destroying buffer failed")
              buffer => this%buffers%next()
           ENDDO
        ENDIF

        end_subroutine()
      END SUBROUTINE integrator_destroy

      template <kind_type:[ppm_kind_single,ppm_kind_double],part_type:[ppm_t_particles_s,ppm_t_particles_d],prop_type:[ppm_t_part_prop_s,ppm_t_part_prop_d],rhsfunc_t:[rhsfunc_s,rhsfunc_d]> nointerface suffixes [s,d]
      SUBROUTINE integrator_step(this,t,dt,istage,info)
        IMPLICIT NONE

        parameter(MK,INTEGER,kind_type)
        !----------------------------------------------------------------------
        ! Arguments
        !----------------------------------------------------------------------
        CLASS(ppm_t_integrator) :: this

        REAL(MK), INTENT(INOUT) :: t
        REAL(MK), INTENT(IN   ) :: dt

        INTEGER,  INTENT(IN   ) :: istage
        INTEGER,  INTENT(  OUT) :: info
        !----------------------------------------------------------------------
        ! Variables
        !----------------------------------------------------------------------
        CLASS(ppm_t_main_abstr), POINTER :: var
        CLASS(ppm_t_main_abstr), POINTER :: change
        CLASS(ppm_t_discr_kind), POINTER :: discr

        REAL(MK) :: rhs_info
        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        start_subroutine("integrator_step")
        ! nothing happening here
        end_subroutine()

      END SUBROUTINE integrator_step

      template <ppm_p_rhsfunc_t:[ppm_p_rhsfunc_s,ppm_p_rhsfunc_d]> nointerface suffixes [s,d]
      SUBROUTINE eulerf_create(this,variables,rhsfunc,rhs_variables,info,options)
        IMPLICIT NONE

        !----------------------------------------------------------------------
        !  Arguments
        !----------------------------------------------------------------------
        CLASS(ppm_t_eulerf)                                  :: this

        CLASS(ppm_v_main_abstr)                              :: variables
        !!! This vector holds all entities to be updated by the integrator the
        !!! elements may be fields or particle discretizations.
        !!!
        !!! fields must either have one discretization or wrapped together with the
        !!! intended discretization in a ppm_t_pair object

        PROCEDURE(ppm_p_rhsfunc_t)                           :: rhsfunc
        !!! The right hand side function to be executed by the eulerf::step function

        CLASS(ppm_v_var_discr_pair)                          :: rhs_variables
        !!! The fields to be passed to the right hand side function.
        !!!
        !!! The elements of the container can be fields, or pairs continaing a field
        !!! and its intended discretization.

        INTEGER,                                INTENT(  OUT) :: info

        CLASS(ppm_t_options), OPTIONAL, TARGET, INTENT(IN   ) :: options

        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        start_subroutine("eulerf_create")

        this%scheme_order   = 1
        this%scheme_memsize = 0
        this%scheme_nstages = 1
        this%scheme_kickoff = ppm_param_ode_scheme_eulerf

        CALL this%ppm_t_integrator%create(variables,rhsfunc,rhs_variables,info)

        end_subroutine()

      END SUBROUTINE eulerf_create

      template <kind_type:[ppm_kind_single,ppm_kind_double],part_type:[ppm_t_particles_s,ppm_t_particles_d],prop_type:[ppm_t_part_prop_s,ppm_t_part_prop_d],rhsfunc_t:[rhsfunc_s,rhsfunc_d]> nointerface suffixes [s,d]
      SUBROUTINE eulerf_step(this,t,dt,istage,info)
        IMPLICIT NONE

        parameter(MK,INTEGER,kind_type)
        !----------------------------------------------------------------------
        ! Arguments
        !----------------------------------------------------------------------
        CLASS(ppm_t_eulerf)     :: this
        REAL(MK), INTENT(INOUT) :: t
        REAL(MK), INTENT(IN   ) :: dt
        INTEGER,  INTENT(IN   ) :: istage
        INTEGER,  INTENT(  OUT) :: info
        !----------------------------------------------------------------------
        ! Variables
        !----------------------------------------------------------------------
        CLASS(ppm_t_main_abstr), POINTER :: var
        CLASS(ppm_t_main_abstr), POINTER :: change
        CLASS(ppm_t_discr_kind), POINTER :: discr

        REAL(MK) :: rhs_info

        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        start_subroutine("eulerf_step")

        rhs_info = this%rhsfunc_t(this%rhs_variables,t,this%changes)

        var    => this%variables%begin()
        change => this%changes%begin()
        discr  => this%discretizations%begin()
        DO WHILE (ASSOCIATED(var))
           check_associated(change,"Fields and changes should have same length")
           check_associated(discr,"Fields and discretizations should have same length")

           foreach e in discrp(discr) with vars([u=>var],[du=>change]) prec(MK) part_type(part_type) prop_type(prop_type)
             for sca
               u_e = u_e + dt*du_e
             for vec
               u_e(:) = u_e(:) + dt*du_e(:)
             for pos
               u_e(:) = u_e(:) + dt*du_e(:)
           end foreach

           var    => this%variables%next()
           change => this%changes%next()
           discr  => this%discretizations%next()
        ENDDO !ASSOCIATED(var)

        check_false("ASSOCIATED(change)","Fields and changes should have same length")
        check_false("ASSOCIATED(discr)","Fields and discretizations should have same length")

        t  = t  + dt
        ! how much to save of this
        ! TODO FIXME
        ! ppm_ode_sent(mid) = 0

        end_subroutine()

      END SUBROUTINE eulerf_step

      template <stsnu:[stsnu_s,stsnu_d],kind_type:[ppm_kind_single,ppm_kind_double],ppm_p_rhsfunc_t:[ppm_p_rhsfunc_s,ppm_p_rhsfunc_d]> nointerface suffixes [s,d]
      SUBROUTINE sts_create(this,variables,rhsfunc,rhs_variables,info,options)
        IMPLICIT NONE

        parameter(MK,INTEGER,kind_type)
        !----------------------------------------------------------------------
        !  Arguments
        !----------------------------------------------------------------------
        CLASS(ppm_t_sts)                                     :: this

        CLASS(ppm_v_main_abstr)                              :: variables
        !!! This vector holds all entities to be updated by the integrator the
        !!! elements may be fields or particle discretizations.
        !!!
        !!! fields must either have one discretization or wrapped together with the
        !!! intended discretization in a ppm_t_pair object

        PROCEDURE(ppm_p_rhsfunc_t)                           :: rhsfunc
        !!! The right hand side function to be executed by the eulerf::step function

        CLASS(ppm_v_var_discr_pair)                          :: rhs_variables
        !!! The fields to be passed to the right hand side function.
        !!!
        !!! The elements of the container can be fields, or pairs continaing a field
        !!! and its intended discretization.

        INTEGER,                                INTENT(  OUT) :: info

        CLASS(ppm_t_options), OPTIONAL, TARGET, INTENT(IN   ) :: options

        !----------------------------------------------------------------------
        ! Local Variables
        !----------------------------------------------------------------------
        CLASS(ppm_t_sts_options_s), POINTER :: sts_options_s
        CLASS(ppm_t_sts_options_d), POINTER :: sts_options_d

        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        start_subroutine("sts_create")

        this%scheme_order   = 1
        this%scheme_memsize = 0
        this%scheme_nstages = 1000
        this%scheme_kickoff = ppm_param_ode_scheme_sts

        this%stsnu(:)  = 0.0_MK
        this%stsnu(1)  = 0.0_MK
        this%stsnu(5)  = 0.04_MK
        this%stsnu(7)  = 0.0015_MK
        this%stsnu(9)  = 0.04_MK
        this%stsnu(20) = 0.006_MK

        SELECT TYPE(options)
        CLASS is (ppm_t_sts_options_s)
           sts_options_s         => options
           this%nsts             = sts_options_s%nsts
           this%stsnu(this%nsts) = sts_options_s%nu

        CLASS is (ppm_t_sts_options_d)
           sts_options_d         => options
           this%nsts             = sts_options_d%nsts
           this%stsnu(this%nsts) = sts_options_d%nu

        END SELECT

        CALL this%ppm_t_integrator%create(variables,rhsfunc,rhs_variables,info)
        or_fail("this%ppm_t_integrator%create")

        end_subroutine()

      END SUBROUTINE sts_create

      template <kind_type:[ppm_kind_single,ppm_kind_double],part_type:[ppm_t_particles_s,ppm_t_particles_d],prop_type:[ppm_t_part_prop_s,ppm_t_part_prop_d],rhsfunc_t:[rhsfunc_s,rhsfunc_d]> nointerface suffixes [s,d]
      SUBROUTINE sts_step(this,t,dt,istage,info)
        IMPLICIT NONE

        parameter(MK,INTEGER,kind_type)
        !----------------------------------------------------------------------
        ! Arguments
        !----------------------------------------------------------------------
        CLASS(ppm_t_sts)        :: this

        REAL(MK), INTENT(INOUT) :: t
        REAL(MK), INTENT(IN   ) :: dt

        INTEGER,  INTENT(IN   ) :: istage
        INTEGER,  INTENT(  OUT) :: info

        !----------------------------------------------------------------------
        ! Local Variables
        !----------------------------------------------------------------------
        CLASS(ppm_t_main_abstr), POINTER :: var
        CLASS(ppm_t_main_abstr), POINTER :: change
        CLASS(ppm_t_discr_kind), POINTER :: discr

        REAL(MK)            :: rhs_info
        REAL(MK)            :: tau
        REAL(MK), PARAMETER :: m_pi=ACOS(-1.0_MK)

        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        start_subroutine("sts_step")

        !-----------------------------------------------------
        !  compute the new dt
        !-----------------------------------------------------
        IF (MK.EQ.ppm_kind_single) THEN
           tau = dt/((this%stsnu_s(this%nsts) - 1.0_MK) * &
           &     COS((2.0_MK * REAL(istage,MK) - 1.0_MK)/REAL(this%nsts,MK) * m_pi * 0.5_MK) + &
           &     1.0_MK + this%stsnu_s(this%nsts))
        ELSE
           tau = dt/((this%stsnu_d(this%nsts) - 1.0_MK) * &
           &     COS((2.0_MK * REAL(istage,MK) - 1.0_MK)/REAL(this%nsts,MK) * m_pi * 0.5_MK) + &
           &     1.0_MK + this%stsnu_d(this%nsts))
        ENDIF

        rhs_info = this%rhsfunc_t(this%rhs_variables,t,this%changes)

        var    => this%variables%begin()
        change => this%changes%begin()
        discr  => this%discretizations%begin()
        DO WHILE (ASSOCIATED(var))
          check_associated(change,"Fields and changes should have same length")
          check_associated(discr,"Fields and discretizations should have same length")

          foreach e in discrp(discr) with vars([u=>var],[du=>change]) prec(MK) part_type(part_type) prop_type(prop_type)
            for sca
              u_e = u_e + tau*du_e
            for vec
              u_e(:) = u_e(:) + tau*du_e(:)
            for pos
              u_e(:) = u_e(:) + tau*du_e(:)
          end foreach

          var    => this%variables%next()
          change => this%changes%next()
          discr  => this%discretizations%next()
        ENDDO !ASSOCIATED(var)

        check_false("ASSOCIATED(change)","Fields and changes should have same length")
        check_false("ASSOCIATED(discr)","Fields and discretizations should have same length")

        t  = t  + tau

        end_subroutine()

      END SUBROUTINE sts_step


      template <ppm_p_rhsfunc_t:[ppm_p_rhsfunc_s,ppm_p_rhsfunc_d]> nointerface suffixes [s,d]
      SUBROUTINE tvdrk2_create(this,variables,rhsfunc,rhs_variables,info,options)
        IMPLICIT NONE

        !----------------------------------------------------------------------
        !  Arguments
        !----------------------------------------------------------------------
        CLASS(ppm_t_tvdrk2)                                   :: this

        CLASS(ppm_v_main_abstr)                               :: variables
        !!! This vector holds all entities to be updated by the integrator the
        !!! elements may be fields or particle discretizations.
        !!!
        !!! fields must either have one discretization or wrapped together with the
        !!! intended discretization in a ppm_t_pair object

        PROCEDURE(ppm_p_rhsfunc_t)                            :: rhsfunc
        !!! The right hand side function to be executed by the eulerf::step function

        CLASS(ppm_v_var_discr_pair)                           :: rhs_variables
        !!! The fields to be passed to the right hand side function.
        !!!
        !!! The elements of the container can be fields, or pairs continaing a field
        !!! and its intended discretization.

        INTEGER,                                INTENT(  OUT) :: info

        CLASS(ppm_t_options), OPTIONAL, TARGET, INTENT(IN   ) :: options

        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        start_subroutine("tvdrk2_create")

        this%scheme_order   = 2
        this%scheme_memsize = 1
        this%scheme_nstages = 2
        this%scheme_kickoff = ppm_param_ode_scheme_tvdrk2

        ! TODO ALLOCATE buffer
        CALL this%ppm_t_integrator%create(variables,rhsfunc,rhs_variables,info)
        or_fail("this%ppm_t_integrator%create")

        end_subroutine()

      END SUBROUTINE tvdrk2_create

      template <kind_type:[ppm_kind_single,ppm_kind_double],part_type:[ppm_t_particles_s,ppm_t_particles_d],prop_type:[ppm_t_part_prop_s,ppm_t_part_prop_d],rhsfunc_t:[rhsfunc_s,rhsfunc_d]> nointerface suffixes [s,d]
      SUBROUTINE tvdrk2_step(this,t,dt,istage,info)
        IMPLICIT NONE

        parameter(MK,INTEGER,kind_type)
        !----------------------------------------------------------------------
        ! Arguments
        !----------------------------------------------------------------------
        CLASS(ppm_t_tvdrk2)     :: this

        REAL(MK), INTENT(INOUT) :: t
        REAL(MK), INTENT(IN   ) :: dt

        INTEGER,  INTENT(IN   ) :: istage
        INTEGER,  INTENT(  OUT) :: info

        !----------------------------------------------------------------------
        ! Local Variables
        !----------------------------------------------------------------------
        CLASS(ppm_t_main_abstr), POINTER :: var
        CLASS(ppm_t_main_abstr), POINTER :: change
        CLASS(ppm_t_field_),     POINTER :: buffer
        CLASS(ppm_t_discr_kind), POINTER :: discr

        REAL(MK) :: rhs_info

        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        start_subroutine("tvdrk2_step")

        ! TODO check if buffers where mapped
        SELECT CASE (istage)
        CASE (1)
          !----------------------------------------------------------------------
          ! STAGE 1: Store the current field in a secondary buffer and update it
          ! with dt*du
          !----------------------------------------------------------------------
          rhs_info = this%rhsfunc_t(this%rhs_variables,t,this%changes)

          var    => this%variables%begin()
          change => this%changes%begin()
          discr  => this%discretizations%begin()
          buffer => this%buffers%begin()
          DO WHILE (ASSOCIATED(var))
            check_associated(change,"Fields and changes should have same length")
            check_associated(discr,"Fields and discretizations should have same length")

            foreach e in discrp(discr) with vars([u=>var],[du=>change]) buffer([bfr=>buffer]) prec(MK) part_type(part_type) prop_type(prop_type)
              for sca
                bfr_e = u_e
                u_e = u_e + dt*du_e
              for vec
                bfr_e(:) = u_e(:)
                u_e(:) = u_e(:) + dt*du_e(:)
              for pos
                bfr_e(:) = u_e(:)
                u_e(:) = u_e(:) + dt*du_e(:)
            end foreach

            var    => this%variables%next()
            change => this%changes%next()
            discr  => this%discretizations%next()
            buffer => this%buffers%next()
          ENDDO !ASSOCIATED(var)

          check_false("ASSOCIATED(change)","Fields and changes should have same length")
          check_false("ASSOCIATED(discr)","Fields and discretizations should have same length")
          check_false("ASSOCIATED(buffer)","Fields and buffers should have same length")

        CASE (2)
          !----------------------------------------------------------------------
          ! STAGE 2: Do another update and interpolate with previously saved data
          !----------------------------------------------------------------------
          rhs_info = this%rhsfunc_t(this%rhs_variables,t+dt,this%changes)

          var    => this%variables%begin()
          change => this%changes%begin()
          discr  => this%discretizations%begin()
          buffer => this%buffers%begin()
          DO WHILE (ASSOCIATED(var))
            check_associated(change,"Fields and changes should have same length")
            check_associated(discr,"Fields and discretizations should have same length")

            foreach e in discrp(discr) with vars([u=>var],[du=>change]) buffer([bfr=>buffer]) prec(MK) part_type(part_type) prop_type(prop_type)
              for sca
                u_e = u_e + dt*du_e
                u_e  = 0.5_MK * (u_e + bfr_e)
              for vec
                u_e(:) = u_e(:) + dt*du_e(:)
                u_e(:)  = 0.5_MK * (u_e(:) + bfr_e(:))
              for pos
                u_e(:) = u_e(:) + dt*du_e(:)
                u_e(:)  = 0.5_MK * (u_e(:) + bfr_e(:))
            end foreach

            var    => this%variables%next()
            change => this%changes%next()
            discr  => this%discretizations%next()
            buffer => this%buffers%next()
          ENDDO !ASSOCIATED(var)

          check_false("ASSOCIATED(change)","Fields and changes should have same length")
          check_false("ASSOCIATED(discr)","Fields and discretizations should have same length")
          check_false("ASSOCIATED(buffer)","Fields and buffers should have same length")

          t  = t  + dt

        END SELECT

        end_subroutine()
      END SUBROUTINE tvdrk2_step


      template <ppm_p_rhsfunc_t:[ppm_p_rhsfunc_s,ppm_p_rhsfunc_d]> nointerface suffixes [s,d]
      SUBROUTINE midrk2_create(this,variables,rhsfunc,rhs_variables,info,options)
        IMPLICIT NONE

        !----------------------------------------------------------------------
        !  Arguments
        !----------------------------------------------------------------------
        CLASS(ppm_t_midrk2)                    :: this

        CLASS(ppm_v_main_abstr)                :: variables
        !!! This vector holds all entities to be updated by the integrator the
        !!! elements may be fields or particle discretizations.
        !!!
        !!! fields must either have one discretization or wrapped together with the
        !!! intended discretization in a ppm_t_pair object

        PROCEDURE(ppm_p_rhsfunc_t)           :: rhsfunc
        !!! The right hand side function to be executed by the eulerf::step function

        CLASS(ppm_v_var_discr_pair)          :: rhs_variables
        !!! The fields to be passed to the right hand side function.
        !!!
        !!! The elements of the container can be fields, or pairs continaing a field
        !!! and its intended discretization.

        INTEGER,                 INTENT(  OUT) :: info

        CLASS(ppm_t_options),TARGET,OPTIONAL,INTENT(IN   ) :: options

        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        start_subroutine("midrk2_create")

        this%scheme_order   = 2
        this%scheme_memsize = 1
        this%scheme_nstages = 2
        this%scheme_kickoff = ppm_param_ode_scheme_midrk2

        ! TODO ALLOCATE buffer
        CALL this%ppm_t_integrator%create(variables,rhsfunc,rhs_variables,info)
        or_fail("this%ppm_t_integrator%create")

        end_subroutine()

      END SUBROUTINE midrk2_create

      template <kind_type:[ppm_kind_single,ppm_kind_double],part_type:[ppm_t_particles_s,ppm_t_particles_d],prop_type:[ppm_t_part_prop_s,ppm_t_part_prop_d],rhsfunc_t:[rhsfunc_s,rhsfunc_d]> nointerface suffixes [s,d]
      SUBROUTINE midrk2_step(this,t,dt,istage,info)
        IMPLICIT NONE

        parameter(MK,INTEGER,kind_type)
        !----------------------------------------------------------------------
        ! Arguments
        !----------------------------------------------------------------------
        CLASS(ppm_t_midrk2)     :: this

        REAL(MK), INTENT(INOUT) :: t
        REAL(MK), INTENT(IN   ) :: dt

        INTEGER , INTENT(IN   ) :: istage
        INTEGER , INTENT(  OUT) :: info

        !----------------------------------------------------------------------
        ! Local Variables
        !----------------------------------------------------------------------
        CLASS(ppm_t_main_abstr), POINTER :: var
        CLASS(ppm_t_main_abstr), POINTER :: change
        CLASS(ppm_t_field_),     POINTER :: buffer
        CLASS(ppm_t_discr_kind), POINTER :: discr

        REAL(MK) :: rhs_info
        REAL(MK) :: hdt

        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        start_subroutine("midrk2_step")

        hdt=0.5_MK*dt

        ! TODO check if buffers where mapped
        SELECT CASE (istage)
        CASE (1)
          !----------------------------------------------------------------------
          ! STAGE 1: Compute midpoint
          ! Store the current field in a secondary buffer and update it
          ! with 0.5*dt*du
          !----------------------------------------------------------------------
          rhs_info = this%rhsfunc_t(this%rhs_variables,t,this%changes)

          var    => this%variables%begin()
          change => this%changes%begin()
          discr  => this%discretizations%begin()
          buffer => this%buffers%begin()
          DO WHILE (ASSOCIATED(var))
            check_associated(change,"Fields and changes should have same length")
            check_associated(discr,"Fields and discretizations should have same length")

            foreach e in discrp(discr) with vars([u=>var],[du=>change]) buffer([bfr=>buffer]) prec(MK) part_type(part_type) prop_type(prop_type)
              for sca
                bfr_e = u_e
                u_e = u_e + hdt*du_e
              for vec
                bfr_e(:) = u_e(:)
                u_e(:) = u_e(:) + hdt*du_e(:)
              for pos
                bfr_e(:) = u_e(:)
                u_e(:) = u_e(:) + hdt*du_e(:)
            end foreach

            var    => this%variables%next()
            change => this%changes%next()
            discr  => this%discretizations%next()
            buffer => this%buffers%next()
          ENDDO !ASSOCIATED(var)

          check_false("ASSOCIATED(change)","Fields and changes should have same length")
          check_false("ASSOCIATED(discr)","Fields and discretizations should have same length")
          check_false("ASSOCIATED(buffer)","Fields and buffers should have same length")

        CASE (2)
          !----------------------------------------------------------------------
          ! STAGE 2: Do another update and interpolate with previously saved data
          !----------------------------------------------------------------------
          rhs_info = this%rhsfunc_t(this%rhs_variables,t+hdt,this%changes)

          var    => this%variables%begin()
          change => this%changes%begin()
          discr  => this%discretizations%begin()
          buffer => this%buffers%begin()
          DO WHILE (ASSOCIATED(var))
            check_associated(change,"Fields and changes should have same length")
            check_associated(discr,"Fields and discretizations should have same length")

            foreach e in discrp(discr) with vars([u=>var],[du=>change]) buffer([bfr=>buffer]) prec(MK) part_type(part_type) prop_type(prop_type)
              for sca
                u_e = bfr_e + dt*du_e
              for vec
                u_e(:) = bfr_e(:) + dt*du_e(:)
              for pos
                u_e(:) = bfr_e(:) + dt*du_e(:)
            end foreach

            var    => this%variables%next()
            change => this%changes%next()
            discr  => this%discretizations%next()
            buffer => this%buffers%next()
          ENDDO !ASSOCIATED(var)

          check_false("ASSOCIATED(change)","Fields and changes should have same length")
          check_false("ASSOCIATED(discr)","Fields and discretizations should have same length")
          check_false("ASSOCIATED(buffer)","Fields and buffers should have same length")

          t  = t  + dt

        END SELECT

        end_subroutine()

      END SUBROUTINE midrk2_step


      template <ppm_p_rhsfunc_t:[ppm_p_rhsfunc_s,ppm_p_rhsfunc_d]> nointerface suffixes [s,d]
      SUBROUTINE rk4_create(this,variables,rhsfunc,rhs_variables,info,options)
        IMPLICIT NONE

        !----------------------------------------------------------------------
        !  Arguments
        !----------------------------------------------------------------------
        CLASS(ppm_t_rk4)                                     :: this

        CLASS(ppm_v_main_abstr)                              :: variables
        !!! This vector holds all entities to be updated by the integrator the
        !!! elements may be fields or particle discretizations.
        !!!
        !!! fields must either have one discretization or wrapped together with the
        !!! intended discretization in a ppm_t_pair object

        PROCEDURE(ppm_p_rhsfunc_t)                            :: rhsfunc
        !!! The right hand side function to be executed by the eulerf::step function

        CLASS(ppm_v_var_discr_pair)                           :: rhs_variables
        !!! The fields to be passed to the right hand side function.
        !!!
        !!! The elements of the container can be fields, or pairs continaing a field
        !!! and its intended discretization.

        INTEGER,                                INTENT(  OUT) :: info

        CLASS(ppm_t_options), OPTIONAL, TARGET, INTENT(IN   ) :: options

        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        start_subroutine("rk4_create")

        this%scheme_order   = 4
        this%scheme_memsize = 4
        this%scheme_nstages = 4
        this%scheme_kickoff = ppm_param_ode_scheme_rk4

        ! TODO allocate buffer
        CALL this%ppm_t_integrator%create(variables,rhsfunc,rhs_variables,info)
        or_fail("this%ppm_t_integrator%create")

        end_subroutine()

      END SUBROUTINE rk4_create


      template <kind_type:[ppm_kind_single,ppm_kind_double],part_type:[ppm_t_particles_s,ppm_t_particles_d],prop_type:[ppm_t_part_prop_s,ppm_t_part_prop_d],rhsfunc_t:[rhsfunc_s,rhsfunc_d]> nointerface suffixes [s,d]
      SUBROUTINE rk4_step(this,t,dt,istage,info)
        IMPLICIT NONE

        parameter(MK,INTEGER,kind_type)
        !----------------------------------------------------------------------
        ! Arguments
        !----------------------------------------------------------------------
        CLASS(ppm_t_rk4)        :: this

        REAL(MK), INTENT(INOUT) :: t
        REAL(MK), INTENT(IN   ) :: dt

        INTEGER , INTENT(IN   ) :: istage
        INTEGER , INTENT(  OUT) :: info

        !----------------------------------------------------------------------
        ! Local Variables
        !----------------------------------------------------------------------
        CLASS(ppm_t_main_abstr), POINTER :: var
        CLASS(ppm_t_main_abstr), POINTER :: change
        CLASS(ppm_t_field_),     POINTER :: buffer
        CLASS(ppm_t_discr_kind), POINTER :: discr

        REAL(MK) :: rhs_info
        REAL(MK) :: hdt,h6dt

        INTEGER :: bs

        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        start_subroutine("rk4_step")

        ! TODO check if buffers where mapped
        SELECT CASE (istage)
        CASE (1)
          hdt=0.5_MK*dt

          !----------------------------------------------------------------------
          ! STAGE 1: x_n + 1/2 dt k1
          !----------------------------------------------------------------------
          rhs_info = this%rhsfunc_t(this%rhs_variables,t,this%changes)

          var    => this%variables%begin()
          change => this%changes%begin()
          discr  => this%discretizations%begin()
          buffer => this%buffers%begin()
          DO WHILE (ASSOCIATED(var))
            check_associated(change,"Fields and changes should have same length")
            check_associated(discr,"Fields and discretizations should have same length")

            foreach e in discrp(discr) with vars([u=>var],[du=>change]) buffer([bfr=>buffer],vecbuf=true) prec(MK) part_type(part_type) prop_type(prop_type)
              for sca
                bfr_e(1) = u_e
                bfr_e(2) = du_e
                u_e = u_e + hdt*du_e
              for vec
                bfr_e(1:ulda) = u_e(:)
                bfr_e(ulda+1:2*ulda) = du_e(:)
                u_e(:) = u_e(:) + hdt*du_e(:)
              for pos
                bfr_e(1:ulda) = u_e(:)
                bfr_e(ulda+1:2*ulda) = du_e(:)
                u_e(:) = u_e(:) + hdt*du_e(:)
            end foreach

            var => this%variables%next()
            change => this%changes%next()
            discr => this%discretizations%next()
            buffer => this%buffers%next()
          ENDDO !ASSOCIATED(var)

          check_false("ASSOCIATED(change)","Fields and changes should have same length")
          check_false("ASSOCIATED(discr)","Fields and discretizations should have same length")
          check_false("ASSOCIATED(buffer)","Fields and buffers should have same length")

        CASE (2)
          hdt=0.5_MK*dt

          !----------------------------------------------------------------------
          ! STAGE 2: x_n + 1/2 dt k2
          !----------------------------------------------------------------------
          rhs_info = this%rhsfunc_t(this%rhs_variables,t+hdt,this%changes)

          var    => this%variables%begin()
          change => this%changes%begin()
          discr  => this%discretizations%begin()
          buffer => this%buffers%begin()
          DO WHILE (ASSOCIATED(var))
            check_associated(change,"Fields and changes should have same length")
            check_associated(discr,"Fields and discretizations should have same length")

            foreach e in discrp(discr) with vars([u=>var],[du=>change]) buffer([bfr=>buffer],vecbuf=true) prec(MK) part_type(part_type) prop_type(prop_type)
              for sca
                bfr_e(3) = du_e
                u_e = bfr_e(1) + hdt*du_e
              for vec
                bfr_e(2*ulda+1:3*ulda) = du_e(:)
                u_e(:) = bfr_e(1:ulda) + hdt*du_e(:)
              for pos
                bfr_e(2*ulda+1:3*ulda) = du_e(:)
                u_e(:) = bfr_e(1:ulda) + hdt*du_e(:)
            end foreach

            var    => this%variables%next()
            change => this%changes%next()
            discr  => this%discretizations%next()
            buffer => this%buffers%next()
          ENDDO !ASSOCIATED(var)

          check_false("ASSOCIATED(change)","Fields and changes should have same length")
          check_false("ASSOCIATED(discr)","Fields and discretizations should have same length")
          check_false("ASSOCIATED(buffer)","Fields and buffers should have same length")

        CASE (3)
          hdt=0.5_MK*dt

          !----------------------------------------------------------------------
          ! STAGE 3: x_n + dt k3
          !----------------------------------------------------------------------
          rhs_info = this%rhsfunc_t(this%rhs_variables,t+hdt,this%changes)

          var    => this%variables%begin()
          change => this%changes%begin()
          discr  => this%discretizations%begin()
          buffer => this%buffers%begin()
          DO WHILE (ASSOCIATED(var))
            check_associated(change,"Fields and changes should have same length")
            check_associated(discr,"Fields and discretizations should have same length")

            foreach e in discrp(discr) with vars([u=>var],[du=>change]) buffer([bfr=>buffer],vecbuf=true) prec(MK) part_type(part_type) prop_type(prop_type)
              for sca
                bfr_e(4) = du_e
                u_e = bfr_e(1) + dt*du_e
              for vec
                bfr_e(3*ulda+1:4*ulda) = du_e(:)
                u_e(:) = bfr_e(1:ulda) + dt*du_e(:)
              for pos
                bfr_e(3*ulda+1:4*ulda) = du_e(:)
                u_e(:) = bfr_e(1:ulda) + dt*du_e(:)
            end foreach

            var    => this%variables%next()
            change => this%changes%next()
            discr  => this%discretizations%next()
            buffer => this%buffers%next()
          ENDDO !ASSOCIATED(var)

          check_false("ASSOCIATED(change)","Fields and changes should have same length")
          check_false("ASSOCIATED(discr)","Fields and discretizations should have same length")
          check_false("ASSOCIATED(buffer)","Fields and buffers should have same length")

        CASE (4)
          h6dt = 1.0_MK/6.0_MK * dt

          !----------------------------------------------------------------------
          ! STAGE 4: Final step, interpolating the previous results
          ! x_n + 1/6 dt (k1 + 2 k2 + 2k3 +k4)
          !----------------------------------------------------------------------
          rhs_info = this%rhsfunc_t(this%rhs_variables,t+dt,this%changes)

          var => this%variables%begin()
          change => this%changes%begin()
          discr => this%discretizations%begin()
          buffer => this%buffers%begin()
          DO WHILE (ASSOCIATED(var))
            check_associated(change,"Fields and changes should have same length")
            check_associated(discr,"Fields and discretizations should have same length")

            foreach e in discrp(discr) with vars([u=>var],[du=>change]) buffer([bfr=>buffer],vecbuf=true) prec(MK) part_type(part_type) prop_type(prop_type)
              for sca
                u_e = bfr_e(1) + h6dt * (bfr_e(2) + 2.0_MK*bfr_e(3) + 2.0_MK*bfr_e(4) + du_e)
              for vec
                u_e(:) = bfr_e(1:ulda) + h6dt * (bfr_e(ulda+1:2*ulda) + 2.0_MK*bfr_e(2*ulda+1:3*ulda) + 2.0_MK*bfr_e(3*ulda+1:4*ulda) + du_e(:))
              for pos
                u_e(:) = bfr_e(1:ulda) + h6dt * (bfr_e(ulda+1:2*ulda) + 2.0_MK*bfr_e(2*ulda+1:3*ulda) + 2.0_MK*bfr_e(3*ulda+1:4*ulda) + du_e(:))
            end foreach

            var    => this%variables%next()
            change => this%changes%next()
            discr  => this%discretizations%next()
            buffer => this%buffers%next()
          ENDDO !ASSOCIATED(var)

          check_false("ASSOCIATED(change)","Fields and changes should have same length")
          check_false("ASSOCIATED(discr)","Fields and discretizations should have same length")
          check_false("ASSOCIATED(buffer)","Fields and buffers should have same length")

          t  = t  + dt

        END SELECT

        end_subroutine()
      END SUBROUTINE rk4_step

