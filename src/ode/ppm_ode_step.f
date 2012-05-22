      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_ode_step
      !-------------------------------------------------------------------------
      !
      !
      !
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License 
      ! as published by the Free Software Foundation, either 
      ! version 3 of the License, or (at your option) any later 
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_ode_step_ss(odeid,rhsfunc,fields_and_parts,discretizaitons,&
      &          bfr,istage,info)
      !!! integrates one stage of a mode
        !-----------------------------------------------------------------------
        !  Includes
        !-----------------------------------------------------------------------
#include "ppm_define.h"
        !-----------------------------------------------------------------------
        !  Modules
        !-----------------------------------------------------------------------
        USE ppm_module_data_ode
        USE ppm_module_data
        USE ppm_module_check_id
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_error
        USE ppm_module_alloc
        USE ppm_module_numerics_data
        USE ppm_module_interfaces
        IMPLICIT NONE
        !-----------------------------------------------------------------------
        !  Arguments
        !-----------------------------------------------------------------------
        INTEGER,                    INTENT(in   ) :: odeid
        !!! id of mode to advance
        PROCEDURE(rhsfunc)                        :: rhsfunc
        !!! function pointer to the rhs
        CLASS(ppm_v_field),       POINTER         :: fields_and_parts
        CLASS(ppm_v_discr_kind),  POINTER         :: discretizations
        REAL(mk), DIMENSION(:,:), POINTER         :: bfr
        !!! buffer
        INTEGER,                    INTENT(in   ) :: istage
        !!! which stage were in
        INTEGER,                    INTENT(  out) :: info
        !!! return status
        !-----------------------------------------------------------------------
        !  Local Variables
        !-----------------------------------------------------------------------
        REAL(mk)                                  :: t, dt
        INTEGER                                   :: scheme, throwaway, i, j
        INTEGER                                   :: umidmin, umidmax
        INTEGER                                   :: mid, ilda
        REAL(mk)                                  :: M_PI
        INTEGER                                   :: stsn
        REAL(mk), DIMENSION(20)                   :: stsnu
        REAL(mk)                                  :: tau
        INTEGER                                   :: topoid
        LOGICAL                                   :: topo_valid
        !-----------------------------------------------------------------------
        !  fill the nu parameters for the sts scheme
        !-----------------------------------------------------------------------
        M_PI = ACOS(-1.0_MK)
        stsnu(1)  = 0.0_mk
        stsnu(5)  = 0.04_mk
        stsnu(7)  = 0.0015_mk
        stsnu(9)  = 0.04_mk
        stsnu(20) = 0.006_mk
        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        CALL substart('ppm_ode_step',t0,info)
        !-----------------------------------------------------------------------
        !  check input arguments
        !-----------------------------------------------------------------------
        IF(Npart.EQ.0) THEN
           !--------------------------------------------------------------------
           ! best case: nothing to do
           !--------------------------------------------------------------------
           time(3) = time(3) + time(4)
           GOTO 9999
        END IF
        IF(ppm_debug.GT.0) THEN
          CALL check_args
          IF (info.NE.0)
              CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                   & 'Error in arguments',__LINE__,info)
              GOTO 9999
          END IF

        END IF ! (ppm_debug.GT.0)
        topoid = ppm_ode_topoid
        mid = ppm_internal_mid(odeid)
        !-----------------------------------------------------------------------
        ! check state if finished, bail out
        !-----------------------------------------------------------------------
        IF(ppm_ode_state(mid).EQ.ppm_ode_state_finished) GOTO 9999
        !-----------------------------------------------------------------------
        ! check istage, if greater that ppm_ode_stages, then
        ! bail out
        !-----------------------------------------------------------------------
        IF(ppm_ode_stages(mid).LT.istage) GOTO 9999
        
        !-----------------------------------------------------------------------
        ! get times and scheme
        !-----------------------------------------------------------------------
        t      = time(3)
        dt     = time(4)
        IF(ppm_ode_state(mid).EQ.ppm_ode_state_kickoff) THEN
           scheme = ppm_ode_kscheme(mid)
        ELSE
           scheme = ppm_ode_ischeme(mid)
        END IF
        IF(ppm_ode_adaptive(mid)) THEN
           !--------------------------------------------------------------------
           ! compute adaptive timestep [TODO]
           !--------------------------------------------------------------------
           CALL ppm_error(ppm_err_argument,'ppm_ode_step',&
                & 'adaptivity not yet implemented',__LINE__,info)

           GOTO 9999
        END IF
        SELECT CASE(scheme)
        CASE(ppm_param_ode_scheme_eulerf)
           !--------------------------------------------------------------------
           !=======
           ! euler:
           !=======
           !--------------------------------------------------------------------

           !--------------------------------------------------------------------
           ! call right hand side and do an euler step
           !--------------------------------------------------------------------
#include "ppm_ode_rhsfunc_macro.h"
#if      __MODE == __SCA           
           DO j=1,npart
              up(j) = up(j) + dt * dup(j)
           END DO
#elif    __MODE == __VEC
           DO i=1,lda
              DO j=1,npart
                 up(i,j) = up(i,j) + dt * dup(i,j)
              END DO
           END DO
#endif 
           t  = t  + dt
           IF(ppm_ode_state(mid).EQ.ppm_ode_state_kickoff) THEN
              ppm_ode_state(mid) = ppm_ode_state_running
           END IF
           ! how much to save of this
           ppm_ode_sent(mid) = 0
         CASE(ppm_param_ode_scheme_sts)
           !-----------------------------------------------------
           !  compute the new dt
           !-----------------------------------------------------
           stsn = ipackdata(1,1)
           if(present(rpackdata)) stsnu(stsn) = rpackdata(1,1)
           tau = dt/((stsnu(stsn)-1.0_mk)*&
           & COS((2.0_mk*REAL(istage,mk)-1.0_mk)/REAL(stsn,mk)*M_PI*0.5_mk)&
           & +1.0_mk+stsnu(stsn))
           !-----------------------------------------------------------------
           !=======
           !  euler:
           !=======
           !-----------------------------------------------------------------
           !  call right hand side and do an euler step
           !-----------------------------------------------------------------
#include "ppm_ode_rhsfunc_macro.h"
#if      __MODE == __SCA           
           DO j=1,npart
              up(j) = up(j) + tau * dup(j)
           END DO
#elif    __MODE == __VEC
           DO i=1,lda
              DO j=1,npart
                 up(i,j) = up(i,j) + tau * dup(i,j)
              END DO
           END DO
#endif
          t  = t  + tau
          IF(ppm_ode_state(mid).EQ.ppm_ode_state_kickoff) THEN
            ppm_ode_state(mid) = ppm_ode_state_running
          END IF
          ! how much to save of this
          ppm_ode_sent(mid) = 0
        CASE(ppm_param_ode_scheme_tvdrk2)
           !--------------------------------------------------------------------
           !============
           ! 2nd tvd rk:
           !============
           !--------------------------------------------------------------------
           SELECT CASE(istage)
           CASE(1)
              !-----------------------------------------------------------------
              ! call rhs, save the old up and do an euler step
              !-----------------------------------------------------------------
#include "ppm_ode_rhsfunc_macro.h"
#if      __MODE == __SCA
              DO i=1,Npart
                 bfr(1,i) = up(i)
                 up(i) = up(i) + dt*dup(i)
              END DO
#elif    __MODE == __VEC
              DO i=1,Npart             
                 DO j=1,lda
                    bfr(j,i) = up(j,i)
                    up(j,i) = up(j,i) + dt*dup(j,i)
                 END DO
              END DO
#endif
              ppm_ode_sent(mid) = 1
           CASE(2) 
              !-----------------------------------------------------------------
              ! call rhs, and do another euler step
              !-----------------------------------------------------------------
#include "ppm_ode_rhsfunc_macro.h"
#if      __MODE == __SCA
              DO i=1,Npart
                 up(i) = up(i) + dt*dup(i)
              END DO
#elif    __MODE == __VEC
              DO i=1,Npart
                 DO j=1,lda
                    up(j,i) = up(j,i) + dt*dup(j,i)
                 END DO
              END DO
#endif
              !-----------------------------------------------------------------
              ! interpolate
              !-----------------------------------------------------------------
#if      __MODE == __SCA
              DO i=1,Npart
                 up(i)   = 0.5_mk * (up(i)   + bfr(1,i)    )
              END DO
#elif    __MODE == __VEC
              DO i=1,Npart                 
                 DO j=1,lda
                 up(j,i) = 0.5_mk * (up(j,i) + bfr(j,i))
              END DO
           END DO
#endif      
              t = t + dt
              ppm_ode_sent(mid) = 0
              IF(ppm_ode_state(mid).EQ.ppm_ode_state_kickoff) THEN
                 ppm_ode_state(mid) = ppm_ode_state_running
              END IF
           END SELECT
        CASE(ppm_param_ode_scheme_midrk2)
           !--------------------------------------------------------------------
           !=============
           ! mid point rk
           !=============
           !--------------------------------------------------------------------
           SELECT CASE(istage)
           CASE(1)
              !-----------------------------------------------------------------
              ! 
              !-----------------------------------------------------------------
#include "ppm_ode_rhsfunc_macro.h"
#if      __MODE == __SCA
              DO i=1,Npart
                 bfr(1,i)     = up(i)
                 up(i) = up(i) + 0.5_mk* dt * dup(i)
              END DO
#elif    __MODE == __VEC
              DO j=1,lda
                 DO i=1,Npart           
                    bfr(j,i) = up(j,i)
                    up(j,i) = up(j,i) + 0.5_mk*dt*dup(j,i)
                 END DO
              END DO
#endif
              ppm_ode_sent(mid) = 1
           CASE(2) 
#include "ppm_ode_rhsfunc_macro.h"
#if      __MODE == __SCA
              DO i=1,Npart
                 up(i)   = bfr(1,i) + dt * dup(i)
              END DO
#elif    __MODE == __VEC
              DO i=1,Npart
                 DO j=1,lda
                    up(j,i) = bfr(j,i) + dt * dup(j,i)
                 END DO
              END DO
#endif                 
              t = t + dt
              ppm_ode_sent(mid) = 0
              IF(ppm_ode_state(mid).EQ.ppm_ode_state_kickoff) THEN
                 ppm_ode_state(mid) = ppm_ode_state_running
              END IF
           END SELECT
        CASE(ppm_param_ode_scheme_rk4)
           !=============
           ! Runge Kutta 4
           !=============
           SELECT CASE(istage)
           CASE(1)
              !-----------------------------------------------------------------
              ! x_n + 1/2 dt k1
              !-----------------------------------------------------------------
#include "ppm_ode_rhsfunc_macro.h"
              DO i=1,Npart
#if      __MODE == __SCA
                 bfr(1,i)     =  up(i)
                 bfr(2,i)     = dup(i) ! k1
#elif    __MODE == __VEC             
                 bfr(1:lda,i) = up(:,i)
                 bfr((lda+1):2*lda,i) = dup(:,i) !k1
#endif
              END DO
              DO i=1,Npart
#if      __MODE == __SCA
                 up(i) = up(i) + 0.5_mk* dt * dup(i)
#elif    __MODE == __VEC
                 DO ilda=1,lda
                    up(ilda,i) = up(ilda,i) + 0.5_mk*dt*dup(ilda,i)
                 END DO
#endif       
              END DO
              ppm_ode_sent(mid) = 2
           CASE(2) 
              !-----------------------------------------------------------------
#include "ppm_ode_rhsfunc_macro.h"
              !-----------------------------------------------------------------
              ! x_n + 1/2 dt k2
              !-----------------------------------------------------------------
              DO i=1,Npart
#if      __MODE == __SCA
                 bfr(3,i)     = dup(i) !k2
#elif    __MODE == __VEC             
                 DO ilda=1,lda
                    bfr((2*lda+ilda),i) = dup(ilda,i) !k2
                 END DO
#endif                 
#if      __MODE == __SCA
                 up(i)   = bfr(1,i) + 0.5_mk* dt * dup(i)
#elif    __MODE == __VEC                 
                 DO ilda=1,lda
                    up(ilda,i) = bfr(ilda,i) + 0.5_mk*dt * dup(ilda,i)
                 END DO
#endif                 
              END DO
              ppm_ode_sent(mid) = 3
           CASE(3) 
              !-----------------------------------------------------------------
              ! 
              !-----------------------------------------------------------------
#include "ppm_ode_rhsfunc_macro.h"
              !-----------------------------------------------------------------
              ! x_n + dt k3
              !-----------------------------------------------------------------
              DO i=1,Npart
#if      __MODE == __SCA
                 bfr(4,i)     = dup(i) !k3
#elif    __MODE == __VEC
                 DO ilda=1,lda
                    bfr((3*lda+ilda),i) = dup(ilda,i) !k3
                 END DO
#endif                 
#if      __MODE == __SCA
                 up(i)   = bfr(1,i) + dt * dup(i)
#elif    __MODE == __VEC                 
                 DO ilda=1,lda
                    up(ilda,i) = bfr(ilda,i) + dt * dup(ilda,i)
                 END DO
#endif                 
              END DO
              ppm_ode_sent(mid) = 4
           CASE(4) 
              !-----------------------------------------------------------------
#include "ppm_ode_rhsfunc_macro.h"
              !-----------------------------------------------------------------
              ! x_n + 1/6 dt (k1 + 2 k2 + 2k3 +k4)
              !-----------------------------------------------------------------
              DO i=1,Npart
#if      __MODE == __SCA
                 up(i) = bfr(1,i) + 1.0_MK/6.0_MK*dt*                    &
     &                   (bfr(2,i) + 2.0_MK*bfr(3,i) + 2.0_MK*bfr(4,i) + up(i))
#elif    __MODE == __VEC
        DO ilda=1,lda
           up(ilda,i) = bfr(ilda,i) + 1.0_MK/6.0_MK*dt*                    &
           &            (bfr((lda+ilda),i) + 2.0_MK*bfr((2*lda+ilda),i) +  &
           &            2.0_MK*bfr((3*lda+ilda),i)+ dup(ilda,i))
        END DO
#endif                 
      END DO
              t = t + dt
              ppm_ode_sent(mid) = 0
              IF(ppm_ode_state(mid).EQ.ppm_ode_state_kickoff) THEN
                 ppm_ode_state(mid) = ppm_ode_state_running
              END IF
           END SELECT
        END SELECT
        !-----------------------------------------------------------------------
        ! pass back dt and time
        !-----------------------------------------------------------------------
        time(3) = t
        time(4) = dt
        !-----------------------------------------------------------------------
        ! stop ode if were ready
        !-----------------------------------------------------------------------
        IF(time(3).GE.time(2)) THEN
           ppm_ode_state(mid) = ppm_ode_state_finished
        END IF
9999    CONTINUE        
        !-----------------------------------------------------------------------
        ! substop
        !-----------------------------------------------------------------------
        CALL substop('ppm_ode_step',t0,info)
        RETURN
           
       SUBROUTINE check_args 
        
        
            IF(.NOT.ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_ode_step',&
                   & 'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
           END IF
           !--------------------------------------------------------------------
           ! check odeid
           !--------------------------------------------------------------------
           umidmin = LBOUND(ppm_internal_mid,1)
           umidmax = UBOUND(ppm_internal_mid,1)
           IF(odeid.LT.umidmin.OR.odeid.GT.umidmax) THEN
              !-----------------------------------------------------------------
              ! user mid does not exist
              !-----------------------------------------------------------------
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                   & 'ODEID does not exist',__LINE__,info)
              GOTO 8888
           ELSE
              IF(ppm_internal_mid(odeid).EQ.-HUGE(odeid)) THEN
                 !--------------------------------------------------------------
                 ! user mid does not exist
                 !--------------------------------------------------------------
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_ode_step',& 
                      & 'ODEID does not exist',__LINE__,info)
                 GOTO 8888
              END IF
           END IF
           !--------------------------------------------------------------------
           ! check dimension
           !--------------------------------------------------------------------
           IF(Npart.LT.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                   & 'Npart cannot be <0',__LINE__,info)
              GOTO 8888
           END IF
           IF(lda.LE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                   & 'LDA must be >00',__LINE__,info)
              GOTO 8888
           END IF
           IF(time(4).LE.0.0_mk) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                   & 'dt must be >=0',__LINE__,info)
              GOTO 8888
           END IF
           IF(time(3).LT.time(1)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                   & 'time must be >= tstart',__LINE__,info)
              GOTO 8888
           END IF
           IF(scheme.EQ.ppm_param_ode_scheme_sts) THEN
                IF(.NOT.PRESENT(ipackdata)) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                         & 'for STS you need to specify N in ipackdata(1,1)',&
                         & __LINE__,info)
                    GOTO 8888
                 ELSE
                    IF(ipackdata(1,1).NE.1.AND.ipackdata(1,1).NE.7.AND.    &
                         & ipackdata(1,1).NE.9.AND.ipackdata(1,1).NE.20) THEN
                       info = ppm_error_error
                       CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                            & 'ipackdata(1,1) must be element of {1,7,9,20}',&
                            & __LINE__,info)
                       GOTO 8888
                    END IF
                END IF
           END IF
           !--------------------------------------------------------------------
           ! check association of up, dup, bfr
           !--------------------------------------------------------------------
           IF(.NOT.ASSOCIATED(up)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                   & 'UP is empty',__LINE__,info)
              GOTO 8888
           END IF
           IF(.NOT.ASSOCIATED(dup)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                   & 'DUP is empty',__LINE__,info)
              GOTO 8888
           END IF
           IF(.NOT.ASSOCIATED(bfr)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                   & 'BFR is empty',__LINE__,info)
              GOTO 8888
           END IF
           CALL ppm_check_topoid(ppm_ode_topoid,topo_valid,info)
           IF (.NOT. topo_valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                   & 'topoid not valid',__LINE__,info)
              GOTO 8888
           ENDIF
8888 CONTINUE
      END SUBROUTINE check_args


#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_ode_step_ss
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_ode_step_ds
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_ode_step_sv
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_ode_step_dv
#endif
#endif
