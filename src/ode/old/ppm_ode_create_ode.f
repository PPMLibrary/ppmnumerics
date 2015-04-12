      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_ode_create_ode
      !-------------------------------------------------------------------------
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

      SUBROUTINE ppm_ode_create_ode(odeid,bfrsize,nstage,ischeme,kscheme,adaptive,info)
      !!! creates a mode using a given schemes
        !-----------------------------------------------------------------------
        !  Modules
        !-----------------------------------------------------------------------
        USE ppm_module_data
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_error
        USE ppm_module_alloc
        USE ppm_module_write
        USE ppm_module_util_invert_list

        USE ppm_module_data_ode
        IMPLICIT NONE

        !-----------------------------------------------------------------------
        !  Arguments
        !-----------------------------------------------------------------------
        INTEGER,           INTENT(INOUT) :: odeid
        !!! user mode/ode id
        !!!
        !!! TIP: creates an odeid if odeid=-1
        INTEGER,           INTENT(  OUT) :: bfrsize
        !!! size of the buffer the schemes need (user must
        !!! allocate lda*bfrsize
        INTEGER,           INTENT(  OUT) :: nstage
        !!! number of stages
        INTEGER,           INTENT(IN   ) :: ischeme
        !!! integration scheme
        INTEGER, OPTIONAL, INTENT(IN   ) :: kscheme
        !!! kickoff scheme
        LOGICAL, OPTIONAL, INTENT(IN   ) :: adaptive
        !!! adaptive dt flag
        INTEGER,           INTENT(  OUT) :: info
        !!! return status

        !-----------------------------------------------------------------------
        !  Local variables
        !-----------------------------------------------------------------------
        INTEGER               :: mid, tkscheme, iopt
        INTEGER, DIMENSION(3) :: tlda

        CHARACTER(LEN=*), PARAMETER :: caller = 'ppm_ode_create_ode'
        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        CALL substart(caller,t0,info)

        !-----------------------------------------------------------------------
        !  check input arguments
        !-----------------------------------------------------------------------
        IF (ppm_debug.GT.0) THEN
           CALL check_args
           IF (info.NE.0) GOTO 9999
        END IF
        !-----------------------------------------------------------------------
        ! create or use a kickoff scheme
        !-----------------------------------------------------------------------
        IF (PRESENT(kscheme)) THEN
           tkscheme = kscheme
        ELSE
           tkscheme = ppm_ode_scheme_kickoff(ischeme)
        END IF
        !-----------------------------------------------------------------------
        ! check if there is enough space in the arrays
        !-----------------------------------------------------------------------
        mid = ppm_max_mid + 1
        IF (mid.GT.ppm_max_mid_allocd) THEN
           ppm_max_mid_allocd = ppm_max_mid_allocd + 1
           iopt = ppm_param_alloc_grow_preserve
           tlda(1) = ppm_max_mid_allocd
#include "ppm_ode_modalloc.h"
        END IF
        !-----------------------------------------------------------------------
        ! if odeid.lt.0 then generate it for the user
        !-----------------------------------------------------------------------
        IF (odeid.LT.0) odeid = mid
        ppm_user_mid(mid) = odeid
        !-----------------------------------------------------------------------
        ! update the inverse list
        !-----------------------------------------------------------------------
        CALL ppm_util_invert_list(ppm_user_mid,ppm_internal_mid,info)

        ppm_max_mid = maxval(ppm_internal_mid)
        !-----------------------------------------------------------------------
        ! store stuff
        !-----------------------------------------------------------------------
        ppm_ode_ischeme(mid) = ischeme
        ppm_ode_kscheme(mid) = tkscheme
        ppm_ode_state(mid)   = ppm_ode_state_inited
        ppm_ode_stages(mid)  = MAX(ppm_ode_scheme_nstages(ischeme),ppm_ode_scheme_nstages(kscheme))
        ppm_ode_bfrsize(mid) = MAX(ppm_ode_scheme_memsize(ischeme),ppm_ode_scheme_memsize(kscheme))

        bfrsize = ppm_ode_bfrsize(mid)
        nstage  = ppm_ode_stages(mid)

        9999 CONTINUE
        !-----------------------------------------------------------------------
        ! substop
        !-----------------------------------------------------------------------
        CALL substop(caller,t0,info)
        RETURN
        CONTAINS
        SUBROUTINE check_args
           !--------------------------------------------------------------------
           ! check if ppm is initialized
           !--------------------------------------------------------------------
           IF (.NOT.ppm_initialized) THEN
              fail('Please call ppm_init first!',ppm_err_ppm_noinit,exit_point=8888,ppm_error=ppm_error_error)
           END IF
           !--------------------------------------------------------------------
           ! check that schemes are valid
           !--------------------------------------------------------------------
           IF (ischeme.GT.0.AND.ischeme.LE.SIZE(ppm_ode_scheme_order)) THEN
              IF (ppm_ode_scheme_order(ischeme).LT.1) THEN
                 !--------------------------------------------------------------
                 ! scheme not yet implemented
                 !--------------------------------------------------------------
                 fail('ISCHEME not implemented',exit_point=8888,ppm_error=ppm_error_error)
              ENDIF
           ELSE
              !-----------------------------------------------------------------
              ! scheme does not exist
              !-----------------------------------------------------------------
              fail('ISCHEME not implemented',exit_point=8888,ppm_error=ppm_error_error)
           ENDIF
           IF (PRESENT(kscheme)) THEN
              IF (kscheme.GT.0.AND.kscheme.LE.SIZE(ppm_ode_scheme_order)) THEN
                 IF (ppm_ode_scheme_order(kscheme).LT.1) THEN
                    !-----------------------------------------------------------
                    ! scheme not yet implemented
                    !-----------------------------------------------------------
                    fail('KSCHEME not implemented',exit_point=8888,ppm_error=ppm_error_error)
                 ENDIF
              ELSE
                 !--------------------------------------------------------------
                 ! scheme does not exist
                 !--------------------------------------------------------------
                 fail('KSCHEME not implemented',exit_point=8888,ppm_error=ppm_error_error)
              ENDIF
           ENDIF
           !--------------------------------------------------------------------
           ! check odeid
           !--------------------------------------------------------------------
           IF (odeid.GE.0) THEN
              IF (ANY(ppm_user_mid(1:ppm_max_mid).EQ.odeid)) THEN
                 !--------------------------------------------------------------
                 ! already exists. sorry
                 !--------------------------------------------------------------
                 fail('ODEID already exists',exit_point=8888,ppm_error=ppm_error_error)
              ENDIF
           ENDIF
        8888 CONTINUE
        END SUBROUTINE check_args
      END SUBROUTINE ppm_ode_create_ode
