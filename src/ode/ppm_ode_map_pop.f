      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_ode_map_pop
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

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_ode_map_pop_s(odeid,bfr,lda,Npart,mpart,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_ode_map_pop_d(odeid,bfr,lda,npart,mpart,info)
#endif
    !!! pops whats needed of the buffer
        !-----------------------------------------------------------------------
        !  Includes
        !-----------------------------------------------------------------------
#include "ppm_define.h"
        !-----------------------------------------------------------------------
        !  Modules
        !-----------------------------------------------------------------------
        USE ppm_module_data_ode
        USE ppm_module_data
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_map
        USE ppm_module_error
        USE ppm_module_alloc
        IMPLICIT NONE
#if     __KIND == __SINGLE_PRECISION
        INTEGER, PARAMETER :: mk = ppm_kind_single
#else
        INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
        !-----------------------------------------------------------------------
        !  Arguments
        !-----------------------------------------------------------------------
        INTEGER,                    INTENT(in   ) :: odeid
        !!! mode to push
        REAL(mk), DIMENSION(:,:),   POINTER       :: bfr
        !!! buffer to pop into
        INTEGER,                    INTENT(in   ) :: lda
        !!! leading dimension
        INTEGER,                    INTENT(in   ) :: Npart
        !!! number of particles
        INTEGER,                    INTENT(inout) :: Mpart
        INTEGER,                    INTENT(  out) :: info
        !!! Return status
        !-----------------------------------------------------------------------
        ! Local Variables
        !-----------------------------------------------------------------------
        INTEGER                                   :: ldasend
        INTEGER                                   :: throwaway
        INTEGER                                   :: mid, iopt
        INTEGER                                   :: umidmin, umidmax
        INTEGER,                    DIMENSION(2)  :: dime
        CHARACTER(LEN=*)  , PARAMETER             :: caller = 'ppm_ode_map_pop'
        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        CALL substart(caller,t0,info)
        IF(Mpart.EQ.0) THEN
           !--------------------------------------------------------------------
           ! just save the number of stages that we would have sent
           !--------------------------------------------------------------------
           ! already happened in ppm_ode_step
           !--------------------------------------------------------------------
           GOTO 9999
        END IF
        !-----------------------------------------------------------------------
        !  check input arguments
        !-----------------------------------------------------------------------
        IF(ppm_debug.GT.0) THEN
          IF (info.NE.0) THEN
              CALL ppm_error(ppm_err_ppm_noinit,caller,&
                   & 'Error in arguments',__LINE__,info)
              GOTO 9999
            ENDIF
        END IF

        mid     = ppm_internal_mid(odeid)
        ldasend = lda*ppm_ode_sent(mid)
        IF(ldasend.EQ.0) GOTO 9999

        !-----------------------------------------------------------------------
        ! get the stuff
        !-----------------------------------------------------------------------
        CALL ppm_map_part_pop(bfr,ldasend,Npart,mpart,info)
        IF(info.NE.0) THEN
           GOTO 9999
        END IF
        !-----------------------------------------------------------------------
        ! have to blow it up to full buffer size again 
        !-----------------------------------------------------------------------
        IF(ppm_ode_sent(mid).LT.ppm_ode_bfrsize(mid)) THEN
           iopt = ppm_param_alloc_fit_preserve
           dime(1) = ppm_ode_bfrsize(mid)*lda
           dime(2) = Mpart
           CALL ppm_alloc(bfr,dime,iopt,info)
           IF(info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,caller, &
                   & 'growing buffer BFR',__LINE__,info)
              GOTO 9999
           END IF
        END IF
9999    CONTINUE        
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
           IF(.NOT.ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,caller,&
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
              CALL ppm_error(ppm_err_argument,caller, &
                   & 'odeid does not exist',__LINE__,info)
              GOTO 8888
           ELSE
              IF(ppm_internal_mid(odeid).EQ.-HUGE(odeid)) THEN
                 !--------------------------------------------------------------
                 ! user mid does not exist
                 !--------------------------------------------------------------
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,caller,& 
                      & 'odeid does not exist',__LINE__,info)
                 GOTO 8888
              END IF
           END IF
           !--------------------------------------------------------------------
           ! check dimension
           !--------------------------------------------------------------------
           IF(Mpart.LT.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller, &
                   & 'Mpart cannot be <0',__LINE__,info)
              GOTO 8888
           END IF
           IF(lda.LE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller, &
                   & 'LDA must be >00',__LINE__,info)
              GOTO 8888
           END IF
8888   CONTINUE
       END SUBROUTINE check_args
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_ode_map_pop_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_ode_map_pop_d
#endif
