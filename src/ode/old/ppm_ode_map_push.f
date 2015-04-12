      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_ode_map_push
      !-------------------------------------------------------------------------
      !
      !
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_ode_map_push.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.8  2004/08/13 15:44:25  michaebe
      !  included mpart as input argument
      !
      !  Revision 1.7  2004/08/12 13:48:22  michaebe
      !  included check for ldasend -> bail out if 0
      !
      !  Revision 1.6  2004/08/12 13:10:39  michaebe
      !  corrected caller specification in substart
      !
      !  Revision 1.5  2004/07/26 13:49:19  ivos
      !  Removed Routines sections from the header comment.
      !
      !  Revision 1.4  2004/07/26 11:33:04  michaebe
      !  inserted the use of the ppm ode data module.
      !
      !  Revision 1.3  2004/07/26 07:50:51  michaebe
      !  Atomized. Otherwise no changes
      !
      !  Revision 1.2  2004/06/10 16:20:04  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.1  2004/02/19 08:33:56  michaebe
      !  initial implementation.
      !
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
      SUBROUTINE ppm_ode_map_push_s(odeid,bfr,lda,Npart,mpart,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_ode_map_push_d(odeid,bfr,lda,Npart,mpart,info)
#endif
      !!! pushes whats needed of the buffer
        !-----------------------------------------------------------------------
        !  Modules
        !-----------------------------------------------------------------------
        USE ppm_module_data
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_error
        USE ppm_module_alloc
        USE ppm_module_write
        USE ppm_module_map
        USE ppm_module_map_part

        USE ppm_module_data_ode
        IMPLICIT NONE

#if     __KIND == __SINGLE_PRECISION
        INTEGER, PARAMETER :: MK = ppm_kind_single
#else
        INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
        !-----------------------------------------------------------------------
        !  Arguments
        !-----------------------------------------------------------------------
        INTEGER,                  INTENT(IN   ) :: odeid
        !!! mode to push
        REAL(MK), DIMENSION(:,:), POINTER       :: bfr
        !!! buffer to push
        INTEGER,                  INTENT(IN   ) :: lda
        !!! leading dimension
        INTEGER,                  INTENT(IN   ) :: Npart
        !!! number of particles
        INTEGER,                  INTENT(INOUT) :: mpart
        INTEGER,                  INTENT(  OUT) :: info
        !!! return status
        !-----------------------------------------------------------------------
        ! Local Variables
        !-----------------------------------------------------------------------
        INTEGER :: ldasend
        INTEGER :: throwaway
        INTEGER :: mid, umidmax, umidmin

        CHARACTER(LEN=*), PARAMETER :: caller = 'ppm_ode_map_push'

        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        CALL substart(caller,t0,info)

        !-----------------------------------------------------------------------
        ! general remark:
        ! were going to push only the stages that are needed to the
        ! other cpu. But that guy will pop the stuff and create a new
        ! array to comprise  only the stages that weve sent? so the
        ! array may shrink? must not happen.
        !-----------------------------------------------------------------------
        IF (Npart.EQ.0) GOTO 9999
        !--------------------------------------------------------------------
        ! just save the number of stages that we would have sent
        !--------------------------------------------------------------------
        ! already happened in ppm_ode_step
        !--------------------------------------------------------------------

        !-----------------------------------------------------------------------
        !  check input arguments
        !-----------------------------------------------------------------------
        IF (ppm_debug.GT.0) THEN
           CALL check
           IF (info.NE.0) GOTO 9999
        ENDIF

        mid     = ppm_internal_mid(odeid)
        ldasend = lda*ppm_ode_sent(mid)

        IF (ldasend.EQ.0) GOTO 9999

        CALL ppm_map_part_push(bfr,ldasend,Npart,info)
        or_fail("ppm_map_part_push")

      9999 CONTINUE
        !-----------------------------------------------------------------------
        ! substop
        !-----------------------------------------------------------------------
        CALL substop(caller,t0,info)
        RETURN
      CONTAINS
        SUBROUTINE check
           !--------------------------------------------------------------------
           ! check if ppm is initialized
           !--------------------------------------------------------------------
           IF (.NOT.ppm_initialized) THEN
              fail('Please call ppm_init first!',ppm_err_ppm_noinit,exit_point=8888)
           ENDIF

           !--------------------------------------------------------------------
           ! check odeid
           !--------------------------------------------------------------------
           umidmin = LBOUND(ppm_internal_mid,1)
           umidmax = UBOUND(ppm_internal_mid,1)
           IF (odeid.LT.umidmin.OR.odeid.GT.umidmax) THEN
              !-----------------------------------------------------------------
              ! user mid does not exist
              !-----------------------------------------------------------------
              fail('odeid does not exist',exit_point=8888)
           ELSE
              IF (ppm_internal_mid(odeid).EQ.-HUGE(odeid)) THEN
                 !--------------------------------------------------------------
                 ! user mid does not exist
                 !--------------------------------------------------------------
                 fail('odeid does not exist',exit_point=8888)
              ENDIF
           ENDIF

           !--------------------------------------------------------------------
           ! check dimension
           !--------------------------------------------------------------------
           IF (Npart.LT.0) THEN
              fail('Npart cannot be <0',exit_point=8888)
           ENDIF

           IF (lda.LE.0) THEN
              fail('LDA must be >0',exit_point=8888)
           ENDIF

           IF (.NOT.ASSOCIATED(bfr)) THEN
              fail('BFR is empty!',exit_point=8888)
           ENDIF

        8888 CONTINUE
        END SUBROUTINE check_args
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_ode_map_push_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_ode_map_push_d
#endif
