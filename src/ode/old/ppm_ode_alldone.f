      !-------------------------------------------------------------------------
      !  Function     :                 ppm_ode_alldone
      !-------------------------------------------------------------------------
      !
      !  Purpose      : will say .true. if all modes are integrated
      !
      !  Input        :
      !
      !  Output       : info                   (I) return status
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_ode_alldone.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.6  2004/07/26 11:59:40  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.5  2004/07/26 11:33:03  michaebe
      !  inserted the use of the ppm ode data module.
      !
      !  Revision 1.4  2004/07/26 11:28:13  michaebe
      !  syntax error... removed.
      !
      !  Revision 1.3  2004/07/26 07:46:51  michaebe
      !  Atomized. Otherwise no changes.
      !
      !  Revision 1.2  2004/06/10 16:20:03  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.1  2004/02/19 08:33:53  michaebe
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
      FUNCTION ppm_ode_alldone(info)
      !-----------------------------------------------------------------------
      !  Modules
      !-----------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write

      USE ppm_module_data_ode, ONLY : ppm_ode_state_finished,ppm_ode_state, &
      &   ppm_max_mid
      IMPLICIT NONE
      !-----------------------------------------------------------------------
      !  Arguments
      !-----------------------------------------------------------------------
      INTEGER, INTENT(  OUT) :: info
      LOGICAL                :: ppm_ode_alldone

      CHARACTER(LEN=ppm_char) :: caller="ppm_ode_alldone"
      !-----------------------------------------------------------------------
      !  call substart
      !-----------------------------------------------------------------------
      CALL substart(caller,t0,info)

      IF (ppm_debug.GT.0) THEN
         CALL check
         IF (info.NE.0) GOTO 9999
      ENDIF

      ppm_ode_alldone = ALL(ppm_ode_state(1:ppm_max_mid).EQ.ppm_ode_state_finished)

      9999 CONTINUE

      !-----------------------------------------------------------------------
      ! substop
      !-----------------------------------------------------------------------
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (.NOT.ppm_initialized) THEN
            fail('Please call ppm_init first!',ppm_err_ppm_noinit,exit_point=8888,ppm_error=ppm_error_error)
         ENDIF
      8888 CONTINUE
      END SUBROUTINE check
      END FUNCTION ppm_ode_alldone
