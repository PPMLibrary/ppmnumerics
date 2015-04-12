      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_ode_init
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

      SUBROUTINE ppm_ode_init(topoid,info)
      !!! This routine initializes the ode solver

        !-----------------------------------------------------------------------
        !  Modules
        !-----------------------------------------------------------------------
        USE ppm_module_data
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_error
        USE ppm_module_write

        USE ppm_module_data_ode
        IMPLICIT NONE

        !-----------------------------------------------------------------------
        !  Arguments
        !-----------------------------------------------------------------------
        INTEGER, INTENT(IN   ) :: topoid
        INTEGER, INTENT(  OUT) :: info
        !-----------------------------------------------------------------------
        ! Local Variables
        !-----------------------------------------------------------------------
        INTEGER               :: iopt, mid
        INTEGER, DIMENSION(3) :: tlda

        CHARACTER(LEN=ppm_char) ::  caller="ppm_ode_init"
        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        CALL substart(caller,t0,info)

        !-----------------------------------------------------------------------
        !  check input arguments
        !-----------------------------------------------------------------------
        IF (ppm_debug.GT.0) THEN
           !--------------------------------------------------------------------
           ! check if ppm is initialized
           !--------------------------------------------------------------------
           IF (.NOT.ppm_initialized) THEN
              fail('Please call ppm_init first!',ppm_err_ppm_noinit)
           ENDIF
        ENDIF

        mid = -1

        !-----------------------------------------------------------------------
        ! nullify some module variables
        !-----------------------------------------------------------------------
        NULLIFY(ppm_ode_ischeme)
        NULLIFY(ppm_ode_adaptive)
        NULLIFY(ppm_ode_stages)
        NULLIFY(ppm_ode_state)
        NULLIFY(ppm_ode_sent)
        NULLIFY(ppm_ode_bfrsize)
        NULLIFY(ppm_internal_mid)
        NULLIFY(ppm_user_mid)
        NULLIFY(ppm_ode_kscheme)

        ppm_max_mid        = 0
        ppm_max_mid_allocd = 0

        !-----------------------------------------------------------------------
        ! register the ID of the topology to be used for the ODE solver
        !-----------------------------------------------------------------------
        ppm_ode_topoid = topoid

      9999 CONTINUE
        !-----------------------------------------------------------------------
        ! substop
        !-----------------------------------------------------------------------
        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE ppm_ode_init
