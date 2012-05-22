      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_ode_finalize
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
      
      SUBROUTINE ppm_ode_finalize(info)
      !!! garbage collection of the ode
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
        USE ppm_module_error
        USE ppm_module_alloc
        IMPLICIT NONE
        
        !-----------------------------------------------------------------------
        !  Arguments
        !-----------------------------------------------------------------------
        INTEGER,                    INTENT(  out) :: info
        INTEGER                                   :: iopt
        INTEGER, DIMENSION(3)                     :: tlda
        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        CALL substart('ppm_ode_finalize',t0,info)
        !-----------------------------------------------------------------------
        !  check input arguments
        !-----------------------------------------------------------------------
        IF(ppm_debug.GT.0) THEN
           !--------------------------------------------------------------------
           ! check if ppm is initialized
           !--------------------------------------------------------------------
           IF(.NOT.ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_ode_finalize',&
                   & 'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
           END IF
        END IF
        iopt = ppm_param_dealloc
        tlda(1) = ppm_max_mid_allocd
#include "ppm_ode_modalloc.h"        
9999    CONTINUE        
        !-----------------------------------------------------------------------
        ! substop
        !-----------------------------------------------------------------------
        CALL substop('ppm_ode_finalize',t0,info)
        RETURN
      END SUBROUTINE ppm_ode_finalize
