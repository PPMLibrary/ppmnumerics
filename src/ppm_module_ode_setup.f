      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_ode_setup
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
      MODULE ppm_module_ode_setup
      !!! Module for ODE setup routines

        !-----------------------------------------------------
        !  Dummy interface
        !-----------------------------------------------------
        INTERFACE ppm_ode_init
           MODULE PROCEDURE ppm_ode_init
        END INTERFACE
        
        INTERFACE ppm_ode_finalize
           MODULE PROCEDURE ppm_ode_finalize
        END INTERFACE
        
        INTERFACE ppm_ode_start
           MODULE PROCEDURE ppm_ode_start
        END INTERFACE
        
        INTERFACE ppm_ode_create_ode
           MODULE PROCEDURE ppm_ode_create_ode
        END INTERFACE

      CONTAINS
#include "ode/ppm_ode_init.f"

#include "ode/ppm_ode_finalize.f"

#include "ode/ppm_ode_start.f"

#include "ode/ppm_ode_create_ode.f"

      END MODULE ppm_module_ode_setup
