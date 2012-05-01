      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_ode_map
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
      
      !-------------------------------------------------------------------------
      !  Define data types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      MODULE ppm_module_ode_map

        !-----------------------------------------------------
        !  Interfaces
        !-----------------------------------------------------

        INTERFACE ppm_ode_map_push
           MODULE PROCEDURE ppm_ode_map_push_s
           MODULE PROCEDURE ppm_ode_map_push_d
        END INTERFACE

        INTERFACE ppm_ode_map_pop
           MODULE PROCEDURE ppm_ode_map_pop_s
           MODULE PROCEDURE ppm_ode_map_pop_d
        END INTERFACE

      CONTAINS
#define __KIND __SINGLE_PRECISION
#include "ode/ppm_ode_map_push.f"
#include "ode/ppm_ode_map_pop.f"
#undef  __KIND
#define __KIND __DOUBLE_PRECISION
#include "ode/ppm_ode_map_push.f"
#include "ode/ppm_ode_map_pop.f"
#undef  __KIND

      END MODULE ppm_module_ode_map
