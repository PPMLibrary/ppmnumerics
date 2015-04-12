      !-------------------------------------------------------------------------
      !  Module       :                 ppm_module_integrator_typedef
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

      MODULE ppm_module_integrator_typedef

      !----------------------------------------------------------------------
      !  Modules
      !----------------------------------------------------------------------
      USE ppm_module_core

      USE ppm_module_numerics_interfaces
      IMPLICIT NONE

      PRIVATE

      !----------------------------------------------------------------------
      !  Types
      !----------------------------------------------------------------------
#include "integrator/integrator_typedef.f"

      PUBLIC :: ppm_t_integrator
      PUBLIC :: ppm_t_sts
      PUBLIC :: ppm_t_sts_options_s
      PUBLIC :: ppm_t_sts_options_d

      PUBLIC :: ppm_t_eulerf
      PUBLIC :: ppm_t_tvdrk2
      PUBLIC :: ppm_t_midrk2
      PUBLIC :: ppm_t_rk4

      !----------------------------------------------------------------------
      !  Type-bound procedures
      !----------------------------------------------------------------------
      CONTAINS

#include "integrator/integrator_typeproc.f"

      END MODULE ppm_module_integrator_typedef
