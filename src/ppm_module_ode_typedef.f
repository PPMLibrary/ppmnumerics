      !-------------------------------------------------------------------------
      !  Module       :                 ppm_module_ode_typedef
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
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !-------------------------------------------------------------------------
      MODULE ppm_module_ode_typedef

      !----------------------------------------------------------------------
      !  Modules
      !----------------------------------------------------------------------
      USE ppm_module_core

      USE ppm_module_numerics_interfaces
      IMPLICIT NONE

      PRIVATE

      !-----------------------------------------------------------------------
      ! ODE states
      !-----------------------------------------------------------------------
      INTEGER, PARAMETER :: ode_state_init      = 0
      INTEGER, PARAMETER :: ode_state_kickoff   = 1 ! currently not used
      INTEGER, PARAMETER :: ode_state_running   = 2
      INTEGER, PARAMETER :: ode_state_finished  = 3

      !----------------------------------------------------------------------
      !  Types
      !----------------------------------------------------------------------
#include "ode/ode_typedef.f"

      PUBLIC :: ode_state_init
      PUBLIC :: ode_state_kickoff
      PUBLIC :: ode_state_running
      PUBLIC :: ode_state_finished

      PUBLIC :: ppm_t_ode

      !----------------------------------------------------------------------
      !  Type-bound procedures
      !----------------------------------------------------------------------
      CONTAINS

#include "ode/ode_typeproc.f"

      END MODULE ppm_module_ode_typedef
