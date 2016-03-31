      !-------------------------------------------------------------------------
      !  Module       :                 ppm_module_data_ode
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

      MODULE ppm_module_data_ode
        !----------------------------------------------------------------------
        !  Modules
        !----------------------------------------------------------------------
        USE ppm_module_core, ONLY : ppm_kind_double
        IMPLICIT NONE

        !-----------------------------------------------------------------------
        !  Includes
        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------
        !  Time
        !-----------------------------------------------------------------------
        REAL(ppm_kind_double) :: t0

        !-----------------------------------------------------------------------
        ! scheme stuff
        ! _o : order
        ! _m : memory needs
        ! _s : number of stages
        ! _k : suitable kick off scheme
        !-----------------------------------------------------------------------
        INTEGER, DIMENSION(7) :: ppm_ode_scheme_order
        INTEGER, DIMENSION(7) :: ppm_ode_scheme_memsize
        INTEGER, DIMENSION(7) :: ppm_ode_scheme_nstages
        INTEGER, DIMENSION(7) :: ppm_ode_scheme_kickoff
        DATA ppm_ode_scheme_order   /1,2,2,4,2,3,1/
        DATA ppm_ode_scheme_memsize /1,2,2,4,2,1,0/
        DATA ppm_ode_scheme_nstages /1,2,2,4,2,3,999999/
        DATA ppm_ode_scheme_kickoff /1,2,2,4,2,6,7/
        ! suggest DATA ppm_ode_scheme_order    /1,2,2,4,1/
        ! suggest: DATA ppm_ode_scheme_memsize /0,1,1,4,0/
        ! suggest: DATA ppm_ode_scheme_nstages /1,2,2,4,999999/
        ! suggest DATA ppm_ode_scheme_kickoff  /1,2,2,4,7/


        !-----------------------------------------------------------------------
        ! what scheme for which mode
        !-----------------------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: ppm_ode_integr_scheme => NULL()
        INTEGER, DIMENSION(:), POINTER :: ppm_ode_kickoff_scheme => NULL()
        !-----------------------------------------------------------------------
        ! use an adaptive timestep
        !-----------------------------------------------------------------------
        LOGICAL, DIMENSION(:), POINTER :: ppm_ode_adaptive => NULL()
        !-----------------------------------------------------------------------
        ! number of stages a mode uses
        !-----------------------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: ppm_ode_stages => NULL()
        !-----------------------------------------------------------------------
        ! state of a mode
        !-----------------------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: ppm_ode_state => NULL()
        !-----------------------------------------------------------------------
        ! number of sent stages
        !-----------------------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: ppm_ode_sent => NULL()
        !-----------------------------------------------------------------------
        ! size of the buffer
        !-----------------------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: ppm_ode_bfrsize => NULL()
        !-----------------------------------------------------------------------
        ! id lists
        !-----------------------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: ppm_ode_mode_id => NULL()

        !-----------------------------------------------------------------------
        ! some stages for ppm_ode_state
        !-----------------------------------------------------------------------
        INTEGER, PARAMETER :: ppm_ode_state_finished  = 3
        INTEGER, PARAMETER :: ppm_ode_state_kickoff   = 2
        INTEGER, PARAMETER :: ppm_ode_state_running   = 1
        INTEGER, PARAMETER :: ppm_ode_state_inited    = 0

        !-----------------------------------------------------------------------
        ! number of modes
        !-----------------------------------------------------------------------
        INTEGER            :: ppm_max_mid
        INTEGER            :: ppm_max_mid_allocd

        !-----------------------------------------------------------------------
        ! topology ID to work on
        !-----------------------------------------------------------------------
        INTEGER            :: ppm_ode_topoid

      END MODULE ppm_module_data_ode
