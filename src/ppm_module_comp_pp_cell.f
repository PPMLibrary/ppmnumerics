      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_comp_pp_cell
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
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION         1
#define __DOUBLE_PRECISION         2
#define __SINGLE_PRECISION_COMPLEX 3
#define __DOUBLE_PRECISION_COMPLEX 4

#define __INTERNAL                 5
#define __USER_FUNCTION            6
#define __LOOKUP_TABLE             7

      MODULE ppm_module_comp_pp_cell
      !!! This module provides the routines concerned
      !!! with computing kernel interactions within a set of particles.
         !----------------------------------------------------------------------
         !  Define interfaces to the cell list versions
         !----------------------------------------------------------------------
         INTERFACE ppm_comp_pp_cell
            MODULE PROCEDURE ppm_comp_pp_cell_si
            MODULE PROCEDURE ppm_comp_pp_cell_di
            MODULE PROCEDURE ppm_comp_pp_cell_sci
            MODULE PROCEDURE ppm_comp_pp_cell_dci
            MODULE PROCEDURE ppm_comp_pp_cell_su
            MODULE PROCEDURE ppm_comp_pp_cell_du
            MODULE PROCEDURE ppm_comp_pp_cell_scu
            MODULE PROCEDURE ppm_comp_pp_cell_dcu
            MODULE PROCEDURE ppm_comp_pp_cell_st
            MODULE PROCEDURE ppm_comp_pp_cell_dt
            MODULE PROCEDURE ppm_comp_pp_cell_sct
            MODULE PROCEDURE ppm_comp_pp_cell_dct
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KERNEL __INTERNAL
#define __KIND __SINGLE_PRECISION
#include "ppm_comp_pp_cell.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_comp_pp_cell.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_comp_pp_cell.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_comp_pp_cell.f"
#undef __KIND
#undef __KERNEL

#define __KERNEL __USER_FUNCTION
#define __KIND __SINGLE_PRECISION
#include "ppm_comp_pp_cell.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_comp_pp_cell.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_comp_pp_cell.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_comp_pp_cell.f"
#undef __KIND
#undef __KERNEL

#define __KERNEL __LOOKUP_TABLE
#define __KIND __SINGLE_PRECISION
#include "ppm_comp_pp_cell.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_comp_pp_cell.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_comp_pp_cell.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_comp_pp_cell.f"
#undef __KIND
#undef __KERNEL

      END MODULE ppm_module_comp_pp_cell
