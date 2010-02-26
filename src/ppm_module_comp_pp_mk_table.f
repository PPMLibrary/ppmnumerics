      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_comp_pp_mk_table
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
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

      MODULE ppm_module_comp_pp_mk_table
      !!! [[ppm_comp_pp_mk_table]]
      !!! This module provides the routines concerned
      !!! with computing kernel interactions within a set of particles.
         !----------------------------------------------------------------------
         !  Define interfaces to the lookup table creation routine
         !----------------------------------------------------------------------
         INTERFACE ppm_comp_pp_mk_table
            MODULE PROCEDURE ppm_comp_pp_mk_table_si
            MODULE PROCEDURE ppm_comp_pp_mk_table_di
            MODULE PROCEDURE ppm_comp_pp_mk_table_sci
            MODULE PROCEDURE ppm_comp_pp_mk_table_dci
            MODULE PROCEDURE ppm_comp_pp_mk_table_su
            MODULE PROCEDURE ppm_comp_pp_mk_table_du
            MODULE PROCEDURE ppm_comp_pp_mk_table_scu
            MODULE PROCEDURE ppm_comp_pp_mk_table_dcu
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KERNEL __INTERNAL
#define __KIND __SINGLE_PRECISION
#include "ppm_comp_pp_mk_table.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_comp_pp_mk_table.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_comp_pp_mk_table.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_comp_pp_mk_table.f"
#undef __KIND
#undef __KERNEL

#define __KERNEL __USER_FUNCTION
#define __KIND __SINGLE_PRECISION
#include "ppm_comp_pp_mk_table.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_comp_pp_mk_table.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_comp_pp_mk_table.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_comp_pp_mk_table.f"
#undef __KIND
#undef __KERNEL

      END MODULE ppm_module_comp_pp_mk_table
