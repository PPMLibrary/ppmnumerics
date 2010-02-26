      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_comp_pp_verlet
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
#define __LOOKUP_TABLE             7

      MODULE ppm_module_comp_pp_verlet
      !!! This module provides the routines concerned
      !!! with computing kernel interactions within a set of particles.
         !----------------------------------------------------------------------
         !  Define interfaces to the Verlet list versions
         !----------------------------------------------------------------------
         INTERFACE ppm_comp_pp_verlet
            MODULE PROCEDURE ppm_comp_pp_verlet_si
            MODULE PROCEDURE ppm_comp_pp_verlet_di
            MODULE PROCEDURE ppm_comp_pp_verlet_sci
            MODULE PROCEDURE ppm_comp_pp_verlet_dci
            MODULE PROCEDURE ppm_comp_pp_verlet_su
            MODULE PROCEDURE ppm_comp_pp_verlet_du
            MODULE PROCEDURE ppm_comp_pp_verlet_scu
            MODULE PROCEDURE ppm_comp_pp_verlet_dcu
            MODULE PROCEDURE ppm_comp_pp_verlet_st
            MODULE PROCEDURE ppm_comp_pp_verlet_dt
            MODULE PROCEDURE ppm_comp_pp_verlet_sct
            MODULE PROCEDURE ppm_comp_pp_verlet_dct
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KERNEL __INTERNAL
#define __KIND __SINGLE_PRECISION
#include "ppm_comp_pp_verlet.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_comp_pp_verlet.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_comp_pp_verlet.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_comp_pp_verlet.f"
#undef __KIND
#undef __KERNEL

#define __KERNEL __USER_FUNCTION
#define __KIND __SINGLE_PRECISION
#include "ppm_comp_pp_verlet.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_comp_pp_verlet.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_comp_pp_verlet.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_comp_pp_verlet.f"
#undef __KIND
#undef __KERNEL

#define __KERNEL __LOOKUP_TABLE
#define __KIND __SINGLE_PRECISION
#include "ppm_comp_pp_verlet.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_comp_pp_verlet.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_comp_pp_verlet.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_comp_pp_verlet.f"
#undef __KIND
#undef __KERNEL

      END MODULE ppm_module_comp_pp_verlet
