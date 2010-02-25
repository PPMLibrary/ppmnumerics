      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_util_gmres
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the GMRES
      !                 solver
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_util_gmres.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2006/05/11 10:27:00  pchatela
      !  Initial insertion
      !
      !  
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      MODULE ppm_module_util_gmres
         INTEGER, PARAMETER   :: ppm_gmres_param_success = 0
         INTEGER, PARAMETER   :: ppm_gmres_param_failure = 1
         INTEGER, PARAMETER   :: ppm_gmres_param_maxiter = 2
               
         !----------------------------------------------------------------------
         !  Define interfaces to the routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_util_gmres
         MODULE PROCEDURE ppm_util_gmres_s
            MODULE PROCEDURE ppm_util_gmres_d
         END INTERFACE
         
         INTERFACE ppm_util_gmres_solveupper
            MODULE PROCEDURE ppm_util_gmres_solveupper_s
            MODULE PROCEDURE ppm_util_gmres_solveupper_d
         END INTERFACE
         
         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_util_gmres_solveupper.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_util_gmres_solveupper.f"
#undef __KIND

#define __KIND __SINGLE_PRECISION
#include "ppm_util_gmres.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_util_gmres.f"
#undef __KIND

      END MODULE ppm_module_util_gmres
