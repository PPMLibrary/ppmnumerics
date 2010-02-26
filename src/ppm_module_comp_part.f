      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_comp_part
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_comp_part
      !!! This module contains all user-callable routines to
      !!! compute kernel-PP interactions on a set of particles.
         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_comp_pp_verlet
         USE ppm_module_comp_pp_cell
         USE ppm_module_comp_pp_ring
         USE ppm_module_comp_pp_correct
         USE ppm_module_comp_pp_mk_table
         
      END MODULE ppm_module_comp_part
