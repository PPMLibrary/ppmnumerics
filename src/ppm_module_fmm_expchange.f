      !-------------------------------------------------------------------------
      ! Module         :            ppm_module_fmm_expchange
      !-------------------------------------------------------------------------
      !
      ! Purpose       :  fast multipole method module, expchange routine
      !               
      !
      ! Remarks       : 
      !
      ! References    : 
      !
      ! Revisions     :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_fmm_expchange.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:58  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.4  2006/06/29 10:28:37  pchatela
      !  Added vector strengths support
      !
      !  Revision 1.3  2005/09/19 13:03:31  polasekb
      !  code cosmetics
      !
      !  Revision 1.2  2005/08/23 14:30:54  polasekb
      !  now making difference between single/double precision
      !
      !  Revision 1.1  2005/05/27 07:58:29  polasekb
      !  initial implementation
      !  
      !  Revision 0  2004/11/11 16:38:33 polasekb
      !  Start
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

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
!#define __INTEGER          3
!#define __LOGICAL          4
!#define __2D               7
!#define __3D               8
#define __SFIELD           9
#define __VFIELD          10

MODULE ppm_module_fmm_expchange   

  !-----------------------------------------------------------------------------
  ! Define Interface
  !-----------------------------------------------------------------------------

  INTERFACE ppm_fmm_expchange
	MODULE PROCEDURE ppm_fmm_expchange_s_sf
	MODULE PROCEDURE ppm_fmm_expchange_d_sf
	MODULE PROCEDURE ppm_fmm_expchange_s_vf
	MODULE PROCEDURE ppm_fmm_expchange_d_vf
  END INTERFACE

  !-----------------------------------------------------------------------------
  ! INCLUDE THE SOURCES
  !-----------------------------------------------------------------------------

CONTAINS

#define __KIND __SINGLE_PRECISION
#define __DIM __SFIELD
#include "ppm_fmm_expchange.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fmm_expchange.f"
#undef __KIND
#undef __DIM

#define __KIND __SINGLE_PRECISION
#define __DIM __VFIELD
#include "ppm_fmm_expchange.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fmm_expchange.f"
#undef __KIND
#undef __DIM

END MODULE ppm_module_fmm_expchange

