      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_ode_step
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
#define __SCA 3
#define __VEC 4

      MODULE ppm_module_ode_step

        !-----------------------------------------------------
        !  Interface
        !-----------------------------------------------------

        INTERFACE ppm_ode_step
           MODULE PROCEDURE ppm_ode_step_ss
           MODULE PROCEDURE ppm_ode_step_ds
           MODULE PROCEDURE ppm_ode_step_sv
           MODULE PROCEDURE ppm_ode_step_dv
        END INTERFACE

#ifdef __F2003
        !-----------------------------------------------------
        !  Abstract RHS Interface
        !-----------------------------------------------------

        ABSTRACT INTERFACE

           ! single precision

           ! scalar
           FUNCTION rhsfunc_ss(topoid,xp,up,dup,lda,npart,ipack,&
                &lpack,rpack,info)
             INTEGER                          , INTENT(IN)  :: topoid
             INTEGER                          , INTENT(IN)  :: lda,npart
             INTEGER                          , INTENT(OUT) :: info
             REAL(KIND(1.0E0)), DIMENSION(:,:), POINTER     :: xp
             REAL(KIND(1.0E0)), DIMENSION(:),   POINTER     :: up
             REAL(KIND(1.0E0)), DIMENSION(:),   POINTER     :: dup
             INTEGER,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: ipack
             LOGICAL,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: lpack
             REAL(kind(1.0E0)), DIMENSION(:,:), INTENT(IN), OPTIONAL :: rpack
             INTEGER                                     :: rhsfunc
           END FUNCTION rhsfunc_ss

           ! vector
           FUNCTION rhsfunc_sv(topoid,xp,up,dup,lda,npart,ipack,&
                &lpack,rpack,info)
             INTEGER                          , INTENT(IN)  :: topoid
             INTEGER                          , INTENT(IN)  :: lda,npart
             INTEGER                          , INTENT(OUT) :: info
             REAL(KIND(1.0E0)), DIMENSION(:,:), POINTER     :: xp,up
             REAL(KIND(1.0E0)), DIMENSION(:,:), POINTER     :: dup
             INTEGER,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: ipack
             LOGICAL,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: lpack
             REAL(kind(1.0E0)), DIMENSION(:,:), INTENT(IN), OPTIONAL :: rpack
             INTEGER                                     :: rhsfunc
           END FUNCTION rhsfunc_sv

           ! double precision

           ! scalar
           FUNCTION rhsfunc_ds(topoid,xp,up,dup,lda,npart,ipack,&
                &lpack,rpack,info)
             INTEGER                          , INTENT(IN)  :: topoid
             INTEGER                          , INTENT(IN)  :: lda,npart
             INTEGER                          , INTENT(OUT) :: info
             REAL(KIND(1.0D0)), DIMENSION(:,:), POINTER     :: xp
             REAL(KIND(1.0D0)), DIMENSION(:),   POINTER     :: up
             REAL(KIND(1.0D0)), DIMENSION(:),   POINTER     :: dup
             INTEGER,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: ipack
             LOGICAL,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: lpack
             REAL(kind(1.0D0)), DIMENSION(:,:), INTENT(IN), OPTIONAL :: rpack
             INTEGER                                     :: rhsfunc
           END FUNCTION rhsfunc_ds

           ! vector
           FUNCTION rhsfunc_dv(topoid,xp,up,dup,lda,npart,ipack,&
                &lpack,rpack,info)
             INTEGER                          , INTENT(IN)  :: topoid
             INTEGER                          , INTENT(IN)  :: lda,npart
             INTEGER                          , INTENT(OUT) :: info
             REAL(KIND(1.0D0)), DIMENSION(:,:), POINTER     :: xp,up
             REAL(KIND(1.0D0)), DIMENSION(:,:), POINTER     :: dup
             INTEGER,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: ipack
             LOGICAL,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: lpack
             REAL(kind(1.0D0)), DIMENSION(:,:), INTENT(IN), OPTIONAL :: rpack
             INTEGER                                     :: rhsfunc
           END FUNCTION rhsfunc_dv

        END INTERFACE

#endif

      CONTAINS
#define __MODE __SCA
#define __KIND __SINGLE_PRECISION
#include "ode/ppm_ode_step.f"
#undef  __KIND
#define __KIND __DOUBLE_PRECISION
#include "ode/ppm_ode_step.f"
#undef  __KIND
#undef  __MODE
#define __MODE __VEC
#define __KIND __SINGLE_PRECISION
#include "ode/ppm_ode_step.f"
#undef  __KIND
#define __KIND __DOUBLE_PRECISION
#include "ode/ppm_ode_step.f"
#undef  __KIND
#undef  __MODE

      END MODULE ppm_module_ode_step
