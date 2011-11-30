      !-----------------------------------------------------------------------
      !  Subroutine   :            ppm_mg_prolong  
      !-----------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
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
      !------------------------------------------------------------------------ 
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_prolong_2d_sca_s(mlev,info)
#elif __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_prolong_2d_sca_d(mlev,info)
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_prolong_3d_sca_s(mlev,info)
#elif __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_prolong_3d_sca_d(mlev,info)
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_prolong_2d_vec_s(mlev,info)
#elif __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_prolong_2d_vec_d(mlev,info)
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_prolong_3d_vec_s(mlev,info)
#elif __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_prolong_3d_vec_d(mlev,info)
#endif
#endif
#endif
         !!! In this routine we prolong the corrections from coarser to 
         !!! finer levels 
         !----------------------------------------------------------------------
         !  Includes
         !----------------------------------------------------------------------
#include "ppm_define.h"
         !----------------------------------------------------------------------
         !  Modules 
         !----------------------------------------------------------------------
         USE ppm_module_data
         USE ppm_module_data_mg
         USE ppm_module_substart
         USE ppm_module_substop
         USE ppm_module_error
         USE ppm_module_alloc
         USE ppm_module_write
         IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION
         INTEGER, PARAMETER :: MK = ppm_kind_single
#else
         INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
         !----------------------------------------------------------------------
         !  Arguments     
         !----------------------------------------------------------------------
         INTEGER,                   INTENT(IN)      ::  mlev
         !!! current level in V-cycle
         INTEGER,                   INTENT(INOUT)   ::  info
         !----------------------------------------------------------------------
         !  Local variables 
         !----------------------------------------------------------------------
         CHARACTER(LEN=256)                         :: cbuf
         INTEGER                                    :: isub,j,j2,i,i2
         INTEGER,DIMENSION(5)                       :: ldl5,ldu5 
         INTEGER,DIMENSION(4)                       :: ldl4,ldu4 
         INTEGER,DIMENSION(4)                       :: ldl3,ldu3 
         INTEGER                                    :: iopt,topoid
         INTEGER                                    :: aa,bb,cc,dd,ee,gg,iface
#if __MESH_DIM == __3D
         INTEGER                                    :: k,k2
#endif
         INTEGER                                    :: mlevp1,ilda
         REAL(MK)                                   :: t0
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
         TYPE(mg_field_2d_sca_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
         TYPE(mg_field_2d_sca_d),DIMENSION(:,:),POINTER :: mgfield
#endif 
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
         TYPE(mg_field_3d_sca_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
         TYPE(mg_field_3d_sca_d),DIMENSION(:,:),POINTER :: mgfield
#endif 
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
         TYPE(mg_field_2d_vec_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
         TYPE(mg_field_2d_vec_d),DIMENSION(:,:),POINTER :: mgfield
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
         TYPE(mg_field_3d_vec_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
         TYPE(mg_field_3d_vec_d),DIMENSION(:,:),POINTER :: mgfield
#endif
#endif
#endif
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
         REAL(MK),DIMENSION(:,:),POINTER :: tuc
         REAL(MK),DIMENSION(:,:),POINTER :: puc
#elif __MESH_DIM == __3D
        REAL(MK),DIMENSION(:,:,:),POINTER :: tuc
        REAL(MK),DIMENSION(:,:,:),POINTER :: puc
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
       REAL(MK),DIMENSION(:,:,:),POINTER :: tuc
       REAL(MK),DIMENSION(:,:,:),POINTER :: puc
#elif __MESH_DIM == __3D
       REAL(MK),DIMENSION(:,:,:,:),POINTER :: tuc
       REAL(MK),DIMENSION(:,:,:,:),POINTER :: puc
#endif
#endif
         !----------------------------------------------------------------------
         !Externals
         !----------------------------------------------------------------------

         !----------------------------------------------------------------------
         !Initialize
         !----------------------------------------------------------------------
         CALL substart('ppm_mg_prolong',t0,info)
         !----------------------------------------------------------------------
         !  Check arguments
         !----------------------------------------------------------------------
         IF (ppm_debug .GT. 0) THEN
             IF (mlev.LT.1) THEN
                   info = ppm_error_error
                   CALL ppm_error(ppm_err_argument,'ppm_mg_prolong',  &
      &                'level must be >0',__LINE__,info)
                   GOTO 9999
             ENDIF
         ENDIF
         !----------------------------------------------------------------------
         !Definition of necessary variables and allocation of arrays
         !----------------------------------------------------------------------
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
         mgfield=>mgfield_2d_sca_s
#elif __KIND == __DOUBLE_PRECISION
         mgfield=>mgfield_2d_sca_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
         mgfield=>mgfield_3d_sca_s
#elif __KIND == __DOUBLE_PRECISION
         mgfield=>mgfield_3d_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
         mgfield=>mgfield_2d_vec_s
#elif __KIND == __DOUBLE_PRECISION
         mgfield=>mgfield_2d_vec_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
         mgfield=>mgfield_3d_vec_s
#elif __KIND == __DOUBLE_PRECISION
         mgfield=>mgfield_3d_vec_d
#endif
#endif
#endif
         !----------------------------------------------------------------------
         !Implementation
         !----------------------------------------------------------------------
         mlevp1 = mlev + 1
         IF (ppm_debug.GT.0) THEN
          WRITE(cbuf,*) 'WELCOME TO THE PROLONG LEVEL:',mlev
          CALL PPM_WRITE(ppm_rank,'mg_prolong',cbuf,info)
         ENDIF
#if __DIM == __SFIELD
#if __MESH_DIM == __2D           
         !----------------------------------------------------------------------
         !prolongation using a 9-point operator
         !----------------------------------------------------------------------
          DO isub=1,nsubs
            tuc=>mgfield(isub,mlevp1)%uc   
            puc=>mgfield(isub,mlev)%uc   
            DO j=1,max_node(2,mlevp1)
               j2=2*j
               DO i=1,max_node(1,mlevp1)
                  i2=2*i
                     puc(i2-1,j2-1) = &
      &                              puc(i2-1,j2-1) + &
      &                              tuc(i,j)
                     puc(i2,j2-1) = & 
      &                              puc(i2,j2-1) + &
      &                    0.5_MK * (tuc(i,j)+&
      &                              tuc(i+1,j))
                     puc(i2-1,j2) = &
      &                              puc(i2-1,j2) + &
      &                   0.5_MK * ( tuc(i,j) + & 
      &                              tuc(i,j+1))
                     puc(i2,j2)  = &
      &                             puc(i2,j2) + &
      &                  0.25_MK * ( tuc(i,j)+&
      &                              tuc(i+1,j) + &
      &                              tuc(i+1,j+1)+&
      &                              tuc(i,j+1))
               ENDDO
            ENDDO
         ENDDO
#elif __MESH_DIM == __3D
         !----------------------------------------------------------------------
         !prolongation using a 27-point operator
         !----------------------------------------------------------------------
          DO isub=1,nsubs
            tuc=>mgfield(isub,mlevp1)%uc 
            puc=>mgfield(isub,mlev)%uc 
            DO k=1,max_node(3,mlevp1)
               k2=2*k
               DO j=1,max_node(2,mlevp1)
                  j2=2*j
                  DO i=1,max_node(1,mlevp1)
                     i2=2*i
                        puc(i2-1,j2-1,k2-1) = &
      &                          puc(i2-1,j2-1,k2-1) + &
      &                          tuc(i,j,k)
                        puc(i2,j2-1,k2-1) = & 
      &                             puc(i2,j2-1,k2-1) + &
      &                   0.5_MK * (tuc(i,j,k)+&
      &                             tuc(i+1,j,k))
                        puc(i2-1,j2,k2-1) = &
      &                             puc(i2-1,j2,k2-1) + &
      &                  0.5_MK * ( tuc(i,j,k) + & 
      &                             tuc(i,j+1,k))
                        puc(i2,j2,k2-1)  = &
      &                            puc(i2,j2,k2-1) + &
      &                 0.25_MK * ( tuc(i,j,k)+&
      &                             tuc(i+1,j,k) + &
      &                             tuc(i+1,j+1,k)+&
      &                             tuc(i,j+1,k))
                        puc(i2-1,j2-1,k2) = &
      &                        puc(i2-1,j2-1,k2) + &
      &                 0.5_MK * (tuc(i,j,k)+&
      &                          tuc(i,j,k+1))
                        puc(i2,j2-1,k2)  = &
      &                             puc(i2,j2-1,k2) + &
      &                 0.25_MK * ( tuc(i,j,k)+&
      &                             tuc(i+1,j,k) + &
      &                             tuc(i+1,j,k+1)+&
      &                             tuc(i,j,k+1))
                        puc(i2-1,j2,k2)  = &
      &                             puc(i2-1,j2,k2) + &
      &                 0.25_MK * ( tuc(i,j,k)+&
      &                             tuc(i,j+1,k) + &
      &                             tuc(i,j,k+1)+&
      &                             tuc(i,j+1,k+1))
                       puc(i2,j2,k2)  = &
      &                             puc(i2,j2,k2) + &
      &                 0.125_MK * (tuc(i,j,k)+&
      &                            tuc(i+1,j,k) + &
      &                            tuc(i+1,j+1,k)+&
      &                            tuc(i,j+1,k)+&
      &                            tuc(i,j,k+1) + &
      &                            tuc(i+1,j,k+1)+&
      &                            tuc(i,j+1,k+1)+&
      &                            tuc(i+1,j+1,k+1))
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D           
         !----------------------------------------------------------------------
         !prolongation using a 9-point operator
         !----------------------------------------------------------------------
          DO isub=1,nsubs
             tuc=>mgfield(isub,mlevp1)%uc 
             puc=>mgfield(isub,mlev)%uc 
             DO j=1,max_node(2,mlevp1)
               j2=2*j
               DO i=1,max_node(1,mlevp1)
                  i2=2*i
                DO ilda=1,vecdim
                     puc(ilda,i2-1,j2-1) = &
      &                              puc(ilda,i2-1,j2-1) + &
      &                              tuc(ilda,i,j)
                     puc(ilda,i2,j2-1) = & 
      &                              puc(ilda,i2,j2-1) + &
      &                    0.5_MK * (tuc(ilda,i,j)+&
      &                              tuc(ilda,i+1,j))
                     puc(ilda,i2-1,j2) = &
      &                              puc(ilda,i2-1,j2) + &
      &                   0.5_MK * ( tuc(ilda,i,j) + & 
      &                              tuc(ilda,i,j+1))
                     puc(ilda,i2,j2)  = &
      &                             puc(ilda,i2,j2) + &
      &                  0.25_MK * ( tuc(ilda,i,j)+&
      &                              tuc(ilda,i+1,j) + &
      &                              tuc(ilda,i+1,j+1)+&
      &                              tuc(ilda,i,j+1))
                ENDDO 
               ENDDO
            ENDDO
         ENDDO
#elif __MESH_DIM == __3D
         !----------------------------------------------------------------------
         !prolongation using a 27-point operator
         !----------------------------------------------------------------------
          DO isub=1,nsubs
             tuc=>mgfield(isub,mlevp1)%uc 
             puc=>mgfield(isub,mlev)%uc 
             DO k=1,max_node(3,mlevp1)
               k2=2*k
               DO j=1,max_node(2,mlevp1)
                  j2=2*j
                  DO i=1,max_node(1,mlevp1)
                     i2=2*i
#ifdef __VECTOR
                        puc(1,i2-1,j2-1,k2-1) = &
      &                          puc(1,i2-1,j2-1,k2-1) + &
      &                          tuc(1,i,j,k)
                        puc(1,i2,j2-1,k2-1) = & 
      &                             puc(1,i2,j2-1,k2-1) + &
      &                   0.5_MK * (tuc(1,i,j,k)+&
      &                             tuc(1,i+1,j,k))
                        puc(1,i2-1,j2,k2-1) = &
      &                             puc(1,i2-1,j2,k2-1) + &
      &                  0.5_MK * ( tuc(1,i,j,k) + & 
      &                             tuc(1,i,j+1,k))
                        puc(1,i2,j2,k2-1)  = &
      &                            puc(1,i2,j2,k2-1) + &
      &                 0.25_MK * ( tuc(1,i,j,k)+&
      &                             tuc(1,i+1,j,k) + &
      &                             tuc(1,i+1,j+1,k)+&
      &                             tuc(1,i,j+1,k))
                        puc(1,i2-1,j2-1,k2) = &
      &                        puc(1,i2-1,j2-1,k2) + &
      &                 0.5_MK * (tuc(1,i,j,k)+&
      &                          tuc(1,i,j,k+1))
                        puc(1,i2,j2-1,k2)  = &
      &                             puc(1,i2,j2-1,k2) + &
      &                 0.25_MK * ( tuc(1,i,j,k)+&
      &                             tuc(1,i+1,j,k) + &
      &                             tuc(1,i+1,j,k+1)+&
      &                             tuc(1,i,j,k+1))
                        puc(1,i2-1,j2,k2)  = &
      &                             puc(1,i2-1,j2,k2) + &
      &                 0.25_MK * ( tuc(1,i,j,k)+&
      &                             tuc(1,i,j+1,k) + &
      &                             tuc(1,i,j,k+1)+&
      &                             tuc(1,i,j+1,k+1))
                       puc(1,i2,j2,k2)  = &
      &                             puc(1,i2,j2,k2) + &
      &                 0.125_MK * (tuc(1,i,j,k)+&
      &                            tuc(1,i+1,j,k) + &
      &                            tuc(1,i+1,j+1,k)+&
      &                            tuc(1,i,j+1,k)+&
      &                            tuc(1,i,j,k+1) + &
      &                            tuc(1,i+1,j,k+1)+&
      &                            tuc(1,i,j+1,k+1)+&
      &                            tuc(1,i+1,j+1,k+1))
                        puc(2,i2-1,j2-1,k2-1) = &
      &                          puc(2,i2-1,j2-1,k2-1) + &
      &                          tuc(2,i,j,k)
                        puc(2,i2,j2-1,k2-1) = & 
      &                             puc(2,i2,j2-1,k2-1) + &
      &                   0.5_MK * (tuc(2,i,j,k)+&
      &                             tuc(2,i+1,j,k))
                        puc(2,i2-1,j2,k2-1) = &
      &                             puc(2,i2-1,j2,k2-1) + &
      &                  0.5_MK * ( tuc(2,i,j,k) + & 
      &                             tuc(2,i,j+1,k))
                        puc(2,i2,j2,k2-1)  = &
      &                            puc(2,i2,j2,k2-1) + &
      &                 0.25_MK * ( tuc(2,i,j,k)+&
      &                             tuc(2,i+1,j,k) + &
      &                             tuc(2,i+1,j+1,k)+&
      &                             tuc(2,i,j+1,k))
                        puc(2,i2-1,j2-1,k2) = &
      &                        puc(2,i2-1,j2-1,k2) + &
      &                 0.5_MK * (tuc(2,i,j,k)+&
      &                          tuc(2,i,j,k+1))
                        puc(2,i2,j2-1,k2)  = &
      &                             puc(2,i2,j2-1,k2) + &
      &                 0.25_MK * ( tuc(2,i,j,k)+&
      &                             tuc(2,i+1,j,k) + &
      &                             tuc(2,i+1,j,k+1)+&
      &                             tuc(2,i,j,k+1))
                        puc(2,i2-1,j2,k2)  = &
      &                             puc(2,i2-1,j2,k2) + &
      &                 0.25_MK * ( tuc(2,i,j,k)+&
      &                             tuc(2,i,j+1,k) + &
      &                             tuc(2,i,j,k+1)+&
      &                             tuc(2,i,j+1,k+1))
                       puc(2,i2,j2,k2)  = &
      &                             puc(2,i2,j2,k2) + &
      &                 0.125_MK * (tuc(2,i,j,k)+&
      &                            tuc(2,i+1,j,k) + &
      &                            tuc(2,i+1,j+1,k)+&
      &                            tuc(2,i,j+1,k)+&
      &                            tuc(2,i,j,k+1) + &
      &                            tuc(2,i+1,j,k+1)+&
      &                            tuc(2,i,j+1,k+1)+&
      &                            tuc(2,i+1,j+1,k+1))
                        puc(3,i2-1,j2-1,k2-1) = &
      &                          puc(3,i2-1,j2-1,k2-1) + &
      &                          tuc(3,i,j,k)
                        puc(3,i2,j2-1,k2-1) = & 
      &                             puc(3,i2,j2-1,k2-1) + &
      &                   0.5_MK * (tuc(3,i,j,k)+&
      &                             tuc(3,i+1,j,k))
                        puc(3,i2-1,j2,k2-1) = &
      &                             puc(3,i2-1,j2,k2-1) + &
      &                  0.5_MK * ( tuc(3,i,j,k) + & 
      &                             tuc(3,i,j+1,k))
                        puc(3,i2,j2,k2-1)  = &
      &                            puc(3,i2,j2,k2-1) + &
      &                 0.25_MK * ( tuc(3,i,j,k)+&
      &                             tuc(3,i+1,j,k) + &
      &                             tuc(3,i+1,j+1,k)+&
      &                             tuc(3,i,j+1,k))
                        puc(3,i2-1,j2-1,k2) = &
      &                        puc(3,i2-1,j2-1,k2) + &
      &                 0.5_MK * (tuc(3,i,j,k)+&
      &                          tuc(3,i,j,k+1))
                        puc(3,i2,j2-1,k2)  = &
      &                             puc(3,i2,j2-1,k2) + &
      &                 0.25_MK * ( tuc(3,i,j,k)+&
      &                             tuc(3,i+1,j,k) + &
      &                             tuc(3,i+1,j,k+1)+&
      &                             tuc(3,i,j,k+1))
                        puc(3,i2-1,j2,k2)  = &
      &                             puc(3,i2-1,j2,k2) + &
      &                 0.25_MK * ( tuc(3,i,j,k)+&
      &                             tuc(3,i,j+1,k) + &
      &                             tuc(3,i,j,k+1)+&
      &                             tuc(3,i,j+1,k+1))
                       puc(3,i2,j2,k2)  = &
      &                             puc(3,i2,j2,k2) + &
      &                 0.125_MK * (tuc(3,i,j,k)+&
      &                            tuc(3,i+1,j,k) + &
      &                            tuc(3,i+1,j+1,k)+&
      &                            tuc(3,i,j+1,k)+&
      &                            tuc(3,i,j,k+1) + &
      &                            tuc(3,i+1,j,k+1)+&
      &                            tuc(3,i,j+1,k+1)+&
      &                            tuc(3,i+1,j+1,k+1))
#else
                   DO ilda=1,vecdim
                        puc(ilda,i2-1,j2-1,k2-1) = &
      &                          puc(ilda,i2-1,j2-1,k2-1) + &
      &                          tuc(ilda,i,j,k)
                        puc(ilda,i2,j2-1,k2-1) = & 
      &                             puc(ilda,i2,j2-1,k2-1) + &
      &                   0.5_MK * (tuc(ilda,i,j,k)+&
      &                             tuc(ilda,i+1,j,k))
                        puc(ilda,i2-1,j2,k2-1) = &
      &                             puc(ilda,i2-1,j2,k2-1) + &
      &                  0.5_MK * ( tuc(ilda,i,j,k) + & 
      &                             tuc(ilda,i,j+1,k))
                        puc(ilda,i2,j2,k2-1)  = &
      &                            puc(ilda,i2,j2,k2-1) + &
      &                 0.25_MK * ( tuc(ilda,i,j,k)+&
      &                             tuc(ilda,i+1,j,k) + &
      &                             tuc(ilda,i+1,j+1,k)+&
      &                             tuc(ilda,i,j+1,k))
                        puc(ilda,i2-1,j2-1,k2) = &
      &                        puc(ilda,i2-1,j2-1,k2) + &
      &                 0.5_MK * (tuc(ilda,i,j,k)+&
      &                          tuc(ilda,i,j,k+1))
                        puc(ilda,i2,j2-1,k2)  = &
      &                             puc(ilda,i2,j2-1,k2) + &
      &                 0.25_MK * ( tuc(ilda,i,j,k)+&
      &                             tuc(ilda,i+1,j,k) + &
      &                             tuc(ilda,i+1,j,k+1)+&
      &                             tuc(ilda,i,j,k+1))
                        puc(ilda,i2-1,j2,k2)  = &
      &                             puc(ilda,i2-1,j2,k2) + &
      &                 0.25_MK * ( tuc(ilda,i,j,k)+&
      &                             tuc(ilda,i,j+1,k) + &
      &                             tuc(ilda,i,j,k+1)+&
      &                             tuc(ilda,i,j+1,k+1))
                       puc(ilda,i2,j2,k2)  = &
      &                             puc(ilda,i2,j2,k2) + &
      &                 0.125_MK * (tuc(ilda,i,j,k)+&
      &                            tuc(ilda,i+1,j,k) + &
      &                            tuc(ilda,i+1,j+1,k)+&
      &                            tuc(ilda,i,j+1,k)+&
      &                            tuc(ilda,i,j,k+1) + &
      &                            tuc(ilda,i+1,j,k+1)+&
      &                            tuc(ilda,i,j+1,k+1)+&
      &                            tuc(ilda,i+1,j+1,k+1))
                    ENDDO
#endif
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
#endif
#endif
         !----------------------------------------------------------------------
         !RETURN
         !----------------------------------------------------------------------
 9999    CONTINUE
         CALL substop('ppm_mg_prolong',t0,info)
         RETURN

#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_prolong_2d_sca_s
#elif __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_prolong_2d_sca_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_prolong_3d_sca_s
#elif __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_prolong_3d_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_prolong_2d_vec_s
#elif __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_prolong_2d_vec_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_prolong_3d_vec_s
#elif __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_prolong_3d_vec_d
#endif
#endif
#endif
