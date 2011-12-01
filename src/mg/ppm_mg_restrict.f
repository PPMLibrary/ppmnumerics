      !----------------------------------------------------------------------
      !  Subroutine   :            ppm_mg_restrict  
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
      SUBROUTINE ppm_mg_restrict_2d_sca_s(topo_id,mlev,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_restrict_2d_sca_d(topo_id,mlev,info)
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_restrict_3d_sca_s(topo_id,mlev,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_restrict_3d_sca_d(topo_id,mlev,info)
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_restrict_2d_vec_s(topo_id,mlev,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_restrict_2d_vec_d(topo_id,mlev,info)
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_restrict_3d_vec_s(topo_id,mlev,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_restrict_3d_vec_d(topo_id,mlev,info)
#endif
#endif
#endif
      !!! In this routine we restrict the error from finer to coarser levels
      !-----------------------------------------------------------------------
      !  Includes
      !-----------------------------------------------------------------------
#include "ppm_define.h"
      !-----------------------------------------------------------------------
      !  Modules 
      !-----------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mg
      USE ppm_module_alloc
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      USE ppm_module_map
      USE ppm_module_map_field
      USE ppm_module_map_field_ghost
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-----------------------------------------------------------------------
      !  Arguments     
      !-----------------------------------------------------------------------
      INTEGER,                   INTENT(IN)      ::  mlev, topo_id
      INTEGER,                   INTENT(INOUT)   ::  info
      !-----------------------------------------------------------------------
      !  Local variables 
      !-----------------------------------------------------------------------
      CHARACTER(LEN=256)                         :: cbuf
      INTEGER                                    :: isub,j,j2,i,i2
      INTEGER                                    :: mlevm1,ilda,iface
      INTEGER,DIMENSION(5)                       :: ldl5,ldu5
      INTEGER,DIMENSION(4)                       :: ldl4,ldu4
      INTEGER,DIMENSION(3)                       :: ldl3,ldu3
      INTEGER                                    :: iopt,topoid
      INTEGER                                    :: a,b,c,d,e,f,g  
#if __MESH_DIM == __3D
      INTEGER                                    :: k,k2
#endif        
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
      REAL(MK),DIMENSION(:,:,:),POINTER :: uc_dummy
#elif __MESH_DIM == __3D
      REAL(MK),DIMENSION(:,:,:,:),POINTER :: uc_dummy
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
      REAL(MK),DIMENSION(:,:,:,:),POINTER :: uc_dummy
#elif __MESH_DIM == __3D
      REAL(MK),DIMENSION(:,:,:,:,:),POINTER :: uc_dummy
#endif
#endif
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
      REAL(MK),DIMENSION(:,:),POINTER :: terr
      REAL(MK),DIMENSION(:,:),POINTER :: pfc
#elif __MESH_DIM == __3D
      REAL(MK),DIMENSION(:,:,:),POINTER :: terr
      REAL(MK),DIMENSION(:,:,:),POINTER :: pfc
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
      REAL(MK),DIMENSION(:,:,:),POINTER :: terr
      REAL(MK),DIMENSION(:,:,:),POINTER :: pfc
#elif __MESH_DIM == __3D
      REAL(MK),DIMENSION(:,:,:,:),POINTER :: terr
      REAL(MK),DIMENSION(:,:,:,:),POINTER :: pfc
#endif
#endif
      !----------------------------------------------------------------------
      !Externals
      !----------------------------------------------------------------------

      !----------------------------------------------------------------------
      !Initialize
      !----------------------------------------------------------------------
      CALL substart('ppm_mg_restrict',t0,info)
      !----------------------------------------------------------------------
      !  Check arguments
      !----------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (mlev.LE.1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_restrict',  &
        &                'level must be >1',__LINE__,info)
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
      mlevm1=mlev-1
      IF (ppm_debug.GT.0) THEN
          WRITE(cbuf,*) 'WELCOME TO THE RESTRICTION LEVEL:',mlev
          CALL PPM_WRITE(ppm_rank,'mg_restrict',cbuf,info)
      ENDIF 
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
      !----------------------------------------------------------------------
      ! Restriction using a 9-point operator(bilinear interpolation)
      ! linear is not accurate enough
      !----------------------------------------------------------------------
      topoid=topo_id
      iopt = ppm_param_alloc_fit
      ldl3(1) = 1-ghostsize(1)
      ldl3(2) = 1-ghostsize(2)
      ldl3(3) = 1
      ldu3(1) = max_node(1,mlevm1)+ghostsize(1)
      ldu3(2) = max_node(2,mlevm1)+ghostsize(2)
      ldu3(3) = nsubs
      CALL ppm_alloc(uc_dummy,ldl3,ldu3,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'restrict',    &
          &                       'uc_dummy',__LINE__,info)
          GOTO 9999
      ENDIF
      DO isub=1,nsubs
          terr=>mgfield(isub,mlevm1)%err 
          DO i=1-ghostsize(1),max_node(1,mlevm1)+ghostsize(1)
              DO j=1-ghostsize(2),max_node(2,mlevm1)+ghostsize(2)
                  uc_dummy(i,j,isub) = terr(i,j)
              ENDDO
          ENDDO   
      ENDDO 
      CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlevm1),&
      &                         ghostsize,info)
      CALL ppm_map_field_push(topoid,mg_meshid(mlevm1),uc_dummy,&
      &                         info)
      CALL ppm_map_field_send(info)
      CALL ppm_map_field_pop(topoid,mg_meshid(mlevm1),uc_dummy,&
      &                          ghostsize,info)
      DO isub=1,nsubs
          terr=>mgfield(isub,mlevm1)%err 
          pfc=>mgfield(isub,mlev)%fc 
          DO i=1-ghostsize(1),max_node(1,mlevm1)+ghostsize(1)
              DO j=1-ghostsize(2),max_node(2,mlevm1)+ghostsize(2)
                  terr(i,j) = uc_dummy(i,j,isub)
              ENDDO
          ENDDO   
          DO j=start(2,isub,mlev),istop(2,isub,mlev)
              j2=2*j 
              DO i=start(1,isub,mlev),istop(1,isub,mlev)
                  i2=2*i
                       pfc(i,j)= &
                       &         0.25_MK * terr(i2-1,j2-1) + &
                       &         0.125_MK * (terr(i2,j2-1) + &
                       &                   terr(i2-2,j2-1) + &
                       &                     terr(i2-1,j2) + &
                       &                  terr(i2-1,j2-2)) + &
                       &        0.0625_MK * (terr(i2,j2-2) + &
                       &                    terr(i2-2,j2) +  &
                       &                     terr(i2-2,j2-2) &
                       &                      + terr(i2,j2)) 
                   ENDDO
               ENDDO
           ENDDO
           iopt = ppm_param_dealloc
           ldl3(1) = 1-ghostsize(1)
           ldl3(2) = 1-ghostsize(2)
           ldl3(3) = 1
           ldu3(1) = max_node(1,mlevm1)+ghostsize(1)
           ldu3(2) = max_node(2,mlevm1)+ghostsize(2)
           ldu3(3) = nsubs
           CALL ppm_alloc(uc_dummy,ldl3,ldu3,iopt,info)
           IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'restrict',    &
              &             'uc_dummy',__LINE__,info)
              GOTO 9999
           ENDIF
#elif __MESH_DIM == __3D
           topoid=topo_id
           iopt = ppm_param_alloc_fit
           ldl4(1) = 1-ghostsize(1)
           ldl4(2) = 1-ghostsize(2)
           ldl4(3) = 1-ghostsize(3)
           ldl4(4) = 1
           ldu4(1) = max_node(1,mlevm1)+ghostsize(1)
           ldu4(2) = max_node(2,mlevm1)+ghostsize(2)
           ldu4(3) = max_node(3,mlevm1)+ghostsize(3)
           ldu4(4) = nsubs
           CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
           IF (info .NE. 0) THEN
               info = ppm_error_fatal
               CALL ppm_error(ppm_err_alloc,'restrict',    &
               &                       'uc_dummy',__LINE__,info)
               GOTO 9999
           ENDIF
           DO isub=1,nsubs
               terr=>mgfield(isub,mlevm1)%err 
               DO i=1-ghostsize(1),max_node(1,mlevm1)+ghostsize(1)
                   DO j=1-ghostsize(2),max_node(2,mlevm1)+ghostsize(2)
                       DO k=1-ghostsize(3),max_node(3,mlevm1)+ghostsize(3)
                           uc_dummy(i,j,k,isub) = terr(i,j,k)
                       ENDDO
                   ENDDO   
               ENDDO
           ENDDO 
           CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlevm1),&
           &                            ghostsize,info)
           CALL ppm_map_field_push(topoid,mg_meshid(mlevm1),uc_dummy,&
           &                       info)
           CALL ppm_map_field_send(info)
           CALL ppm_map_field_pop(topoid,mg_meshid(mlevm1),uc_dummy,&
           &                      ghostsize,info)

           DO isub=1,nsubs
               terr=>mgfield(isub,mlevm1)%err 
               pfc=>mgfield(isub,mlev)%fc 
               DO i=1-ghostsize(1),max_node(1,mlevm1)+ghostsize(1)
                   DO j=1-ghostsize(2),max_node(2,mlevm1)+ghostsize(2)
                       DO k=1-ghostsize(3),max_node(3,mlevm1)+ghostsize(3)
                           terr(i,j,k) = uc_dummy(i,j,k,isub)
                       ENDDO
                   ENDDO
               ENDDO
               DO k=start(3,isub,mlev),istop(3,isub,mlev)
                   k2=2*k
                   DO j=start(2,isub,mlev),istop(2,isub,mlev)
                       j2=2*j 
                       DO i=start(1,isub,mlev),istop(1,isub,mlev)
                           i2=2*i
                           pfc(i,j,k) = &
                           & 0.125_MK * terr(i2-1,j2-1,k2-1) + &
                           & 0.0625_MK * (terr(i2,j2-1,k2-1) + &
                           &             terr(i2-2,j2-1,k2-1)+ &
                           &              terr(i2-1,j2,k2-1) + &
                           &             terr(i2-1,j2-2,k2-1))+&
                           & 0.03125_MK * (&
                           &               terr(i2,j2-2,k2-1)+ &
                           &             terr(i2-2,j2,k2-1) +  &
                           &             terr(i2-2,j2-2,k2-1) +&
                           &                terr(i2,j2,k2-1)) 

                          pfc(i,j,k)=              pfc(i,j,k) +&
                          & 0.0625_MK * terr(i2-1,j2-1,k2-2) + &
                          & 0.03125_MK *(terr(i2,j2-1,k2-2)  + &
                          &              terr(i2-2,j2-1,k2-2)+ &
                          &                 terr(i2-1,j2,k2-2)+&
                          &              terr(i2-1,j2-2,k2-2))+&
                          & 0.015625_MK*(&
                          &                 terr(i2,j2-2,k2-2)+&
                          &              terr(i2-2,j2,k2-2) +  &
                          &               terr(i2-2,j2-2,k2-2)+&
                          &                  terr(i2,j2,k2-2)) 
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
          iopt = ppm_param_dealloc
          ldl4(1) = 1-ghostsize(1)
          ldl4(2) = 1-ghostsize(2)
          ldl4(3) = 1-ghostsize(3)
          ldl4(4) = 1
          ldu4(1) = max_node(1,mlevm1)+ghostsize(1)
          ldu4(2) = max_node(2,mlevm1)+ghostsize(2)
          ldu4(3) = max_node(3,mlevm1)+ghostsize(3)
          ldu4(4) = nsubs
          CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'restrict',    &
              &             'uc_dummy',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
          !----------------------------------------------------------------------
          ! Restriction using a 9-point operator(bilinear interpolation)
          ! linear is not accurate enough
          !---------------------------------------------------------------------
          topoid=topo_id
          iopt = ppm_param_alloc_fit
          ldl4(1) = 1
          ldl4(2) = 1-ghostsize(1)
          ldl4(3) = 1-ghostsize(2)
          ldl4(4) = 1
          ldu4(1) = vecdim
          ldu4(2) = max_node(1,mlevm1)+ghostsize(1)
          ldu4(3) = max_node(2,mlevm1)+ghostsize(2)
          ldu4(4) = nsubs
          CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'restrict',    &
              &             'uc_dummy',__LINE__,info)
              GOTO 9999
          ENDIF
          DO isub=1,nsubs
              terr=>mgfield(isub,mlevm1)%err
              DO i=1-ghostsize(1),max_node(1,mlevm1)+ghostsize(1)
                  DO j=1-ghostsize(2),max_node(2,mlevm1)+ghostsize(2)
                      DO ilda=1,vecdim
                          uc_dummy(ilda,i,j,isub) = terr(ilda,i,j)
                      ENDDO
                  ENDDO   
              ENDDO
          ENDDO 

          CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlevm1),&
          &                         ghostsize,info)
          CALL ppm_map_field_push(topoid,mg_meshid(mlevm1),uc_dummy,&
          &                         vecdim,info)
          CALL ppm_map_field_send(info)
          CALL ppm_map_field_pop(topoid,mg_meshid(mlevm1),uc_dummy,&
          &                          vecdim,ghostsize,info)
          DO isub=1,nsubs
              terr=>mgfield(isub,mlevm1)%err
              pfc=>mgfield(isub,mlev)%fc
              DO i=1-ghostsize(1),max_node(1,mlevm1)+ghostsize(1)
                  DO j=1-ghostsize(2),max_node(2,mlevm1)+ghostsize(2)
                      DO ilda=1,vecdim
                          terr(ilda,i,j) = uc_dummy(ilda,i,j,isub)
                      ENDDO
                  ENDDO   
              ENDDO
              DO j=start(2,isub,mlev),istop(2,isub,mlev)
                  j2=2*j 
                  DO i=start(1,isub,mlev),istop(1,isub,mlev)
                      i2=2*i
                      DO ilda=1,vecdim
                          pfc(ilda,i,j)= &
                          & 0.25_MK * terr(ilda,i2-1,j2-1) + &
                          & 0.125_MK * (terr(ilda,i2,j2-1) + &
                          &            terr(ilda,i2-2,j2-1)+ &
                          &             terr(ilda,i2-1,j2) + &
                          &            terr(ilda,i2-1,j2-2))+&
                          &  0.0625_MK * (terr(ilda,i2,j2-2)+&
                          &            terr(ilda,i2-2,j2) +  &
                          &              terr(ilda,i2-2,j2-2)&
                          &               + terr(ilda,i2,j2)) 
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
          iopt = ppm_param_dealloc
          ldl4(1) = 1
          ldl4(2) = 1-ghostsize(1)
          ldl4(3) = 1-ghostsize(2)
          ldl4(4) = 1
          ldu4(1) = vecdim
          ldu4(2) = max_node(1,mlevm1)+ghostsize(1)
          ldu4(3) = max_node(2,mlevm1)+ghostsize(2)
          ldu4(4) = nsubs
          CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'restrict',    &
              &                       'uc_dummy',__LINE__,info)
              GOTO 9999
          ENDIF
#elif __MESH_DIM == __3D
          topoid=topo_id
          iopt = ppm_param_alloc_fit
          ldl5(1) = 1
          ldl5(2) = 1-ghostsize(1)
          ldl5(3) = 1-ghostsize(2)
          ldl5(4) = 1-ghostsize(3)
          ldl5(5) = 1
          ldu5(1) = vecdim
          ldu5(2) = max_node(1,mlevm1)+ghostsize(1)
          ldu5(3) = max_node(2,mlevm1)+ghostsize(2)
          ldu5(4) = max_node(3,mlevm1)+ghostsize(3)
          ldu5(5) = nsubs
          CALL ppm_alloc(uc_dummy,ldl5,ldu5,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'restrict',    &
              &                       'uc_dummy',__LINE__,info)
              GOTO 9999
          ENDIF
          DO isub=1,nsubs
              terr=>mgfield(isub,mlevm1)%err
              DO k=1-ghostsize(3),max_node(3,mlevm1)+ghostsize(3)
                  DO j=1-ghostsize(2),max_node(2,mlevm1)+ghostsize(2)
                      DO i=1-ghostsize(1),max_node(1,mlevm1)+ghostsize(1)
#ifdef __VECTOR
                          uc_dummy(1,i,j,k,isub) = terr(1,i,j,k)
                          uc_dummy(2,i,j,k,isub) = terr(2,i,j,k)
                          uc_dummy(3,i,j,k,isub) = terr(3,i,j,k)
#else
                          DO ilda=1,vecdim
                              uc_dummy(ilda,i,j,k,isub) = terr(ilda,i,j,k)
#endif
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO  
      ENDDO 
      CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlevm1),&
      &                         ghostsize,info)
      CALL ppm_map_field_push(topoid,mg_meshid(mlevm1),uc_dummy,&
      &                         vecdim,info)
      CALL ppm_map_field_send(info)
      CALL ppm_map_field_pop(topoid,mg_meshid(mlevm1),uc_dummy,&
      &                          vecdim,ghostsize,info)
      DO isub=1,nsubs
          terr=>mgfield(isub,mlevm1)%err
          pfc=>mgfield(isub,mlev)%fc

          DO k=1-ghostsize(3),max_node(3,mlevm1)+ghostsize(3)
              DO j=1-ghostsize(2),max_node(2,mlevm1)+ghostsize(2)
                  DO i=1-ghostsize(1),max_node(1,mlevm1)+ghostsize(1)
#ifdef __VECTOR
                      terr(1,i,j,k) = uc_dummy(1,i,j,k,isub)
                      terr(2,i,j,k) = uc_dummy(2,i,j,k,isub)
                      terr(3,i,j,k) = uc_dummy(3,i,j,k,isub)
#else
                      DO ilda=1,vecdim
                          terr(ilda,i,j,k) = uc_dummy(ilda,i,j,k,isub)
#endif       
                     ENDDO
                 ENDDO
             ENDDO
         ENDDO
         DO k=start(3,isub,mlev),istop(3,isub,mlev)
             k2=2*k
             DO j=start(2,isub,mlev),istop(2,isub,mlev)
                 j2=2*j 
                 DO i=start(1,isub,mlev),istop(1,isub,mlev)
                     i2=2*i
#ifdef __VECTOR
                     pfc(1,i,j,k)= &
                     &                  0.125_MK * &
                     &    terr(1,i2-1,j2-1,k2-1) + &
                     &                0.0625_MK * (&
                     &       terr(1,i2,j2-1,k2-1) +&
                     &     terr(1,i2-2,j2-1,k2-1)+ &
                     &      terr(1,i2-1,j2,k2-1) + &
                     &     terr(1,i2-1,j2-2,k2-1))+&
                     &               0.03125_MK * (&
                     &       terr(1,i2,j2-2,k2-1)+ &
                     &     terr(1,i2-2,j2,k2-1) +  &
                     &     terr(1,i2-2,j2-2,k2-1) +&
                     &        terr(1,i2,j2,k2-1)) 
                     pfc(1,i,j,k)= &
                     &                pfc(1,i,j,k)+&
                     &                 0.0625_MK * &
                     &      terr(1,i2-1,j2-1,k2) + &
                     &               0.03125_MK * (&
                     &         terr(1,i2,j2-1,k2) +&
                     &       terr(1,i2-2,j2-1,k2)+ &
                     &         terr(1,i2-1,j2,k2) +&
                     &       terr(1,i2-1,j2-2,k2))+&
                     &              0.015625_MK * (&
                     &         terr(1,i2,j2-2,k2)+ &
                     &       terr(1,i2-2,j2,k2) +  &
                     &      terr(1,i2-2,j2-2,k2) + &
                     &         terr(1,i2,j2,k2)) 
                     pfc(1,i,j,k)= &
                     &               pfc(1,i,j,k) +&
                     &                 0.0625_MK * &
                     &    terr(1,i2-1,j2-1,k2-2) + &
                     &                0.03125_MK *(&
                     &       terr(1,i2,j2-1,k2-2) +&
                     &     terr(1,i2-2,j2-1,k2-2)+ &
                     &        terr(1,i2-1,j2,k2-2)+&
                     &     terr(1,i2-1,j2-2,k2-2))+&
                     &                0.015625_MK*(&
                     &        terr(1,i2,j2-2,k2-2)+&
                     &     terr(1,i2-2,j2,k2-2) +  &
                     &      terr(1,i2-2,j2-2,k2-2)+&
                     &         terr(1,i2,j2,k2-2)) 
                     pfc(2,i,j,k)= &
                     &                  0.125_MK * &
                     &    terr(2,i2-1,j2-1,k2-1) + &
                     &                0.0625_MK * (&
                     &       terr(2,i2,j2-1,k2-1) +&
                     &     terr(2,i2-2,j2-1,k2-1)+ &
                     &      terr(2,i2-1,j2,k2-1) + &
                     &     terr(2,i2-1,j2-2,k2-1))+&
                     &               0.03125_MK * (&
                     &       terr(2,i2,j2-2,k2-1)+ &
                     &     terr(2,i2-2,j2,k2-1) +  &
                     &     terr(2,i2-2,j2-2,k2-1) +&
                     &        terr(2,i2,j2,k2-1)) 
                     pfc(2,i,j,k)= &
                     &                pfc(2,i,j,k)+&
                     &                 0.0625_MK * &
                     &      terr(2,i2-1,j2-1,k2) + &
                     &               0.03125_MK * (&
                     &         terr(2,i2,j2-1,k2) +&
                     &       terr(2,i2-2,j2-1,k2)+ &
                     &         terr(2,i2-1,j2,k2) +&
                     &       terr(2,i2-1,j2-2,k2))+&
                     &              0.015625_MK * (&
                     &         terr(2,i2,j2-2,k2)+ &
                     &       terr(2,i2-2,j2,k2) +  &
                     &      terr(2,i2-2,j2-2,k2) + &
                     &         terr(2,i2,j2,k2)) 
                     pfc(2,i,j,k)= &
                     &               pfc(2,i,j,k) +&
                     &                 0.0625_MK * &
                     &    terr(2,i2-1,j2-1,k2-2) + &
                     &                0.03125_MK *(&
                     &       terr(2,i2,j2-1,k2-2) +&
                     &     terr(2,i2-2,j2-1,k2-2)+ &
                     &        terr(2,i2-1,j2,k2-2)+&
                     &     terr(2,i2-1,j2-2,k2-2))+&
                     &                0.015625_MK*(&
                     &        terr(2,i2,j2-2,k2-2)+&
                     &     terr(2,i2-2,j2,k2-2) +  &
                     &      terr(2,i2-2,j2-2,k2-2)+&
                     &         terr(2,i2,j2,k2-2)) 
                     pfc(3,i,j,k)= &
                     &                  0.125_MK * &
                     &    terr(3,i2-1,j2-1,k2-1) + &
                     &                0.0625_MK * (&
                     &       terr(3,i2,j2-1,k2-1) +&
                     &     terr(3,i2-2,j2-1,k2-1)+ &
                     &      terr(3,i2-1,j2,k2-1) + &
                     &     terr(3,i2-1,j2-2,k2-1))+&
                     &               0.03125_MK * (&
                     &       terr(3,i2,j2-2,k2-1)+ &
                     &     terr(3,i2-2,j2,k2-1) +  &
                     &     terr(3,i2-2,j2-2,k2-1) +&
                     &        terr(3,i2,j2,k2-1)) 
                     pfc(3,i,j,k)= &
                     &                pfc(3,i,j,k)+&
                     &                 0.0625_MK * &
                     &      terr(3,i2-1,j2-1,k2) + &
                     &               0.03125_MK * (&
                     &         terr(3,i2,j2-1,k2) +&
                     &       terr(3,i2-2,j2-1,k2)+ &
                     &         terr(3,i2-1,j2,k2) +&
                     &       terr(3,i2-1,j2-2,k2))+&
                     &              0.015625_MK * (&
                     &         terr(3,i2,j2-2,k2)+ &
                     &       terr(3,i2-2,j2,k2) +  &
                     &      terr(3,i2-2,j2-2,k2) + &
                     &         terr(3,i2,j2,k2)) 
                     pfc(3,i,j,k)= &
                     &               pfc(3,i,j,k) +&
                     &                 0.0625_MK * &
                     &    terr(3,i2-1,j2-1,k2-2) + &
                     &                0.03125_MK *(&
                     &       terr(3,i2,j2-1,k2-2) +&
                     &     terr(3,i2-2,j2-1,k2-2)+ &
                     &        terr(3,i2-1,j2,k2-2)+&
                     &     terr(3,i2-1,j2-2,k2-2))+&
                     &                0.015625_MK*(&
                     &        terr(3,i2,j2-2,k2-2)+&
                     &     terr(3,i2-2,j2,k2-2) +  &
                     &      terr(3,i2-2,j2-2,k2-2)+&
                     &         terr(3,i2,j2,k2-2)) 
#else
                     DO ilda=1,vecdim
                         pfc(ilda,i,j,k)= &
                         &                   0.125_MK * &
                         &  terr(ilda,i2-1,j2-1,k2-1) + &
                         &                 0.0625_MK * (&
                         &     terr(ilda,i2,j2-1,k2-1) +&
                         &   terr(ilda,i2-2,j2-1,k2-1)+ &
                         &    terr(ilda,i2-1,j2,k2-1) + &
                         &   terr(ilda,i2-1,j2-2,k2-1))+&
                         &                0.03125_MK * (&
                         &     terr(ilda,i2,j2-2,k2-1)+ &
                         &   terr(ilda,i2-2,j2,k2-1) +  &
                         &   terr(ilda,i2-2,j2-2,k2-1) +&
                         &      terr(ilda,i2,j2,k2-1)) 
                         pfc(ilda,i,j,k)= &
                         &              pfc(ilda,i,j,k)+&
                         &                  0.0625_MK * &
                         &    terr(ilda,i2-1,j2-1,k2) + &
                         &                0.03125_MK * (&
                         &       terr(ilda,i2,j2-1,k2) +&
                         &     terr(ilda,i2-2,j2-1,k2)+ &
                         &       terr(ilda,i2-1,j2,k2) +&
                         &     terr(ilda,i2-1,j2-2,k2))+&
                         &               0.015625_MK * (&
                         &       terr(ilda,i2,j2-2,k2)+ &
                         &     terr(ilda,i2-2,j2,k2) +  &
                         &    terr(ilda,i2-2,j2-2,k2) + &
                         &       terr(ilda,i2,j2,k2)) 
                         pfc(ilda,i,j,k)= &
                         &             pfc(ilda,i,j,k) +&
                         &                  0.0625_MK * &
                         &  terr(ilda,i2-1,j2-1,k2-2) + &
                         &                 0.03125_MK *(&
                         &     terr(ilda,i2,j2-1,k2-2) +&
                         &   terr(ilda,i2-2,j2-1,k2-2)+ &
                         &      terr(ilda,i2-1,j2,k2-2)+&
                         &   terr(ilda,i2-1,j2-2,k2-2))+&
                         &                 0.015625_MK*(&
                         &      terr(ilda,i2,j2-2,k2-2)+&
                         &   terr(ilda,i2-2,j2,k2-2) +  &
                         &    terr(ilda,i2-2,j2-2,k2-2)+&
                         &       terr(ilda,i2,j2,k2-2)) 
                     ENDDO
#endif
                 ENDDO
             ENDDO
         ENDDO
     ENDDO
     iopt = ppm_param_dealloc
     ldl5(1) = 1
     ldl5(2) = 1-ghostsize(1)
     ldl5(3) = 1-ghostsize(2)
     ldl5(4) = 1-ghostsize(3)
     ldl5(5) = 1
     ldu5(1) = vecdim
     ldu5(2) = max_node(1,mlevm1)+ghostsize(1)
     ldu5(3) = max_node(2,mlevm1)+ghostsize(2)
     ldu5(4) = max_node(3,mlevm1)+ghostsize(3)
     ldu5(5) = nsubs
     CALL ppm_alloc(uc_dummy,ldl5,ldu5,iopt,info)
     IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'restrict',    &
         &                       'uc_dummy',__LINE__,info)
         GOTO 9999
     ENDIF
#endif
#endif
        !----------------------------------------------------------------------
        ! Return
        !----------------------------------------------------------------------
   9999    CONTINUE
           CALL substop('ppm_mg_restrict',t0,info)
           RETURN

#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_restrict_2d_sca_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_restrict_2d_sca_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_restrict_3d_sca_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_restrict_3d_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_restrict_2d_vec_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_restrict_2d_vec_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_restrict_3d_vec_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_restrict_3d_vec_d
#endif
#endif
#endif
