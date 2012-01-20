
      !-------------------------------------------------------------------------
      !  Subroutine   :            ppm_mg_smooth_coarse    
      !-------------------------------------------------------------------------
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
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_smooth_coarse_2D_sca_s(topo_id,nsweep,mlev,&
     &                                         c1,c2,c3,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_smooth_coarse_2D_sca_d(topo_id,nsweep,mlev,&
     &                                         c1,c2,c3,info)
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_smooth_coarse_3D_sca_s(topo_id,nsweep,mlev,&
     &                                         c1,c2,c3,c4,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_smooth_coarse_3D_sca_d(topo_id,nsweep,mlev,&
     &                                         c1,c2,c3,c4,info)
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_smooth_coarse_2D_vec_s(topo_id,nsweep,mlev,&
     &                                         c1,c2,c3,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_smooth_coarse_2D_vec_d(topo_id,nsweep,mlev,&
     &                                         c1,c2,c3,info)
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_smooth_coarse_3D_vec_s(topo_id,nsweep,mlev,&
     &                                         c1,c2,c3,c4,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_smooth_coarse_3D_vec_d(topo_id,nsweep,mlev,&
     &                                         c1,c2,c3,c4,info)
#endif
#endif
#endif
         !!! In this routine we compute the corrections for the function 
         !!! based on the Gauss-Seidel iteration
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
         USE ppm_module_map
         USE ppm_module_data_mesh
         USE ppm_module_write
         USE ppm_module_map_field
         USE ppm_module_map_field_ghost
         IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION
         INTEGER, PARAMETER :: MK = ppm_kind_single
#else
         INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
         !------------------------------------------------------------------
         !  Arguments     
         !----------------------------------------------------------------------
         INTEGER,                   INTENT(IN)      ::  nsweep
         INTEGER,                   INTENT(IN)      ::  mlev, topo_id
#if  __MESH_DIM == __2D
         REAL(MK),                  INTENT(IN)      ::  c1,c2,c3 
#elif __MESH_DIM == __3D
         REAL(MK),                  INTENT(IN)      ::  c1,c2,c3,c4 
#endif
         INTEGER,                   INTENT(INOUT)   ::  info
         !--------------------------------------------------------------------
         !  Local variables 
         !----------------------------------------------------------------------
         CHARACTER(LEN=256) :: cbuf
         INTEGER                                    ::  i,j,isub,color,colos
         INTEGER,DIMENSION(:,:),POINTER             ::  lorig => NULL()
         INTEGER,DIMENSION(:,:),POINTER             ::  lext => NULL()
         INTEGER,DIMENSION(:),POINTER               ::  a => NULL()
         INTEGER,DIMENSION(:),POINTER               ::  b => NULL()
         INTEGER,DIMENSION(:),POINTER               ::  c => NULL()
         INTEGER,DIMENSION(:),POINTER               ::  d => NULL()
         INTEGER,DIMENSION(:),POINTER               ::  e => NULL()
         INTEGER,DIMENSION(:),POINTER               ::  g => NULL()
         REAL(MK)                                   ::  c11,c22,c33,c44 
         INTEGER                                    ::  ilda,isweep,count
         INTEGER                                    ::  k,idom
         REAL(MK)                                   ::  x,y,dx,dy
         REAL(MK)                                   ::  omega
         INTEGER,DIMENSION(1)                       ::  ldu1,ldl1
         INTEGER,DIMENSION(2)                       ::  ldu2
#if __MESH_DIM == __2D
         INTEGER,DIMENSION(4)                       ::  ldl4,ldu4
         INTEGER,DIMENSION(3)                       ::  ldl3,ldu3
#endif
#if __MESH_DIM == __3D
         INTEGER,DIMENSION(5)                       ::  ldl5,ldu5
         INTEGER,DIMENSION(4)                       ::  ldl4,ldu4
         REAL(MK)                                   ::  dz
#endif
         INTEGER                                    ::  iopt,iface,topoid
         REAL(MK)                                   ::  t0
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
         TYPE(mg_field_2d_sca_s),DIMENSION(:,:),POINTER :: mgfield => NULL()
#elif __KIND == __DOUBLE_PRECISION
         TYPE(mg_field_2d_sca_d),DIMENSION(:,:),POINTER :: mgfield => NULL()
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
        TYPE(mg_field_3d_sca_s),DIMENSION(:,:),POINTER :: mgfield => NULL()
#elif __KIND == __DOUBLE_PRECISION
         TYPE(mg_field_3d_sca_d),DIMENSION(:,:),POINTER :: mgfield => NULL()
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
         TYPE(mg_field_2d_vec_s),DIMENSION(:,:),POINTER :: mgfield => NULL()
#elif __KIND == __DOUBLE_PRECISION
         TYPE(mg_field_2d_vec_d),DIMENSION(:,:),POINTER :: mgfield => NULL()
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
         TYPE(mg_field_3d_vec_s),DIMENSION(:,:),POINTER :: mgfield => NULL()
#elif __KIND == __DOUBLE_PRECISION
         TYPE(mg_field_3d_vec_d),DIMENSION(:,:),POINTER :: mgfield => NULL()
#endif
#endif
#endif
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
         REAL(MK),DIMENSION(:,:,:),POINTER :: uc   => NULL()
#elif __MESH_DIM == __3D
         REAL(MK),DIMENSION(:,:,:,:),POINTER :: uc   => NULL()
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
         REAL(MK),DIMENSION(:,:,:,:),POINTER :: uc   => NULL()
#elif __MESH_DIM == __3D
         REAL(MK),DIMENSION(:,:,:,:,:),POINTER :: uc   => NULL()
#endif
#endif
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
         REAL(MK),DIMENSION(:,:,:),POINTER :: oldu => NULL()
#elif __MESH_DIM == __3D
         REAL(MK),DIMENSION(:,:,:,:),POINTER :: oldu => NULL()
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
         REAL(MK),DIMENSION(:,:,:,:),POINTER :: oldu   => NULL()
#elif __MESH_DIM == __3D
         REAL(MK),DIMENSION(:,:,:,:,:),POINTER :: oldu => NULL()
#endif
#endif
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
      REAL(MK) :: moldu
#elif __MESH_DIM == __3D
      REAL(MK) :: moldu
#endif
#elif  __DIM == __VFIELD
#if __MESH_DIM == __2D
      REAL(MK),DIMENSION(:),POINTER :: moldu => NULL()
#elif __MESH_DIM == __3D
      REAL(MK),DIMENSION(:),POINTER :: moldu => NULL()
#endif
#endif
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
         REAL(MK),DIMENSION(:,:),POINTER :: tuc => NULL()
#elif __MESH_DIM == __3D
        REAL(MK),DIMENSION(:,:,:),POINTER :: tuc => NULL()
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
       REAL(MK),DIMENSION(:,:,:),POINTER :: tuc => NULL()
#elif __MESH_DIM == __3D
       REAL(MK),DIMENSION(:,:,:,:),POINTER :: tuc => NULL()
#endif
#endif
#if __KIND == __SINGLE_PRECISION
       omega=omega_s
       dx=dx_s
       dy=dy_s
#if __MESH_DIM == __3D
       dz=dz_s
#endif
#elif __KIND == __DOUBLE_PRECISION
       omega=omega_d
       dx=dx_d
       dy=dy_d
#if __MESH_DIM == __3D
       dz=dz_d
#endif
#endif
         !----------------------------------------------------------------------
         !Externals
         !----------------------------------------------------------------------

         !----------------------------------------------------------------------
         !Initialize
         !----------------------------------------------------------------------
         CALL substart('ppm_mg_smooth_coarse',t0,info)
         IF (ppm_debug.GT.0) THEN 
          WRITE (cbuf,*) 'SMOOTHER entering ','mlev:',mlev
          CALL PPM_WRITE(ppm_rank,'mg_smooth',cbuf,info)
         ENDIF
         !----------------------------------------------------------------------
         !  Check arguments
         !----------------------------------------------------------------------
         IF (ppm_debug .GT. 0) THEN
           IF (nsweep.LT.1) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_coarse',  &
      &            'nsweep must be >=1',__LINE__,info)
               GOTO 9999
           ENDIF
           IF (mlev.LE.1) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_coarse',  &
      &            'level must be >1',__LINE__,info)
               GOTO 9999
           ENDIF
           IF (c1.LE.0.0_MK) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_coarse',  &
      &            'Factor c1 must be >0',__LINE__,info)
               GOTO 9999
           ENDIF
           IF (c2.LE.0.0_MK) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_coarse',  &
      &            'Factor c2 must be >0',__LINE__,info)
               GOTO 9999
           ENDIF
           IF (c3.LE.0.0_MK) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_coarse',  &
      &            'Factor c3 must be >0',__LINE__,info)
               GOTO 9999
           ENDIF
#if __MESH_DIM == __3D
           IF (c4.LE.0.0_MK) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_coarse',  &
      &            'Factor c4 must be >0',__LINE__,info)
               GOTO 9999
          ENDIF
#endif
         ENDIF
         !----------------------------------------------------------------------
         !Definition of necessary variables and allocation of arrays
         !----------------------------------------------------------------------
         topoid=topo_id
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
            iopt = ppm_param_alloc_fit
            ldu2(1) = ppm_dim
            ldu2(2) = nsubs
            CALL ppm_alloc(lorig,ldu2,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_fine',    &
      &                       'origi',__LINE__,info)
            GOTO 9999
            ENDIF
            CALL ppm_alloc(lext,ldu2,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_fine',    &
      &                       'origi',__LINE__,info)
            GOTO 9999
            ENDIF
             iopt = ppm_param_alloc_fit
             ldl1(1) = 1
             ldu1(1) = nsubs
             CALL ppm_alloc(a,ldl1,ldu1,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_coarse',    &
       &                       'a',__LINE__,info)
             GOTO 9999
             ENDIF
             CALL ppm_alloc(b,ldl1,ldu1,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_coarse',    &
       &                       'b',__LINE__,info)
             GOTO 9999
             ENDIF
             CALL ppm_alloc(c,ldl1,ldu1,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_coarse',    &
       &                       'c',__LINE__,info)
             GOTO 9999
             ENDIF
             CALL ppm_alloc(d,ldl1,ldu1,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_coarse',    &
       &                       'd',__LINE__,info)
             GOTO 9999
             ENDIF
             CALL ppm_alloc(e,ldl1,ldu1,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_coarse',    &
       &                       'e',__LINE__,info)
             GOTO 9999
             ENDIF
             CALL ppm_alloc(g,ldl1,ldu1,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_coarse',    &
       &                       'g',__LINE__,info)
             GOTO 9999
             ENDIF
            lorig=0
            lext=0
#if  __DIM == __SFIELD
#if  __MESH_DIM == __2D
         !----------------------------------------------------------------------
         !Implementation
         !---------------------------------------------------------------------
             iopt = ppm_param_alloc_fit
             ldl3(1) = 1-ghostsize(1)
             ldl3(2) = 1-ghostsize(2)
             ldl3(3) = 1
             ldu3(1) = max_node(1,mlev)+ghostsize(1)
             ldu3(2) = max_node(2,mlev)+ghostsize(2)
             ldu3(3) = nsubs
             CALL ppm_alloc(uc,ldl3,ldu3,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_coarse',    &
       &                       'uc',__LINE__,info)
             GOTO 9999
             ENDIF
             ! write data from mgfield DS to temporary uc field
             DO isub=1,nsubs
               tuc=>mgfield(isub,mlev)%uc
               uc(:,:,isub)=tuc(:,:)
             ENDDO
             DO isweep=1,nsweep
               DO color=0,1
                DO isub=1,nsubs
                DO iface=1,4
                  IF (bcdef_sca(isub,iface).EQ.ppm_param_bcdef_dirichlet) THEN
                      IF (iface.EQ.1) THEN
                          i=1
                          DO j=1,max_node(2,mlev)
                            uc(i,j,isub)=0.0_MK
                          ENDDO
                      ELSEIF (iface.EQ.2) THEN
                          i=max_node(1,mlev)
                          DO j=1,max_node(2,mlev)
                            uc(i,j,isub)=0.0_MK
                          enddo
                      ELSEIF (iface.EQ.3) THEN
                          j=1
                          DO i=1,max_node(1,mlev)
                            uc(i,j,isub)=0.0_MK
                          ENDDO
                      ELSEIF (iface.EQ.4) THEN
                          j=max_node(2,mlev)
                          DO i=1,max_node(1,mlev)
                            uc(i,j,isub)=0.0_MK
                          ENDDO
                      ENDIF !iface
                  ENDIF !bckind
                ENDDO!iface
                ENDDO!DO isub 
               !----------------------------------------------------------------
               !Communicate
               !----------------------------------------------------------------

               CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlev),&
               &                         ghostsize,info)
               CALL ppm_map_field_push(topoid,mg_meshid(mlev),uc,&
               &                         info)
               CALL ppm_map_field_send(info)
               CALL ppm_map_field_pop(topoid,mg_meshid(mlev),uc,&
               &                          ghostsize,info)
               DO isub=1,nsubs
              DO j=lorig(2,isub),lext(2,isub)
                 DO i=lorig(1,isub)+mod(j+color,2), &
                    & lext(1,isub)-mod(j+color,2),2
                     IF ((i.GE.1.AND.i.LE.max_node(1,mlev)).AND.&
                     (j.GE.1.AND.j.LE.max_node(2,mlev))) THEN
                     uc(i,j,isub) = uc(i,j,isub)+omega*(c1*( &
                     &  (uc(i-1,j,isub)+uc(i+1,j,isub))*c2 + &
                     &  (uc(i,j-1,isub)+uc(i,j+1,isub))*c3 + &
                     &           mgfield(isub,mlev)%fc(i,j)) &
                     &                        -uc(i,j,isub)) 

                   ENDIF
                 ENDDO
               ENDDO
             ENDDO!isub
           ENDDO!DO color   
           IF (isweep.EQ.nsweep) THEN
             CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlev),&
             &                         ghostsize,info)
             CALL ppm_map_field_push(topoid,mg_meshid(mlev),uc,&
             &                         info)
             CALL ppm_map_field_send(info)
             CALL ppm_map_field_pop(topoid,mg_meshid(mlev),uc,&
             &                          ghostsize,info)
           ENDIF
         ENDDO!DO nsweep
         DO isub=1,nsubs
           tuc=>mgfield(isub,mlev)%uc
           tuc(:,:)=uc(:,:,isub)
         ENDDO  
             iopt = ppm_param_dealloc
             CALL ppm_alloc(uc,ldl3,ldu3,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_coarse',    &
       &                       'uc',__LINE__,info)
             GOTO 9999
             ENDIF
#elif __MESH_DIM == __3D
         !----------------------------------------------------------------------
         !Implementation
         !---------------------------------------------------------------------
             iopt = ppm_param_alloc_fit
             ldl4(1) = 1-ghostsize(1)
             ldl4(2) = 1-ghostsize(2)
             ldl4(3) = 1-ghostsize(3)
             ldl4(4) = 1
             ldu4(1) = max_node(1,mlev)+ghostsize(1)
             ldu4(2) = max_node(2,mlev)+ghostsize(2)
             ldu4(3) = max_node(3,mlev)+ghostsize(3)
             ldu4(4) = nsubs
             CALL ppm_alloc(uc,ldl4,ldu4,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_coarse',    &
       &                       'uc',__LINE__,info)
             GOTO 9999
             ENDIF
             ! write data from mgfield DS to temporary uc field
             DO isub=1,nsubs
               tuc=>mgfield(isub,mlev)%uc
               uc(:,:,:,isub)=tuc(:,:,:)
             ENDDO
         DO isweep=1,nsweep 
            DO color=0,1
                DO isub=1,nsubs
                DO iface=1,6
                  IF (bcdef_sca(isub,iface).EQ.ppm_param_bcdef_dirichlet) THEN
                      IF (iface.EQ.1) THEN
                          i=1
                          DO k=1,max_node(3,mlev)
                            DO j=1,max_node(2,mlev)
                              uc(i,j,k,isub)=0.0_MK
                            ENDDO
                          ENDDO
                      ELSEIF (iface.EQ.2) THEN
                          i=max_node(1,mlev)
                          DO k=1,max_node(3,mlev)
                            DO j=1,max_node(2,mlev)
                              uc(i,j,k,isub)=0.0_MK
                            ENDDO
                          ENDDO
                      ELSEIF (iface.EQ.3) THEN
                          j=1
                          DO k=1,max_node(3,mlev)
                            DO i=1,max_node(1,mlev)
                              uc(i,j,k,isub)=0.0_MK
                            ENDDO
                          ENDDO
                      ELSEIF (iface.EQ.4) THEN
                          j=max_node(2,mlev)
                          DO k=1,max_node(3,mlev)
                            DO i=1,max_node(1,mlev)
                              uc(i,j,k,isub)=0.0_MK
                            ENDDO
                          ENDDO
                      ELSEIF (iface.EQ.5) THEN
                          k=1
                          DO j=1,max_node(2,mlev)
                            DO i=1,max_node(1,mlev)
                              uc(i,j,k,isub)=0.0_MK
                            ENDDO
                          ENDDO
                      ELSEIF (iface.EQ.6) THEN
                          k=max_node(3,mlev)
                          DO j=1,max_node(2,mlev)
                            DO i=1,max_node(1,mlev)
                              uc(i,j,k,isub)=0.0_MK
                            ENDDO
                          ENDDO
                      ENDIF !iface
                  ENDIF !bckind
                ENDDO!iface
                ENDDO!DO isub 
               !----------------------------------------------------------------
               !Communicate
               !----------------------------------------------------------------
              CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlev),&
        &                         ghostsize,info)
              CALL ppm_map_field_push(topoid,mg_meshid(mlev),uc,&
        &                         info)
              CALL ppm_map_field_send(info)
              CALL ppm_map_field_pop(topoid,mg_meshid(mlev),uc,&
        &                          ghostsize,info)
                  
              DO isub=1,nsubs
                DO k=lorig(3,isub),lext(3,isub)
                  DO j=lorig(2,isub),lext(2,isub)
                    DO i=lorig(1,isub)+mod(j+color,2), &
                      & lext(1,isub)-mod(j+color,2),2
                      IF ((i.GE.1.AND.i.LE.max_node(1,mlev)).AND.&
                      &   (j.GE.1.AND.j.LE.max_node(2,mlev)).AND.&
                      &   (k.GE.1.AND.k.LE.max_node(3,mlev))) THEN
                      uc(i,j,k,isub) = uc(i,j,k,isub)+omega*(c1*( &
                      &  (uc(i-1,j,k,isub)+uc(i+1,j,k,isub))*c2 + &
                      &  (uc(i,j-1,k,isub)+uc(i,j+1,k,isub))*c3 + &
                      &  (uc(i,j,k-1,isub)+uc(i,j,k+1,isub))*c4 + &
                      &           mgfield(isub,mlev)%fc(i,j,k)) &
                      &                        -uc(i,j,k,isub)) 

                      ENDIF
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO!isub
            ENDDO!DO color
               IF (isweep.EQ.nsweep) THEN
                CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlev),&
        &                         ghostsize,info)
                CALL ppm_map_field_push(topoid,mg_meshid(mlev),uc,&
        &                         info)
                CALL ppm_map_field_send(info)

                CALL ppm_map_field_pop(topoid,mg_meshid(mlev),uc,&
        &                          ghostsize,info)
               ENDIF
           ENDDO!Do isweep
           DO isub=1,nsubs
             tuc=>mgfield(isub,mlev)%uc
             tuc(:,:,:)=uc(:,:,:,isub)
           ENDDO  
           iopt = ppm_param_dealloc
           CALL ppm_alloc(uc,ldl4,ldu4,iopt,info)
           IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_coarse',    &
             &                       'uc',__LINE__,info)
             GOTO 9999
           ENDIF
#endif
#elif __DIM == __VFIELD
#if  __MESH_DIM == __2D
         !----------------------------------------------------------------------
         !Implementation
         !---------------------------------------------------------------------
             iopt = ppm_param_alloc_fit
             ldl4(1) = 1
             ldl4(2) = 1-ghostsize(1)
             ldl4(3) = 1-ghostsize(2)
             ldl4(4) = 1
             ldu4(1) = vecdim
             ldu4(2) = max_node(1,mlev)+ghostsize(1)
             ldu4(3) = max_node(2,mlev)+ghostsize(2)
             ldu4(4) = nsubs
             CALL ppm_alloc(uc,ldl4,ldu4,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_coarse',    &
       &                       'uc',__LINE__,info)
             GOTO 9999
             ENDIF
         DO isweep=1,nsweep
            DO color=0,1
               DO isub=1,nsubs
                  tuc=>mgfield(isub,mlev)%uc
                  uc(:,:,:,isub)=tuc(:,:,:)
               ENDDO!DO isub 
               !----------------------------------------------------------------
               !Communicate 
               !----------------------------------------------------------------
               CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlev),&
        &                         ghostsize,info)
              CALL ppm_map_field_push(topoid,mg_meshid(mlev),uc,&
        &                         vecdim,info)
              CALL ppm_map_field_send(info)
              CALL ppm_map_field_pop(topoid,mg_meshid(mlev),uc,&
        &                          vecdim,ghostsize,info)
               DO isub=1,nsubs
                  tuc=>mgfield(isub,mlev)%uc
                  tuc(:,:,:)=uc(&
      &                         :,:,:,isub)
                  DO j=start(2,isub,mlev),istop(2,isub,mlev)
                     DO i=start(1,isub,mlev)+mod(j+color,2),istop(1,isub,mlev),2
                      DO ilda=1,vecdim
                           tuc(ilda,i,j) = c1*(&
      &                                   (tuc(ilda,i-1,j)+ &
      &                                tuc(ilda,i+1,j))*c2 + &
      &                                 (tuc(ilda,i,j-1)+&
      &                                  tuc(ilda,i,j+1))*c3-&
      &                                         mgfield(isub,mlev)%fc(ilda,i,j))
                      ENDDO  
                     ENDDO
                  ENDDO
               ENDDO
                    IF (isweep.EQ.nsweep) THEN
                     IF (color.EQ.1) THEN
                      DO isub=1,nsubs
                       tuc=>mgfield(isub,mlev)%uc
                       uc(:,:,:,isub)=tuc(:,:,:)
                      ENDDO
                     ENDIF
                    ENDIF
            ENDDO!DO color   
              IF (isweep.EQ.nsweep) THEN
              CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlev),&
        &                         ghostsize,info)
              CALL ppm_map_field_push(topoid,mg_meshid(mlev),uc,&
        &                         vecdim,info)
              CALL ppm_map_field_send(info)
              CALL ppm_map_field_pop(topoid,mg_meshid(mlev),uc,&
        &                          vecdim,ghostsize,info)
               DO isub=1,nsubs
                  tuc=>mgfield(isub,mlev)%uc
                  tuc(:,:,:)=uc(:,:,:,isub)
               ENDDO
              ENDIF 
         ENDDO!DO nsweep
             iopt = ppm_param_dealloc
             ldl4(1) = 1
             ldl4(2) = 1-ghostsize(1)
             ldl4(3) = 1-ghostsize(2)
             ldl4(4) = 1
             ldu4(1) = vecdim
             ldu4(2) = max_node(1,mlev)+ghostsize(1)
             ldu4(3) = max_node(2,mlev)+ghostsize(2)
             ldu4(4) = nsubs
             CALL ppm_alloc(uc,ldl4,ldu4,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_coarse',    &
       &                       'uc',__LINE__,info)
             GOTO 9999
             ENDIF
#elif __MESH_DIM == __3D
         !----------------------------------------------------------------------
         !Implementation
         !---------------------------------------------------------------------
             iopt = ppm_param_alloc_fit
             ldl5(1) = 1
             ldl5(2) = 1-ghostsize(1)
             ldl5(3) = 1-ghostsize(2)
             ldl5(4) = 1-ghostsize(3)
             ldl5(5) = 1
             ldu5(1) = vecdim
             ldu5(2) = max_node(1,mlev)+ghostsize(1)
             ldu5(3) = max_node(2,mlev)+ghostsize(2)
             ldu5(4) = max_node(3,mlev)+ghostsize(3)
             ldu5(5) = nsubs
             CALL ppm_alloc(uc,ldl5,ldu5,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_coarse',    &
       &                       'uc',__LINE__,info)
             GOTO 9999
             ENDIF
            iopt = ppm_param_alloc_fit
            ldu1(1)=vecdim
            CALL ppm_alloc(moldu,ldu1,iopt,info)
            IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_coarse',    &
       &                       'moldu',__LINE__,info)
             GOTO 9999
            ENDIF
         DO isweep=1,nsweep 
            DO color=0,1
              DO isub=1,nsubs
                  !-------------------------------------------------------------
                  !Impose boundaries 
                  !-------------------------------------------------------------
                  tuc=>mgfield(isub,mlev)%uc
                  DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                    DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                      DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
#ifdef __VECTOR
                        uc(1,i,j,k,isub)=tuc(1,i,j,k)
                        uc(2,i,j,k,isub)=tuc(2,i,j,k)
                        uc(3,i,j,k,isub)=tuc(3,i,j,k)
#else
                       DO ilda=1,vecdim 
                        uc(ilda,i,j,k,isub)=tuc(ilda,i,j,k)
                       ENDDO 
#endif
                      ENDDO
                     ENDDO
                    ENDDO 
               ENDDO!DO isub 
               !----------------------------------------------------------------
               !Communicate 
               !----------------------------------------------------------------
               CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlev),&
        &                         ghostsize,info)
               CALL ppm_map_field_push(topoid,mg_meshid(mlev),uc,&
        &                         vecdim,info)
               CALL ppm_map_field_send(info)
               CALL ppm_map_field_pop(topoid,mg_meshid(mlev),uc,&
        &                          vecdim,ghostsize,info)
                 a=0
                 b=0
                 c=0
                 d=0
                 e=0
                 g=0
               DO isub=1,nsubs
                  tuc=>mgfield(isub,mlev)%uc
                    DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                     DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                      DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
#ifdef __VECTOR
                          tuc(1,i,j,k)=uc(1,i,j,k,isub)
                          tuc(2,i,j,k)=uc(2,i,j,k,isub)
                          tuc(3,i,j,k)=uc(3,i,j,k,isub)
#else
                       DO ilda=1,vecdim 
                          tuc(ilda,i,j,k)=uc(ilda,i,j,k,isub)
                       ENDDO
#endif
                      ENDDO
                     ENDDO
                   ENDDO
                 DO  ilda=1,vecdim
                  IF (.NOT.lperiodic) THEN
                    DO iface=1,6
                    IF (bcdef_vec(ilda,isub,iface).EQ.ppm_param_bcdef_periodic) THEN
                     !DO NOTHING
                    ELSEIF (bcdef_vec(ilda,isub,iface).EQ.ppm_param_bcdef_dirichlet) THEN
                      IF (iface.EQ.1) THEN
                        a(isub)=1
                        IF (bcdef_vec(ilda,isub,2).EQ.0) THEN
                          b(isub)=-1
                        ENDIF
                        i=1
                        DO j=1,max_node(2,mlev)
                          DO k=1,max_node(3,mlev)
                             tuc(ilda,i,j,k)=0.0_MK
                          ENDDO
                        ENDDO
                     ELSEIF (iface.EQ.2) THEN
                        b(isub)=1
                        IF (bcdef_vec(ilda,isub,1).EQ.0) THEN
                          a(isub)=-1
                        ENDIF
                       i=max_node(1,mlev)
                        DO j=1,max_node(2,mlev)
                          DO k=1,max_node(3,mlev)
                            tuc(ilda,i,j,k)=0.0_MK
                          ENDDO
                        ENDDO
                     ELSEIF (iface.EQ.3) THEN
                       c(isub)=1
                        IF (bcdef_vec(ilda,isub,4).EQ.0) THEN
                         d(isub)=-1
                        ENDIF
                       j=1
                        DO i=1,max_node(1,mlev)
                         Do k=1,max_node(3,mlev)
                              tuc(ilda,i,j,k)=0.0_MK

                         enddo
                        ENDDO
                     ELSEIF (iface.EQ.4) THEN
                       d(isub)=1
                        IF (bcdef_vec(ilda,isub,3).EQ.0) THEN
                         c(isub)=-1
                        ENDIF
                       j=max_node(2,mlev)
                        DO i=1,max_node(1,mlev)
                         Do k=1,max_node(3,mlev)
                              tuc(ilda,i,j,k)=0.0_MK
                         enddo
                        ENDDO
                     ELSEIF (iface.EQ.5) Then
                       e(isub)=1
                       IF (bcdef_vec(ilda,isub,6).EQ.0) THEN
                         g(isub)=-1
                       ENDIF
                       k=1
                        DO i=1,max_node(1,mlev)
                          DO j=1,max_node(2,mlev)
                             tuc(ilda,i,j,k)=0.0_MK
                          ENDDO
                        ENDDO
                      ELSEIF (iface.EQ.6) THEN
                        g(isub)=1
                        IF (bcdef_vec(ilda,isub,5).EQ.0) THEN
                         e(isub)=-1
                        ENDIF
                        k=max_node(3,mlev)
                        DO i=1,max_node(1,mlev)
                         Do j=1,max_node(2,mlev)
                             tuc(ilda,i,j,k)=0.0_MK
                         ENDDO
                        ENDDO
                      ENDIF
                  ENDIF
                 ENDDO!face
                ENDIF
                ENDDO!ilda
                  DO k=start(3,isub,mlev)+e(isub),istop(3,isub,mlev)-g(isub)  
                     DO j=start(2,isub,mlev)+c(isub),istop(2,isub,mlev)-d(isub)
                        DO i=start(1,isub,mlev)+mod(j+k+color,2)+a(isub), &
      &                 istop(1,isub,mlev)-b(isub)-mod(j+k+color,2),2
                         IF ((i.GE.1.AND.i.LE.max_node(1,mlev)).AND.(j.GE.1.AND.j.LE.max_node(2,mlev)) &
      &                    .AND.(k.GE.1.AND.k.LE.max_node(3,mlev))) THEN
#ifdef __VECTOR
                         moldu(1) = tuc(1,i,j,k)
                         moldu(2) = tuc(2,i,j,k)
                         moldu(3) = tuc(3,i,j,k)
#else
                      do ilda=1,vecdim
                         moldu(ilda) = tuc(ilda,i,j,k)
                      end do
#endif
#ifdef __VECTOR
                              tuc(1,i,j,k) = moldu(1)+&
      &                             omega*(& 
      &                             c1*((tuc(1,i-1,j,k)+ &
      &                            tuc(1,i+1,j,k))*c2 + &
      &                                 (tuc(1,i,j-1,k)+&
      &                            tuc(1,i,j+1,k))*c3 + &
      &                           (tuc(1,i,j,k-1)+&
      &                            tuc(1,i,j,k+1))*c4 - &
      &                            mgfield(isub,mlev)%fc(1,i,j,k))&
      &                            -moldu(1))
                              tuc(2,i,j,k) = moldu(2)+&
      &                             omega*(& 
      &                             c1*((tuc(2,i-1,j,k)+ &
      &                            tuc(2,i+1,j,k))*c2 + &
      &                                 (tuc(2,i,j-1,k)+&
      &                            tuc(2,i,j+1,k))*c3 + &
      &                           (tuc(2,i,j,k-1)+&
      &                            tuc(2,i,j,k+1))*c4 - &
      &                            mgfield(isub,mlev)%fc(2,i,j,k))&
      &                            -moldu(2))
                              tuc(3,i,j,k) = moldu(3)+&
      &                             omega*(& 
      &                             c1*((tuc(3,i-1,j,k)+ &
      &                            tuc(3,i+1,j,k))*c2 + &
      &                                 (tuc(3,i,j-1,k)+&
      &                            tuc(3,i,j+1,k))*c3 + &
      &                           (tuc(3,i,j,k-1)+&
      &                            tuc(3,i,j,k+1))*c4 - &
      &                            mgfield(isub,mlev)%fc(3,i,j,k))&
      &                            -moldu(3))
#else
                      DO ilda=1,vecdim
                              tuc(ilda,i,j,k) = moldu(ilda)+&
      &                             omega*(& 
      &                             c1*((tuc(ilda,i-1,j,k)+ &
      &                            tuc(ilda,i+1,j,k))*c2 + &
      &                                 (tuc(ilda,i,j-1,k)+&
      &                            tuc(ilda,i,j+1,k))*c3 + &
      &                           (tuc(ilda,i,j,k-1)+&
      &                            tuc(ilda,i,j,k+1))*c4 - &
      &                            mgfield(isub,mlev)%fc(ilda,i,j,k))&
      &                            -moldu(ilda))
                         ENDDO!ilda
#endif
                        ENDIF
                        ENDDO!i
                     ENDDO!j
                  ENDDO!k
               ENDDO!isubs   
                   IF (isweep.EQ.nsweep) THEN
                    IF (color.EQ.1) THEN
                     DO isub=1,nsubs
                       tuc=>mgfield(isub,mlev)%uc
                       DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                         DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                           DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
                             DO ilda=1,vecdim 
                                 uc(ilda,i,j,k,isub)=tuc(ilda,i,j,k)
                             ENDDO 
                           ENDDO
                         ENDDO
                       ENDDO 
                     ENDDO!isub   

                    ENDIF
                   ENDIF 
           ENDDO!DO color
          IF (isweep.EQ.nsweep) THEN
              CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlev),&
        &                         ghostsize,info)
              CALL ppm_map_field_push(topoid,mg_meshid(mlev),uc,&
        &                         vecdim,info)
              CALL ppm_map_field_send(info)
              CALL ppm_map_field_pop(topoid,mg_meshid(mlev),uc,&
        &                          vecdim,ghostsize,info)
                   DO isub=1,nsubs 
                    tuc=>mgfield(isub,mlev)%uc
                    DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                     DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                      DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
                       DO ilda=1,vecdim 
                          tuc(ilda,i,j,k)=uc(ilda,i,j,k,isub)
                       ENDDO
                      ENDDO
                     ENDDO
                    ENDDO
                   ENDDO
          ENDIF
         ENDDO!Do isweep
             iopt = ppm_param_dealloc
             ldl5(1) = 1
             ldl5(2) = 1-ghostsize(1)
             ldl5(3) = 1-ghostsize(2)
             ldl5(4) = 1-ghostsize(3)
             ldl5(5) = 1
             ldu5(1) = vecdim
             ldu5(2) = max_node(1,mlev)+ghostsize(1)
             ldu5(4) = max_node(2,mlev)+ghostsize(2)
             ldu5(4) = max_node(3,mlev)+ghostsize(3)
             ldu5(5) = nsubs
             CALL ppm_alloc(uc,ldl5,ldu5,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_coarse',    &
       &                       'uc',__LINE__,info)
             GOTO 9999
             ENDIF
            iopt = ppm_param_dealloc
            ldu1(1)=vecdim
            CALL ppm_alloc(moldu,ldu1,iopt,info)
            IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_coarse',    &
       &                       'moldu',__LINE__,info)
             GOTO 9999
            ENDIF
#endif
#endif
        !---------------------------------------------------------------------
        !  Deallocate local work arrays
        !----------------------------------------------------------------------
            iopt = ppm_param_dealloc
            CALL ppm_alloc(a,ldl1,ldu1,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_fine',    &
      &                       'a',__LINE__,info)
            GOTO 9999
            ENDIF
            CALL ppm_alloc(b,ldl1,ldu1,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_fine',    &
      &                       'b',__LINE__,info)
            GOTO 9999
            ENDIF
            CALL ppm_alloc(c,ldl1,ldu1,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_fine',    &
      &                       'c',__LINE__,info)
            GOTO 9999
            ENDIF
            CALL ppm_alloc(d,ldl1,ldu1,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_fine',    &
      &                       'd',__LINE__,info)
            GOTO 9999
            ENDIF
            CALL ppm_alloc(e,ldl1,ldu1,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_fine',    &
      &                       'e',__LINE__,info)
            GOTO 9999
            ENDIF
            CALL ppm_alloc(g,ldl1,ldu1,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_fine',    &
      &                       'g',__LINE__,info)
            GOTO 9999
            ENDIF
         !---------------------------------------------------------------------
         !  Return 
         !----------------------------------------------------------------------
 9999    CONTINUE
         CALL substop('ppm_mg_smooth_coarse',t0,info)
         RETURN

#if __DIM == __SFIELD
#if   __MESH_DIM   == __2D
#if    __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_smooth_coarse_2D_sca_s
#elif  __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_smooth_coarse_2D_sca_d
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_smooth_coarse_3D_sca_s
#elif  __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_smooth_coarse_3D_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if   __MESH_DIM   == __2D
#if    __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_smooth_coarse_2D_vec_s
#elif  __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_smooth_coarse_2D_vec_d
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_smooth_coarse_3D_vec_s
#elif  __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_smooth_coarse_3D_vec_d
#endif
#endif
#endif
