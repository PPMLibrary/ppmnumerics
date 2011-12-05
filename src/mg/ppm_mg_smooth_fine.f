      !-----------------------------------------------------------------------
      !  Subroutine   :            ppm_mg_smooth_fine
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
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_smooth_fine_2D_sca_s(topo_id,u,f,nsweep,mlev,&
     &                                       c1,c2,c3,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_smooth_fine_2D_sca_d(topo_id,u,f,nsweep,mlev,&
     &                                       c1,c2,c3,info)
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_smooth_fine_3D_sca_s(topo_id,u,f,nsweep,mlev,&
     &                                       c1,c2,c3,c4,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_smooth_fine_3D_sca_d(topo_id,u,f,nsweep,mlev,&
     &                                       c1,c2,c3,c4,info)
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_smooth_fine_2D_vec_s(topo_id,u,f,nsweep,mlev,&
     &                                       c1,c2,c3,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_smooth_fine_2D_vec_d(topo_id,u,f,nsweep,mlev,&
     &                                       c1,c2,c3,info)
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_smooth_fine_3D_vec_s(topo_id,u,f,nsweep,mlev,&
     &                                       c1,c2,c3,c4,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_smooth_fine_3D_vec_d(topo_id,u,f,nsweep,mlev,&
     &                                       c1,c2,c3,c4,info)
#endif
#endif
#endif
         !!! In this routine we compute the corrections for the function based 
         !!! on the Gauss-Seidel iteration
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
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
         REAL(MK),DIMENSION(:,:,:),POINTER     ::  u
         REAL(MK),DIMENSION(:,:,:),POINTER     ::  f
#elif __MESH_DIM == __3D
         REAL(MK),DIMENSION(:,:,:,:),POINTER   ::  u
         REAL(MK),DIMENSION(:,:,:,:),POINTER   ::  f
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
         REAL(MK),DIMENSION(:,:,:,:),POINTER     ::  u
         REAL(MK),DIMENSION(:,:,:,:),POINTER     ::  f
#elif __MESH_DIM == __3D
         REAL(MK),DIMENSION(:,:,:,:,:),POINTER   ::  u
         REAL(MK),DIMENSION(:,:,:,:,:),POINTER   ::  f
#endif
#endif
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
         INTEGER                                    ::  i,j,isub,color
         INTEGER                                    ::  ilda,isweep,count
         REAL(MK)                                   ::  c11,c22,c33,c44
         REAL(MK)                                   ::  dx,dy
         INTEGER,DIMENSION(:),POINTER               ::  a => NULL()
         INTEGER,DIMENSION(:),POINTER               ::  b => NULL()
         INTEGER,DIMENSION(:),POINTER               ::  c => NULL()
         INTEGER,DIMENSION(:),POINTER               ::  d => NULL()
         INTEGER,DIMENSION(:),POINTER               ::  e => NULL()
         INTEGER,DIMENSION(:),POINTER               ::  g => NULL()
         INTEGER                                    ::  k,idom
         REAL(MK)                                   ::  x,y
         REAL(MK)                                   ::  omega
         INTEGER,DIMENSION(1)                       ::  ldl1,ldu1
#if __MESH_DIM == __2D
         INTEGER,DIMENSION(4)                       ::  ldl4,ldu4
         INTEGER,DIMENSION(3)                       ::  ldl3,ldu3
#endif
#if __MESH_DIM == __3D
        REAL(MK)                                   ::  dz
        INTEGER,DIMENSION(5)                       ::  ldl5,ldu5
        INTEGER,DIMENSION(4)                       ::  ldl4,ldu4
#endif
        INTEGER                                    ::  iopt,iface,topoid
        REAL(MK)                                   ::  t0
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
     REAL(MK),DIMENSION(:,:,:),POINTER :: oldu
#elif __MESH_DIM == __3D
     REAL(MK),DIMENSION(:,:,:,:),POINTER :: oldu
#endif
#elif  __DIM == __VFIELD
#if __MESH_DIM == __2D
     REAL(MK),DIMENSION(:,:,:,:),POINTER :: oldu
#elif __MESH_DIM == __3D
     REAL(MK),DIMENSION(:,:,:,:,:),POINTER :: oldu
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
     REAL(MK),DIMENSION(:),POINTER :: moldu
#elif __MESH_DIM == __3D
     REAL(MK),DIMENSION(:),POINTER :: moldu
#endif
#endif
        !----------------------------------------------------------------------
        !Externals
        !----------------------------------------------------------------------

        !----------------------------------------------------------------------
        !Initialize
        !----------------------------------------------------------------------
        CALL substart('ppm_mg_smooth_fine',t0,info)
        !----------------------------------------------------------------------
        !  Check arguments
        !----------------------------------------------------------------------
        IF (ppm_debug .GT. 0) THEN
          IF (nsweep.LT.1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_fine',  &
     &            'nsweep must be >=1',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (mlev.LT.1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_fine',  &
     &            'level must be >1',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (c1.LE.0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_fine',  &
     &            'Factor c1 must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (c2.LE.0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_fine',  &
     &            'Factor c2 must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (c3.LE.0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_fine',  &
     &            'Factor c3 must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
#if __MESH_DIM == __3D
          IF (c4.LE.0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_fine',  &
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
            iopt = ppm_param_alloc_fit
            ldl1(1) = 1
            ldu1(1) = nsubs
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
#if  __DIM == __SFIELD
#if  __MESH_DIM == __2D
        !----------------------------------------------------------------------
        !Implementation
        !---------------------------------------------------------------------
        count = 0
        DO isweep=1,nsweep
          DO color=0,1
            a=0
            b=0
            c=0
            d=0
            DO isub=1,nsubs
              !-------------------------------------------------------------
              !Impose boundaries on even if color=0 or odd if color=1
              !-------------------------------------------------------------
              IF (.NOT.lperiodic) THEN
                DO iface=1,4
                  IF (bcdef_sca(isub,iface).EQ.ppm_param_bcdef_periodic) THEN
                    !DO NOTHING
                  ELSEIF (bcdef_sca(isub,iface).EQ.ppm_param_bcdef_dirichlet) THEN
                      IF (iface.EQ.1) THEN
                          a(isub)=1
                          IF (bcdef_sca(isub,2).EQ.0) THEN
                            b(isub)=-1
                          ENDIF
                          i=1
                          DO j=1,max_node(2,mlev)
                            u(i,j,isub)=mgfield(isub,1)%bcvalue(iface)%pbcvalue(j)
                          ENDDO
                      ELSEIF (iface.EQ.2) THEN
                          b(isub)=1
                          IF (bcdef_sca(isub,1).EQ.0) THEN
                            a(isub)=-1
                          ENDIF
                          i=max_node(1,mlev)
                          DO j=1,max_node(2,mlev)
                            u(i,j,isub)=mgfield(isub,1)%bcvalue(iface)%pbcvalue(j)
                          enddo
                      ELSEIF (iface.EQ.3) THEN
                          c(isub)= 1
                          IF (bcdef_sca(isub,4).EQ.0) THEN
                            d(isub)=-1
                          ENDIF
                          j=1
                          DO i=1,max_node(1,mlev)
                            u(i,j,isub)=mgfield(isub,1)%bcvalue(iface)%pbcvalue(i)
                          ENDDO
                      ELSEIF (iface.EQ.4) THEN
                          d(isub)=1
                          IF (bcdef_sca(isub,3).EQ.0) THEN
                            c(isub)=-1
                          ENDIF
                          j=max_node(2,mlev)
                          DO i=1,max_node(1,mlev)
                            u(i,j,isub)=mgfield(isub,1)%bcvalue(iface)%pbcvalue(i)
                          ENDDO
                      ENDIF !iface
                    ENDIF !bckind
                  ENDDO!iface
                ENDIF !periodic
              ENDDO!DO isub 
              !----------------------------------------------------------------
              !Communicate
              !----------------------------------------------------------------

              CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlev),&
        &                         ghostsize,info)
              CALL ppm_map_field_push(topoid,mg_meshid(mlev),u,info)
              CALL ppm_map_field_send(info)

              CALL ppm_map_field_pop(topoid,mg_meshid(mlev),u,ghostsize,info)
            DO isub=1,nsubs
              DO j=start(2,isub,1)+c(isub),istop(2,isub,1)-d(isub)
                 DO i=start(1,isub,1)+mod(j+k+color,2)+a(isub), &
     &                istop(1,isub,1)-b(isub)-mod(j+k+color,2),2
                     IF ((i.GE.1.AND.i.LE.max_node(1,mlev)).AND. &
     &                   (j.GE.1.AND.j.LE.max_node(2,mlev))) THEN
                       u(i,j,isub)=u(i,j,isub) + omega*(c1*(  &
                       &   (u(i-1,j,isub)+u(i+1,j,isub))*c2 + &
                       &   (u(i,j-1,isub)+u(i,j+1,isub))*c3 - &
                       &        f(i,j,isub)) - u(i,j,isub))
                     ENDIF
                 ENDDO
              ENDDO
             ENDDO !isub
        ENDDO!DO color
        IF (isweep.EQ.nsweep) THEN
              CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlev),&
        &                         ghostsize,info)
              CALL ppm_map_field_push(topoid,mg_meshid(mlev),u,info)
              CALL ppm_map_field_send(info)

              CALL ppm_map_field_pop(topoid,mg_meshid(mlev),u,ghostsize,info)

        ENDIF
      ENDDO
#elif __MESH_DIM == __3D
        !----------------------------------------------------------------------
        !Implementation
        !---------------------------------------------------------------------
        DO isweep=1,nsweep
           DO color=0,1
              a=0
              b=0
              c=0
              d=0
              e=0
              g=0
              DO isub=1,nsubs
                 !-------------------------------------------------------------
                 !Impose boundaries on even if color=0 or odd if color=1
                 !-------------------------------------------------------------
                IF (.NOT.lperiodic) THEN
                 DO iface=1,6
                  IF (bcdef_sca(isub,iface).EQ.ppm_param_bcdef_periodic) THEN
                   !DO NOTHING
                  ELSEIF (bcdef_sca(isub,iface).EQ.ppm_param_bcdef_dirichlet) THEN
                    IF (iface.EQ.1) THEN
                        a(isub)=1
                       IF (bcdef_sca(isub,2).EQ.0) THEN
                        b(isub)=-1
                       ENDIF
                      i=1
                       DO j=1,max_node(2,mlev)
                        DO k=1,max_node(3,mlev)
                           u(i,j,k,isub)=mgfield(isub,1)%bcvalue(iface)%pbcvalue(j,k)
                        ENDDO
                       ENDDO
                    ELSEIF (iface.EQ.2) THEN
                       b(isub)=1
                       IF (bcdef_sca(isub,1).EQ.0) THEN
                        a(isub)=-1
                       ENDIF
                      i=max_node(1,mlev)
                       DO j=1,max_node(2,mlev)
                        DO k=1,max_node(3,mlev)
                          u(i,j,k,isub)=mgfield(isub,1)%bcvalue(iface)%pbcvalue(j,k)
                        ENDDO
                       enddo
                    ELSEIF (iface.EQ.3) THEN
                        c(isub)= 1
                       IF (bcdef_sca(isub,4).EQ.0) THEN
                        d(isub)=-1
                       ENDIF
                      j=1
                       DO i=1,max_node(1,mlev)
                        Do k=1,max_node(3,mlev)
                          u(i,j,k,isub)=mgfield(isub,1)%bcvalue(iface)%pbcvalue(i,k)
                        enddo
                       ENDDO
                    ELSEIF (iface.EQ.4) THEN
                      d(isub)=1
                       IF (bcdef_sca(isub,3).EQ.0) THEN
                        c(isub)=-1
                       ENDIF
                      j=max_node(2,mlev)
                       DO i=1,max_node(1,mlev)
                        Do k=1,max_node(3,mlev)
                          u(i,j,k,isub)=mgfield(isub,1)%bcvalue(iface)%pbcvalue(i,k)
                        enddo
                       ENDDO
                    ELSEIF (iface.EQ.5) Then
                        e(isub)=1
                       IF (bcdef_sca(isub,6).EQ.0) THEN
                        g(isub)=-1
                       ENDIF
                       k=1
                       DO i=1,max_node(1,mlev)
                        Do j=1,max_node(2,mlev)
                          u(i,j,k,isub)=mgfield(isub,1)%bcvalue(iface)%pbcvalue(i,j)
                        enddo
                       ENDDO
                     ELSEIF (iface.EQ.6) Then
                       g(isub)=1
                       IF (bcdef_sca(isub,5).EQ.0) THEN
                        e(isub)=-1
                       ENDIF
                       k=max_node(3,mlev)
                       DO i=1,max_node(1,mlev)
                        Do j=1,max_node(2,mlev)
                          u(i,j,k,isub)=mgfield(isub,1)%bcvalue(iface)%pbcvalue(i,j)

                        ENDDO
                       ENDDO
                     ENDIF
                 ENDIF
              ENDDO!iface
           ENDIF
         ENDDO!DO isub
              !----------------------------------------------------------------
              !Communicate red(even) if color==0 or communicate black(odd)
              !if color==1
              !----------------------------------------------------------------
              CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlev),&
        &                         ghostsize,info)
              CALL ppm_map_field_push(topoid,mg_meshid(mlev),u,info)
              CALL ppm_map_field_send(info)
              CALL ppm_map_field_pop(topoid,mg_meshid(mlev),u,ghostsize,info)
              DO isub=1,nsubs
              DO k=start(3,isub,1)+e(isub),istop(3,isub,1)-g(isub)
                 DO j=start(2,isub,1)+c(isub),istop(2,isub,1)-d(isub)
                    DO i=start(1,isub,1)+mod(j+k+color,2)+a(isub), &
     &                istop(1,isub,1)-b(isub)-mod(j+k+color,2),2
                        IF ((i.GE.1.AND.i.LE.max_node(1,mlev)).AND. &
     &                      (j.GE.1.AND.j.LE.max_node(2,mlev)).AND. &
     &                      (k.GE.1.AND.k.LE.max_node(3,mlev))) THEN
                          u(i,j,k,isub) = u(i,j,k,isub)+omega*(&
     &                       c1*((u(i-1,j,k,isub)+u(i+1,j,k,isub))*c2+ &
     &                           (u(i,j-1,k,isub)+u(i,j+1,k,isub))*c3+ &
     &                           (u(i,j,k-1,isub)+u(i,j,k+1,isub))*c4- &
     &                               f(i,j,k,isub))-u(i,j,k,isub))
                      ENDIF
                    ENDDO
                ENDDO
              ENDDO
           ENDDO!subs
        ENDDO!DO color
        IF (isweep.EQ.nsweep) THEN
            CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlev),&
            &                         ghostsize,info)
            CALL ppm_map_field_push(topoid,mg_meshid(mlev),u,info)
            CALL ppm_map_field_send(info)
            CALL ppm_map_field_pop(topoid,mg_meshid(mlev),u,ghostsize,info)
        ENDIF
      ENDDO
#endif
#elif __DIM == __VFIELD
#if  __MESH_DIM == __2D
        !----------------------------------------------------------------------
        !Implementation
        !---------------------------------------------------------------------
        count = 0
        DO isweep=1,nsweep
           DO color=0,1
              !----------------------------------------------------------------
              !Communicate red(even) if color==0 or communicate black(odd)
              !if color==1
              !----------------------------------------------------------------
              CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlev),&
        &                         ghostsize,info)
              CALL ppm_map_field_push(topoid,mg_meshid(mlev),u,vecdim,info)
              CALL ppm_map_field_send(info)
              CALL ppm_map_field_pop(topoid,mg_meshid(mlev),u,&
        &                          vecdim,ghostsize,info)
           DO isub=1,nsubs
              DO j=start(2,isub,1),istop(2,isub,1)
                 DO i=start(1,isub,1)+mod(j+color,2),istop(1,isub,1),2
                  DO ilda=1,vecdim
                       u(ilda,i,j,isub)=c1*((u(ilda,i-1,j,isub)+&
     &                                       u(ilda,i+1,j,isub))*c2 &
     &                  +(u(ilda,i,j-1,isub)+u(ilda,i,j+1,isub))*c3 -  &
     &                                             f(ilda,i,j,isub))
                  ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO!DO color
        IF (isweep.EQ.nsweep) THEN
              CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlev),&
        &                         ghostsize,info)
              CALL ppm_map_field_push(topoid,mg_meshid(mlev),u,vecdim,info)
              CALL ppm_map_field_send(info)
              CALL ppm_map_field_pop(topoid,mg_meshid(mlev),u,&
        &                          vecdim,ghostsize,info)
         ENDIF
       ENDDO
#elif __MESH_DIM == __3D
        !----------------------------------------------------------------------
        !Implementation
        !---------------------------------------------------------------------
            iopt = ppm_param_alloc_fit
            ldu1(1)=vecdim
            CALL ppm_alloc(moldu,ldu1,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_fine',    &
      &                       'moldu',__LINE__,info)
            GOTO 9999
            ENDIF
        DO isweep=1,nsweep
           DO color=0,1
              a=0
              b=0
              c=0
              d=0
              e=0
              g=0
              DO isub=1,nsubs
                DO ilda=1,vecdim
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
                            u(ilda,i,j,k,isub)=mgfield(isub,1)%bcvalue(iface)%pbcvalue(ilda,j,k)
                        enddo
                       ENDDO
                    ELSEIF (iface.EQ.2) THEN
                       b(isub)=1
                       IF (bcdef_vec(ilda,isub,1).EQ.0) THEN
                        a(isub)=-1
                       ENDIF
                       i=max_node(1,mlev)
                       DO j=1,max_node(2,mlev)
                        DO k=1,max_node(3,mlev)
                            u(ilda,i,j,k,isub)=mgfield(isub,1)%bcvalue(iface)%pbcvalue(ilda,j,k)
                        ENDDO
                       enddo
                    ELSEIF (iface.EQ.3) THEN
                      c(isub)=1
                       IF (bcdef_vec(ilda,isub,4).EQ.0) THEN
                        d(isub)=-1
                       ENDIF
                      j=1
                       DO i=1,max_node(1,mlev)
                        Do k=1,max_node(3,mlev)
                          u(ilda,i,j,k,isub)=mgfield(isub,1)%bcvalue(iface)%pbcvalue(ilda,i,k)
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
                           u(ilda,i,j,k,isub)=mgfield(isub,1)%bcvalue(iface)%pbcvalue(ilda,i,k)
                        enddo
                       ENDDO
                    ELSEIF (iface.EQ.5) Then
                       e(isub)=1
                       IF (bcdef_vec(ilda,isub,6).EQ.0) THEN
                        g(isub)=-1
                       ENDIF
                       k=1
                       DO i=1,max_node(1,mlev)
                        Do j=1,max_node(2,mlev)
                            u(ilda,i,j,k,isub)=mgfield(isub,1)%bcvalue(iface)%pbcvalue(ilda,i,j)
                        enddo
                       ENDDO
                     ELSEIF (iface.EQ.6) Then
                       g(isub)=1
                       IF (bcdef_vec(ilda,isub,5).EQ.0) THEN
                        e(isub)=-1
                       ENDIF
                       k=max_node(3,mlev)
                       DO i=1,max_node(1,mlev)
                        Do j=1,max_node(2,mlev)
                             u(ilda,i,j,k,isub)=mgfield(isub,1)%bcvalue(iface)%pbcvalue(ilda,i,j)
                        ENDDO
                       ENDDO
                    ENDIF
           ENDIF
                     ENddo !iface
                endif !periodic
              Enddo !ilda
         ENDDO!DO isub
              !----------------------------------------------------------------
              !Communicate red(even) if color==0 or communicate black(odd)
              !if color==1
              !----------------------------------------------------------------
              CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlev),&
        &                         ghostsize,info)
              CALL ppm_map_field_push(topoid,mg_meshid(mlev),u,vecdim,info)
              CALL ppm_map_field_send(info)
              CALL ppm_map_field_pop(topoid,mg_meshid(mlev),u,&
        &                          vecdim,ghostsize,info)
#ifdef  __VECTOR
             DO isub=1,nsubs
              DO k=start(3,isub,1)+e(isub),istop(3,isub,1)-g(isub)
                 DO j=start(2,isub,1)+c(isub),istop(2,isub,1)-d(isub)
                    DO i=start(1,isub,1)+mod(j+k+color,2)+a(isub), &
     &                istop(1,isub,1)-b(isub)-mod(j+k+color,2),2
                        IF ((i.GE.1.AND.i.LE.max_node(1,mlev)).AND. &
     & (j.GE.1.AND.j.LE.max_node(2,mlev)).AND.(k.GE.1.AND.k.LE.max_node(3,mlev))) THEN
                        moldu(1) = u(1,i,j,k,isub)
                        moldu(2) = u(2,i,j,k,isub)
                        moldu(3) = u(3,i,j,k,isub)
                       u(1,i,j,k,isub)=moldu(1)+omega*&
     &                           (&
                       &c1*((u(1,i-1,j,k,isub)+ &
     &                                          u(1,i+1,j,k,isub))*c2 &
     &                  +(u(1,i,j-1,k,isub)+u(1,i,j+1,k,isub))*c3 &
     &                  +(u(1,i,j,k-1,isub)+u(1,i,j,k+1,isub))*c4- &
     &                                                 f(1,i,j,k,isub))&
&-moldu(1))
                       u(2,i,j,k,isub)=moldu(2)+omega*&
     &                           (&
                       &c1*((u(2,i-1,j,k,isub)+ &
     &                                          u(2,i+1,j,k,isub))*c2 &
     &                  +(u(2,i,j-1,k,isub)+u(2,i,j+1,k,isub))*c3 &
     &                  +(u(2,i,j,k-1,isub)+u(2,i,j,k+1,isub))*c4- &
     &                                                 f(2,i,j,k,isub))&
&-moldu(2))
                       u(3,i,j,k,isub)=moldu(3)+omega*&
     &                           (&
                       &c1*((u(3,i-1,j,k,isub)+ &
     &                                          u(3,i+1,j,k,isub))*c2 &
     &                  +(u(3,i,j-1,k,isub)+u(3,i,j+1,k,isub))*c3 &
     &                  +(u(3,i,j,k-1,isub)+u(3,i,j,k+1,isub))*c4- &
     &                                                 f(3,i,j,k,isub))&
&-moldu(3))
                    ENDDO
                ENDDO
              ENDDO
            ENDDO!subs
#else
             DO isub=1,nsubs
              DO k=start(3,isub,1)+e(isub),istop(3,isub,1)-g(isub)
                 DO j=start(2,isub,1)+c(isub),istop(2,isub,1)-d(isub)
                    DO i=start(1,isub,1)+mod(j+k+color,2)+a(isub), &
     &                istop(1,isub,1)-b(isub)-mod(j+k+color,2),2
                        IF ((i.GE.1.AND.i.LE.max_node(1,mlev)).AND. &
     &                    (j.GE.1.AND.j.LE.max_node(2,mlev)).AND.(k.GE.1.AND.k.LE.max_node(3,mlev))) THEN
                     do ilda=1,vecdim
                        moldu(ilda) = u(ilda,i,j,k,isub)
                     end do
                     DO ilda=1,vecdim
                       u(ilda,i,j,k,isub)=moldu(ilda)+omega*&
     &                           (&
                       &c1*((u(ilda,i-1,j,k,isub)+ &
     &                                          u(ilda,i+1,j,k,isub))*c2 &
     &                  +(u(ilda,i,j-1,k,isub)+u(ilda,i,j+1,k,isub))*c3 &
     &                  +(u(ilda,i,j,k-1,isub)+u(ilda,i,j,k+1,isub))*c4- &
     &                                                 f(ilda,i,j,k,isub))&
&-moldu(ilda))
                     ENDDO
                 ENDIF
                    ENDDO
                ENDDO
              ENDDO
            ENDDO!subs
#endif
          ENDDO!DO color
            IF (isweep.EQ.nsweep) THEN
              CALL ppm_map_field_ghost_get(topoid,mg_meshid(mlev),&
        &                         ghostsize,info)
              CALL ppm_map_field_push(topoid,mg_meshid(mlev),u,vecdim,info)
              CALL ppm_map_field_send(info)
              CALL ppm_map_field_pop(topoid,mg_meshid(mlev),u,&
        &                          vecdim,ghostsize,info)
           ENDIF
       ENDDO
            iopt = ppm_param_dealloc
            ldu1(1)=vecdim
            CALL ppm_alloc(moldu,ldu1,iopt,info)
            IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_mg_smooth_fine',    &
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
        CALL substop('ppm_mg_smooth_fine',t0,info)
        RETURN

#if __DIM == __SFIELD
#if   __MESH_DIM   == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_smooth_fine_2D_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_smooth_fine_2D_sca_d
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_smooth_fine_3D_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_smooth_fine_3D_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if   __MESH_DIM   == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_smooth_fine_2D_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_smooth_fine_2D_vec_d
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_smooth_fine_3D_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_smooth_fine_3D_vec_d
#endif
#endif
#endif
