      !------------------------------------------------------------------------
      !  Subroutine   :                  ppm_mg_solv 
      !------------------------------------------------------------------------
      ! Copyright (c) 2011 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
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

       
#if   __DIM   == __SFIELD
#if   __MESH_DIM   == __2D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_solv_2d_sca_s(topo_id,u,f,&
      &                               initsweep,finsweep,restrsweep,prolsweep,&
      &                               Eu,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_solv_2d_sca_d(topo_id,u,f,&
      &                               initsweep,finsweep,restrsweep,prolsweep,&
      &                               Eu,info)
#endif 
#elif __MESH_DIM   == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_solv_3d_sca_s(topo_id,u,f,&
      &                               initsweep,finsweep,restrsweep,prolsweep,&
      &                               Eu,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_solv_3d_sca_d(topo_id,u,f,&
      &                               initsweep,finsweep,restrsweep,prolsweep,&
      &                               Eu,info)
#endif 
#endif 
#elif __DIM == __VFIELD
#if   __MESH_DIM   == __2D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_solv_2d_vec_s(topo_id,u,f,lda,&
      &                               initsweep,finsweep,restrsweep,prolsweep,&
      &                               Eu,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_solv_2d_vec_d(topo_id,u,f,lda,&
      &                               initsweep,finsweep,restrsweep,prolsweep,&
      &                               Eu,info)
#endif
#elif __MESH_DIM   == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_solv_3d_vec_s(topo_id,u,f,lda,&
      &                               initsweep,finsweep,restrsweep,prolsweep,&
      &                               Eu,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_solv_3d_vec_d(topo_id,u,f,lda,&
      &                               initsweep,finsweep,restrsweep,prolsweep,&
      &                               Eu,info)
#endif
#endif
#endif
    !!! Solves the given equation using the multigrid method


#include "ppm_define.h"
         !---------------------------------------------------------------------
         !  Modules 
         !----------------------------------------------------------------------
        USE ppm_module_data
        USE ppm_module_data_mg
        USE ppm_module_data_mesh
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_error
        USE ppm_module_alloc
        USE ppm_module_map
        USE ppm_module_mg_core
        USE ppm_module_mg_res
        USE ppm_module_mg_prolong
        USE ppm_module_mg_smooth
        USE ppm_module_write
         IMPLICIT NONE
#ifdef __MPI
       INCLUDE  'mpif.h'
#endif
#if    __KIND == __SINGLE_PRECISION
         INTEGER, PARAMETER :: MK = ppm_kind_single
#else
         INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
         !----------------------------------------------------------------------
         !  Arguments (for u and f index: local mesh locations and isub) 
         !----------------------------------------------------------------------
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
         REAL(MK),DIMENSION(:,:,:),POINTER     ::  u
         !!! the field of the solution (including the ghost layer)
         REAL(MK),DIMENSION(:,:,:),POINTER     ::  f
         !!! the field of the right hand side (without the ghost layer)
#elif __MESH_DIM == __3D
         REAL(MK),DIMENSION(:,:,:,:),POINTER   ::  u
         !!! the field of the solution (including the ghost layer)
         REAL(MK),DIMENSION(:,:,:,:),POINTER   ::  f
         !!! the field of the right hand side (without the ghost layer)
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
         REAL(MK),DIMENSION(:,:,:,:),POINTER     ::  u
         !!! the field of the solution (including the ghost layer)
         REAL(MK),DIMENSION(:,:,:,:),POINTER     ::  f
         !!! the field of the right hand side (without the ghost layer)
#elif __MESH_DIM == __3D
         REAL(MK),DIMENSION(:,:,:,:,:),POINTER   ::  u
         !!! the field of the solution (including the ghost layer)
         REAL(MK),DIMENSION(:,:,:,:,:),POINTER   ::  f
         !!! the field of the right hand side (without the ghost layer)
#endif
#endif
#if __DIM == __VFIELD
         INTEGER,INTENT(IN)                      :: lda
         !!! leading dimension
#endif
         INTEGER,                   INTENT(IN)   ::  initsweep
         !!! initial smoothing sweeps in the finest level
         INTEGER,                   INTENT(IN)   ::  finsweep
         !!! final smoothing sweeps in the finest level
         INTEGER,                   INTENT(IN)   ::  restrsweep
         !!! Number of smoothing sweeps after each restriction step
         !!! IMPORTANT: This parameter is considere important
         INTEGER,                   INTENT(IN)   ::  prolsweep
         !!! Number of smoothing sweeps after each prolongation step
         REAL(MK),                  INTENT(OUT)  ::  Eu  
         INTEGER,                   INTENT(INOUT)   ::  info
         INTEGER,                   INTENT(IN   )   ::  topo_id
         !----------------------------------------------------------------------
         !  Local variables 
         !----------------------------------------------------------------------
         REAL(MK)                             :: t0
         REAL(MK)                             :: E,res
         INTEGER                              :: iface,count,k
         INTEGER                              :: ix,iy  
         CHARACTER(LEN=256)                   :: cbuf
         INTEGER                              :: mlev,color,it
         INTEGER                              :: ncalls=0
         REAL(MK)                             :: c1,c2,c3,c4  
         INTEGER                              :: isub,i,j
         REAL(MK)                             :: x,y
         REAL(MK)                             :: gEu 
         INTEGER                              :: MPI_PREC
         TYPE(ppm_t_topo),      POINTER       :: topo => NULL()
         TYPE(ppm_t_equi_mesh), POINTER       :: mesh => NULL()
#if __MESH_DIM == __3D
         REAL(MK)                             :: c5,dz,rdz2
         INTEGER,DIMENSION(4)                 :: ldl4,ldu4
         INTEGER,DIMENSION(5)                 :: ldl5,ldu5
#endif
         INTEGER                              :: ilda
         REAL(MK)                             :: rdx2,rdy2
         REAL(MK)                             :: dx,dy
#if __MESH_DIM == __2D
         INTEGER,DIMENSION(3)                 :: ldl3,ldu3
         INTEGER,DIMENSION(4)                 :: ldl4,ldu4
#endif
         INTEGER                              :: topoid,iopt,idom
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
     REAL(MK),DIMENSION(:,:),POINTER :: tuc
#elif __MESH_DIM == __3D
     REAL(MK),DIMENSION(:,:,:),POINTER :: tuc
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
     REAL(MK),DIMENSION(:,:,:),POINTER :: tuc
#elif __MESH_DIM == __3D
     REAL(MK),DIMENSION(:,:,:,:),POINTER :: tuc
#endif
#endif
         !----------------------------------------------------------------------
         !  Externals 
         !----------------------------------------------------------------------

         !----------------------------------------------------------------------
         !  Initialize 
         !----------------------------------------------------------------------
         CALL substart('ppm_mg_solv',t0,info)
#ifdef __MPI
        IF (ppm_kind.EQ.ppm_kind_single) THEN
           MPI_PREC = MPI_REAL
        ELSE
           MPI_PREC = MPI_DOUBLE_PRECISION
        ENDIF
#endif
        topoid=topo_id
        topo => ppm_topo(topo_id)%t
        mesh => topo%mesh(mg_meshid(1))
        !----------------------------------------------------------------------
        !  Check arguments
        !----------------------------------------------------------------------
        IF (ppm_debug .GT. 0) THEN
            CALL check()
        ENDIF
        !----------------------------------------------------------------------
        ! Definition of necessary variables and allocation of arrays
        !----------------------------------------------------------------------
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
#if __DIM == __SFIELD
        mgfield=>mgfield_2d_sca_s
#elif __DIM == __VFIELD
        mgfield=>mgfield_2d_vec_s
#endif
        rdx2=rdx2_s
        rdy2=rdy2_s
        dx=dx_s
        dy=dy_s
#elif __KIND == __DOUBLE_PRECISION
#if __DIM == __SFIELD
        mgfield=>mgfield_2d_sca_d
#elif __DIM == __VFIELD
        mgfield=>mgfield_2d_vec_d
#endif
        rdx2=rdx2_d
        rdy2=rdy2_d
        dx=dx_d
        dy=dy_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
#if __DIM == __SFIELD
        mgfield=>mgfield_3d_sca_s
#elif __DIM == __VFIELD 
        mgfield=>mgfield_3d_vec_s
#endif
        rdx2=rdx2_s
        rdy2=rdy2_s
        rdz2=rdz2_s
        dx=dx_s
        dy=dy_s
        dz=dz_s
#elif __KIND == __DOUBLE_PRECISION
#if __DIM == __SFIELD
        mgfield=>mgfield_3d_sca_d
#elif __DIM == __VFIELD
        mgfield=>mgfield_3d_vec_d
#endif
        rdx2=rdx2_d
        rdy2=rdy2_d
        rdz2=rdz2_d
        dx=dx_d
        dy=dy_d
        dz=dz_d
#endif
#endif
     topoid=topo_id
     ncalls=ncalls+1
     IF (ncalls.EQ.1) THEN
        DO i=1,maxlev
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
             iopt = ppm_param_alloc_fit
             ldl3(1) = 1-ghostsize(1)
             ldl3(2) = 1-ghostsize(2)
             ldl3(3) = 1
             ldu3(1) = max_node(1,i)+ghostsize(1)
             ldu3(2) = max_node(2,i)+ghostsize(2)
             ldu3(3) = nsubs
             CALL ppm_alloc(uc_dummy,ldl3,ldu3,iopt,info)
             IF (info .NE. 0) THEN
               info = ppm_error_fatal
               CALL ppm_error(ppm_err_alloc,'ppm_mg_solv',    &
     &                       'uc_dummy',__LINE__,info)
               GOTO 9999
              ENDIF
             uc_dummy(:,:,:)=0.0_MK
#elif __MESH_DIM ==__3D
             iopt = ppm_param_alloc_fit
             ldl4(1) = 1-ghostsize(1)
             ldl4(2) = 1-ghostsize(2)
             ldl4(3) = 1-ghostsize(3)
             ldl4(4) = 1
             ldu4(1) = max_node(1,i)+ghostsize(1)
             ldu4(2) = max_node(2,i)+ghostsize(2)
             ldu4(3) = max_node(3,i)+ghostsize(3)
             ldu4(4) = nsubs
             CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
             IF (info .NE. 0) THEN
               info = ppm_error_fatal
               CALL ppm_error(ppm_err_alloc,'ppm_mg_solv',    &
      &                       'uc_dummy',__LINE__,info)
               GOTO 9999
              ENDIF
             uc_dummy(:,:,:,:)=0.0_MK
#endif
#if __MESH_DIM == __2D
             iopt = ppm_param_dealloc
             ldl3(1) = 1-ghostsize(1)
             ldl3(2) = 1-ghostsize(2)
             ldl3(3) = 1
             ldu3(1) = max_node(1,i)+ghostsize(1)
             ldu3(2) = max_node(2,i)+ghostsize(2)
             ldu3(3) = nsubs
             CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
             IF (info .NE. 0) THEN
               info = ppm_error_fatal
               CALL ppm_error(ppm_err_alloc,'ppm_mg_solv',    &
     &                       'uc_dummy',__LINE__,info)
               GOTO 9999
              ENDIF
#elif __MESH_DIM ==__3D
             iopt = ppm_param_dealloc
             ldl4(1) = 1-ghostsize(1)
             ldl4(2) = 1-ghostsize(2)
             ldl4(3) = 1-ghostsize(3)
             ldl4(4) = 1
             ldu4(1) = max_node(1,i)+ghostsize(1)
             ldu4(2) = max_node(2,i)+ghostsize(2)
             ldu4(3) = max_node(3,i)+ghostsize(3)
             ldu4(4) = nsubs
             CALL ppm_alloc(uc_dummy,ldl5,ldu5,iopt,info)
             IF (info .NE. 0) THEN
               info = ppm_error_fatal
               CALL ppm_error(ppm_err_alloc,'ppm_mg_solv',    &
      &                       'uc_dummy',__LINE__,info)
               GOTO 9999
              ENDIF
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
             iopt = ppm_param_alloc_fit
             ldl4(1) = 1
             ldl4(2) = 1-ghostsize(1)
             ldl4(3) = 1-ghostsize(2)
             ldl4(4) = 1
             ldu4(1) = vecdim
             ldu4(2) = max_node(1,i)+ghostsize(1)
             ldu4(3) = max_node(2,i)+ghostsize(2)
             ldu4(4) = nsubs
             CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
             IF (info .NE. 0) THEN
               info = ppm_error_fatal
               CALL ppm_error(ppm_err_alloc,'ppm_mg_solv',    &
     &                       'uc_dummy',__LINE__,info)
               GOTO 9999
              ENDIF
             uc_dummy(:,:,:,:)=0.0_MK
#elif __MESH_DIM ==__3D
             iopt = ppm_param_alloc_fit
             ldl5(1) = 1
             ldl5(2) = 1-ghostsize(1)
             ldl5(3) = 1-ghostsize(2)
             ldl5(4) = 1-ghostsize(3)
             ldl5(5) = 1
             ldu5(1) = vecdim
             ldu5(2) = max_node(1,i)+ghostsize(1)
             ldu5(3) = max_node(2,i)+ghostsize(2)
             ldu5(4) = max_node(3,i)+ghostsize(3)
             ldu5(5) = nsubs
             CALL ppm_alloc(uc_dummy,ldl5,ldu5,iopt,info)
             IF (info .NE. 0) THEN
               info = ppm_error_fatal
               CALL ppm_error(ppm_err_alloc,'ppm_mg_solv',    &
      &                       'uc_dummy',__LINE__,info)
               GOTO 9999
              ENDIF
             uc_dummy(:,:,:,:,:)=0.0_MK
#endif
#if __MESH_DIM == __2D
             iopt = ppm_param_dealloc
             ldl4(1) = 1-ghostsize(1)
             ldl4(1) = 1-ghostsize(2)
             ldl4(1) = 1
             ldu4(1) = max_node(1,i)+ghostsize(1)
             ldu4(2) = max_node(2,i)+ghostsize(2)
             ldu4(3) = nsubs
             CALL ppm_alloc(uc_dummy,ldl3,ldu3,iopt,info)
             IF (info .NE. 0) THEN
               info = ppm_error_fatal
               CALL ppm_error(ppm_err_alloc,'ppm_mg_solv',    &
     &                       'uc_dummy',__LINE__,info)
               GOTO 9999
              ENDIF
#elif __MESH_DIM ==__3D
             iopt = ppm_param_dealloc
             ldl5(1) = 1-ghostsize(1)
             ldl5(2) = 1-ghostsize(2)
             ldl5(3) = 1-ghostsize(3)
             ldl5(4) = 1
             ldu5(1) = max_node(1,i)+ghostsize(1)
             ldu5(2) = max_node(2,i)+ghostsize(2)
             ldu5(3) = max_node(3,i)+ghostsize(3)
             ldu5(4) = nsubs
             CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
             IF (info .NE. 0) THEN
               info = ppm_error_fatal
               CALL ppm_error(ppm_err_alloc,'ppm_mg_solv',    &
      &                       'uc_dummy',__LINE__,info)
               GOTO 9999
              ENDIF
#endif
#endif
        ENDDO
        ncalls=ncalls+1
      ENDIF
        !----------------------------------------------------------------------
        ! DO initial sweeps in the finest mesh with the smoother to get the 
        ! initial solution 
        !----------------------------------------------------------------------
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
        c1 = 1.0_MK/(2.0_MK*(rdx2+rdy2))  
        c2 = rdx2
        c3 = rdy2     
        c4 = 2.0_MK*c2+2.0_MK*c3
        count = 0
        CALL ppm_mg_smooth_sca(topo_id,u,f,initsweep,1,c1,c2,c3,info)
        !----------------------------------------------------------------------
        ! Compute residual
        !----------------------------------------------------------------------
        CALL ppm_mg_res_sca(topo_id,u,f,c1,c2,c3,c4,E,info)
#ifdef __MPI
        CALL MPI_AllReduce(E,gEu,1,MPI_PREC,MPI_MAX,ppm_comm,info)
        E=gEu
#endif
      IF (info .NE. 0) THEN 
         GOTO 9999
      ENDIF 
      IF (ppm_debug.GT.0) THEN 
        WRITE(cbuf,*) 'Eu:',E
        CALL PPM_Write(ppm_rank,'mg_solv',cbuf,info)
      ENDIF
        !---------------------------------------------------------------------
        !Initiation of the function correction. (We start on purpose with lev=2)
        !----------------------------------------------------------------------
        DO mlev=2,maxlev
           DO isub=1,nsubs
              tuc=>mgfield(isub,mlev)%uc
              DO j=start(2,isub,mlev),istop(2,isub,mlev)
                 DO i=start(1,isub,mlev),istop(1,isub,mlev)
                       tuc(i,j)=0.0_MK
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        !----------------------------------------------------------------------
        ! Run Multigrid core routine (execute smoothing, restriction,
        ! prologongation cycles)
        !----------------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
        CALL ppm_mg_core_2d_sca_s(topo_id,2,restrsweep,prolsweep,info)  
#elif __KIND == __DOUBLE_PRECISION
        CALL ppm_mg_core_2d_sca_d(topo_id,2,restrsweep,prolsweep,info)  
#endif   
        !----------------------------------------------------------------------
        !PROLONG the solution to the finest grid
        !----------------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
         CALL ppm_mg_prolong_2d_sca_s(1,info)
#elif __KIND == __DOUBLE_PRECISION
         CALL ppm_mg_prolong_2d_sca_d(1,info)
#endif
        !----------------------------------------------------------------------
        !UPDATE THE FUNCTION
        !----------------------------------------------------------------------
        DO isub=1,nsubs
           !tuc=>mgfield(isub,mlev)%uc
           tuc=>mgfield(isub,1)%uc
           DO j=start(2,isub,1),istop(2,isub,1)   
              DO i=start(1,isub,1),istop(1,isub,1)
                    u(i,j,isub)=tuc(i,j) 
              ENDDO
           ENDDO
        ENDDO
        !----------------------------------------------------------------------
        !DO the final sweeps
        !--------------------------------------------------------------------
        CALL ppm_mg_smooth_sca(topo_id,u,f,finsweep,1,c1,c2,c3,info)
        CALL ppm_mg_res_sca(topo_id,u,f,c1,c2,c3,c4,E,info)
#ifdef __MPI
        CALL MPI_AllReduce(E,gEu,1,MPI_PREC,MPI_MAX,ppm_comm,info)
        Eu=gEu
#else
        Eu=E
#endif
#elif __MESH_DIM == __3D
        c1 = 1.0_MK/(2.0_MK*(rdx2+rdy2+rdz2))
        c2 = rdx2
        c3 = rdy2
        c4 = rdz2 
        c5 = 2.0_MK*c2+2.0_MK*c3+2.0_MK*c4
        CALL ppm_mg_smooth_sca(topo_id,u,f,initsweep,1,c1,c2,c3,c4,info)
        !----------------------------------------------------------------------
        ! Compute residual
        !----------------------------------------------------------------------
        CALL ppm_mg_res_sca(topo_id,u,f,c1,c2,c3,c4,c5,E,info)
#ifdef __MPI
        CALL MPI_AllReduce(E,gEu,1,MPI_PREC,MPI_MAX,ppm_comm,info)
        E=gEu
#endif
      IF (ppm_debug.GT.0) THEN 
        WRITE(cbuf,*) 'Eu:',E
        CALL PPM_WRITE(ppm_rank,'mg_solv',cbuf,info)
      ENDIF
         !---------------------------------------------------------------------
        !Initiation of the function correction. (We start on purpose with lev=2)
        !----------------------------------------------------------------------
        DO mlev=2,maxlev
           DO isub=1,nsubs
              tuc=>mgfield(isub,mlev)%uc
              DO k=start(3,isub,mlev),istop(3,isub,mlev) 
                 DO j=start(2,isub,mlev),istop(2,isub,mlev)
                    DO i=start(1,isub,mlev),istop(1,isub,mlev)
                          tuc(i,j,k)=0.0_MK
                    ENDDO
                ENDDO
              ENDDO
           ENDDO
        ENDDO
        !----------------------------------------------------------------------
        !CALL THE MULTIGRID TO DO NICE STUFF TO OUR FUNCTION
        !----------------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
        CALL ppm_mg_core_3d_sca_s(topo_id,2,restrsweep,prolsweep,info)
#elif __KIND == __DOUBLE_PRECISION
        CALL ppm_mg_core_3d_sca_d(topo_id,2,restrsweep,prolsweep,info)
#endif
        !----------------------------------------------------------------------
        !PROLONG the solution to the finest grid
        !----------------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
        CALL ppm_mg_prolong_3d_sca_s(1,info)
#elif __KIND == __DOUBLE_PRECISION
        CALL ppm_mg_prolong_3d_sca_d(1,info)
#endif
        !----------------------------------------------------------------------
        !UPDATE THE FUNCTION
        !----------------------------------------------------------------------
        DO isub=1,nsubs
              !tuc=>mgfield(isub,mlev)%uc
              tuc=>mgfield(isub,1)%uc
           DO k=start(3,isub,1),istop(3,isub,1)
              DO j=start(2,isub,1),istop(2,isub,1)
                 DO i=start(1,isub,1),istop(1,isub,1)
                       u(i,j,k,isub)=tuc(i,j,k)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        !----------------------------------------------------------------------
        !DO the final sweeps
          !--------------------------------------------------------------------
        CALL ppm_mg_smooth_sca(topo_id,u,f,finsweep,1,c1,c2,c3,c4,info)
        CALL ppm_mg_res_sca(topo_id,u,f,c1,c2,c3,c4,c5,E,info)
#ifdef __MPI
        CALL MPI_AllReduce(E,gEu,1,MPI_PREC,MPI_MAX,ppm_comm,info)
        Eu=gEu
#else
        Eu=E
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
        c1 = 1.0_MK/(2.0_MK*(rdx2+rdy2))  
        c2 = rdx2
        c3 = rdy2     
        c4 = 2.0_MK*c2+2.0_MK*c3
        count = 0
        CALL ppm_mg_smooth_vec(topo_id,u,f,initsweep,1,c1,c2,c3,info) 
        !----------------------------------------------------------------------
        ! Compute residual
        !----------------------------------------------------------------------
        CALL ppm_mg_res_vec(topo_id,u,f,c1,c2,c3,c4,E,info)   
#ifdef __MPI
        CALL MPI_AllReduce(E,gEu,1,MPI_PREC,MPI_MAX,ppm_comm,info)
        E=gEu
#endif
        IF (ppm_debug.GT.0) THEN 
         WRITE(cbuf,*) 'Eu:',E
         CALL PPM_WRITE(ppm_rank,'mg_solv',cbuf,info)
        ENDIF

         !---------------------------------------------------------------------
        !Initiation of the function correction. (We start on purpose with lev=2)
        !----------------------------------------------------------------------
        DO mlev=2,maxlev
           DO isub=1,nsubs
              tuc=>mgfield(isub,mlev)%uc
              DO j=start(2,isub,mlev),istop(2,isub,mlev)
                 DO i=start(1,isub,mlev),istop(1,isub,mlev)
                  DO ilda=1,vecdim
                       tuc(ilda,i,j)=0.0_MK
                  ENDDO 
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        !----------------------------------------------------------------------
        !CALL THE MULTIGRID TO DO NICE STUFF TO OUR FUNCTION
        !----------------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
        CALL ppm_mg_core_2d_vec_s(topo_id,2,restrsweep,prolsweep,info)  
#elif __KIND == __DOUBLE_PRECISION
        CALL ppm_mg_core_2d_vec_d(topo_id,2,restrsweep,prolsweep,info)  
#endif   
        !----------------------------------------------------------------------
        !PROLONG the solution to the finest grid
        !----------------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
        CALL ppm_mg_prolong_2d_vec_s(1,info)
#elif __KIND == __DOUBLE_PRECISION
        CALL ppm_mg_prolong_2d_vec_d(1,info)
#endif   
        !----------------------------------------------------------------------
        !UPDATE THE FUNCTION
        !----------------------------------------------------------------------
        DO isub=1,nsubs
           !tuc=>mgfield(isub,mlev)%uc
           tuc=>mgfield(isub,1)%uc
           DO j=start(2,isub,1),istop(2,isub,1)   
              DO i=start(1,isub,1),istop(1,isub,1)
               DO ilda=1,vecdim
                    u(ilda,i,j,isub)=tuc(ilda,i,j)
               ENDDO 
              ENDDO
           ENDDO
        ENDDO
        !----------------------------------------------------------------------
        !DO the final sweeps
        !--------------------------------------------------------------------
        CALL ppm_mg_smooth_vec(topo_id,u,f,finsweep,1,c1,c2,c3,info)
        CALL ppm_mg_res_vec(topo_id,u,f,c1,c2,c3,c4,E,info)   

#ifdef __MPI
        CALL MPI_AllReduce(E,gEu,1,MPI_PREC,MPI_MAX,ppm_comm,info)
        Eu=gEu        
#else
        Eu=E
#endif
#elif __MESH_DIM == __3D
        c1 = 1.0_MK/(2.0_MK*(rdx2+rdy2+rdz2))
        c2 = rdx2
        c3 = rdy2
        c4 = rdz2 
        c5 = 2.0_MK*c2+2.0_MK*c3+2.0_MK*c4
        CALL ppm_mg_smooth_vec(topo_id,u,f,initsweep,1,c1,c2,c3,c4,info)
        !-----------------------------------------------------------------
        ! Compute residual
        !-----------------------------------------------------------------

        CALL ppm_mg_res_vec(topo_id,u,f,c1,c2,c3,c4,c5,E,info)
#ifdef __MPI
        CALL MPI_AllReduce(E,gEu,1,MPI_PREC,MPI_MAX,ppm_comm,info)
        E=gEu        
#endif
        IF (ppm_debug.GT.0) THEN 
         WRITE(cbuf,*) 'Eu:',E
         CALL PPM_WRITE(ppm_rank,'mg_solv',cbuf,info)
        ENDIF
        !---------------------------------------------------------------------
        !Initiation of the function correction. (We start on purpose with lev=2)
        !----------------------------------------------------------------------
        DO mlev=2,maxlev
           DO isub=1,nsubs
              tuc=>mgfield(isub,mlev)%uc
              DO k=start(3,isub,mlev),istop(3,isub,mlev) 
                 DO j=start(2,isub,mlev),istop(2,isub,mlev)
                    DO i=start(1,isub,mlev),istop(1,isub,mlev)
#ifdef __VECTOR
                          tuc(1,i,j,k)=0.0_MK
                          tuc(2,i,j,k)=0.0_MK
                          tuc(3,i,j,k)=0.0_MK
#else
                     DO ilda=1,vecdim 
                          tuc(ilda,i,j,k)=0.0_MK
                     ENDDO
#endif
                    ENDDO
                ENDDO
              ENDDO
           ENDDO
        ENDDO
        !----------------------------------------------------------------------
        !CALL THE MULTIGRID TO DO NICE STUFF TO OUR FUNCTION
        !----------------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
        CALL ppm_mg_core_3d_vec_s(topo_id,2,restrsweep,prolsweep,info)
#elif __KIND == __DOUBLE_PRECISION
        CALL ppm_mg_core_3d_vec_d(topo_id,2,restrsweep,prolsweep,info)
#endif
        !----------------------------------------------------------------------
        !PROLONG the solution to the finest grid
        !----------------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
        CALL ppm_mg_prolong_3d_vec_s(1,info)
#elif __KIND == __DOUBLE_PRECISION
        CALL ppm_mg_prolong_3d_vec_d(1,info)
#endif
        !----------------------------------------------------------------------
        !UPDATE THE FUNCTION
        !----------------------------------------------------------------------
        DO isub=1,nsubs
           !tuc=>mgfield(isub,mlev)%uc
           tuc=>mgfield(isub,1)%uc
           DO k=start(3,isub,1),istop(3,isub,1)
              DO j=start(2,isub,1),istop(2,isub,1)
                 DO i=start(1,isub,1),istop(1,isub,1)
#ifdef __VECTOR
                       u(1,i,j,k,isub)=tuc(1,i,j,k)
                       u(2,i,j,k,isub)=tuc(2,i,j,k)
                       u(3,i,j,k,isub)=tuc(3,i,j,k)
#else
                  DO ilda=1,vecdim
                       u(ilda,i,j,k,isub)=tuc(ilda,i,j,k)
                  ENDDO
#endif
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        !----------------------------------------------------------------------
        !DO the final sweeps
        !--------------------------------------------------------------------
        CALL ppm_mg_smooth_vec(topo_id,u,f,finsweep,1,c1,c2,c3,c4,info)
        CALL ppm_mg_res_vec(topo_id,u,f,c1,c2,c3,c4,c5,E,info)   
#ifdef __MPI
        CALL MPI_AllReduce(E,gEu,1,MPI_PREC,MPI_MAX,ppm_comm,info)
        Eu=gEu
#else
        Eu=E
#endif
#endif
#endif
        !----------------------------------------------------------------------
        !  Return 
        !----------------------------------------------------------------------
9999    CONTINUE
        CALL substop('ppm_mg_solv',t0,info)
        RETURN
        CONTAINS

        SUBROUTINE check
        
#if __DIM == __SFIELD        
#if __MESH_DIM == __2D        
           IF (SIZE(u,3) .LT. nsubs) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'solution exist on nsubs subdomains',__LINE__,info)        
              GOTO 8888
           ENDIF
           DO i=1,nsubs
              idom=topo%isublist(i)
              IF (SIZE(u(:,:,i),1).LT.mesh%nnodes(1,idom)+ &
     &                                              2*ghostsize(1)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'solution mess with mesh points in x-dir!',__LINE__,info)  
                 GOTO 8888    
              ENDIF
              IF (SIZE(u(:,:,i),2).LT.mesh%nnodes(2,idom) &
     &                                             +2*ghostsize(2)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'solution mess with mesh points in y-dir!',__LINE__,info)
                 GOTO 8888
              ENDIF
           ENDDO
           IF (SIZE(f,3) .LT. nsubs) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'rhs exist on nsubs subdomains!',__LINE__,info)  
              GOTO 8888
           ENDIF
           DO i=1,nsubs
              idom=topo%isublist(i)
              IF (SIZE(f(:,:,i),1).LT. mesh%nnodes(1,idom)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'rhs mess with mesh points in x-dir!',__LINE__,info)
                 GOTO 8888
              ENDIF
              IF (SIZE(f(:,:,i),2).LT. mesh%nnodes(2,idom)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'rhs mess with mesh points in y-dir!',__LINE__,info)
                 GOTO 8888
              ENDIF
           ENDDO
#elif __MESH_DIM == __3D
           IF (SIZE(u,4) .LT. nsubs) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'solution exist on nsubs subdomains!',__LINE__,info)        
              GOTO 8888
           ENDIF
           DO i=1,nsubs
              idom=topo%isublist(i)
              IF (SIZE(u(:,:,:,i),1).LT.mesh%nnodes(1,idom)+ &
     &                                              2*ghostsize(1)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'solution mess with mesh points in x-dir!',__LINE__,info)  
                 GOTO 8888    
              ENDIF
              IF (SIZE(u(:,:,:,i),2).LT.mesh%nnodes(2,idom)+ &
     &                                              2*ghostsize(1)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'solution mess with mesh points in y-dir!',__LINE__,info)
                 GOTO 8888
              ENDIF
              IF (SIZE(u(:,:,:,i),3).LT.mesh%nnodes(3,idom)+ &
     &                                              2*ghostsize(1)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'solution mess with mesh points in z-dir!',__LINE__,info)
                 GOTO 8888
              ENDIF
          ENDDO
           IF (SIZE(f,4) .LT. nsubs) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'rhs exist on nsubs subdomains!',__LINE__,info)  
              GOTO 8888
           ENDIF
           DO i=1,nsubs
              idom=topo%isublist(i)
              IF (SIZE(f(:,:,:,i),1).LT.mesh%nnodes(1,idom)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &                   'rhs mess with mesh points in x-dir!',__LINE__,info)
                 GOTO 8888
              ENDIF
              IF (SIZE(f(:,:,:,i),2).LT.mesh%nnodes(2,idom)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &            'rhs mess with mesh points in y-dir!',__LINE__,info)
                 GOTO 8888
              ENDIF
              IF (SIZE(f(:,:,:,i),3).LT.mesh%nnodes(3,idom)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'rhs mess with mesh points in z-dir!',__LINE__,info)
                 GOTO 8888
              ENDIF
           ENDDO
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D        
           IF (SIZE(u,4) .LT. nsubs) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'solution exist on nsubs subdomains',__LINE__,info)        
              GOTO 8888
           ENDIF
           DO i=1,nsubs
              idom=topo%isublist(i)
              IF (SIZE(u(:,:,:,i),2).LT.mesh%nnodes(1,idom)+2) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'solution mess with mesh points in x-dir!',__LINE__,info)  
                 GOTO 8888    
              ENDIF
              IF (SIZE(u(:,:,:,i),3).LT.mesh%nnodes(2,idom)+2) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'solution mess with mesh points in y-dir!',__LINE__,info)
                 GOTO 8888
              ENDIF
           ENDDO
           IF (SIZE(f,4) .LT. nsubs) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'rhs exist on nsubs subdomains!',__LINE__,info)  
              GOTO 8888
           ENDIF
           DO i=1,nsubs
              idom=topo%isublist(i)
              IF (SIZE(f(:,:,:,i),2).LT.mesh%nnodes(1,idom)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'rhs mess with mesh points in x-dir!',__LINE__,info)
                 GOTO 8888
              ENDIF
              IF (SIZE(f(:,:,:,i),3).LT.mesh%nnodes(2,idom)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'rhs mess with mesh points in y-dir!',__LINE__,info)
                 GOTO 8888
              ENDIF
           ENDDO
#elif __MESH_DIM == __3D
           IF (SIZE(u,5) .LT. nsubs) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'solution exist on nsubs subdomains!',__LINE__,info)        
              GOTO 8888
           ENDIF
           DO i=1,nsubs
              idom=topo%isublist(i)
              IF (SIZE(u(:,:,:,:,i),2).LT.mesh%nnodes(1,idom)+ &
     &                                              2*ghostsize(1)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'solution mess with mesh points in x-dir!',__LINE__,info)  
                 GOTO 8888    
              ENDIF
              IF (SIZE(u(:,:,:,:,i),3).LT.mesh%nnodes(2,idom)+ &
     &                                              2*ghostsize(1)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'solution mess with mesh points in y-dir!',__LINE__,info)
                 GOTO 8888
              ENDIF
              IF (SIZE(u(:,:,:,:,i),4).LT.mesh%nnodes(3,idom)+ &
     &                                              2*ghostsize(1)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'solution mess with mesh points in z-dir!',__LINE__,info)
                 GOTO 8888
              ENDIF
          ENDDO
           IF (SIZE(f,5) .LT. nsubs) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'rhs exist on nsubs subdomains!',__LINE__,info)  
              GOTO 8888
           ENDIF
           DO i=1,nsubs
              idom=topo%isublist(i)
              IF (SIZE(f(:,:,:,:,i),2).LT.mesh%nnodes(1,idom)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &                   'rhs mess with mesh points in x-dir!',__LINE__,info)
                 GOTO 8888
              ENDIF
              IF (SIZE(f(:,:,:,:,i),3).LT.mesh%nnodes(2,idom)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &            'rhs mess with mesh points in y-dir!',__LINE__,info)
                 GOTO 8888
              ENDIF
              IF (SIZE(f(:,:,:,:,i),4).LT.mesh%nnodes(3,idom)) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
     &             'rhs mess with mesh points in z-dir!',__LINE__,info)
                 GOTO 8888
              ENDIF
           ENDDO
#endif
#endif

8888    CONTINUE
        RETURN

        END SUBROUTINE check

#if    __DIM == __SFIELD
#if    __MESH_DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_solv_2d_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_solv_2d_sca_d
#endif
#elif  __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_solv_3d_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_solv_3d_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if    __MESH_DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_solv_2d_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_solv_2d_vec_d
#endif
#elif  __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_solv_3d_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_solv_3d_vec_d
#endif
#endif
#endif
