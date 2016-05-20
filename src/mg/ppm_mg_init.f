      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_mg_init
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


#if __DIM == __SFIELD
#if __MESH_DIM  == __2D
#if __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_init_2d_sca_s(topo_id,mesh_id,equation,ighostsize,&
       &          smoother,ibcdef,bcvalue,limlev,wcycle,omega,info)
#elif  __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_init_2d_sca_d(topo_id,mesh_id,equation,ighostsize,&
       &          smoother,ibcdef,bcvalue,limlev,wcycle,omega,info)
#endif
#elif  __MESH_DIM  == __3D
#if    __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_init_3d_sca_s(topo_id,mesh_id,equation,ighostsize,&
       &          smoother,ibcdef,bcvalue,limlev,wcycle,omega,info)
#elif  __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_init_3d_sca_d(topo_id,mesh_id,equation,ighostsize,&
       &          smoother,ibcdef,bcvalue,limlev,wcycle,omega,info)
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM  == __2D
#if __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_init_2d_vec_s(topo_id,mesh_id,equation,ighostsize,&
       &          smoother,lda,ibcdef,bcvalue,limlev,wcycle,omega,info)
#elif  __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_init_2d_vec_d(topo_id,mesh_id,equation,ighostsize,&
       &          smoother,lda,ibcdef,bcvalue,limlev,wcycle,omega,info)
#endif
#elif  __MESH_DIM  == __3D
#if    __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_init_3d_vec_s(topo_id,mesh_id,equation,ighostsize,&
       &          smoother,lda,ibcdef,bcvalue,limlev,wcycle,omega,info)
#elif  __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_init_3d_vec_d(topo_id,mesh_id,equation,ighostsize,&
       &          smoother,lda,ibcdef,bcvalue,limlev,wcycle,omega,info)
#endif
#endif
#endif
         !!!  This routine initializes the multigrid solver for 2D and 3D
         !!!  problems
         !!!
         !!! [NOTE]
         !!! Please pay attention that in order to be able to coarsen the mesh
         !!! it should be divisible with 2.
         !!! If you want to solve different equations the whole machinery should
         !!! be called twice. Also the solver is currently programmed for the
         !!! Poisson problem. A future improvement woudl be to use a general
         !!! stencil.
         !----------------------------------------------------------------------
         !  Includes
         !----------------------------------------------------------------------
#include "ppm_define.h"
         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
         USE ppm_module_data
         USE ppm_module_data_mg
         USE ppm_module_mg_alloc
         USE ppm_module_alloc
         USE ppm_module_error
         USE ppm_module_write
         USE ppm_module_mesh
         USE ppm_module_mesh_derive
         USE ppm_module_substart
         USE ppm_module_substop
         USE ppm_module_typedef
         IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION
         INTEGER, PARAMETER :: MK = ppm_kind_single
#else
         INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
         !----------------------------------------------------------------------
         !  Arguments
         !----------------------------------------------------------------------
         INTEGER,  INTENT(IN)                               :: topo_id
         !!! ID of current topology
         INTEGER,  INTENT(IN)                               :: mesh_id
         !!! ID of mesh of data fields
         INTEGER, INTENT(IN)                                :: equation
         !!! Kind of equation to be solved.
         !!!
         !!! Currently only ppm_param_eq_poisson supported
         INTEGER,DIMENSION(:),INTENT(IN)                    :: ighostsize
         !!! Ghostlayer size
         INTEGER, INTENT(IN)                                :: smoother
         !!! Smoother to be used.
         !!!
         !!! Currently only ppm_param_smooth_rbsor supported
#if __DIM == __VFIELD
         INTEGER,              INTENT(IN)                   ::  lda
         !!! Leading dimension
#endif
#if __DIM == __SFIELD
         INTEGER,DIMENSION(:)                               ::  ibcdef
         !!! Boundary condition types. Any of
         !!!
         !!! * ppm_param_bcdef_periodic
         !!! * ppm_param_bcdef_dirichlet
#if __MESH_DIM == __2D
         REAL(MK),DIMENSION(:,:,:)                          ::  bcvalue
#elif __MESH_DIM == __3D
         REAL(MK),DIMENSION(:,:,:,:)                        ::  bcvalue
#endif
         !!! Boundary condition values to be used.
         !!!
         !!! In the case of periodic BC, the content of bcvalue is ignored
         !!! The indeces (and their sizes) are (4. only in 3D):
         !!!
         !!! 1. isub (nsublist)
         !!! 2. dim (2*ppm_dim) (west,east,south,north,bottom,top)
         !!! 3. index1 (maximum extent of field in any direction)
         !!! 4. (index2 (maximum extent of field in any direction))
#elif __DIM == __VFIELD
         INTEGER,DIMENSION(:,:)                               ::  ibcdef
         !!! Boundary condition types. Any of
         !!!
         !!! * ppm_param_bcdef_periodic
         !!! * ppm_param_bcdef_dirichlet
#if __MESH_DIM == __2D
         REAL(MK),DIMENSION(:,:,:,:)                          ::  bcvalue
#elif __MESH_DIM == __3D
         REAL(MK),DIMENSION(:,:,:,:,:)                        ::  bcvalue
#endif
         !!! Boundary condition values to be used.
         !!!
         !!! In the case of periodic BC, the content of bcvalue is ignored
         !!! The indeces (and their sizes) are (5. only in 3D):
         !!!
         !!! 1. vector index
         !!! 2. isub (nsublist)
         !!! 3. dim (2*ppm_dim) (west,east,south,north,bottom,top)
         !!! 4. index1 (maximum extent of field in any direction)
         !!! 5. (index2 (maximum extent of field in any direction))
#endif
         INTEGER,INTENT(IN)                                 :: limlev
         !!! Number of levels to coarsen.
         LOGICAL,INTENT(IN)                                 :: wcycle
         !!! TRUE if the user wants W-cycle otherwise FALSE
         REAL(MK),INTENT(IN)                                :: omega
         !!! relaxation parameter for SOR
         INTEGER, INTENT(OUT)                               :: info
         !!! return status. 0 upon success.
         !----------------------------------------------------------------------
         !  Local variables
         !----------------------------------------------------------------------
         REAL(MK)                             :: t0
         REAL(MK)                             :: lmyeps
         INTEGER                              :: mlev,isub
         INTEGER                              :: idom
         INTEGER                              :: count,ilda,iface
         INTEGER                              :: i,j,k
         INTEGER                              :: kk
         TYPE(ppm_t_topo),      POINTER       :: topo
         TYPE(ppm_t_equi_mesh), POINTER       :: mesh
#if __MESH_DIM == __2D
         INTEGER                              :: dir
#endif
         INTEGER                              :: iter1,iter2,ix,iy
         INTEGER                              :: ipoint,jpoint
         INTEGER                              :: meshid,newmeshid
         INTEGER , DIMENSION(1)               :: ldu1
         INTEGER , DIMENSION(2)               :: ldu2,ldl2 ,direc
         INTEGER , DIMENSION(3)               :: ldu3,ldl3
#if __MESH_DIM == __3D
         INTEGER                              :: dir1,dir2,jj,iz
         INTEGER , DIMENSION(4)               :: ldu4,ldl4
#endif
         INTEGER , DIMENSION(ppm_dim)         :: Nml
         REAL(MK), DIMENSION(ppm_dim)         :: min_phys,max_phys
         REAL(MK), DIMENSION(:,:), POINTER    :: min_sub
         REAL(MK), DIMENSION(:,:), POINTER    :: max_sub
         INTEGER                              :: iopt,topoid
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
      REAL(MK),DIMENSION(:,:),POINTER :: terr
#elif __MESH_DIM == __3D
      REAL(MK),DIMENSION(:,:,:),POINTER :: tuc
      REAL(MK),DIMENSION(:,:,:),POINTER :: terr
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
      REAL(MK),DIMENSION(:,:,:),POINTER :: tuc
      REAL(MK),DIMENSION(:,:,:),POINTER :: terr
#elif __MESH_DIM == __3D
      REAL(MK),DIMENSION(:,:,:,:),POINTER :: tuc
      REAL(MK),DIMENSION(:,:,:,:),POINTER :: terr
#endif
#endif
         !----------------------------------------------------------------------
         !  Externals
         !----------------------------------------------------------------------

         !----------------------------------------------------------------------
         !  Initialize
         !----------------------------------------------------------------------
         CALL substart('ppm_mg_init',t0,info)
         !----------------------------------------------------------------------
         !  Check arguments
         !----------------------------------------------------------------------
           IF (ppm_debug.GT.0) THEN
#if __DIM == __VFIELD
           IF (lda.LE.0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_poiss_mg_init',  &
      &            'lda must be >0',__LINE__,info)
               GOTO 9999
           ENDIF
#endif
         ENDIF
         !----------------------------------------------------------------------
         ! Definition of necessary variables and allocation of arrays
         !----------------------------------------------------------------------
#if __DIM == __SFIELD
         vecdim = 1
#elif __DIM == __VFIELD
         vecdim = lda
#endif
         w_cycle=wcycle

         topoid = topo_id
         meshid = mesh_id
         topo => ppm_topo(topo_id)%t
         mesh => topo%mesh(mesh_id)
         nsubs  = topo%nsublist
#if    __KIND == __SINGLE_PRECISION
         min_phys(:)=topo%min_physs(:)
         max_phys(:)=topo%max_physs(:)
         min_sub => topo%min_subs(:,:)
         max_sub => topo%max_subs(:,:)
         omega_s=omega
         lmyeps=ppm_myepss
#elif  __KIND == __DOUBLE_PRECISION
         min_phys(:)=topo%min_physd(:)
         max_phys(:)=topo%max_physd(:)
         min_sub => topo%min_subd(:,:)
         max_sub => topo%max_subd(:,:)
         omega_d=omega
         lmyeps=ppm_myepsd
#endif
#if __MESH_DIM == __2D
         Nml(1) = mesh%Nm(1)
         Nml(2) = mesh%Nm(2)
         maxlev = INT(log10(Nml(1)*Nml(2)*REAL(ppm_nproc,MK))/log10(2.0_MK))
         IF (maxlev.GT.limlev) THEN
            maxlev=limlev
         ENDIF
#if __KIND == __SINGLE_PRECISION
         dx_s = (max_phys(1)-min_phys(1))/REAL((Nml(1)-1),MK)
         dy_s = (max_phys(2)-min_phys(2))/REAL((Nml(2)-1),MK)
         rdx2_s  = 1.0_MK/(dx_s*dx_s)
         rdy2_s  = 1.0_MK/(dy_s*dy_s)
#elif __KIND == __DOUBLE_PRECISION
         dx_d = (max_phys(1)-min_phys(1))/REAL((Nml(1)-1),MK)
         dy_d = (max_phys(2)-min_phys(2))/REAL((Nml(2)-1),MK)

         rdx2_d  = 1.0_MK/(dx_d*dx_d)
         rdy2_d  = 1.0_MK/(dy_d*dy_d)

#endif
#elif __MESH_DIM == __3D
         Nml(1) = mesh%Nm(1)
         Nml(2) = mesh%Nm(2)
         Nml(3) = mesh%Nm(3)
         maxlev = INT(log10(Nml(1)*Nml(2)*Nml(3)* &
      &           REAL(ppm_nproc,MK))/log10(2.0_MK))

         IF (maxlev.GT.limlev) THEN
          maxlev=limlev
         ENDIF
#if __KIND == __SINGLE_PRECISION
         dx_s = (max_phys(1)-min_phys(1))/REAL((Nml(1)-1),MK)
         dy_s = (max_phys(2)-min_phys(2))/REAL((Nml(2)-1),MK)
         dz_s = (max_phys(3)-min_phys(3))/REAL((Nml(3)-1),MK)
         rdx2_s = 1.0_MK/(dx_s*dx_s)
         rdy2_s = 1.0_MK/(dy_s*dy_s)
         rdz2_s = 1.0_MK/(dz_s*dz_s)
#elif __KIND == __DOUBLE_PRECISION
         dx_d = (max_phys(1)-min_phys(1))/REAL((Nml(1)-1),MK)
         dy_d = (max_phys(2)-min_phys(2))/REAL((Nml(2)-1),MK)
         dz_d = (max_phys(3)-min_phys(3))/REAL((Nml(3)-1),MK)
         rdx2_d = 1.0_MK/(dx_d*dx_d)
         rdy2_d = 1.0_MK/(dy_d*dy_d)
         rdz2_d = 1.0_MK/(dz_d*dz_d)
#endif
#endif
#if __DIM == __SFIELD
         iopt = ppm_param_alloc_fit
         ldu2(1) = nsubs
         ldu2(2) = 2*ppm_dim
         CALL ppm_alloc(bcdef_sca,ldu2,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
      &                   'Boundary condiotions',__LINE__,info)
            GOTO 9999
         ENDIF
         bcdef_sca(:,:)=0
         DO isub=1,nsubs
          idom=topo%isublist(isub)
          !---------------------------------------------------------------------
          !  compare the west boundary
          !---------------------------------------------------------------------
          IF (ABS(min_sub(1,idom)-min_phys(1)) .LT. &
      &       lmyeps*(max_sub(1,idom)-min_sub(1,idom))) THEN
             bcdef_sca(isub,1)=ibcdef(1)
          ENDIF
          !---------------------------------------------------------------------
          !  compare the east boundary
          !---------------------------------------------------------------------
          IF (ABS(max_sub(1,idom)-max_phys(1)) .LT. &
      &       lmyeps*(max_sub(1,idom)-min_sub(1,idom))) THEN
             bcdef_sca(isub,2)=ibcdef(2)
          ENDIF
          !---------------------------------------------------------------------
          !  compare the south boundary
          !---------------------------------------------------------------------
          IF (ABS(min_sub(2,idom)-min_phys(2)) .LT. &
      &       lmyeps*(max_sub(2,idom)-min_sub(2,idom))) THEN
             bcdef_sca(isub,3)=ibcdef(3)
          ENDIF
          !---------------------------------------------------------------------
          !  compare the north boundary
          !---------------------------------------------------------------------
          IF (ABS(max_sub(2,idom)-max_phys(2)) .LT. &
      &       lmyeps*(max_sub(2,idom)-min_sub(2,idom))) THEN
             bcdef_sca(isub,4)=ibcdef(4)
          ENDIF
#if __MESH_DIM == __3D
          !-----------------------------------------------------------------
          !  compare the bottom boundary
          !---------------------------------------------------------------------
          IF (ABS(min_sub(3,idom)-min_phys(3)) .LT. &
      &       lmyeps*(max_sub(3,idom)-min_sub(3,idom))) THEN
             bcdef_sca(isub,5)=ibcdef(5)
          ENDIF
          !---------------------------------------------------------------------
          !  compare the top boundary
          !---------------------------------------------------------------------
          IF (ABS(max_sub(3,idom)-max_phys(3)) .LT. &
      &       lmyeps*(max_sub(3,idom)-min_sub(3,idom))) THEN
             bcdef_sca(isub,6)=ibcdef(6)
          ENDIF
#endif
         ENDDO
         lperiodic=.TRUE.
         DO isub=1,nsubs
           DO i=1,2*ppm_dim
             IF (bcdef_sca(isub,i).NE.ppm_param_bcdef_periodic) THEN
               lperiodic=.FALSE.
               EXIT
             ENDIF
           ENDDO
         ENDDO
#elif __DIM == __VFIELD
         iopt = ppm_param_alloc_fit
         ldu3(1) = vecdim
         ldu3(2) = nsubs
         ldu3(3) = 2*ppm_dim
         CALL ppm_alloc(bcdef_vec,ldu3,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
      &                   'Boundary condiotions',__LINE__,info)
            GOTO 9999
         ENDIF
         bcdef_vec(:,:,:)=0
         DO isub=1,nsubs
           idom=topo%isublist(isub)
           DO ilda=1,vecdim
           !------------------------------------------------------------------
           !  compare the west boundary
           !---------------------------------------------------------------------
           IF (ABS(min_sub(1,idom)-min_phys(1)) .LT. &
      &       lmyeps*(max_sub(1,idom)-min_sub(1,idom))) THEN
             bcdef_vec(ilda,isub,1)=ibcdef(ilda,1)
           ENDIF
           !---------------------------------------------------------------------
           !  compare the east boundary
           !---------------------------------------------------------------------
           IF (ABS(max_sub(1,idom)-max_phys(1)) .LT. &
       &       lmyeps*(max_sub(1,idom)-min_sub(1,idom))) THEN
             bcdef_vec(ilda,isub,2)=ibcdef(ilda,2)
          ENDIF
          !---------------------------------------------------------------------
          !  compare the south boundary
          !---------------------------------------------------------------------
          IF (ABS(min_sub(2,idom)-min_phys(2)) .LT. &
      &       lmyeps*(max_sub(2,idom)-min_sub(2,idom))) THEN
             bcdef_vec(ilda,isub,3)=ibcdef(ilda,3)
          ENDIF
          !---------------------------------------------------------------------
          !  compare the north boundary
          !---------------------------------------------------------------------
          IF (ABS(max_sub(2,idom)-max_phys(2)) .LT. &
      &       lmyeps*(max_sub(2,idom)-min_sub(2,idom))) THEN
             bcdef_vec(ilda,isub,4)=ibcdef(ilda,4)
          ENDIF
#if __MESH_DIM == __3D
          !-----------------------------------------------------------------
          !  compare the bottom boundary
          !---------------------------------------------------------------------
          IF (ABS(min_sub(3,idom)-min_phys(3)) .LT. &
      &       lmyeps*(max_sub(3,idom)-min_sub(3,idom))) THEN
             bcdef_vec(ilda,isub,5)=ibcdef(ilda,5)
          ENDIF
          !---------------------------------------------------------------------
          !  compare the top boundary
          !---------------------------------------------------------------------
          IF (ABS(max_sub(3,idom)-max_phys(3)) .LT. &
      &       lmyeps*(max_sub(3,idom)-min_sub(3,idom))) THEN
             bcdef_vec(ilda,isub,6)=ibcdef(ilda,6)
          ENDIF
#endif
         enddo
         enddo
         lperiodic=.TRUE.
         Do isub=1,nsubs
           DO i=1,2*ppm_dim
            DO ilda=1,vecdim
             IF (bcdef_vec(ilda,isub,i).NE.ppm_param_bcdef_periodic) Then
                 lperiodic=.FALSE.
                 EXIT
             ENDIF
            ENDDO
           ENDDO
         ENDDO
#endif
        !-----------------------------------------------------------------------
        ! Allocation of the ghostsize
        !-----------------------------------------------------------------------
         iopt = ppm_param_alloc_fit
         ldu1(1) = ppm_dim
         CALL ppm_alloc(ghostsize,ldu1,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
      &                   'ghostsize',__LINE__,info)
            GOTO 9999
         ENDIF
         ghostsize=ighostsize
         !----------------------------------------------------------------------
         ! Allocation of the factor for coarsening (later set to 2)
         !----------------------------------------------------------------------
         iopt = ppm_param_alloc_fit
         ldu1(1) = ppm_dim
         CALL ppm_alloc(factor,ldu1,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
      &                   'factor',__LINE__,info)
            GOTO 9999
         ENDIF
         factor(:) = 2
         !----------------------------------------------------------------------
         ! IDs for the meshes on the different levels
         !----------------------------------------------------------------------
         iopt = ppm_param_alloc_fit
         ldu1(1) = maxlev
         CALL ppm_alloc(mg_meshid,ldu1,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
       &                  'mg_meshid',__LINE__,info)
            GOTO 9999
         ENDIF

         !----------------------------------------------------------------------
         ! Allocating the start index for the iteration through the mesh points.
         !----------------------------------------------------------------------
         iopt = ppm_param_alloc_fit
         ldu3(1) = ppm_dim
         ldu3(2) = nsubs
         ldu3(3) = maxlev
         CALL ppm_alloc(start,ldu3,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
      &             'starting indices when updating the field',__LINE__,info)
            GOTO 9999
         ENDIF
         !----------------------------------------------------------------------
         ! Allocating the stop index for the iteration through the mesh points.
         !----------------------------------------------------------------------
         iopt = ppm_param_alloc_fit
         ldu3(1) = ppm_dim
         ldu3(2) = nsubs
         ldu3(3) = maxlev
         CALL ppm_alloc(istop,ldu3,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
      &        'istopping indices when updating the field',__LINE__,info)
            GOTO 9999
         ENDIF

         !----------------------------------------------------------------------
         ! Allocating the multigrid fields used on the different levels
         !----------------------------------------------------------------------
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
         iopt = ppm_param_alloc_fit
         ldu2(1) = nsubs
         ldu2(2) = maxlev
         CALL ppm_mg_alloc(mgfield_2d_sca_s,ldu2,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
      &        'Multigrid fields used on the different levels',__LINE__,info)
            GOTO 9999
         ENDIF
         mgfield => mgfield_2d_sca_s
#elif __KIND == __DOUBLE_PRECISION
         iopt = ppm_param_alloc_fit
         ldu2(1) = nsubs
         ldu2(2) = maxlev
         CALL ppm_mg_alloc(mgfield_2d_sca_d,ldu2,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
      &        'Multigrid fields used on the different levels',__LINE__,info)
            GOTO 9999
         ENDIF
         mgfield => mgfield_2d_sca_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
         iopt = ppm_param_alloc_fit
         ldu2(1) = nsubs
         ldu2(2) = maxlev
         CALL ppm_mg_alloc(mgfield_3d_sca_s,ldu2,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
                 &        'Multigrid fields used on the different levels',__LINE__,info)
            GOTO 9999
         ENDIF
         mgfield => mgfield_3d_sca_s
#elif __KIND == __DOUBLE_PRECISION
         iopt = ppm_param_alloc_fit
         ldu2(1) = nsubs
         ldu2(2) = maxlev
         CALL ppm_mg_alloc(mgfield_3d_sca_d,ldu2,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
      &        'Multigrid fields used on the different levels',__LINE__,info)
            GOTO 9999
         ENDIF
         mgfield => mgfield_3d_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
         iopt = ppm_param_alloc_fit
         ldu2(1) = nsubs
         ldu2(2) = maxlev
         CALL ppm_mg_alloc(mgfield_2d_vec_s,ldu2,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
      &        'Multigrid fields used on the different levels',__LINE__,info)
            GOTO 9999
         ENDIF
         mgfield => mgfield_2d_vec_s
#elif __KIND == __DOUBLE_PRECISION
         iopt = ppm_param_alloc_fit
         ldu2(1) = nsubs
         ldu2(2) = maxlev
         CALL ppm_mg_alloc(mgfield_2d_vec_d,ldu2,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
      &        'Multigrid fields used on the different levels',__LINE__,info)
            GOTO 9999
         ENDIF
         mgfield => mgfield_2d_vec_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
         iopt = ppm_param_alloc_fit
         ldu2(1) = nsubs
         ldu2(2) = maxlev
         CALL ppm_mg_alloc(mgfield_3d_vec_s,ldu2,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
         &      'Multigrid fields used on the different levels',__LINE__,info)
            GOTO 9999
         ENDIF
         mgfield => mgfield_3d_vec_s
#elif __KIND == __DOUBLE_PRECISION
         iopt = ppm_param_alloc_fit
         ldu2(1) = nsubs
         ldu2(2) = maxlev
         CALL ppm_mg_alloc(mgfield_3d_vec_d,ldu2,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
      &        'Multigrid fields used on the different levels',__LINE__,info)
            GOTO 9999
         ENDIF
         mgfield => mgfield_3d_vec_d
#endif
#endif
#endif
         iopt = ppm_param_alloc_fit
         ldu2(1) = 2*ppm_dim
         ldu2(2) = nsubs
         CALL ppm_alloc(lboundary,ldu2,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
      &        'Problem with the boundary alloc.',__LINE__,info)
            GOTO 9999
         ENDIF

         iopt = ppm_param_alloc_fit
         ldu2(1) = ppm_dim
         ldu2(2) = maxlev
         CALL ppm_alloc(max_node,ldu2,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
      &        'Problem with a maximum number alloc.',__LINE__,info)
            GOTO 9999
         ENDIF
         ldu3(1) = ppm_dim
         ldu3(2) = nsubs
         ldu3(3) = maxlev
         CALL ppm_alloc(istart,ldu3,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
      &        'Problem with istart alloc.',__LINE__,info)
            GOTO 9999
         ENDIF
         max_node(:,:) = 0

         lboundary(:,:) = .FALSE.
         start(:,:,:) = 1
         !----------------------------------------------------------------------
         ! Derive coarser meshes
         !----------------------------------------------------------------------
        DO mlev=1,maxlev
#if __MESH_DIM == __2D
           !-------------------------------------------------------------------
           ! Go through the subs, define the istopping indices on each mesh,
           ! check and store if it is on the boundary, allocate the
           ! multigrid fields and pass the boundary values.
           !-------------------------------------------------------------------
           DO i=1,nsubs
              idom=topo%isublist(i)
              istop(:,i,mlev)= mesh%nnodes(:,idom)
              istart(:,i,mlev) = mesh%istart(:,isub)
              DO j=1,ppm_dim
                 IF (max_node(j,mlev).LT.istop(j,i,mlev)) THEN
                    max_node(j,mlev)=istop(j,i,mlev)
                 ENDIF
              ENDDO
              !----------------------------------------------------------------
              ! Allocate the function correction, the restricted errors,
              ! the residuals and the values on the boundary on each level.
              !----------------------------------------------------------------
#if __DIM == __SFIELD
              iopt = ppm_param_alloc_fit
              ldl2(1) = 1-ghostsize(1)
              ldl2(2) = 1-ghostsize(2)
              ldu2(1) = mesh%nnodes(1,idom)+ghostsize(1)
              ldu2(2) = mesh%nnodes(2,idom)+ghostsize(2)
              CALL ppm_alloc(mgfield(i,mlev)%uc,ldl2,ldu2,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the function corr. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF
              tuc => mgfield(i,mlev)%uc
              tuc = 0.0_MK
              iopt = ppm_param_alloc_fit
              ldu2(1) = mesh%nnodes(1,idom)
              ldu2(2) = mesh%nnodes(2,idom)
              CALL ppm_alloc(mgfield(i,mlev)%fc,ldu2,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
         &        'Problem with the restricted err. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF
              mgfield(i,mlev)%fc(:,:)=0.0_MK
              iopt = ppm_param_alloc_fit
              ldl2(1) = 1-ghostsize(1)
              ldl2(2) = 1-ghostsize(2)
              ldu2(1) = mesh%nnodes(1,idom)+ghostsize(1)
              ldu2(2) = mesh%nnodes(2,idom)+ghostsize(2)
              CALL ppm_alloc(mgfield(i,mlev)%err,ldl2,ldu2,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &           'Problem with the residual alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF
              terr => mgfield(i,mlev)%err
              terr(:,:)=0.0_MK
              ! ALLOCATE THE BCVALUE(IT IS A TYPE!!)
              IF (.NOT.lperiodic) THEN
                 iopt = ppm_param_alloc_fit
                 ldu1(1) = 2*ppm_dim
                 CALL ppm_mg_alloc(mgfield(i,mlev)%bcvalue,ldu1,iopt,info)
                 IF (info .NE. 0) THEN
                    info = ppm_error_fatal
                    CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &              'Problem with the BOUNDARY alloc.',__LINE__,info)
                    GOTO 9999
                 ENDIF
                 !ALLOCATE THE PBCVALUE
                 DO iface = 1,2*ppm_dim
                    iopt = ppm_param_alloc_fit
                    IF (iface.EQ.1.OR.iface.EQ.2) THEN
                       ldu1(1) = max_node(2,mlev)
                    ELSE
                       ldu1(1) = max_node(1,mlev)
                    ENDIF
                    CALL ppm_alloc(mgfield(i,mlev)%bcvalue(iface)%pbcvalue,&
                    &              ldu1,iopt,info)
                    IF (info .NE. 0) THEN
                       info = ppm_error_fatal
                       CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
      &                'Problem with the BOUNDARY alloc.',__LINE__,info)
                       GOTO 9999
                    ENDIF
                 ENDDO !iface
                 DO iface=1,2*ppm_dim
                 IF (iface.EQ.1.OR.iface.EQ.2) THEN
                     direc(1)=2
                 ELSEIF (iface.EQ.3.OR.iface.EQ.4) then
                     direc(1)=1
                 ENDIF
                 DO ipoint=1,max_node(direc(1),mlev)
                    IF (mlev.EQ.1) THEN
                       mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ipoint) = &
                    &                  bcvalue(i,iface,ipoint)
                    ELSE
                       IF(bcdef_sca(i,iface).EQ.ppm_param_bcdef_neumann) THEN
                           mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ipoint) = &
                       &   mgfield(i,mlev-1)%bcvalue(iface)%pbcvalue(2*ipoint-1)
                       ELSE
                       ! NO CORRECTIONS FOR THE DIRICHLET
                          mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ipoint)=0.0_MK
                       ENDIF
                    ENDIF
                 ENDDO !ipoint
                 ENDDO !faces
              ENDIF!lperiodic
#elif __DIM == __VFIELD
              iopt = ppm_param_alloc_fit
              ldl3(1) = 1
              ldl3(2) = 1-ghostsize(1)
              ldl3(3) = 1-ghostsize(2)
              ldu3(1) = vecdim
              ldu3(2) = mesh%nnodes(1,idom)+ghostsize(1)
              ldu3(3) = mesh%nnodes(2,idom)+ghostsize(2)
              CALL ppm_alloc(mgfield(i,mlev)%uc,ldl3,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
              &       'Problem with the function corr. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF
              tuc => mgfield(i,mlev)%uc
              tuc = 0.0_MK
              iopt = ppm_param_alloc_fit
              ldu3(1) = vecdim
              ldu3(2) = mesh%nnodes(1,idom)
              ldu3(3) = mesh%nnodes(2,idom)
              CALL ppm_alloc(mgfield(i,mlev)%fc,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
              &       'Problem with the restricted err. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF
              mgfield(i,mlev)%fc(:,:,:)=0.0_MK
              iopt = ppm_param_alloc_fit
              ldl3(1) = 1
              ldl3(2) = 1-ghostsize(1)
              ldl3(3) = 1-ghostsize(2)
              ldu3(1) = vecdim
              ldu3(2) = mesh%nnodes(1,idom)+ghostsize(1)
              ldu3(3) = mesh%nnodes(2,idom)+ghostsize(2)
              CALL ppm_alloc(mgfield(i,mlev)%err,ldl3,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
              & 'Problem with the residual alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF
              terr=>mgfield(i,mlev)%err
              terr(:,:,:)=0.0_MK
#endif
           ENDDO!DO 1,nsubs
#elif __MESH_DIM == __3D
           DO i=1,nsubs
              idom=topo%isublist(i)
              istop(:,i,mlev) = mesh%nnodes(:,idom)
              istart(:,i,mlev) = mesh%istart(:,isub)
              DO j=1,ppm_dim
                 IF (max_node(j,mlev).LT.istop(j,i,mlev)) THEN
                    max_node(j,mlev)=istop(j,i,mlev)
                 ENDIF
              ENDDO
              IF (topo%subs_bc(1,idom).EQ.1) THEN
                 lboundary(1,i)=.TRUE.
              ELSEIF (topo%subs_bc(3,idom).EQ.1) THEN
                 lboundary(3,i)=.TRUE.
              ELSEIF (topo%subs_bc(2,idom).EQ.1) THEN
                 lboundary(2,i)=.TRUE.
              ELSEIF (topo%subs_bc(4,idom).EQ.1) THEN
                 lboundary(4,i)=.TRUE.
              ELSEIF (topo%subs_bc(5,idom).EQ.1) THEN
                 lboundary(5,i)=.TRUE.
              ELSEIF (topo%subs_bc(6,idom).EQ.1) THEN
                 lboundary(6,i)=.TRUE.
              ENDIF
              !----------------------------------------------------------------
              ! Allocate the function correction, the restricted errors and the
              !residuals on each level.
              !----------------------------------------------------------------
#if __DIM == __SFIELD
              iopt = ppm_param_alloc_fit
              ldl3(1) = 1-ghostsize(1)
              ldl3(2) = 1-ghostsize(2)
              ldl3(3) = 1-ghostsize(3)
              ldu3(1) = mesh%nnodes(1,idom)+ghostsize(1)
              ldu3(2) = mesh%nnodes(2,idom)+ghostsize(2)
              ldu3(3) = mesh%nnodes(3,idom)+ghostsize(3)
              CALL ppm_alloc(mgfield(i,mlev)%uc,ldl3,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
      &          'Problem with the function corr. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF
              tuc=>mgfield(i,mlev)%uc
              tuc=0.0_MK
              iopt = ppm_param_alloc_fit
              ldu3(1) = mesh%nnodes(1,idom)
              ldu3(2) = mesh%nnodes(2,idom)
              ldu3(3) = mesh%nnodes(3,idom)
              CALL ppm_alloc(mgfield(i,mlev)%fc,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
      &          'Problem with the restricted err. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF
              mgfield(i,mlev)%fc=0.0_MK
              iopt = ppm_param_alloc_fit
              ldl3(1) = 1-ghostsize(1)
              ldl3(2) = 1-ghostsize(2)
              ldl3(3) = 1-ghostsize(3)
              ldu3(1) = mesh%nnodes(1,idom)+ghostsize(1)
              ldu3(2) = mesh%nnodes(2,idom)+ghostsize(2)
              ldu3(3) = mesh%nnodes(3,idom)+ghostsize(3)
              CALL ppm_alloc(mgfield(i,mlev)%err,ldl3,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
      &          'Problem with the residual alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF
              terr=>mgfield(i,mlev)%err
              terr=0.0_MK
              !ALLOCATE THE BCVALUE(IT IS A TYPE!!)
              IF (.NOT.lperiodic) THEN
                 iopt = ppm_param_alloc_fit
                 ldu1(1) = 2*ppm_dim
                 CALL ppm_mg_alloc(mgfield(i,mlev)%bcvalue,ldu1,iopt,info)
                 IF (info .NE. 0) THEN
                    info = ppm_error_fatal
                    CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
      &        'Problem with the BOUNDARY alloc.',__LINE__,info)
                    GOTO 9999
                 ENDIF
                 !ALLOCATE THE PBCVALUE
                 DO iface=1,2*ppm_dim
                    iopt = ppm_param_alloc_fit
                    IF (iface.EQ.1.OR.iface.EQ.2) THEN
                       ldu2(1) = max_node(2,mlev)
                       ldu2(2)= max_node(3,mlev)
                    ELSEIF (iface.EQ.3.OR. iface.EQ.4) then
                       ldu2(1) = max_node(1,mlev)
                       ldu2(2)=max_node(3,mlev)
                    ELSE
                       ldu2(1)=max_node(1,mlev)
                       ldu2(2)=max_node(2,mlev)
                    ENDIF
                    CALL ppm_alloc(mgfield(i,mlev)%bcvalue(iface)%pbcvalue,&
                 &                 ldu2,iopt,info)
                    IF (info .NE. 0) THEN
                       info = ppm_error_fatal
                       CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
      &                 'Problem with the BOUNDARY alloc.',__LINE__,info)
                       GOTO 9999
                    ENDIF
                 ENDDO
                 DO iface=1,2*ppm_dim
                    IF (iface.EQ.1.OR.iface.EQ.2) THEN
                       direc(1)=2
                       direc(2)=3
                    ELSEIF (iface.EQ.3.OR.iface.EQ.4) THEN
                       direc(1)=1
                       direc(2)=3
                    ELSE
                       direc(1)=1
                       direc(2)=2
                    ENDIF
                    DO ipoint=1,max_node(direc(1),mlev)
                    DO jpoint=1,max_node(direc(2),mlev)
                       IF (mlev.EQ.1) THEN
                          mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ipoint,jpoint) =&
                          &               bcvalue(i,iface,ipoint,jpoint)
                       ELSE
                          IF(bcdef_sca(i,iface).EQ.ppm_param_bcdef_neumann) THEN
                             mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ipoint,jpoint)=&
                   &         mgfield(i,mlev-1)%bcvalue(iface)%pbcvalue(2*ipoint-1,2*jpoint-1)
                          ELSE
                       !NO CORRECTIONS FOR THE DIRICHLET
                             mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ipoint,jpoint)=&
                             &            0.0_MK
                          ENDIF
                       ENDIF
                    ENDDO
                    ENDDO
                 ENDDO!faces
              ENDIF !lperiodic
#elif __DIM == __VFIELD
              iopt = ppm_param_alloc_fit
              ldl4(1) = 1
              ldl4(2) = 1-ghostsize(1)
              ldl4(3) = 1-ghostsize(2)
              ldl4(4) = 1-ghostsize(3)
              ldu4(1) = vecdim
              ldu4(2) = mesh%nnodes(1,idom)+ghostsize(1)
              ldu4(3) = mesh%nnodes(2,idom)+ghostsize(2)
              ldu4(4) = mesh%nnodes(3,idom)+ghostsize(3)
              CALL ppm_alloc(mgfield(i,mlev)%uc,ldl4,ldu4,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
         &            'Problem with the function corr. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF
              tuc=>mgfield(i,mlev)%uc
              tuc=0.0_MK
              iopt = ppm_param_alloc_fit
              ldu4(1) = vecdim
              ldu4(2) = mesh%nnodes(1,idom)
              ldu4(3) = mesh%nnodes(2,idom)
              ldu4(4) = mesh%nnodes(3,idom)
              CALL ppm_alloc(mgfield(i,mlev)%fc,ldu4,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
                 &        'Problem with the restricted err. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF
              mgfield(i,mlev)%fc=0.0_MK
              iopt = ppm_param_alloc_fit
              ldl4(1) = 1
              ldl4(2) = 1-ghostsize(1)
              ldl4(3) = 1-ghostsize(2)
              ldl4(4) = 1-ghostsize(3)
              ldu4(1) = vecdim
              ldu4(2) = mesh%nnodes(1,idom)+ghostsize(1)
              ldu4(3) = mesh%nnodes(2,idom)+ghostsize(2)
              ldu4(4) = mesh%nnodes(3,idom)+ghostsize(3)
              CALL ppm_alloc(mgfield(i,mlev)%err,ldl4,ldu4,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
                 &        'Problem with the residual alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF
              terr=>mgfield(i,mlev)%err
              terr=0.0_MK
              !ALLOCATE THE BCVALUE(IT IS A TYPE!!)
              IF (.NOT.lperiodic) THEN
                 iopt = ppm_param_alloc_fit
                 ldu1=2*ppm_dim
                 CALL ppm_mg_alloc(mgfield(i,mlev)%bcvalue,ldu1,iopt,info)
                 IF (info .NE. 0) THEN
                    info = ppm_error_fatal
                     CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
                     &        'Problem with the BOUNDARY alloc.',__LINE__,info)
                     GOTO 9999
                 ENDIF
                 !ALLOCATE THE PBCVALUE
                 DO iface=1,2*ppm_dim
                    iopt = ppm_param_alloc_fit
                    ldu3(1)=vecdim
                    IF (iface.EQ.1.OR.iface.EQ.2) THEN
                       ldu3(2) = max_node(2,mlev)
                       ldu3(3)= max_node(3,mlev)
                    ELSEIF (iface.EQ.3.OR. iface.EQ.4) then
                       ldu3(2) = max_node(1,mlev)
                       ldu3(3)=max_node(3,mlev)
                    ELSE
                       ldu3(2)=max_node(1,mlev)
                       ldu3(3)=max_node(2,mlev)
                    ENDIF
                    CALL ppm_alloc(mgfield(i,mlev)%bcvalue(iface)%pbcvalue,ldu3,iopt,info)
                    IF (info .NE. 0) THEN
                       info = ppm_error_fatal
                       CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
                       &        'Problem with the BOUNDARY alloc.',__LINE__,info)
                       GOTO 9999
                    ENDIF
                 ENDDO
                 DO iface=1,2*ppm_dim
                    IF (iface.EQ.1.OR.iface.EQ.2) THEN
                       direc(1)=2
                       direc(2)=3
                    ELSEIF (iface.EQ.3.OR.iface.EQ.4) THEN
                       direc(1)=1
                       direc(2)=3
                    ELSE
                       direc(1)=1
                       direc(2)=2
                    ENDIF
                    DO ipoint=1,max_node(direc(1),mlev)
                    DO jpoint=1,max_node(direc(2),mlev)
                    DO ilda=1,vecdim
                       IF (mlev.EQ.1) THEN
                          mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ilda,ipoint,jpoint) &
        &                                =bcvalue(ilda,i,iface,ipoint,jpoint)
                       ELSE
                       IF(bcdef_vec(ilda,i,iface).EQ.ppm_param_bcdef_neumann) THEN
                          mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ilda,ipoint,jpoint)=&
        &                 mgfield(i,mlev-1)%bcvalue(iface)%pbcvalue(ilda,2*ipoint-1,2*jpoint-1)
                       ELSE
                              !NO CORRECTIONS FOR THE DIRICHLET
                          mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ilda,ipoint,jpoint)=0.0_MK
                       ENDIF
                       ENDIF
                    ENDDO
                    ENDDO
                    ENDDO
                 ENDDO
              ENDIF !lperiodic
#endif
           ENDDO!DO i=1,nsubs
#endif
           mg_meshid(mlev)=meshid
           newmeshid=-1
           IF (mlev.LT.maxlev) THEN
             CALL ppm_mesh_derive(topoid,meshid,newmeshid,&
     &                          ppm_param_mesh_coarsen,factor,info)
             meshid = newmeshid
             mesh => topo%mesh(meshid)
           ENDIF
         ENDDO!DO mlev=1,maxlev
         !----------------------------------------------------------------------
         !  Return
         !----------------------------------------------------------------------
 9999    CONTINUE
         CALL substop('ppm_mg_init',t0,info)
         RETURN
#if    __DIM       == __SFIELD
#if    __MESH_DIM  == __2D
#if    __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_init_2d_sca_s
#elif  __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_init_2d_sca_d
#endif
#elif  __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_init_3d_sca_s
#elif  __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_init_3d_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if    __MESH_DIM  == __2D
#if    __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_init_2d_vec_s
#elif  __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_init_2d_vec_d
#endif
#elif  __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_init_3d_vec_s
#elif  __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_init_3d_vec_d
#endif
#endif
#endif
