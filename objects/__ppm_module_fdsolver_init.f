




      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_fdsolver_solve
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for field solver 
      !                        routines.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
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












   MODULE ppm_module_fdsolver_init

         !----------------------------------------------------------------------
         !  Define interface to ppm_fdsolver_init
         !----------------------------------------------------------------------
         INTERFACE ppm_fdsolver_init
            MODULE PROCEDURE ppm_fdsolver_init_2d_sca_s
            MODULE PROCEDURE ppm_fdsolver_init_2d_sca_d
            MODULE PROCEDURE ppm_fdsolver_init_2d_vec_s
            MODULE PROCEDURE ppm_fdsolver_init_2d_vec_d
            MODULE PROCEDURE ppm_fdsolver_init_3d_sca_s
            MODULE PROCEDURE ppm_fdsolver_init_3d_sca_d
            MODULE PROCEDURE ppm_fdsolver_init_3d_vec_s
            MODULE PROCEDURE ppm_fdsolver_init_3d_vec_d
         END INTERFACE
         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_fdsolver_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : 
      !                initializes fieldsolver by creating all FFT-plans apriori
      !
      !  Input        :  
      !                  DATA_fv(:,:,:,:) (F) field data
      !                  lda_fv      (I) size of leading dimension in vector
      !                                  case
      !                  mesh_id_data(I) mesh ID of the current data field mesh
      !                  topo_ids(2) (I) topology IDs on which the FFTs are 
      !                                  performed
      !                                  topo_ids(1) initial topology(xpencils)
      !                                  topo_ids(2) second  topology(ypencils)
      !                                  topo_ids(3) third   topology(zpencils)
      !                                         (3D only!!)
      !                                 
      !                  mesh_ids(3) (I) mesh IDs where the FFTs are performed
      !                                  mesh_ids(1) first mesh 
      !                                         (xpencils,real
      !                                  mesh_ids(2) second mesh 
      !                                         (xpencils,complex)
      !                                  mesh_ids(3) third mesh 
      !                                         (ypencils,complex)
      !                                  mesh_ids(4) forth mesh 
      !                                         (zpencils,complex) (3D only!!)
      !                  ghostsize(3) (I)ghostsize
      !
      !  Input/output : 
      !                                
      !                                
      !
      !  Output       :   
      !                 info          (I) return status. =0 if no error.
      !
      !  Remarks      : 
      !                                                  
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fdsolver_init.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.13  2006/09/04 18:34:43  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.11  2006/04/10 08:54:29  pchatela
      !  Made xp a REAL, DIMENSION(:,:), POINTER to get rid of warnings
      !
      !  Revision 1.10  2006/04/07 17:41:04  hiebers
      !  Changed type of variable xp to POINTER
      !
      !  Revision 1.9  2005/11/29 10:56:58  hiebers
      !  Changed arguments for fftwRoutines from INTEGER to INTEGER,DIMENSION(1)
      !  for  NAG compiler
      !
      !  Revision 1.8  2005/06/04 00:40:23  michaebe
      !  Cosmetics of cosmetics?
      !
      !  Revision 1.7  2005/06/04 00:37:49  michaebe
      !  cosmetics
      !
      !  Revision 1.6  2005/06/04 00:36:12  michaebe
      !  __ppm_module_fdsolver_init.f, line 491.25: 1516-023 (E) Subscript is
      !  out of bounds.
      !  This compiler warning eventually kinda got on my nerves so I acted.
      !
      !  Revision 1.5  2005/02/19 07:32:31  ivos
      !  Resolved CVS conflicts.
      !
      !  Revision 1.4  2005/02/18 08:01:55  hiebers
      !  minor changes in error messages
      !
      !  Revision 1.3  2005/02/16 22:22:59  ivos
      !  Bugfix: replaced non-existing ppm_module_data_fdsolver with
      !  ppm_module_data_fieldsolver.
      !
      !  Revision 1.2  2005/02/16 12:41:13  hiebers
      !  exchange FFTW_ESTIMATE by FFTW_MEASURE
      !
      !  Revision 1.1  2005/02/16 10:22:34  hiebers
      !  initial implementation
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fdsolver_init_2d_sca_s(DATA_fv,mesh_id_data,topo_ids, &
     &  field_topoid,mesh_ids,ghostsize, info)
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_mktopo
      USE ppm_module_typedef
      USE ppm_module_data_fieldsolver
      IMPLICIT NONE
      INTEGER, PARAMETER :: MK = ppm_kind_single
      !-------------------------------------------------------------------------
      !  FFTW include
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      ! POINTER to data
      REAL(MK), DIMENSION(:,:,:),        POINTER   :: data_fv
      ! mesh ID of the data
      INTEGER                    , INTENT(IN)      :: mesh_id_data
      ! topo ID of the field
      INTEGER                    , INTENT(IN   )   :: field_topoid
      ! topo / mesh ID for the mapping

      INTEGER, DIMENSION(2)      , INTENT(IN   )   :: topo_ids
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: mesh_ids
      INTEGER, DIMENSION(2)      , INTENT(IN   )   :: ghostsize





      INTEGER                    , INTENT(  OUT)   :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(ppm_kind_double)                   :: t0
      ! counters
      INTEGER                                 :: k, i, j
      CHARACTER(LEN=ppm_char)                 :: mesg

      REAL(MK), DIMENSION(:,:),     POINTER   :: data_real
      COMPLEX(MK), DIMENSION(:,:),  POINTER   :: data_comp, data_compl
      INTEGER, DIMENSION(2)                   :: lda





      ! size of the data_in 
      INTEGER,DIMENSION(1)                    :: MB_in
      INTEGER                                 :: Nx_in, Ny_in, Nz_in
      ! size of the data_out 
      INTEGER                                 :: Nx_out, Ny_out, Nz_out
      ! variables
      REAL(MK), DIMENSION(:,:),POINTER        :: xp
      INTEGER                                 :: Npart
      INTEGER                                 :: decomp, assign
      REAL(MK), DIMENSION(3  )                :: min_phys, max_phys
      REAL(MK), DIMENSION(3  )                :: length
      REAL(MK), DIMENSION(3  )                :: length_phys
      INTEGER , DIMENSION(6  )                :: bcdef 
      INTEGER                                 :: nsubs,topo_id, mesh_id
      INTEGER                                 :: yhmax, zhmax
      INTEGER                                 :: mesh_id_xpen, mesh_id_ypen
      INTEGER                                 :: mesh_id_xpen_complex
      INTEGER                                 :: mesh_id_zpen
      INTEGER                                 :: mesh_id_slab
      INTEGER                                 :: mesh_id_slab_complex
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      INTEGER , DIMENSION(:,:), POINTER       :: istart, istart_xpen_complex  
      INTEGER , DIMENSION(:,:), POINTER       :: istart_ypen, istart_trans
      INTEGER , DIMENSION(:,:), POINTER       :: istart_zpen
      INTEGER , DIMENSION(:,:), POINTER       :: ndata, ndata_xpen_complex
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_ypen, ndata_trans
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_zpen, ndata_slab
      INTEGER                                 :: dim, n,idom
      INTEGER                                 :: iopt
      INTEGER                                 :: topo_id_xpen, topo_id_ypen
      INTEGER                                 :: topo_id_zpen
      INTEGER                                 :: topo_id_slab
      INTEGER, DIMENSION(3)                   :: Nm, Nm_com, Nm_poisson
      INTEGER, DIMENSION(2)                   :: Nm_slab
      LOGICAL                                 :: Its_xpencil_topo
      LOGICAL                                 :: Its_xyslab_topo
      TYPE(ppm_t_topo)        , POINTER       :: f_topo
      TYPE(ppm_t_equi_mesh)   , POINTER       :: f_mesh
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_fdsolver_init',t0,info)
      f_topo => ppm_topo(field_topoid)%t
      f_mesh => f_topo%mesh(mesh_id_data)
      
      !-------------------------------------------------------------------------
      !  Error if FFTW library or NEC is not available
      !-------------------------------------------------------------------------
      info = ppm_error_error
      CALL ppm_error(ppm_err_nofftw,'ppm_fdsolver_init',  &
     &    'PPM was compiled without fftw support',279,info)
      CALL ppm_error(ppm_err_noMathKeisan,'ppm_fdsolver_init',  &
     &    'PPM was compiled without MATHKEISAN support',283,info)
      GOTO 9999   
     IF (info .NE. 0) THEN
          WRITE(mesg,'(A)') 'could not deallocate memory'
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_fdsolver_init',mesg,974,&
     &                                                                 info)
          GOTO 9999
      ENDIF
 9999 CONTINUE
      CALL substop('ppm_fdsolver_init',t0,info)
      RETURN

      END SUBROUTINE ppm_fdsolver_init_2d_sca_s 

      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_fdsolver_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : 
      !                initializes fieldsolver by creating all FFT-plans apriori
      !
      !  Input        :  
      !                  DATA_fv(:,:,:,:) (F) field data
      !                  lda_fv      (I) size of leading dimension in vector
      !                                  case
      !                  mesh_id_data(I) mesh ID of the current data field mesh
      !                  topo_ids(2) (I) topology IDs on which the FFTs are 
      !                                  performed
      !                                  topo_ids(1) initial topology(xpencils)
      !                                  topo_ids(2) second  topology(ypencils)
      !                                  topo_ids(3) third   topology(zpencils)
      !                                         (3D only!!)
      !                                 
      !                  mesh_ids(3) (I) mesh IDs where the FFTs are performed
      !                                  mesh_ids(1) first mesh 
      !                                         (xpencils,real
      !                                  mesh_ids(2) second mesh 
      !                                         (xpencils,complex)
      !                                  mesh_ids(3) third mesh 
      !                                         (ypencils,complex)
      !                                  mesh_ids(4) forth mesh 
      !                                         (zpencils,complex) (3D only!!)
      !                  ghostsize(3) (I)ghostsize
      !
      !  Input/output : 
      !                                
      !                                
      !
      !  Output       :   
      !                 info          (I) return status. =0 if no error.
      !
      !  Remarks      : 
      !                                                  
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fdsolver_init.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.13  2006/09/04 18:34:43  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.11  2006/04/10 08:54:29  pchatela
      !  Made xp a REAL, DIMENSION(:,:), POINTER to get rid of warnings
      !
      !  Revision 1.10  2006/04/07 17:41:04  hiebers
      !  Changed type of variable xp to POINTER
      !
      !  Revision 1.9  2005/11/29 10:56:58  hiebers
      !  Changed arguments for fftwRoutines from INTEGER to INTEGER,DIMENSION(1)
      !  for  NAG compiler
      !
      !  Revision 1.8  2005/06/04 00:40:23  michaebe
      !  Cosmetics of cosmetics?
      !
      !  Revision 1.7  2005/06/04 00:37:49  michaebe
      !  cosmetics
      !
      !  Revision 1.6  2005/06/04 00:36:12  michaebe
      !  __ppm_module_fdsolver_init.f, line 491.25: 1516-023 (E) Subscript is
      !  out of bounds.
      !  This compiler warning eventually kinda got on my nerves so I acted.
      !
      !  Revision 1.5  2005/02/19 07:32:31  ivos
      !  Resolved CVS conflicts.
      !
      !  Revision 1.4  2005/02/18 08:01:55  hiebers
      !  minor changes in error messages
      !
      !  Revision 1.3  2005/02/16 22:22:59  ivos
      !  Bugfix: replaced non-existing ppm_module_data_fdsolver with
      !  ppm_module_data_fieldsolver.
      !
      !  Revision 1.2  2005/02/16 12:41:13  hiebers
      !  exchange FFTW_ESTIMATE by FFTW_MEASURE
      !
      !  Revision 1.1  2005/02/16 10:22:34  hiebers
      !  initial implementation
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fdsolver_init_2d_sca_d(DATA_fv,mesh_id_data,topo_ids, &
     &  field_topoid,mesh_ids,ghostsize, info)
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_mktopo
      USE ppm_module_typedef
      USE ppm_module_data_fieldsolver
      IMPLICIT NONE
      INTEGER, PARAMETER :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  FFTW include
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      ! POINTER to data
      REAL(MK), DIMENSION(:,:,:),        POINTER   :: data_fv
      ! mesh ID of the data
      INTEGER                    , INTENT(IN)      :: mesh_id_data
      ! topo ID of the field
      INTEGER                    , INTENT(IN   )   :: field_topoid
      ! topo / mesh ID for the mapping
      INTEGER, DIMENSION(2)      , INTENT(IN   )   :: topo_ids
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: mesh_ids
      INTEGER, DIMENSION(2)      , INTENT(IN   )   :: ghostsize
      INTEGER                    , INTENT(  OUT)   :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(ppm_kind_double)                   :: t0
      ! counters
      INTEGER                                 :: k, i, j
      CHARACTER(LEN=ppm_char)                 :: mesg
      REAL(MK), DIMENSION(:,:),     POINTER   :: data_real
      COMPLEX(MK), DIMENSION(:,:),  POINTER   :: data_comp, data_compl
      INTEGER, DIMENSION(2)                   :: lda
      ! size of the data_in 
      INTEGER,DIMENSION(1)                    :: MB_in
      INTEGER                                 :: Nx_in, Ny_in, Nz_in
      ! size of the data_out 
      INTEGER                                 :: Nx_out, Ny_out, Nz_out
      ! variables
      REAL(MK), DIMENSION(:,:),POINTER        :: xp
      INTEGER                                 :: Npart
      INTEGER                                 :: decomp, assign
      REAL(MK), DIMENSION(3  )                :: min_phys, max_phys
      REAL(MK), DIMENSION(3  )                :: length
      REAL(MK), DIMENSION(3  )                :: length_phys
      INTEGER , DIMENSION(6  )                :: bcdef 
      INTEGER                                 :: nsubs,topo_id, mesh_id
      INTEGER                                 :: yhmax, zhmax
      INTEGER                                 :: mesh_id_xpen, mesh_id_ypen
      INTEGER                                 :: mesh_id_xpen_complex
      INTEGER                                 :: mesh_id_zpen
      INTEGER                                 :: mesh_id_slab
      INTEGER                                 :: mesh_id_slab_complex
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      INTEGER , DIMENSION(:,:), POINTER       :: istart, istart_xpen_complex  
      INTEGER , DIMENSION(:,:), POINTER       :: istart_ypen, istart_trans
      INTEGER , DIMENSION(:,:), POINTER       :: istart_zpen
      INTEGER , DIMENSION(:,:), POINTER       :: ndata, ndata_xpen_complex
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_ypen, ndata_trans
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_zpen, ndata_slab
      INTEGER                                 :: dim, n,idom
      INTEGER                                 :: iopt
      INTEGER                                 :: topo_id_xpen, topo_id_ypen
      INTEGER                                 :: topo_id_zpen
      INTEGER                                 :: topo_id_slab
      INTEGER, DIMENSION(3)                   :: Nm, Nm_com, Nm_poisson
      INTEGER, DIMENSION(2)                   :: Nm_slab
      LOGICAL                                 :: Its_xpencil_topo
      LOGICAL                                 :: Its_xyslab_topo
      TYPE(ppm_t_topo)        , POINTER       :: f_topo
      TYPE(ppm_t_equi_mesh)   , POINTER       :: f_mesh
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_fdsolver_init',t0,info)
      f_topo => ppm_topo(field_topoid)%t
      f_mesh => f_topo%mesh(mesh_id_data)
      
      !-------------------------------------------------------------------------
      !  Error if FFTW library or NEC is not available
      !-------------------------------------------------------------------------
      info = ppm_error_error
      CALL ppm_error(ppm_err_nofftw,'ppm_fdsolver_init',  &
     &    'PPM was compiled without fftw support',279,info)
      CALL ppm_error(ppm_err_noMathKeisan,'ppm_fdsolver_init',  &
     &    'PPM was compiled without MATHKEISAN support',283,info)
      GOTO 9999   
     IF (info .NE. 0) THEN
          WRITE(mesg,'(A)') 'could not deallocate memory'
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_fdsolver_init',mesg,974,&
     &                                                                 info)
          GOTO 9999
      ENDIF
 9999 CONTINUE
      CALL substop('ppm_fdsolver_init',t0,info)
      RETURN

      END SUBROUTINE ppm_fdsolver_init_2d_sca_d 


      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_fdsolver_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : 
      !                initializes fieldsolver by creating all FFT-plans apriori
      !
      !  Input        :  
      !                  DATA_fv(:,:,:,:) (F) field data
      !                  lda_fv      (I) size of leading dimension in vector
      !                                  case
      !                  mesh_id_data(I) mesh ID of the current data field mesh
      !                  topo_ids(2) (I) topology IDs on which the FFTs are 
      !                                  performed
      !                                  topo_ids(1) initial topology(xpencils)
      !                                  topo_ids(2) second  topology(ypencils)
      !                                  topo_ids(3) third   topology(zpencils)
      !                                         (3D only!!)
      !                                 
      !                  mesh_ids(3) (I) mesh IDs where the FFTs are performed
      !                                  mesh_ids(1) first mesh 
      !                                         (xpencils,real
      !                                  mesh_ids(2) second mesh 
      !                                         (xpencils,complex)
      !                                  mesh_ids(3) third mesh 
      !                                         (ypencils,complex)
      !                                  mesh_ids(4) forth mesh 
      !                                         (zpencils,complex) (3D only!!)
      !                  ghostsize(3) (I)ghostsize
      !
      !  Input/output : 
      !                                
      !                                
      !
      !  Output       :   
      !                 info          (I) return status. =0 if no error.
      !
      !  Remarks      : 
      !                                                  
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fdsolver_init.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.13  2006/09/04 18:34:43  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.11  2006/04/10 08:54:29  pchatela
      !  Made xp a REAL, DIMENSION(:,:), POINTER to get rid of warnings
      !
      !  Revision 1.10  2006/04/07 17:41:04  hiebers
      !  Changed type of variable xp to POINTER
      !
      !  Revision 1.9  2005/11/29 10:56:58  hiebers
      !  Changed arguments for fftwRoutines from INTEGER to INTEGER,DIMENSION(1)
      !  for  NAG compiler
      !
      !  Revision 1.8  2005/06/04 00:40:23  michaebe
      !  Cosmetics of cosmetics?
      !
      !  Revision 1.7  2005/06/04 00:37:49  michaebe
      !  cosmetics
      !
      !  Revision 1.6  2005/06/04 00:36:12  michaebe
      !  __ppm_module_fdsolver_init.f, line 491.25: 1516-023 (E) Subscript is
      !  out of bounds.
      !  This compiler warning eventually kinda got on my nerves so I acted.
      !
      !  Revision 1.5  2005/02/19 07:32:31  ivos
      !  Resolved CVS conflicts.
      !
      !  Revision 1.4  2005/02/18 08:01:55  hiebers
      !  minor changes in error messages
      !
      !  Revision 1.3  2005/02/16 22:22:59  ivos
      !  Bugfix: replaced non-existing ppm_module_data_fdsolver with
      !  ppm_module_data_fieldsolver.
      !
      !  Revision 1.2  2005/02/16 12:41:13  hiebers
      !  exchange FFTW_ESTIMATE by FFTW_MEASURE
      !
      !  Revision 1.1  2005/02/16 10:22:34  hiebers
      !  initial implementation
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fdsolver_init_2d_vec_s(DATA_fv,lda_fv,mesh_id_data, &
     &   field_topoid,topo_ids, field_topoid, mesh_ids,ghostsize, info)
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_mktopo
      USE ppm_module_typedef
      USE ppm_module_data_fieldsolver
      IMPLICIT NONE
      INTEGER, PARAMETER :: MK = ppm_kind_single
      !-------------------------------------------------------------------------
      !  FFTW include
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      ! POINTER to data
      REAL(MK), DIMENSION(:,:,:,:),      POINTER   :: data_fv
      INTEGER,                           INTENT(IN):: lda_fv
      ! mesh ID of the data
      INTEGER                    , INTENT(IN)      :: mesh_id_data
      ! topo ID of the field
      INTEGER                    , INTENT(IN   )   :: field_topoid
      ! topo / mesh ID for the mapping
      INTEGER, DIMENSION(2)      , INTENT(IN   )   :: topo_ids
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: mesh_ids
      INTEGER, DIMENSION(2)      , INTENT(IN   )   :: ghostsize
      INTEGER                    , INTENT(  OUT)   :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(ppm_kind_double)                   :: t0
      ! counters
      INTEGER                                 :: k, i, j
      CHARACTER(LEN=ppm_char)                 :: mesg
      REAL(MK), DIMENSION(:,:),     POINTER   :: data_real
      COMPLEX(MK), DIMENSION(:,:),  POINTER   :: data_comp, data_compl
      INTEGER, DIMENSION(2)                   :: lda
      ! size of the data_in 
      INTEGER,DIMENSION(1)                    :: MB_in
      INTEGER                                 :: Nx_in, Ny_in, Nz_in
      ! size of the data_out 
      INTEGER                                 :: Nx_out, Ny_out, Nz_out
      ! variables
      REAL(MK), DIMENSION(:,:),POINTER        :: xp
      INTEGER                                 :: Npart
      INTEGER                                 :: decomp, assign
      REAL(MK), DIMENSION(3  )                :: min_phys, max_phys
      REAL(MK), DIMENSION(3  )                :: length
      REAL(MK), DIMENSION(3  )                :: length_phys
      INTEGER , DIMENSION(6  )                :: bcdef 
      INTEGER                                 :: nsubs,topo_id, mesh_id
      INTEGER                                 :: yhmax, zhmax
      INTEGER                                 :: mesh_id_xpen, mesh_id_ypen
      INTEGER                                 :: mesh_id_xpen_complex
      INTEGER                                 :: mesh_id_zpen
      INTEGER                                 :: mesh_id_slab
      INTEGER                                 :: mesh_id_slab_complex
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      INTEGER , DIMENSION(:,:), POINTER       :: istart, istart_xpen_complex  
      INTEGER , DIMENSION(:,:), POINTER       :: istart_ypen, istart_trans
      INTEGER , DIMENSION(:,:), POINTER       :: istart_zpen
      INTEGER , DIMENSION(:,:), POINTER       :: ndata, ndata_xpen_complex
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_ypen, ndata_trans
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_zpen, ndata_slab
      INTEGER                                 :: dim, n,idom
      INTEGER                                 :: iopt
      INTEGER                                 :: topo_id_xpen, topo_id_ypen
      INTEGER                                 :: topo_id_zpen
      INTEGER                                 :: topo_id_slab
      INTEGER, DIMENSION(3)                   :: Nm, Nm_com, Nm_poisson
      INTEGER, DIMENSION(2)                   :: Nm_slab
      LOGICAL                                 :: Its_xpencil_topo
      LOGICAL                                 :: Its_xyslab_topo
      TYPE(ppm_t_topo)        , POINTER       :: f_topo
      TYPE(ppm_t_equi_mesh)   , POINTER       :: f_mesh
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_fdsolver_init',t0,info)
      f_topo => ppm_topo(field_topoid)%t
      f_mesh => f_topo%mesh(mesh_id_data)
      
      !-------------------------------------------------------------------------
      !  Error if FFTW library or NEC is not available
      !-------------------------------------------------------------------------
      info = ppm_error_error
      CALL ppm_error(ppm_err_nofftw,'ppm_fdsolver_init',  &
     &    'PPM was compiled without fftw support',279,info)
      CALL ppm_error(ppm_err_noMathKeisan,'ppm_fdsolver_init',  &
     &    'PPM was compiled without MATHKEISAN support',283,info)
      GOTO 9999   
     IF (info .NE. 0) THEN
          WRITE(mesg,'(A)') 'could not deallocate memory'
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_fdsolver_init',mesg,974,&
     &                                                                 info)
          GOTO 9999
      ENDIF
 9999 CONTINUE
      CALL substop('ppm_fdsolver_init',t0,info)
      RETURN

      END SUBROUTINE ppm_fdsolver_init_2d_vec_s 

      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_fdsolver_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : 
      !                initializes fieldsolver by creating all FFT-plans apriori
      !
      !  Input        :  
      !                  DATA_fv(:,:,:,:) (F) field data
      !                  lda_fv      (I) size of leading dimension in vector
      !                                  case
      !                  mesh_id_data(I) mesh ID of the current data field mesh
      !                  topo_ids(2) (I) topology IDs on which the FFTs are 
      !                                  performed
      !                                  topo_ids(1) initial topology(xpencils)
      !                                  topo_ids(2) second  topology(ypencils)
      !                                  topo_ids(3) third   topology(zpencils)
      !                                         (3D only!!)
      !                                 
      !                  mesh_ids(3) (I) mesh IDs where the FFTs are performed
      !                                  mesh_ids(1) first mesh 
      !                                         (xpencils,real
      !                                  mesh_ids(2) second mesh 
      !                                         (xpencils,complex)
      !                                  mesh_ids(3) third mesh 
      !                                         (ypencils,complex)
      !                                  mesh_ids(4) forth mesh 
      !                                         (zpencils,complex) (3D only!!)
      !                  ghostsize(3) (I)ghostsize
      !
      !  Input/output : 
      !                                
      !                                
      !
      !  Output       :   
      !                 info          (I) return status. =0 if no error.
      !
      !  Remarks      : 
      !                                                  
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fdsolver_init.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.13  2006/09/04 18:34:43  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.11  2006/04/10 08:54:29  pchatela
      !  Made xp a REAL, DIMENSION(:,:), POINTER to get rid of warnings
      !
      !  Revision 1.10  2006/04/07 17:41:04  hiebers
      !  Changed type of variable xp to POINTER
      !
      !  Revision 1.9  2005/11/29 10:56:58  hiebers
      !  Changed arguments for fftwRoutines from INTEGER to INTEGER,DIMENSION(1)
      !  for  NAG compiler
      !
      !  Revision 1.8  2005/06/04 00:40:23  michaebe
      !  Cosmetics of cosmetics?
      !
      !  Revision 1.7  2005/06/04 00:37:49  michaebe
      !  cosmetics
      !
      !  Revision 1.6  2005/06/04 00:36:12  michaebe
      !  __ppm_module_fdsolver_init.f, line 491.25: 1516-023 (E) Subscript is
      !  out of bounds.
      !  This compiler warning eventually kinda got on my nerves so I acted.
      !
      !  Revision 1.5  2005/02/19 07:32:31  ivos
      !  Resolved CVS conflicts.
      !
      !  Revision 1.4  2005/02/18 08:01:55  hiebers
      !  minor changes in error messages
      !
      !  Revision 1.3  2005/02/16 22:22:59  ivos
      !  Bugfix: replaced non-existing ppm_module_data_fdsolver with
      !  ppm_module_data_fieldsolver.
      !
      !  Revision 1.2  2005/02/16 12:41:13  hiebers
      !  exchange FFTW_ESTIMATE by FFTW_MEASURE
      !
      !  Revision 1.1  2005/02/16 10:22:34  hiebers
      !  initial implementation
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fdsolver_init_2d_vec_d(DATA_fv,lda_fv,mesh_id_data, &
     &   field_topoid,topo_ids, field_topoid, mesh_ids,ghostsize, info)
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_mktopo
      USE ppm_module_typedef
      USE ppm_module_data_fieldsolver
      IMPLICIT NONE
      INTEGER, PARAMETER :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  FFTW include
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      ! POINTER to data
      REAL(MK), DIMENSION(:,:,:,:),      POINTER   :: data_fv
      INTEGER,                           INTENT(IN):: lda_fv
      ! mesh ID of the data
      INTEGER                    , INTENT(IN)      :: mesh_id_data
      ! topo ID of the field
      INTEGER                    , INTENT(IN   )   :: field_topoid
      ! topo / mesh ID for the mapping
      INTEGER, DIMENSION(2)      , INTENT(IN   )   :: topo_ids
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: mesh_ids
      INTEGER, DIMENSION(2)      , INTENT(IN   )   :: ghostsize
      INTEGER                    , INTENT(  OUT)   :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(ppm_kind_double)                   :: t0
      ! counters
      INTEGER                                 :: k, i, j
      CHARACTER(LEN=ppm_char)                 :: mesg
      REAL(MK), DIMENSION(:,:),     POINTER   :: data_real
      COMPLEX(MK), DIMENSION(:,:),  POINTER   :: data_comp, data_compl
      INTEGER, DIMENSION(2)                   :: lda
      ! size of the data_in 
      INTEGER,DIMENSION(1)                    :: MB_in
      INTEGER                                 :: Nx_in, Ny_in, Nz_in
      ! size of the data_out 
      INTEGER                                 :: Nx_out, Ny_out, Nz_out
      ! variables
      REAL(MK), DIMENSION(:,:),POINTER        :: xp
      INTEGER                                 :: Npart
      INTEGER                                 :: decomp, assign
      REAL(MK), DIMENSION(3  )                :: min_phys, max_phys
      REAL(MK), DIMENSION(3  )                :: length
      REAL(MK), DIMENSION(3  )                :: length_phys
      INTEGER , DIMENSION(6  )                :: bcdef 
      INTEGER                                 :: nsubs,topo_id, mesh_id
      INTEGER                                 :: yhmax, zhmax
      INTEGER                                 :: mesh_id_xpen, mesh_id_ypen
      INTEGER                                 :: mesh_id_xpen_complex
      INTEGER                                 :: mesh_id_zpen
      INTEGER                                 :: mesh_id_slab
      INTEGER                                 :: mesh_id_slab_complex
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      INTEGER , DIMENSION(:,:), POINTER       :: istart, istart_xpen_complex  
      INTEGER , DIMENSION(:,:), POINTER       :: istart_ypen, istart_trans
      INTEGER , DIMENSION(:,:), POINTER       :: istart_zpen
      INTEGER , DIMENSION(:,:), POINTER       :: ndata, ndata_xpen_complex
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_ypen, ndata_trans
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_zpen, ndata_slab
      INTEGER                                 :: dim, n,idom
      INTEGER                                 :: iopt
      INTEGER                                 :: topo_id_xpen, topo_id_ypen
      INTEGER                                 :: topo_id_zpen
      INTEGER                                 :: topo_id_slab
      INTEGER, DIMENSION(3)                   :: Nm, Nm_com, Nm_poisson
      INTEGER, DIMENSION(2)                   :: Nm_slab
      LOGICAL                                 :: Its_xpencil_topo
      LOGICAL                                 :: Its_xyslab_topo
      TYPE(ppm_t_topo)        , POINTER       :: f_topo
      TYPE(ppm_t_equi_mesh)   , POINTER       :: f_mesh
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_fdsolver_init',t0,info)
      f_topo => ppm_topo(field_topoid)%t
      f_mesh => f_topo%mesh(mesh_id_data)
      
      !-------------------------------------------------------------------------
      !  Error if FFTW library or NEC is not available
      !-------------------------------------------------------------------------
      info = ppm_error_error
      CALL ppm_error(ppm_err_nofftw,'ppm_fdsolver_init',  &
     &    'PPM was compiled without fftw support',279,info)
      CALL ppm_error(ppm_err_noMathKeisan,'ppm_fdsolver_init',  &
     &    'PPM was compiled without MATHKEISAN support',283,info)
      GOTO 9999   
     IF (info .NE. 0) THEN
          WRITE(mesg,'(A)') 'could not deallocate memory'
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_fdsolver_init',mesg,974,&
     &                                                                 info)
          GOTO 9999
      ENDIF
 9999 CONTINUE
      CALL substop('ppm_fdsolver_init',t0,info)
      RETURN

      END SUBROUTINE ppm_fdsolver_init_2d_vec_d




      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_fdsolver_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : 
      !                initializes fieldsolver by creating all FFT-plans apriori
      !
      !  Input        :  
      !                  DATA_fv(:,:,:,:) (F) field data
      !                  lda_fv      (I) size of leading dimension in vector
      !                                  case
      !                  mesh_id_data(I) mesh ID of the current data field mesh
      !                  topo_ids(2) (I) topology IDs on which the FFTs are 
      !                                  performed
      !                                  topo_ids(1) initial topology(xpencils)
      !                                  topo_ids(2) second  topology(ypencils)
      !                                  topo_ids(3) third   topology(zpencils)
      !                                         (3D only!!)
      !                                 
      !                  mesh_ids(3) (I) mesh IDs where the FFTs are performed
      !                                  mesh_ids(1) first mesh 
      !                                         (xpencils,real
      !                                  mesh_ids(2) second mesh 
      !                                         (xpencils,complex)
      !                                  mesh_ids(3) third mesh 
      !                                         (ypencils,complex)
      !                                  mesh_ids(4) forth mesh 
      !                                         (zpencils,complex) (3D only!!)
      !                  ghostsize(3) (I)ghostsize
      !
      !  Input/output : 
      !                                
      !                                
      !
      !  Output       :   
      !                 info          (I) return status. =0 if no error.
      !
      !  Remarks      : 
      !                                                  
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fdsolver_init.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.13  2006/09/04 18:34:43  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.11  2006/04/10 08:54:29  pchatela
      !  Made xp a REAL, DIMENSION(:,:), POINTER to get rid of warnings
      !
      !  Revision 1.10  2006/04/07 17:41:04  hiebers
      !  Changed type of variable xp to POINTER
      !
      !  Revision 1.9  2005/11/29 10:56:58  hiebers
      !  Changed arguments for fftwRoutines from INTEGER to INTEGER,DIMENSION(1)
      !  for  NAG compiler
      !
      !  Revision 1.8  2005/06/04 00:40:23  michaebe
      !  Cosmetics of cosmetics?
      !
      !  Revision 1.7  2005/06/04 00:37:49  michaebe
      !  cosmetics
      !
      !  Revision 1.6  2005/06/04 00:36:12  michaebe
      !  __ppm_module_fdsolver_init.f, line 491.25: 1516-023 (E) Subscript is
      !  out of bounds.
      !  This compiler warning eventually kinda got on my nerves so I acted.
      !
      !  Revision 1.5  2005/02/19 07:32:31  ivos
      !  Resolved CVS conflicts.
      !
      !  Revision 1.4  2005/02/18 08:01:55  hiebers
      !  minor changes in error messages
      !
      !  Revision 1.3  2005/02/16 22:22:59  ivos
      !  Bugfix: replaced non-existing ppm_module_data_fdsolver with
      !  ppm_module_data_fieldsolver.
      !
      !  Revision 1.2  2005/02/16 12:41:13  hiebers
      !  exchange FFTW_ESTIMATE by FFTW_MEASURE
      !
      !  Revision 1.1  2005/02/16 10:22:34  hiebers
      !  initial implementation
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fdsolver_init_3d_sca_s(DATA_fv,mesh_id_data,topo_ids, &
     &  field_topoid,mesh_ids,ghostsize, info)
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_mktopo
      USE ppm_module_typedef
      USE ppm_module_data_fieldsolver
      IMPLICIT NONE
      INTEGER, PARAMETER :: MK = ppm_kind_single
      !-------------------------------------------------------------------------
      !  FFTW include
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      ! POINTER to data
      REAL(MK), DIMENSION(:,:,:,:),      POINTER   :: data_fv
      ! mesh ID of the data
      INTEGER                    , INTENT(IN)      :: mesh_id_data
      ! topo ID of the field
      INTEGER                    , INTENT(IN   )   :: field_topoid
      ! topo / mesh ID for the mapping
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: topo_ids
      INTEGER, DIMENSION(4)      , INTENT(IN   )   :: mesh_ids
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: ghostsize
      INTEGER                    , INTENT(  OUT)   :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(ppm_kind_double)                   :: t0
      ! counters
      INTEGER                                 :: k, i, j
      CHARACTER(LEN=ppm_char)                 :: mesg
      REAL(MK), DIMENSION(:,:,:),   POINTER   :: data_real
      COMPLEX(MK), DIMENSION(:,:,:),POINTER   :: data_comp, data_compl
      INTEGER, DIMENSION(3)                   :: lda
      ! size of the data_in 
      INTEGER,DIMENSION(1)                    :: MB_in
      INTEGER                                 :: Nx_in, Ny_in, Nz_in
      ! size of the data_out 
      INTEGER                                 :: Nx_out, Ny_out, Nz_out
      ! variables
      REAL(MK), DIMENSION(:,:),POINTER        :: xp
      INTEGER                                 :: Npart
      INTEGER                                 :: decomp, assign
      REAL(MK), DIMENSION(3  )                :: min_phys, max_phys
      REAL(MK), DIMENSION(3  )                :: length
      REAL(MK), DIMENSION(3  )                :: length_phys
      INTEGER , DIMENSION(6  )                :: bcdef 
      INTEGER                                 :: nsubs,topo_id, mesh_id
      INTEGER                                 :: yhmax, zhmax
      INTEGER                                 :: mesh_id_xpen, mesh_id_ypen
      INTEGER                                 :: mesh_id_xpen_complex
      INTEGER                                 :: mesh_id_zpen
      INTEGER                                 :: mesh_id_slab
      INTEGER                                 :: mesh_id_slab_complex
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      INTEGER , DIMENSION(:,:), POINTER       :: istart, istart_xpen_complex  
      INTEGER , DIMENSION(:,:), POINTER       :: istart_ypen, istart_trans
      INTEGER , DIMENSION(:,:), POINTER       :: istart_zpen
      INTEGER , DIMENSION(:,:), POINTER       :: ndata, ndata_xpen_complex
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_ypen, ndata_trans
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_zpen, ndata_slab
      INTEGER                                 :: dim, n,idom
      INTEGER                                 :: iopt
      INTEGER                                 :: topo_id_xpen, topo_id_ypen
      INTEGER                                 :: topo_id_zpen
      INTEGER                                 :: topo_id_slab
      INTEGER, DIMENSION(3)                   :: Nm, Nm_com, Nm_poisson
      INTEGER, DIMENSION(2)                   :: Nm_slab
      LOGICAL                                 :: Its_xpencil_topo
      LOGICAL                                 :: Its_xyslab_topo
      TYPE(ppm_t_topo)        , POINTER       :: f_topo
      TYPE(ppm_t_equi_mesh)   , POINTER       :: f_mesh
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_fdsolver_init',t0,info)
      f_topo => ppm_topo(field_topoid)%t
      f_mesh => f_topo%mesh(mesh_id_data)
      
      !-------------------------------------------------------------------------
      !  Error if FFTW library or NEC is not available
      !-------------------------------------------------------------------------
      info = ppm_error_error
      CALL ppm_error(ppm_err_nofftw,'ppm_fdsolver_init',  &
     &    'PPM was compiled without fftw support',279,info)
      CALL ppm_error(ppm_err_noMathKeisan,'ppm_fdsolver_init',  &
     &    'PPM was compiled without MATHKEISAN support',283,info)
      GOTO 9999   
     IF (info .NE. 0) THEN
          WRITE(mesg,'(A)') 'could not deallocate memory'
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_fdsolver_init',mesg,974,&
     &                                                                 info)
          GOTO 9999
      ENDIF
 9999 CONTINUE
      CALL substop('ppm_fdsolver_init',t0,info)
      RETURN

      END SUBROUTINE ppm_fdsolver_init_3d_sca_s 

      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_fdsolver_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : 
      !                initializes fieldsolver by creating all FFT-plans apriori
      !
      !  Input        :  
      !                  DATA_fv(:,:,:,:) (F) field data
      !                  lda_fv      (I) size of leading dimension in vector
      !                                  case
      !                  mesh_id_data(I) mesh ID of the current data field mesh
      !                  topo_ids(2) (I) topology IDs on which the FFTs are 
      !                                  performed
      !                                  topo_ids(1) initial topology(xpencils)
      !                                  topo_ids(2) second  topology(ypencils)
      !                                  topo_ids(3) third   topology(zpencils)
      !                                         (3D only!!)
      !                                 
      !                  mesh_ids(3) (I) mesh IDs where the FFTs are performed
      !                                  mesh_ids(1) first mesh 
      !                                         (xpencils,real
      !                                  mesh_ids(2) second mesh 
      !                                         (xpencils,complex)
      !                                  mesh_ids(3) third mesh 
      !                                         (ypencils,complex)
      !                                  mesh_ids(4) forth mesh 
      !                                         (zpencils,complex) (3D only!!)
      !                  ghostsize(3) (I)ghostsize
      !
      !  Input/output : 
      !                                
      !                                
      !
      !  Output       :   
      !                 info          (I) return status. =0 if no error.
      !
      !  Remarks      : 
      !                                                  
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fdsolver_init.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.13  2006/09/04 18:34:43  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.11  2006/04/10 08:54:29  pchatela
      !  Made xp a REAL, DIMENSION(:,:), POINTER to get rid of warnings
      !
      !  Revision 1.10  2006/04/07 17:41:04  hiebers
      !  Changed type of variable xp to POINTER
      !
      !  Revision 1.9  2005/11/29 10:56:58  hiebers
      !  Changed arguments for fftwRoutines from INTEGER to INTEGER,DIMENSION(1)
      !  for  NAG compiler
      !
      !  Revision 1.8  2005/06/04 00:40:23  michaebe
      !  Cosmetics of cosmetics?
      !
      !  Revision 1.7  2005/06/04 00:37:49  michaebe
      !  cosmetics
      !
      !  Revision 1.6  2005/06/04 00:36:12  michaebe
      !  __ppm_module_fdsolver_init.f, line 491.25: 1516-023 (E) Subscript is
      !  out of bounds.
      !  This compiler warning eventually kinda got on my nerves so I acted.
      !
      !  Revision 1.5  2005/02/19 07:32:31  ivos
      !  Resolved CVS conflicts.
      !
      !  Revision 1.4  2005/02/18 08:01:55  hiebers
      !  minor changes in error messages
      !
      !  Revision 1.3  2005/02/16 22:22:59  ivos
      !  Bugfix: replaced non-existing ppm_module_data_fdsolver with
      !  ppm_module_data_fieldsolver.
      !
      !  Revision 1.2  2005/02/16 12:41:13  hiebers
      !  exchange FFTW_ESTIMATE by FFTW_MEASURE
      !
      !  Revision 1.1  2005/02/16 10:22:34  hiebers
      !  initial implementation
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fdsolver_init_3d_sca_d(DATA_fv,mesh_id_data,topo_ids, &
     &  field_topoid,mesh_ids,ghostsize, info)
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_mktopo
      USE ppm_module_typedef
      USE ppm_module_data_fieldsolver
      IMPLICIT NONE
      INTEGER, PARAMETER :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  FFTW include
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      ! POINTER to data
      REAL(MK), DIMENSION(:,:,:,:),      POINTER   :: data_fv
      ! mesh ID of the data
      INTEGER                    , INTENT(IN)      :: mesh_id_data
      ! topo ID of the field
      INTEGER                    , INTENT(IN   )   :: field_topoid
      ! topo / mesh ID for the mapping
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: topo_ids
      INTEGER, DIMENSION(4)      , INTENT(IN   )   :: mesh_ids
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: ghostsize
      INTEGER                    , INTENT(  OUT)   :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(ppm_kind_double)                   :: t0
      ! counters
      INTEGER                                 :: k, i, j
      CHARACTER(LEN=ppm_char)                 :: mesg
      REAL(MK), DIMENSION(:,:,:),   POINTER   :: data_real
      COMPLEX(MK), DIMENSION(:,:,:),POINTER   :: data_comp, data_compl
      INTEGER, DIMENSION(3)                   :: lda
      ! size of the data_in 
      INTEGER,DIMENSION(1)                    :: MB_in
      INTEGER                                 :: Nx_in, Ny_in, Nz_in
      ! size of the data_out 
      INTEGER                                 :: Nx_out, Ny_out, Nz_out
      ! variables
      REAL(MK), DIMENSION(:,:),POINTER        :: xp
      INTEGER                                 :: Npart
      INTEGER                                 :: decomp, assign
      REAL(MK), DIMENSION(3  )                :: min_phys, max_phys
      REAL(MK), DIMENSION(3  )                :: length
      REAL(MK), DIMENSION(3  )                :: length_phys
      INTEGER , DIMENSION(6  )                :: bcdef 
      INTEGER                                 :: nsubs,topo_id, mesh_id
      INTEGER                                 :: yhmax, zhmax
      INTEGER                                 :: mesh_id_xpen, mesh_id_ypen
      INTEGER                                 :: mesh_id_xpen_complex
      INTEGER                                 :: mesh_id_zpen
      INTEGER                                 :: mesh_id_slab
      INTEGER                                 :: mesh_id_slab_complex
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      INTEGER , DIMENSION(:,:), POINTER       :: istart, istart_xpen_complex  
      INTEGER , DIMENSION(:,:), POINTER       :: istart_ypen, istart_trans
      INTEGER , DIMENSION(:,:), POINTER       :: istart_zpen
      INTEGER , DIMENSION(:,:), POINTER       :: ndata, ndata_xpen_complex
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_ypen, ndata_trans
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_zpen, ndata_slab
      INTEGER                                 :: dim, n,idom
      INTEGER                                 :: iopt
      INTEGER                                 :: topo_id_xpen, topo_id_ypen
      INTEGER                                 :: topo_id_zpen
      INTEGER                                 :: topo_id_slab
      INTEGER, DIMENSION(3)                   :: Nm, Nm_com, Nm_poisson
      INTEGER, DIMENSION(2)                   :: Nm_slab
      LOGICAL                                 :: Its_xpencil_topo
      LOGICAL                                 :: Its_xyslab_topo
      TYPE(ppm_t_topo)        , POINTER       :: f_topo
      TYPE(ppm_t_equi_mesh)   , POINTER       :: f_mesh
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_fdsolver_init',t0,info)
      f_topo => ppm_topo(field_topoid)%t
      f_mesh => f_topo%mesh(mesh_id_data)
      
      !-------------------------------------------------------------------------
      !  Error if FFTW library or NEC is not available
      !-------------------------------------------------------------------------
      info = ppm_error_error
      CALL ppm_error(ppm_err_nofftw,'ppm_fdsolver_init',  &
     &    'PPM was compiled without fftw support',279,info)
      CALL ppm_error(ppm_err_noMathKeisan,'ppm_fdsolver_init',  &
     &    'PPM was compiled without MATHKEISAN support',283,info)
      GOTO 9999   
     IF (info .NE. 0) THEN
          WRITE(mesg,'(A)') 'could not deallocate memory'
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_fdsolver_init',mesg,974,&
     &                                                                 info)
          GOTO 9999
      ENDIF
 9999 CONTINUE
      CALL substop('ppm_fdsolver_init',t0,info)
      RETURN

      END SUBROUTINE ppm_fdsolver_init_3d_sca_d 


      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_fdsolver_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : 
      !                initializes fieldsolver by creating all FFT-plans apriori
      !
      !  Input        :  
      !                  DATA_fv(:,:,:,:) (F) field data
      !                  lda_fv      (I) size of leading dimension in vector
      !                                  case
      !                  mesh_id_data(I) mesh ID of the current data field mesh
      !                  topo_ids(2) (I) topology IDs on which the FFTs are 
      !                                  performed
      !                                  topo_ids(1) initial topology(xpencils)
      !                                  topo_ids(2) second  topology(ypencils)
      !                                  topo_ids(3) third   topology(zpencils)
      !                                         (3D only!!)
      !                                 
      !                  mesh_ids(3) (I) mesh IDs where the FFTs are performed
      !                                  mesh_ids(1) first mesh 
      !                                         (xpencils,real
      !                                  mesh_ids(2) second mesh 
      !                                         (xpencils,complex)
      !                                  mesh_ids(3) third mesh 
      !                                         (ypencils,complex)
      !                                  mesh_ids(4) forth mesh 
      !                                         (zpencils,complex) (3D only!!)
      !                  ghostsize(3) (I)ghostsize
      !
      !  Input/output : 
      !                                
      !                                
      !
      !  Output       :   
      !                 info          (I) return status. =0 if no error.
      !
      !  Remarks      : 
      !                                                  
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fdsolver_init.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.13  2006/09/04 18:34:43  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.11  2006/04/10 08:54:29  pchatela
      !  Made xp a REAL, DIMENSION(:,:), POINTER to get rid of warnings
      !
      !  Revision 1.10  2006/04/07 17:41:04  hiebers
      !  Changed type of variable xp to POINTER
      !
      !  Revision 1.9  2005/11/29 10:56:58  hiebers
      !  Changed arguments for fftwRoutines from INTEGER to INTEGER,DIMENSION(1)
      !  for  NAG compiler
      !
      !  Revision 1.8  2005/06/04 00:40:23  michaebe
      !  Cosmetics of cosmetics?
      !
      !  Revision 1.7  2005/06/04 00:37:49  michaebe
      !  cosmetics
      !
      !  Revision 1.6  2005/06/04 00:36:12  michaebe
      !  __ppm_module_fdsolver_init.f, line 491.25: 1516-023 (E) Subscript is
      !  out of bounds.
      !  This compiler warning eventually kinda got on my nerves so I acted.
      !
      !  Revision 1.5  2005/02/19 07:32:31  ivos
      !  Resolved CVS conflicts.
      !
      !  Revision 1.4  2005/02/18 08:01:55  hiebers
      !  minor changes in error messages
      !
      !  Revision 1.3  2005/02/16 22:22:59  ivos
      !  Bugfix: replaced non-existing ppm_module_data_fdsolver with
      !  ppm_module_data_fieldsolver.
      !
      !  Revision 1.2  2005/02/16 12:41:13  hiebers
      !  exchange FFTW_ESTIMATE by FFTW_MEASURE
      !
      !  Revision 1.1  2005/02/16 10:22:34  hiebers
      !  initial implementation
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fdsolver_init_3d_vec_s(DATA_fv,lda_fv,mesh_id_data, &
     &   field_topoid,topo_ids, field_topoid, mesh_ids,ghostsize, info)
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_mktopo
      USE ppm_module_typedef
      USE ppm_module_data_fieldsolver
      IMPLICIT NONE
      INTEGER, PARAMETER :: MK = ppm_kind_single
      !-------------------------------------------------------------------------
      !  FFTW include
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      ! POINTER to data
      REAL(MK), DIMENSION(:,:,:,:,:),    POINTER   :: data_fv
      INTEGER,                           INTENT(IN):: lda_fv
      ! mesh ID of the data
      INTEGER                    , INTENT(IN)      :: mesh_id_data
      ! topo ID of the field
      INTEGER                    , INTENT(IN   )   :: field_topoid
      ! topo / mesh ID for the mapping
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: topo_ids
      INTEGER, DIMENSION(4)      , INTENT(IN   )   :: mesh_ids
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: ghostsize
      INTEGER                    , INTENT(  OUT)   :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(ppm_kind_double)                   :: t0
      ! counters
      INTEGER                                 :: k, i, j
      CHARACTER(LEN=ppm_char)                 :: mesg
      REAL(MK), DIMENSION(:,:,:),   POINTER   :: data_real
      COMPLEX(MK), DIMENSION(:,:,:),POINTER   :: data_comp, data_compl
      INTEGER, DIMENSION(3)                   :: lda
      ! size of the data_in 
      INTEGER,DIMENSION(1)                    :: MB_in
      INTEGER                                 :: Nx_in, Ny_in, Nz_in
      ! size of the data_out 
      INTEGER                                 :: Nx_out, Ny_out, Nz_out
      ! variables
      REAL(MK), DIMENSION(:,:),POINTER        :: xp
      INTEGER                                 :: Npart
      INTEGER                                 :: decomp, assign
      REAL(MK), DIMENSION(3  )                :: min_phys, max_phys
      REAL(MK), DIMENSION(3  )                :: length
      REAL(MK), DIMENSION(3  )                :: length_phys
      INTEGER , DIMENSION(6  )                :: bcdef 
      INTEGER                                 :: nsubs,topo_id, mesh_id
      INTEGER                                 :: yhmax, zhmax
      INTEGER                                 :: mesh_id_xpen, mesh_id_ypen
      INTEGER                                 :: mesh_id_xpen_complex
      INTEGER                                 :: mesh_id_zpen
      INTEGER                                 :: mesh_id_slab
      INTEGER                                 :: mesh_id_slab_complex
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      INTEGER , DIMENSION(:,:), POINTER       :: istart, istart_xpen_complex  
      INTEGER , DIMENSION(:,:), POINTER       :: istart_ypen, istart_trans
      INTEGER , DIMENSION(:,:), POINTER       :: istart_zpen
      INTEGER , DIMENSION(:,:), POINTER       :: ndata, ndata_xpen_complex
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_ypen, ndata_trans
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_zpen, ndata_slab
      INTEGER                                 :: dim, n,idom
      INTEGER                                 :: iopt
      INTEGER                                 :: topo_id_xpen, topo_id_ypen
      INTEGER                                 :: topo_id_zpen
      INTEGER                                 :: topo_id_slab
      INTEGER, DIMENSION(3)                   :: Nm, Nm_com, Nm_poisson
      INTEGER, DIMENSION(2)                   :: Nm_slab
      LOGICAL                                 :: Its_xpencil_topo
      LOGICAL                                 :: Its_xyslab_topo
      TYPE(ppm_t_topo)        , POINTER       :: f_topo
      TYPE(ppm_t_equi_mesh)   , POINTER       :: f_mesh
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_fdsolver_init',t0,info)
      f_topo => ppm_topo(field_topoid)%t
      f_mesh => f_topo%mesh(mesh_id_data)
      
      !-------------------------------------------------------------------------
      !  Error if FFTW library or NEC is not available
      !-------------------------------------------------------------------------
      info = ppm_error_error
      CALL ppm_error(ppm_err_nofftw,'ppm_fdsolver_init',  &
     &    'PPM was compiled without fftw support',279,info)
      CALL ppm_error(ppm_err_noMathKeisan,'ppm_fdsolver_init',  &
     &    'PPM was compiled without MATHKEISAN support',283,info)
      GOTO 9999   
     IF (info .NE. 0) THEN
          WRITE(mesg,'(A)') 'could not deallocate memory'
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_fdsolver_init',mesg,974,&
     &                                                                 info)
          GOTO 9999
      ENDIF
 9999 CONTINUE
      CALL substop('ppm_fdsolver_init',t0,info)
      RETURN

      END SUBROUTINE ppm_fdsolver_init_3d_vec_s

      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_fdsolver_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : 
      !                initializes fieldsolver by creating all FFT-plans apriori
      !
      !  Input        :  
      !                  DATA_fv(:,:,:,:) (F) field data
      !                  lda_fv      (I) size of leading dimension in vector
      !                                  case
      !                  mesh_id_data(I) mesh ID of the current data field mesh
      !                  topo_ids(2) (I) topology IDs on which the FFTs are 
      !                                  performed
      !                                  topo_ids(1) initial topology(xpencils)
      !                                  topo_ids(2) second  topology(ypencils)
      !                                  topo_ids(3) third   topology(zpencils)
      !                                         (3D only!!)
      !                                 
      !                  mesh_ids(3) (I) mesh IDs where the FFTs are performed
      !                                  mesh_ids(1) first mesh 
      !                                         (xpencils,real
      !                                  mesh_ids(2) second mesh 
      !                                         (xpencils,complex)
      !                                  mesh_ids(3) third mesh 
      !                                         (ypencils,complex)
      !                                  mesh_ids(4) forth mesh 
      !                                         (zpencils,complex) (3D only!!)
      !                  ghostsize(3) (I)ghostsize
      !
      !  Input/output : 
      !                                
      !                                
      !
      !  Output       :   
      !                 info          (I) return status. =0 if no error.
      !
      !  Remarks      : 
      !                                                  
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fdsolver_init.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.13  2006/09/04 18:34:43  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.11  2006/04/10 08:54:29  pchatela
      !  Made xp a REAL, DIMENSION(:,:), POINTER to get rid of warnings
      !
      !  Revision 1.10  2006/04/07 17:41:04  hiebers
      !  Changed type of variable xp to POINTER
      !
      !  Revision 1.9  2005/11/29 10:56:58  hiebers
      !  Changed arguments for fftwRoutines from INTEGER to INTEGER,DIMENSION(1)
      !  for  NAG compiler
      !
      !  Revision 1.8  2005/06/04 00:40:23  michaebe
      !  Cosmetics of cosmetics?
      !
      !  Revision 1.7  2005/06/04 00:37:49  michaebe
      !  cosmetics
      !
      !  Revision 1.6  2005/06/04 00:36:12  michaebe
      !  __ppm_module_fdsolver_init.f, line 491.25: 1516-023 (E) Subscript is
      !  out of bounds.
      !  This compiler warning eventually kinda got on my nerves so I acted.
      !
      !  Revision 1.5  2005/02/19 07:32:31  ivos
      !  Resolved CVS conflicts.
      !
      !  Revision 1.4  2005/02/18 08:01:55  hiebers
      !  minor changes in error messages
      !
      !  Revision 1.3  2005/02/16 22:22:59  ivos
      !  Bugfix: replaced non-existing ppm_module_data_fdsolver with
      !  ppm_module_data_fieldsolver.
      !
      !  Revision 1.2  2005/02/16 12:41:13  hiebers
      !  exchange FFTW_ESTIMATE by FFTW_MEASURE
      !
      !  Revision 1.1  2005/02/16 10:22:34  hiebers
      !  initial implementation
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fdsolver_init_3d_vec_d(DATA_fv,lda_fv,mesh_id_data, &
     &   field_topoid,topo_ids, field_topoid, mesh_ids,ghostsize, info)
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_mktopo
      USE ppm_module_typedef
      USE ppm_module_data_fieldsolver
      IMPLICIT NONE
      INTEGER, PARAMETER :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  FFTW include
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      ! POINTER to data
      REAL(MK), DIMENSION(:,:,:,:,:),    POINTER   :: data_fv
      INTEGER,                           INTENT(IN):: lda_fv
      ! mesh ID of the data
      INTEGER                    , INTENT(IN)      :: mesh_id_data
      ! topo ID of the field
      INTEGER                    , INTENT(IN   )   :: field_topoid
      ! topo / mesh ID for the mapping
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: topo_ids
      INTEGER, DIMENSION(4)      , INTENT(IN   )   :: mesh_ids
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: ghostsize
      INTEGER                    , INTENT(  OUT)   :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(ppm_kind_double)                   :: t0
      ! counters
      INTEGER                                 :: k, i, j
      CHARACTER(LEN=ppm_char)                 :: mesg
      REAL(MK), DIMENSION(:,:,:),   POINTER   :: data_real
      COMPLEX(MK), DIMENSION(:,:,:),POINTER   :: data_comp, data_compl
      INTEGER, DIMENSION(3)                   :: lda
      ! size of the data_in 
      INTEGER,DIMENSION(1)                    :: MB_in
      INTEGER                                 :: Nx_in, Ny_in, Nz_in
      ! size of the data_out 
      INTEGER                                 :: Nx_out, Ny_out, Nz_out
      ! variables
      REAL(MK), DIMENSION(:,:),POINTER        :: xp
      INTEGER                                 :: Npart
      INTEGER                                 :: decomp, assign
      REAL(MK), DIMENSION(3  )                :: min_phys, max_phys
      REAL(MK), DIMENSION(3  )                :: length
      REAL(MK), DIMENSION(3  )                :: length_phys
      INTEGER , DIMENSION(6  )                :: bcdef 
      INTEGER                                 :: nsubs,topo_id, mesh_id
      INTEGER                                 :: yhmax, zhmax
      INTEGER                                 :: mesh_id_xpen, mesh_id_ypen
      INTEGER                                 :: mesh_id_xpen_complex
      INTEGER                                 :: mesh_id_zpen
      INTEGER                                 :: mesh_id_slab
      INTEGER                                 :: mesh_id_slab_complex
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      INTEGER , DIMENSION(:,:), POINTER       :: istart, istart_xpen_complex  
      INTEGER , DIMENSION(:,:), POINTER       :: istart_ypen, istart_trans
      INTEGER , DIMENSION(:,:), POINTER       :: istart_zpen
      INTEGER , DIMENSION(:,:), POINTER       :: ndata, ndata_xpen_complex
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_ypen, ndata_trans
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_zpen, ndata_slab
      INTEGER                                 :: dim, n,idom
      INTEGER                                 :: iopt
      INTEGER                                 :: topo_id_xpen, topo_id_ypen
      INTEGER                                 :: topo_id_zpen
      INTEGER                                 :: topo_id_slab
      INTEGER, DIMENSION(3)                   :: Nm, Nm_com, Nm_poisson
      INTEGER, DIMENSION(2)                   :: Nm_slab
      LOGICAL                                 :: Its_xpencil_topo
      LOGICAL                                 :: Its_xyslab_topo
      TYPE(ppm_t_topo)        , POINTER       :: f_topo
      TYPE(ppm_t_equi_mesh)   , POINTER       :: f_mesh
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_fdsolver_init',t0,info)
      f_topo => ppm_topo(field_topoid)%t
      f_mesh => f_topo%mesh(mesh_id_data)
      
      !-------------------------------------------------------------------------
      !  Error if FFTW library or NEC is not available
      !-------------------------------------------------------------------------
      info = ppm_error_error
      CALL ppm_error(ppm_err_nofftw,'ppm_fdsolver_init',  &
     &    'PPM was compiled without fftw support',279,info)
      CALL ppm_error(ppm_err_noMathKeisan,'ppm_fdsolver_init',  &
     &    'PPM was compiled without MATHKEISAN support',283,info)
      GOTO 9999   
     IF (info .NE. 0) THEN
          WRITE(mesg,'(A)') 'could not deallocate memory'
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_fdsolver_init',mesg,974,&
     &                                                                 info)
          GOTO 9999
      ENDIF
 9999 CONTINUE
      CALL substop('ppm_fdsolver_init',t0,info)
      RETURN

      END SUBROUTINE ppm_fdsolver_init_3d_vec_d


      END MODULE ppm_module_fdsolver_init
