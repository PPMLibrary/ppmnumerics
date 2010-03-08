      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_gmm_init
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_gmm_init(meshid,Nest,prec,info)
      !!! This routine initializes the ppm_gmm module and allocates all data
      !!! structures.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_gmm
      USE ppm_module_data_mesh
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN   ) :: Nest
      !!! Estimated number of grid points adjacent to the zero level. This is
      !!! only used to estimate the size of the memory needed and will be grown
      !!! automatically if too small.
      INTEGER, INTENT(IN   ) :: prec
      !!! Precision of the field data to be reinitialized. One of:
      !!!
      !!! *ppm_kind_single
      !!! *ppm_kind_double
      INTEGER, INTENT(IN   ) :: meshid
      !!! Mesh ID (user numbering) for which a GMM should be initialized.
      INTEGER, INTENT(  OUT) :: info
      !!! Return status. 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)  :: ldu
      INTEGER                :: iopt,i,isub
      LOGICAL                :: lok
      REAL(ppm_kind_double)  :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_gmm_init',t0,info)
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_init',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_check_meshid(ppm_param_id_user,meshid,ppm_field_topoid, &
     &        lok,info)
          IF (.NOT. lok) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_init',  &
     &            'meshid out of range',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (Nest .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_init',  &
     &            'Nest must be >0!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF ((prec.NE.ppm_kind_single).AND.(prec.NE.ppm_kind_double)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_init',  &
     &            'Illegal precision specifiec!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      !  Nullify the pointers to allow proper inquiry of ASSOCIATED status.
      !-------------------------------------------------------------------------
      NULLIFY(gmm_ipos)
      NULLIFY(gmm_phis)
      NULLIFY(gmm_phid)
      NULLIFY(gmm_state2d)
      NULLIFY(gmm_state3d)
      !-------------------------------------------------------------------------
      !  Allocate sparse work space structure
      !-------------------------------------------------------------------------
      gmm_lsiz = Nest
      iopt = ppm_param_alloc_fit
      ldu(1) = gmm_lsiz
      IF (prec .EQ. ppm_kind_double) THEN
          CALL ppm_alloc(gmm_phid,ldu,iopt,info)
      ELSE
          CALL ppm_alloc(gmm_phis,ldu,iopt,info)
      ENDIF
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_gmm_init',     &
     &        'sparse data values GMM_PHI',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = ppm_dim + 1
      ldu(2) = gmm_lsiz
      CALL ppm_alloc(gmm_ipos,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_gmm_init',     &
     &        'sparse data locations GMM_IPOS',__LINE__,info)
          GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      !  Set topoid for all GMM operations
      !-------------------------------------------------------------------------
      gmm_topoid = ppm_field_topoid
      !-------------------------------------------------------------------------
      !  Translate meshid to internal numbering and store it
      !-------------------------------------------------------------------------
      gmm_meshid = ppm_meshid(gmm_topoid)%internal(meshid)
      !-------------------------------------------------------------------------
      !  Determine max extent of mesh in any sub
      !-------------------------------------------------------------------------
      maxxhi = 0
      maxyhi = 0
      maxzhi = 0
      DO i=1,ppm_nsublist(gmm_topoid)
          isub = ppm_isublist(i,gmm_topoid)
          IF (ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(1,isub).GT.maxxhi) &
     &        maxxhi = ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(1,isub)
          IF (ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(2,isub).GT.maxyhi) &
     &        maxyhi = ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(2,isub)
          IF (ppm_dim .GT. 2) THEN 
             IF (ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(3,isub).GT.maxzhi)&
     &           maxzhi = ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(3,isub)
          ENDIF
      ENDDO
      !-------------------------------------------------------------------------
      !  Memory increment step size 
      !-------------------------------------------------------------------------
      incr = MAX(maxxhi,maxyhi,maxzhi)
      incr = 10*incr*incr
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_gmm_init',t0,info)
      RETURN

      END SUBROUTINE ppm_gmm_init
