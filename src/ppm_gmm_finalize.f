      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_gmm_finalize
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_gmm_finalize(info)
      !!! This routine finalizes the ppm_gmm module and deallocates all
      !!! data structures.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_gmm
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER, INTENT(  OUT) :: info
      !!! Return status. 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)  :: ldu
      INTEGER                :: iopt
      REAL(ppm_kind_double)  :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_gmm_finalize',t0,info)
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_finalize',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      !  Deallocate work space structures
      !-------------------------------------------------------------------------
      gmm_lsiz = -1
      iopt = ppm_param_dealloc
      CALL ppm_alloc(gmm_phis,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_finalize',     &
     &        'sparse data values GMM_PHIS',__LINE__,info)
      ENDIF
      CALL ppm_alloc(gmm_phid,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_finalize',     &
     &        'sparse data values GMM_PHID',__LINE__,info)
      ENDIF
      CALL ppm_alloc(gmm_ipos,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_finalize',     &
     &        'sparse data locations GMM_IPOS',__LINE__,info)
      ENDIF
      CALL ppm_alloc(gmm_clod,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_finalize',     &
     &        'close point locations GMM_CLOD',__LINE__,info)
      ENDIF
      CALL ppm_alloc(gmm_clos,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_finalize',     &
     &        'close point locations GMM_CLOS',__LINE__,info)
      ENDIF
      CALL ppm_alloc(gmm_clod2,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_finalize',     &
     &        'close point locations GMM_CLOD2',__LINE__,info)
      ENDIF
      CALL ppm_alloc(gmm_clos2,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_finalize',     &
     &        'close point locations GMM_CLOS2',__LINE__,info)
      ENDIF
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_gmm_finalize',t0,info)
      RETURN

      END SUBROUTINE ppm_gmm_finalize
