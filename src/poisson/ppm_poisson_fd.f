      !-------------------------------------------------------------------------
      ! ppm_poisson_fd.f90
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      SUBROUTINE __ROUTINE(topoid,meshid,fieldin,fieldout,dtype,info)

      USE ppm_module_topo_get

      IMPLICIT NONE
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN)                                         :: topoid
      INTEGER, INTENT(IN)                                         :: meshid
      REAL(__PREC),DIMENSION(:,:,:,:,:),POINTER                   :: fieldin
      REAL(__PREC),DIMENSION(:,:,:,:,:),POINTER                   :: fieldout
      !INTEGER,DIMENSION(__DIM),INTENT(IN)                         :: gstw
      INTEGER, INTENT(IN)                                         :: dtype
      INTEGER, INTENT(OUT)                                        :: info

      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      INTEGER,PARAMETER                 :: MK = __PREC
      REAL(__PREC)                      :: t0
      TYPE(ppm_t_topo),POINTER          :: topology
      TYPE(ppm_t_equi_mesh)             :: mesh
      REAL(__PREC)                      :: dx,dy,dz
      REAL(__PREC)                      :: facx,facy,facz
      INTEGER                           :: isub,isubl
      INTEGER                           :: i,j,k

      !-------------------------------------------------------------------------
      ! Initialise routine
      !-------------------------------------------------------------------------
      CALL substart('ppm_poisson_fd',t0,info)
      !-------------------------------------------------------------------------
      ! Get topology and mesh values
      !-------------------------------------------------------------------------
      CALL ppm_topo_get(topoid,topology,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_poisson_init_predef','Failed to get topology.',isub)
        GOTO 9999
      ENDIF
      mesh  = topology%mesh(meshid)

      dx = (topology%max_physd(1)-topology%min_physd(1))/REAL(mesh%nm(1)-1)
      dy = (topology%max_physd(2)-topology%min_physd(2))/REAL(mesh%nm(2)-1)
      dz = (topology%max_physd(3)-topology%min_physd(3))/REAL(mesh%nm(3)-1)

      !-----------------------------------------------------------------------
      ! Do the finite difference calculation
      !-----------------------------------------------------------------------
      !-----------------------------------------------------------------------
      ! Curl, 2nd order FD
      !-----------------------------------------------------------------------
      IF (dtype .EQ. ppm_poisson_drv_curl_fd2) THEN
        facx = 1.0_MK/(2.0_MK*dx)
        facy = 1.0_MK/(2.0_MK*dy)
        facz = 1.0_MK/(2.0_MK*dz)

        DO isub=1,topology%nsublist
          isubl=topology%isublist(isub)
          DO k=1,mesh%nnodes(3,isubl)
            DO j=1,mesh%nnodes(2,isubl)
              DO i=1,mesh%nnodes(1,isubl)
                fieldout(1,i,j,k,isub) = &
                &  facz*(fieldin(2,i  ,j  ,k+1,isub)- &
                       & fieldin(2,i  ,j  ,k-1,isub)) &
                & -facy*(fieldin(3,i  ,j+1,k  ,isub)- &
                       & fieldin(3,i  ,j-1,k  ,isub))
                fieldout(2,i,j,k,isub) = &
                &  facx*(fieldin(3,i+1,j  ,k  ,isub)- &
                       & fieldin(3,i-1,j  ,k  ,isub)) &
                & -facz*(fieldin(1,i  ,j  ,k+1,isub)- &
                       & fieldin(1,i  ,j  ,k-1,isub))
                fieldout(3,i,j,k,isub) = &
                &  facy*(fieldin(1,i  ,j+1,k  ,isub)- &
                       & fieldin(1,i  ,j-1,k  ,isub)) &
                & -facx*(fieldin(2,i+1,j  ,k  ,isub)- &
                       & fieldin(2,i-1,j  ,k  ,isub))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      !-----------------------------------------------------------------------
      ! Curl, 4th order FD - untested
      !-----------------------------------------------------------------------
      ELSE IF (dtype .EQ. ppm_poisson_drv_curl_fd4) THEN
        facx = 1.0_MK/(12.0_MK*dx)
        facy = 1.0_MK/(12.0_MK*dy)
        facz = 1.0_MK/(12.0_MK*dz)

        DO isub=1,topology%nsublist
          isubl=topology%isublist(isub)
          DO k=1,mesh%nnodes(3,isubl)
            DO j=1,mesh%nnodes(2,isubl)
              DO i=1,mesh%nnodes(1,isubl)
                fieldout(1,i,j,k,isub) = &
                &  facz*(         -fieldin(2,i  ,j  ,k+2,isub)  &
                       &   +8.0_MK*fieldin(2,i  ,j  ,k+1,isub)  &
                       &   -8.0_MK*fieldin(2,i  ,j  ,k-1,isub)  &
                       &          +fieldin(2,i  ,j  ,k-2,isub)) &
                & -facy*(         -fieldin(3,i  ,j+2,k  ,isub)  &
                       &   +8.0_MK*fieldin(3,i  ,j+1,k  ,isub)  &
                       &   -8.0_MK*fieldin(3,i  ,j-1,k  ,isub)  &
                       &          +fieldin(3,i  ,j-2,k  ,isub))
                fieldout(2,i,j,k,isub) = &
                &  facx*(         -fieldin(3,i+2,j  ,k  ,isub)  &
                       &   +8.0_MK*fieldin(3,i+1,j  ,k  ,isub)  &
                       &   -8.0_MK*fieldin(3,i-1,j  ,k  ,isub)  &
                       &          +fieldin(3,i-2,j  ,k  ,isub)) &
                & -facz*(         -fieldin(1,i  ,j  ,k+2,isub)  &
                       &   +8.0_MK*fieldin(1,i  ,j  ,k+1,isub)  &
                       &   -8.0_MK*fieldin(1,i  ,j  ,k-1,isub)  &
                       &          +fieldin(1,i  ,j  ,k-2,isub))
                fieldout(3,i,j,k,isub) = &
                &  facy*(         -fieldin(1,i  ,j+2,k  ,isub)  &
                       &   +8.0_MK*fieldin(1,i  ,j+1,k  ,isub)  &
                       &   -8.0_MK*fieldin(1,i  ,j-1,k  ,isub)  &
                       &          +fieldin(1,i  ,j-2,k  ,isub)) &
                & -facx*(         -fieldin(2,i+2,j  ,k  ,isub)  &
                       &   +8.0_MK*fieldin(2,i+1,j  ,k  ,isub)  &
                       &   -8.0_MK*fieldin(2,i-1,j  ,k  ,isub)  &
                       &          +fieldin(2,i-2,j  ,k  ,isub))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_poisson_fd',t0,info)
      RETURN

      END SUBROUTINE __ROUTINE
