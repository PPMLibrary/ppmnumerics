      !-------------------------------------------------------------------------
      ! ppm_poisson_solve.f90
      !-------------------------------------------------------------------------
      !@mapping calls can be cleaned up
      !-------------------------------------------------------------------------
!!#define __ROUTINE ppm_poisson_solve
!!#define __DIM  3
!!#define __NCOM  3
!!#define __ZEROSI (/0,0,0/)
!#define __NOPE
      SUBROUTINE __ROUTINE(topoid,meshid,ppmpoisson,fieldin,fieldout,gstw,info,&
                         & tmpcase)

      USE ppm_module_map_field
      USE ppm_module_map_field_global
      USE ppm_module_map


      IMPLICIT NONE
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN)                                         :: topoid
      INTEGER, INTENT(IN)                                         :: meshid
      TYPE(ppm_poisson_plan),INTENT(INOUT)                        :: ppmpoisson
      !REAL(__PREC),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT)     :: fieldin
      REAL(__PREC),DIMENSION(:,:,:,:,:),POINTER                   :: fieldin
      !REAL(__PREC),DIMENSION(:,:,:,:,:),POINTER,INTENT(INOUT)     :: fieldout
      REAL(__PREC),DIMENSION(:,:,:,:,:),POINTER                   :: fieldout
      INTEGER,DIMENSION(__DIM),INTENT(IN)                         :: gstw
      INTEGER, INTENT(OUT)                                        :: info
      INTEGER,OPTIONAL,INTENT(IN)                                 :: tmpcase

      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      INTEGER,PARAMETER                 :: MK = __PREC
      REAL(__PREC)                      :: t0
      INTEGER                           :: isub,isubl
      INTEGER                           :: i,j,k
      INTEGER                           :: info2
      INTEGER                           :: presentcase
      REAL(__PREC)                      :: wdotk
      INTEGER                           :: gi,gj,gk
      REAL(__PREC)                      :: kx,ky,kz
      REAL(__PREC)                      :: phix,phiy,phiz
      REAL(__PREC)                      :: normfac
#ifndef __NOPE
INTEGER                           :: trank !@
trank =0
#endif
      !-------------------------------------------------------------------------
      ! Initialise routine
      !-------------------------------------------------------------------------
      CALL substart('ppm_poisson_solve',t0,info)

      !-------------------------------------------------------------------------
      ! Check if we run a different/temporary case
      !-------------------------------------------------------------------------
      IF (PRESENT(tmpcase)) THEN
         presentcase = tmpcase
      ELSE
         presentcase = ppmpoisson%case
      ENDIF
!@write(*,*) 'solving ....  1', ppm_rank

      !-------------------------------------------------------------------------
      !@ Perhaps check if ghostlayer suffices for a given fd stencil
      !-------------------------------------------------------------------------


      !-------------------------------------------------------------------------
      ! Perhaps allocate (and deallocate) arrays !@
      !-------------------------------------------------------------------------
#ifdef __VARMESH
  !@tmp only works for one subdomain (not parallel)
        write(*,*) 'fieldin', ppm_rank
        DO isub=1,ppmpoisson%nsublistxy
          isubl=ppmpoisson%isublistxy(isub)
          DO k=1,ppmpoisson%ndataxy(3,isubl)/2-1
            write(*,*) 'z',k
            DO i=1,ppmpoisson%ndataxy(1,isubl)/2-1
              DO j=1,ppmpoisson%ndataxy(2,isubl)/2-1
                write(*,'(A,E12.4,A,$)') '  (',fieldin(1,i,j,k,isub),')'
              END DO
              write(*,*)
            END DO
          END DO
        END DO
#endif
        !-----------------------------------------------------------------------
        ! Set the real xy slabs 0 (part of the 0 padding) for free-space
        !@ free-space calculations and reprojection may cause problems !why?
        !-----------------------------------------------------------------------
        IF (presentcase .EQ. ppm_poisson_grn_pois_fre) THEN
          ppmpoisson%fldxyr = 0.0_MK
        ENDIF
        !-----------------------------------------------------------------------
        ! Map data globally to the slabs (XY)
        ! This is where the vorticity is extended and padded with 0 for free-space
        !-----------------------------------------------------------------------
        !Initialise
        CALL ppm_map_field_global(&
        & topoid, &
        & ppmpoisson%topoidxy, &
        & meshid, &
        & ppmpoisson%meshidxy,info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_poisson_solve','Failed to initialise field mapping.',info2)
          GOTO 9999
        ENDIF

        !Push the data
        CALL ppm_map_field_push(&
        & topoid, &
        & meshid,fieldin,__NCOM,info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to push vector field.',info2)
          GOTO 9999
        ENDIF

        !Send
        CALL ppm_map_field_send(info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to send field.',info2)
          GOTO 9999
        ENDIF

        !Retrieve
        CALL ppm_map_field_pop(&
        & ppmpoisson%topoidxy, &
        & ppmpoisson%meshidxy,ppmpoisson%fldxyr, &
        & __NCOM,__ZEROSI,info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',info2)
          GOTO 9999
        ENDIF
!@write(*,*) 'solving ....  2', ppm_rank
#ifdef __VARMESH
!@tmp
      write(*,*) 'fldxyr', ppm_rank
      DO isub=1,ppmpoisson%nsublistxy
        isubl=ppmpoisson%isublistxy(isub)
        DO k=1,ppmpoisson%ndataxy(3,isubl)-1
          write(*,*) 'z',k
          DO i=1,ppmpoisson%ndataxy(1,isubl)-1
            DO j=1,ppmpoisson%ndataxy(2,isubl)-1
              write(*,'(A,E12.4,A,$)') '  (',ppmpoisson%fldxyr(1,i,j,k,isub),')'
            END DO
            write(*,*)
          END DO
        END DO
      END DO
#endif

      !-----------------------------------------------------------------------
      ! Do slab FFT (XY) - use the non-reduced topology !@what does this mean
      !-----------------------------------------------------------------------
      CALL ppm_fft_execute_2d(ppmpoisson%topoidxy,&
      & ppmpoisson%meshidxy, ppmpoisson%planfxy, &
      & ppmpoisson%fldxyr, ppmpoisson%fldxyc, &
      & info)
!@write(*,*) 'solving ....  3', ppm_rank

#ifdef __NOPE
      if (ppm_rank .EQ. trank)   THEN
      write(*,*) 'Rfftx12', ppm_rank
      DO isub=1,ppmpoisson%nsublistxy
        isubl=ppmpoisson%isublistxy(isub)
        DO k=1,ppmpoisson%ndataxy(3,isubl)-1
          write(*,*) 'z',k
            DO i=1,ppmpoisson%ndataxy(1,isubl)-1
          DO j=1,ppmpoisson%ndataxy(2,isubl)-1
              write(*,'(A,E12.4,A,$)') '  (',ppmpoisson%fldxyr(1,i,j,k,isub),')'
            END DO
            write(*,*)
          END DO
        END DO
      END DO
      write(*,*) 'Rffty12', ppm_rank
      DO isub=1,ppmpoisson%nsublistxy
        isubl=ppmpoisson%isublistxy(isub)
        DO k=1,ppmpoisson%ndataxy(3,isubl)-1
          write(*,*) 'z',k
            DO i=1,ppmpoisson%ndataxy(1,isubl)-1
          DO j=1,ppmpoisson%ndataxy(2,isubl)-1
              write(*,'(A,E12.4,A,$)') '  (',ppmpoisson%fldxyr(2,i,j,k,isub),')'
            END DO
            write(*,*)
          END DO
        END DO
      END DO
      write(*,*) 'Rfftz12', ppm_rank
      DO isub=1,ppmpoisson%nsublistxy
        isubl=ppmpoisson%isublistxy(isub)
        DO k=1,ppmpoisson%ndataxy(3,isubl)-1
          write(*,*) 'z',k
            DO i=1,ppmpoisson%ndataxy(1,isubl)-1
          DO j=1,ppmpoisson%ndataxy(2,isubl)-1
              write(*,'(A,E12.4,A,$)') '  (',ppmpoisson%fldxyr(3,i,j,k,isub),')'
            END DO
            write(*,*)
          END DO
        END DO
      END DO
      ENDIF
#endif


#ifdef __NOPE
      if (ppm_rank .EQ. trank)   THEN
      write(*,*) 'fftx12', ppm_rank
      DO isub=1,ppmpoisson%nsublistxy
        isubl=ppmpoisson%isublistxy(isub)
        DO k=1,ppmpoisson%ndataxy(3,isubl)-1
          write(*,*) 'z',k
            DO i=1,ppmpoisson%ndataxy(1,isubl)-1
          DO j=1,ppmpoisson%ndataxy(2,isubl)-1
              write(*,'(A,E12.4,E12.4,A,$)') '  (',ppmpoisson%fldxyc(1,i,j,k,isub),')'
            END DO
            write(*,*)
          END DO
        END DO
      END DO
      write(*,*) 'ffty12', ppm_rank
      DO isub=1,ppmpoisson%nsublistxy
        isubl=ppmpoisson%isublistxy(isub)
        DO k=1,ppmpoisson%ndataxy(3,isubl)-1
          write(*,*) 'z',k
            DO i=1,ppmpoisson%ndataxy(1,isubl)-1
          DO j=1,ppmpoisson%ndataxy(2,isubl)-1
              write(*,'(A,E12.4,E12.4,A,$)') '  (',ppmpoisson%fldxyc(2,i,j,k,isub),')'
            END DO
            write(*,*)
          END DO
        END DO
      END DO
      write(*,*) 'fftz12', ppm_rank
      DO isub=1,ppmpoisson%nsublistxy
        isubl=ppmpoisson%isublistxy(isub)
        DO k=1,ppmpoisson%ndataxy(3,isubl)-1
          write(*,*) 'z',k
            DO i=1,ppmpoisson%ndataxy(1,isubl)-1
          DO j=1,ppmpoisson%ndataxy(2,isubl)-1
              write(*,'(A,E12.4,E12.4,A,$)') '  (',ppmpoisson%fldxyc(3,i,j,k,isub),')'
            END DO
            write(*,*)
          END DO
        END DO
      END DO
      ENDIF
#endif
!@write(*,*) 'solving ....  4', ppm_rank

      !-----------------------------------------------------------------------
      ! Map to the pencils (Z)
      !-----------------------------------------------------------------------
      !Initialise
      CALL ppm_map_field_global(&
      & ppmpoisson%topoidxy, &
      & ppmpoisson%topoidz, &
      !& ppmpoisson%meshidxyc, & !@to be used when meshid>1 works
      & ppmpoisson%meshidxy, &
      & ppmpoisson%meshidz,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to initialise field mapping.',info2)
        GOTO 9999
      ENDIF

      !Push the data
      CALL ppm_map_field_push(&
      & ppmpoisson%topoidxy, &
      & ppmpoisson%meshidxyc,ppmpoisson%fldxyc,__NCOM,info)
      !& ppmpoisson%meshidxyc,ppmpoisson%fldxyc,__NCOM,info)!@to be used when meshid>1 works
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to push vector field.',info2)
        GOTO 9999
      ENDIF

      !Send
      CALL ppm_map_field_send(info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to send field.',info2)
        GOTO 9999
      ENDIF

      !Retrieve
      CALL ppm_map_field_pop(&
      & ppmpoisson%topoidz, &
      & ppmpoisson%meshidz,ppmpoisson%fldzc1, &
      & __NCOM,__ZEROSI,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',info2)
        GOTO 9999
      ENDIF
!@write(*,*) 'solving ....  5', ppm_rank

      !-----------------------------------------------------------------------
      ! Do pencil FFT (Z)
      !-----------------------------------------------------------------------
      CALL ppm_fft_execute_1d(ppmpoisson%topoidz,&
      & ppmpoisson%meshidz, ppmpoisson%planfz, &
      & ppmpoisson%fldzc1, ppmpoisson%fldzc2, &
      & info)

!@write(*,*) 'solving ....  6', ppm_rank
#ifdef __NOPE
      if (ppm_rank .EQ. trank)   THEN
            write(*,*)
            write(*,*) !ppmpoisson%fldzc1 = 0.0_MK
      write(*,*) 'fftx123', ppm_rank
      DO isub=1,ppmpoisson%nsublistz
        isubl=ppmpoisson%isublistz(isub)
        DO k=1,ppmpoisson%ndataz(3,isubl)-1
          write(*,*) 'z',k
            DO i=1,ppmpoisson%ndataz(1,isubl)-1
          DO j=1,ppmpoisson%ndataz(2,isubl)-1
              write(*,'(A,E12.4,E12.4,A,$)') '  (',ppmpoisson%fldzc2(1,i,j,k,isub),')'
            END DO
            write(*,*)
          END DO
        END DO
      END DO
      write(*,*) 'ffty123', ppm_rank
      DO isub=1,ppmpoisson%nsublistz
        isubl=ppmpoisson%isublistz(isub)
        DO k=1,ppmpoisson%ndataz(3,isubl)-1
          write(*,*) 'z',k
            DO i=1,ppmpoisson%ndataz(1,isubl)-1
          DO j=1,ppmpoisson%ndataz(2,isubl)-1
              write(*,'(A,E12.4,E12.4,A,$)') '  (',ppmpoisson%fldzc2(2,i,j,k,isub),')'
            END DO
            write(*,*)
          END DO
        END DO
      END DO
      write(*,*) 'fftz123', ppm_rank
      DO isub=1,ppmpoisson%nsublistz
        isubl=ppmpoisson%isublistz(isub)
        DO k=1,ppmpoisson%ndataz(3,isubl)-1
          write(*,*) 'z',k
            DO i=1,ppmpoisson%ndataz(1,isubl)-1
          DO j=1,ppmpoisson%ndataz(2,isubl)-1
              write(*,'(A,E12.4,E12.4,A,$)') '  (',ppmpoisson%fldzc2(3,i,j,k,isub),')'
            END DO
            write(*,*)
          END DO
        END DO
      END DO

      write(*,*)
      write(*,*)
      write(*,*) 'green', ppm_rank
      DO isub=1,ppmpoisson%nsublistz
        isubl=ppmpoisson%isublistz(isub)
        DO k=1,ppmpoisson%ndataz(3,isubl)-1
          write(*,*) 'z',k
            DO i=1,ppmpoisson%ndataz(1,isubl)-1
          DO j=1,ppmpoisson%ndataz(2,isubl)-1
              write(*,'(E12.4,$)') REAL(ppmpoisson%fldgrnr(i,j,k,isub))
            END DO
            write(*,*)
          END DO
        END DO
      END DO
      ENDIF
#endif
!@write(*,*) 'solving ....  7', ppm_rank

      !-----------------------------------------------------------------------
      ! Apply the periodic Green's function
      !-----------------------------------------------------------------------
      IF (presentcase .EQ. ppm_poisson_grn_pois_per) THEN
        DO isub=1,ppmpoisson%nsublistz
          isubl=ppmpoisson%isublistz(isub)
          DO k=1,ppmpoisson%ndataz(3,isubl)
            DO j=1,ppmpoisson%ndataz(2,isubl)
              DO i=1,ppmpoisson%ndataz(1,isubl)
                ppmpoisson%fldzc2(1,i,j,k,isub) = ppmpoisson%fldgrnr( i,j,k,isub)*&
                                                & ppmpoisson%fldzc2(1,i,j,k,isub)
                ppmpoisson%fldzc2(2,i,j,k,isub) = ppmpoisson%fldgrnr( i,j,k,isub)*&
                                                & ppmpoisson%fldzc2(2,i,j,k,isub)
                ppmpoisson%fldzc2(3,i,j,k,isub) = ppmpoisson%fldgrnr( i,j,k,isub)*&
                                                & ppmpoisson%fldzc2(3,i,j,k,isub)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      !-----------------------------------------------------------------------
      ! Apply the free-space Green's function
      !-----------------------------------------------------------------------
      ELSE IF (presentcase .EQ. ppm_poisson_grn_pois_fre) THEN
        DO isub=1,ppmpoisson%nsublistz
          isubl=ppmpoisson%isublistz(isub)
          DO k=1,ppmpoisson%ndataz(3,isubl)
            DO j=1,ppmpoisson%ndataz(2,isubl)
              DO i=1,ppmpoisson%ndataz(1,isubl)
                ppmpoisson%fldzc2(1,i,j,k,isub) = ppmpoisson%fldgrnc( i,j,k,isub)*&
                                                & ppmpoisson%fldzc2(1,i,j,k,isub)
                ppmpoisson%fldzc2(2,i,j,k,isub) = ppmpoisson%fldgrnc( i,j,k,isub)*&
                                                & ppmpoisson%fldzc2(2,i,j,k,isub)
                ppmpoisson%fldzc2(3,i,j,k,isub) = ppmpoisson%fldgrnc( i,j,k,isub)*&
                                                & ppmpoisson%fldzc2(3,i,j,k,isub)
                !!ppmpoisson%fldzc2(1,i,j,k,isub) = ppmpoisson%fldgrnc( i,j,k,isub)
                !!ppmpoisson%fldzc2(2,i,j,k,isub) = ppmpoisson%fldgrnc( i,j,k,isub)
                !!ppmpoisson%fldzc2(3,i,j,k,isub) = ppmpoisson%fldgrnc( i,j,k,isub)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!@write(*,*) 'solving ....  8', ppm_rank
      !-----------------------------------------------------------------------
      ! Spectral derivatives
      ! normkx, etc contains 2pi/Lx
      !-----------------------------------------------------------------------
      IF (ppmpoisson%derivatives .EQ. ppm_poisson_drv_curl_sp .AND.&
         & (presentcase .EQ. ppm_poisson_grn_pois_per .OR. &
            presentcase .EQ. ppm_poisson_grn_pois_fre)) THEN
        write(*,*) 'curl fft style! not working yet though!'
        !remembering to normalize the FFT
        normfac = 1.0_MK/ REAL((ppmpoisson%nmz(1)-1)* & !vertex
                             & (ppmpoisson%nmz(2)-1)* &
                             & (ppmpoisson%nmz(3)-1),MK)
        DO isub=1,ppmpoisson%nsublistz
          isubl=ppmpoisson%isublistz(isub)
          DO k=1,ppmpoisson%ndataz(3,isubl)
            gk = k - 1 + (ppmpoisson%istartz(3,isubl)-1)
            IF (gk .GT. (ppmpoisson%nmz(3)-1)/2) gk = gk-(ppmpoisson%nmz(3)-1)
            kz = CMPLX(0.0_MK,gk*ppmpoisson%normkz,MK)
            DO j=1,ppmpoisson%ndataz(2,isubl)
              gj = j - 1 + (ppmpoisson%istartz(2,isubl)-1)
              IF (gj .GT. (ppmpoisson%nmz(2)-1)/2) gj = gj-(ppmpoisson%nmz(2)-1)
              ky = CMPLX(0.0_MK,gj*ppmpoisson%normky,MK)
              DO i=1,ppmpoisson%ndataz(1,isubl)
                gi = i - 1 + (ppmpoisson%istartz(1,isubl)-1)
                IF (gi .GT. (ppmpoisson%nmz(1)-1)/2) gi = gi-(ppmpoisson%nmz(1)-1)
                kx = CMPLX(0.0_MK,gi*ppmpoisson%normkx,MK)

                phix = ppmpoisson%fldzc2(1,i,j,k,isub)
                phiy = ppmpoisson%fldzc2(2,i,j,k,isub)
                phiz = ppmpoisson%fldzc2(3,i,j,k,isub)

                !maybe normfac on kx,ky,kz?
                ppmpoisson%fldzc2(1,i,j,k,isub) = normfac*(ky*phiz-kz*phiy)
                ppmpoisson%fldzc2(2,i,j,k,isub) = normfac*(kz*phix-kx*phiz)
                ppmpoisson%fldzc2(3,i,j,k,isub) = normfac*(kx*phiy-ky*phix)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      !-----------------------------------------------------------------------
      ! Vorticity re-projection
      !@ has the domain length really been included in these wave numbers
      !-----------------------------------------------------------------------
      ELSE IF (presentcase .EQ. ppm_poisson_grn_reprojec) THEN
        !remembering to normalize the FFT
        normfac = 1.0_MK/ REAL((ppmpoisson%nmz(1)-1)* & !vertex
                              & (ppmpoisson%nmz(2)-1)* &
                              & (ppmpoisson%nmz(3)-1),MK)
        write(*,*) 'reprojection ppm style' !@
        DO isub=1,ppmpoisson%nsublistz
          isubl=ppmpoisson%isublistz(isub)
          DO k=1,ppmpoisson%ndataz(3,isubl)
            gk = k - 1 + (ppmpoisson%istartz(3,isubl)-1)
            IF (gk .GT. (ppmpoisson%nmz(3)-1)/2) gk = gk-(ppmpoisson%nmz(3)-1)
            kz = REAL(gk,MK)
            DO j=1,ppmpoisson%ndataz(2,isubl)
              gj = j - 1 + (ppmpoisson%istartz(2,isubl)-1)
              IF (gj .GT. (ppmpoisson%nmz(2)-1)/2) gj = gj-(ppmpoisson%nmz(2)-1)
              ky = REAL(gj,MK)
              DO i=1,ppmpoisson%ndataz(1,isubl)
                gi = i - 1 + (ppmpoisson%istartz(1,isubl)-1)
                IF (gi .GT. (ppmpoisson%nmz(1)-1)/2) gi = gi-(ppmpoisson%nmz(1)-1)
                kx = REAL(gi,MK)

                IF (gi .EQ. 0 .AND. gj .EQ. 0 .AND. gk .EQ. 0) THEN
                  wdotk = 0.0_mk
                ELSE
                  wdotk = (ppmpoisson%fldzc2(1,i,j,k,isub) * kx +  &
                        &  ppmpoisson%fldzc2(2,i,j,k,isub) * ky +  &
                        &  ppmpoisson%fldzc2(3,i,j,k,isub) * kz) / &
                        &  (kx*kx+ky*ky+kz*kz)
                ENDIF

                ppmpoisson%fldzc2(1,i,j,k,isub) = &
                  & (ppmpoisson%fldzc2(1,i,j,k,isub) - wdotk*kx)*normfac
                ppmpoisson%fldzc2(2,i,j,k,isub) = &
                  & (ppmpoisson%fldzc2(2,i,j,k,isub) - wdotk*ky)*normfac
                ppmpoisson%fldzc2(3,i,j,k,isub) = &
                  & (ppmpoisson%fldzc2(3,i,j,k,isub) - wdotk*kz)*normfac
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!@write(*,*) 'solving ....  9', ppm_rank
      !-----------------------------------------------------------------------
      ! IFFT pencil (Z)
      !-----------------------------------------------------------------------
      CALL ppm_fft_execute_1d(ppmpoisson%topoidz,&
      & ppmpoisson%meshidz, ppmpoisson%planbz, &
      & ppmpoisson%fldzc2, ppmpoisson%fldzc1, &
      & info)

!@write(*,*) 'solving ....  11', ppm_rank
#ifdef __NOPE
      if (ppm_rank .EQ. trank)   THEN
            write(*,*)
            write(*,*)
      write(*,*) 'ifftx12', ppm_rank
      DO isub=1,ppmpoisson%nsublistz
        isubl=ppmpoisson%isublistz(isub)
        DO k=1,ppmpoisson%ndataz(3,isubl)-1
          write(*,*) 'z',k
            DO i=1,ppmpoisson%ndataz(1,isubl)-1
          DO j=1,ppmpoisson%ndataz(2,isubl)-1
              write(*,'(A,E12.4,E12.4,A,$)') '  (',ppmpoisson%fldzc1(1,i,j,k,isub),')'
            END DO
            write(*,*)
          END DO
        END DO
      END DO
      write(*,*) 'iffty12', ppm_rank
      DO isub=1,ppmpoisson%nsublistz
        isubl=ppmpoisson%isublistz(isub)
        DO k=1,ppmpoisson%ndataz(3,isubl)-1
          write(*,*) 'z',k
            DO i=1,ppmpoisson%ndataz(1,isubl)-1
          DO j=1,ppmpoisson%ndataz(2,isubl)-1
              write(*,'(A,E12.4,E12.4,A,$)') '  (',ppmpoisson%fldzc1(2,i,j,k,isub),')'
            END DO
            write(*,*)
          END DO
        END DO
      END DO
      write(*,*) 'ifftz12', ppm_rank
      DO isub=1,ppmpoisson%nsublistz
        isubl=ppmpoisson%isublistz(isub)
        DO k=1,ppmpoisson%ndataz(3,isubl)-1
          write(*,*) 'z',k
            DO i=1,ppmpoisson%ndataz(1,isubl)-1
          DO j=1,ppmpoisson%ndataz(2,isubl)-1
              write(*,'(A,E12.4,E12.4,A,$)') '  (',ppmpoisson%fldzc1(3,i,j,k,isub),')'
            END DO
            write(*,*)
          END DO
        END DO
      END DO
      ENDIF
#endif

      !-----------------------------------------------------------------------
      ! Map back to slabs (XY)
      !-----------------------------------------------------------------------
      !Initialise
      CALL ppm_map_field_global(&
      & ppmpoisson%topoidz, &
      & ppmpoisson%topoidxy, &
      & ppmpoisson%meshidz, &
      & ppmpoisson%meshidxyc,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to initialise field mapping.',info2)
        GOTO 9999
      ENDIF

      !Push the data
      CALL ppm_map_field_push(&
      & ppmpoisson%topoidz, &
      & ppmpoisson%meshidz,ppmpoisson%fldzc1,__NCOM,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to push vector field.',info2)
        GOTO 9999
      ENDIF

      !Send
      CALL ppm_map_field_send(info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to send field.',info2)
        GOTO 9999
      ENDIF

      !Retrieve
      CALL ppm_map_field_pop(&
      & ppmpoisson%topoidxy, &
      & ppmpoisson%meshidxyc,ppmpoisson%fldxyc, &
      & __NCOM,__ZEROSI,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',info2)
        GOTO 9999
      ENDIF

!@write(*,*) 'solving ....  12', ppm_rank
      !-----------------------------------------------------------------------
      ! IFFT (XY) use the non-reduced topology
      !-----------------------------------------------------------------------
      CALL ppm_fft_execute_2d(ppmpoisson%topoidxy,&
      & ppmpoisson%meshidxy, ppmpoisson%planbxy, &
      & ppmpoisson%fldxyc, ppmpoisson%fldxyr, &
      & info)

!@write(*,*) 'solving ....  13', ppm_rank
#ifdef __NOPE
      if (ppm_rank .EQ. trank)   THEN
            write(*,*)
            write(*,*)
      write(*,*) 'ifftx', ppm_rank
      DO isub=1,ppmpoisson%nsublistxy
        isubl=ppmpoisson%isublistxy(isub)
        DO k=1,ppmpoisson%ndataxy(3,isubl)-1
          write(*,*) 'z',k
            DO i=1,ppmpoisson%ndataxy(1,isubl)-1
          DO j=1,ppmpoisson%ndataxy(2,isubl)-1
              write(*,'(E12.4,$)') ppmpoisson%fldxyr(1,i,j,k,isub)
            END DO
            write(*,*)
          END DO
        END DO
      END DO
      write(*,*) 'iffty', ppm_rank
      DO isub=1,ppmpoisson%nsublistxy
        isubl=ppmpoisson%isublistxy(isub)
        DO k=1,ppmpoisson%ndataxy(3,isubl)-1
          write(*,*) 'z',k
            DO i=1,ppmpoisson%ndataxy(1,isubl)-1
          DO j=1,ppmpoisson%ndataxy(2,isubl)-1
              write(*,'(E12.4,$)') ppmpoisson%fldxyr(2,i,j,k,isub)
            END DO
            write(*,*)
          END DO
        END DO
      END DO
      write(*,*) 'ifftz', ppm_rank
      DO isub=1,ppmpoisson%nsublistxy
        isubl=ppmpoisson%isublistxy(isub)
        DO k=1,ppmpoisson%ndataxy(3,isubl)-1
          write(*,*) 'z',k
            DO i=1,ppmpoisson%ndataxy(1,isubl)-1
          DO j=1,ppmpoisson%ndataxy(2,isubl)-1
              write(*,'(E12.4,$)') ppmpoisson%fldxyr(3,i,j,k,isub)
            END DO
            write(*,*)
          END DO
        END DO
      END DO
      ENDIF
#endif

#ifdef __VARMESH
!@tmp
      write(*,*) 'fldxyr2', ppm_rank
      DO isub=1,ppmpoisson%nsublistxy
        isubl=ppmpoisson%isublistxy(isub)
        DO k=1,ppmpoisson%ndataxy(3,isubl)
          write(*,*) 'z',k
          DO i=1,ppmpoisson%ndataxy(1,isubl)
            DO j=1,ppmpoisson%ndataxy(2,isubl)
              write(*,'(A,E12.4,A,$)') '  (',ppmpoisson%fldxyr(1,i,j,k,isub),')'
            END DO
            write(*,*)
          END DO
        END DO
      END DO
#endif
      !-----------------------------------------------------------------------
      ! Map back to standard topology (XYZ)
      !-----------------------------------------------------------------------
      !Initialise
      CALL ppm_map_field_global(&
      & ppmpoisson%topoidxy, &
      & topoid, &
      & ppmpoisson%meshidxy, &
      & meshid,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to initialise field mapping.',info2)
        GOTO 9999
      ENDIF
!@write(*,*) 'solving ....  13a', ppm_rank

      !Push the data
      CALL ppm_map_field_push(&
      & ppmpoisson%topoidxy, &
      & ppmpoisson%meshidxy,ppmpoisson%fldxyr,__NCOM,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to push vector field.',info2)
        GOTO 9999
      ENDIF
!@write(*,*) 'solving ....  13b', ppm_rank, ppmpoisson%topoidxy, topoid, ppmpoisson%meshidxy, meshid

      !Send
      CALL ppm_map_field_send(info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to send field.',info2)
        GOTO 9999
      ENDIF

!@write(*,*) 'solving ....  14', ppm_rank
      !-------------------------------------------------------------------------
      ! FINAL RETRIEVE - Here we do different things depending on the task
      ! i.e. the receiver varies
      !-------------------------------------------------------------------------
      IF ((ppmpoisson%derivatives .EQ. ppm_poisson_drv_curl_fd2 .OR. &
           ppmpoisson%derivatives .EQ. ppm_poisson_drv_curl_fd4) .AND. &
         & (presentcase .EQ. ppm_poisson_grn_pois_per .OR. &  !@these may be unnecessary - perhaps just the derive value. Or maybe not: in case of vorticity reprojection we could get lost
            presentcase .EQ. ppm_poisson_grn_pois_fre)) THEN
!@write(*,*) 'solving ....  14a', ppm_rank
        CALL ppm_map_field_pop(&
        & topoid, &
        & meshid,ppmpoisson%drv_vr, &
        & __NCOM,gstw,info) !!! gst
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',info2)
          GOTO 9999
        ENDIF
!@write(*,*) 'solving ....  14b', ppm_rank
        !-------------------------------------------------------------------------
        ! Ghost the temporary array for derivatives (drv_vr)
        !-------------------------------------------------------------------------
!@write(*,*) 'solving ....  14c', ppm_rank
        CALL ppm_map_field_ghost_get(topoid,meshid,gstw,info)
        IF (info .NE. 0) THEN
           CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to initialise ghosts.',info2)
           GOTO 9999
        ENDIF
!@write(*,*) 'solving ....  14d', ppm_rank
        CALL ppm_map_field_push(topoid,meshid,ppmpoisson%drv_vr,__NCOM,info)
        IF (info .NE. 0) THEN
           CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to push ghosts.',info2)
           GOTO 9999
        ENDIF
!@write(*,*) 'solving ....  14e', ppm_rank
        CALL ppm_map_field_send(info)
        IF (info .NE. 0) THEN
           CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to send ghosts.',info2)
           GOTO 9999
        ENDIF
!@write(*,*) 'solving ....  14f', ppm_rank
        CALL ppm_map_field_pop(topoid,meshid,ppmpoisson%drv_vr,__NCOM,gstw,info)
        IF (info .NE. 0) THEN
           CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop ghosts.',info2)
           GOTO 9999
        ENDIF

      ELSE
!@write(*,*) 'solving ....  14g', ppm_rank
        CALL ppm_map_field_pop(&
        & topoid, &
        & meshid,fieldout, &
        & __NCOM,gstw,info) !!! gst
!@write(*,*) 'solving ....  14h', ppm_rank
      ENDIF
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',info2)
        GOTO 9999
      ENDIF

!@write(*,*) 'solving ....  15', ppm_rank
#ifdef __VARMESH
!@tmp only works for one subdomain (not parallel)
      write(*,*) 'fieldout', ppm_rank
      DO isub=1,ppmpoisson%nsublistxy
        isubl=ppmpoisson%isublistxy(isub)
        DO k=1,ppmpoisson%ndataxy(3,isubl)/2-1
          write(*,*) 'z',k
          DO i=1,ppmpoisson%ndataxy(1,isubl)/2-1
            DO j=1,ppmpoisson%ndataxy(2,isubl)/2-1
              write(*,'(A,E12.4,A,$)') '  (',fieldout(1,i,j,k,isub),')'
            END DO
            write(*,*)
          END DO
        END DO
      END DO
#endif

      !-------------------------------------------------------------------------
      !@To come: treat ghost layer to make FD stencils work
      !-------------------------------------------------------------------------
      IF (ppmpoisson%derivatives .EQ. ppm_poisson_drv_curl_fd2 .AND.&
         & (presentcase .EQ. ppm_poisson_grn_pois_fre)) THEN
         CALL ppm_poisson_extrapolateghost(topoid,meshid,ppmpoisson%drv_vr,&
                                     & 2,4,gstw,info)
      ENDIF
      IF (ppmpoisson%derivatives .EQ. ppm_poisson_drv_curl_fd4 .AND.&
         & (presentcase .EQ. ppm_poisson_grn_pois_fre)) THEN
         CALL ppm_poisson_extrapolateghost(topoid,meshid,ppmpoisson%drv_vr,&
                                     & 2,4,gstw,info)
      ENDIF

!@write(*,*) 'solving ....  16', ppm_rank
      !-------------------------------------------------------------------------
      ! Optionally do derivatives
      ! Perhaps make ppm_poisson_fd take _none as argument. Then maybe no
      ! if-statement is required
      !-------------------------------------------------------------------------
      IF (ppmpoisson%derivatives .EQ. ppm_poisson_drv_curl_fd2 .AND.&
         & (presentcase .EQ. ppm_poisson_grn_pois_per .OR. &  !@these may be unnecessary - perhaps just the derive value
            presentcase .EQ. ppm_poisson_grn_pois_fre)) THEN
         CALL ppm_poisson_fd(topoid,meshid,ppmpoisson%drv_vr,fieldout,&
                           & ppm_poisson_drv_curl_fd2,info)
      ENDIF
      IF (ppmpoisson%derivatives .EQ. ppm_poisson_drv_curl_fd4 .AND.&
         & (presentcase .EQ. ppm_poisson_grn_pois_per .OR. &  !@these may be unnecessary - perhaps just the derive value
            presentcase .EQ. ppm_poisson_grn_pois_fre)) THEN
         CALL ppm_poisson_fd(topoid,meshid,ppmpoisson%drv_vr,fieldout,&
                           & ppm_poisson_drv_curl_fd4,info)
      ENDIF

!@write(*,*) 'solving ....  17', ppm_rank
      !-------------------------------------------------------------------------
      ! Finally ghost the velocity/stream function field before returning it
      ! Also extrapolate if freespace
      !-------------------------------------------------------------------------
      CALL ppm_map_field_ghost_get(topoid,meshid,gstw,info)
      CALL ppm_map_field_push(topoid,meshid,fieldout,__NCOM,info)
      CALL ppm_map_field_send(info)
      CALL ppm_map_field_pop(topoid,meshid,fieldout,__NCOM,gstw,info)
      IF (presentcase .EQ. ppm_poisson_grn_pois_fre) THEN
         CALL ppm_poisson_extrapolateghost(topoid,meshid,fieldout,&
                                     & 2,4,gstw,info)
      ENDIF

!@write(*,*) 'solving ....  18', ppm_rank
      !-------------------------------------------------------------------------
      ! Perhaps allocate (and deallocate) arrays !@
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_poisson_solve',t0,info)
      RETURN

      END SUBROUTINE __ROUTINE

