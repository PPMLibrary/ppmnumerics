      !-------------------------------------------------------------------------
      ! ppm_poisson_solve.f90
      !-------------------------------------------------------------------------
      !@mapping calls can be cleaned up
      !-------------------------------------------------------------------------
      #define __ROUTINE ppm_poisson_solve
      #define __DIM  3
      #define __NCOM  3
      #define __ZEROSI (/0,0,0/)
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


      !-------------------------------------------------------------------------
      ! Perhaps allocate (and deallocate) arrays !@
      !-------------------------------------------------------------------------


      !------------------------------------------------------------------------
      ! Map data globally to the slabs (XY)
      !------------------------------------------------------------------------
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

      !!ppmpoisson%fldxyc= 0.0_MK
      !!ppmpoisson%fldzc1= 0.0_MK
      !!ppmpoisson%fldzc2= 0.0_MK
      !-----------------------------------------------------------------------
      ! Do slab FFT (XY) - use the non-reduced topology
      !-----------------------------------------------------------------------
      CALL ppm_fft_execute_2d(ppmpoisson%topoidxy,&
      & ppmpoisson%meshidxy, ppmpoisson%planfxy, &
      & ppmpoisson%fldxyr, ppmpoisson%fldxyc, &
      & info)

      !DO isub=1,ppmpoisson%nsublistxy
        !isubl=ppmpoisson%isublistxy(isub)
        !!Copy x/y boundaries
        !DO k=1,ppmpoisson%ndataxy(3,isubl)-1
          !DO i=1,ppmpoisson%ndataxy(1,isubl)
            !ppmpoisson%fldxyc(1,i,1,k,isub) = ppmpoisson%fldxyc(1,i,ppmpoisson%ndataxy(2,isubl),k,isub)
          !END DO
          !DO j=1,ppmpoisson%ndataxy(2,isubl)
            !ppmpoisson%fldxyc(1,1,j,k,isub) = ppmpoisson%fldxyc(1,ppmpoisson%ndataxy(1,isubl),j,k,isub)
          !END DO
        !END DO
      !END DO

      !!ppmpoisson%fldxyr = 0.0_MK
      !!fieldin           = 0.0_MK
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

      !#ifdef __NOPE
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

      !-----------------------------------------------------------------------
      ! Do pencil FFT (Z)
      !-----------------------------------------------------------------------
      CALL ppm_fft_execute_1d(ppmpoisson%topoidz,&
      & ppmpoisson%meshidz, ppmpoisson%planfz, &
      & ppmpoisson%fldzc1, ppmpoisson%fldzc2, &
      & info)

      #ifdef __NOPE
      if (ppm_rank .EQ. trank)   THEN
            write(*,*)
            write(*,*)
      !ppmpoisson%fldzc1 = 0.0_MK
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

      !-----------------------------------------------------------------------
      ! Apply the Green's function
      ! ... depending on the case
      !-----------------------------------------------------------------------
      !!kill?: ppmpoisson%fldgrnc = -10000000.0_MK
      !!kill?: DO isub=1,ppmpoisson%nsublistz
        !!kill?: isubl=ppmpoisson%isublistz(isub)
        !!kill?: DO k=1,ppmpoisson%ndataz(3,isubl)
          !!kill?: DO j=1,ppmpoisson%ndataz(2,isubl)
            !!kill?: DO i=1,ppmpoisson%ndataz(1,isubl)
              !!kill?: ppmpoisson%fldgrnc(1,i,j,k,isub) = ppmpoisson%fldzc2(1,i,j,k,isub)
              !!kill?: ppmpoisson%fldgrnc(2,i,j,k,isub) = ppmpoisson%fldzc2(2,i,j,k,isub)
              !!kill?: ppmpoisson%fldgrnc(3,i,j,k,isub) = ppmpoisson%fldzc2(3,i,j,k,isub)
              !!kill?: !ppmpoisson%fldzc2(1,i,j,k,isub) = ppmpoisson%fldgrnr(i,j,k,isub)
              !!kill?: !ppmpoisson%fldzc2(2,i,j,k,isub) = ppmpoisson%fldgrnr(i,j,k,isub)
              !!kill?: !ppmpoisson%fldzc2(3,i,j,k,isub) = ppmpoisson%fldgrnr(i,j,k,isub)
              !!kill?: !!ppmpoisson%fldzc2(1,i,j,k,isub) = CMPLX((REAL(i*j*k)),0.0)
              !!kill?: !!ppmpoisson%fldzc2(2,i,j,k,isub) = CMPLX((REAL(i*j*k)),0.0)
              !!kill?: !!ppmpoisson%fldzc2(3,i,j,k,isub) = CMPLX((REAL(i*j*k)),0.0)
              !!kill?: !fieldout(1,i,j,k,isub) = REAL(ppmpoisson%fldgrnr(i,j,k,isub))
              !!kill?: !fieldout(2,i,j,k,isub) = REAL(ppmpoisson%fldgrnr(i,j,k,isub))
              !!kill?: !fieldout(3,i,j,k,isub) = REAL(ppmpoisson%fldgrnr(i,j,k,isub))
            !!kill?: ENDDO
          !!kill?: ENDDO
        !!kill?: ENDDO
      !!kill?: ENDDO

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
      ELSE IF (presentcase .EQ. ppm_poisson_grn_reprojec) THEN
        !remember to normalize the FFT
        normfac = 1.0_MK/ REAL((ppmpoisson%nmz(1)-1)* &
                              & (ppmpoisson%nmz(2)-1)* &
                              & (ppmpoisson%nmz(3)-1),MK)
        write(*,*) 'reprojection ppm style'
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
  
                !!write(*,*) ppmpoisson%fldzc2(:,i,j,k,isub),wdotk

                ppmpoisson%fldzc2(1,i,j,k,isub) = &
                  & (ppmpoisson%fldzc2(1,i,j,k,isub) - wdotk*kx)*normfac
                ppmpoisson%fldzc2(2,i,j,k,isub) = &
                  & (ppmpoisson%fldzc2(2,i,j,k,isub) - wdotk*ky)*normfac
                ppmpoisson%fldzc2(3,i,j,k,isub) = &
                  & (ppmpoisson%fldzc2(3,i,j,k,isub) - wdotk*kz)*normfac
                !!write(*,*) ppmpoisson%fldzc2(:,i,j,k,isub)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      !write(*,*) ppmpoisson%fldzc2
      !-----------------------------------------------------------------------
      ! IFFT pencil (Z)
      !-----------------------------------------------------------------------
      CALL ppm_fft_execute_1d(ppmpoisson%topoidz,&
      & ppmpoisson%meshidz, ppmpoisson%planbz, &
      & ppmpoisson%fldzc2, ppmpoisson%fldzc1, &
      & info)
      !!CALL ppm_fft_normalize(ppmpoisson%topoidz,&
      !!& ppmpoisson%meshidz, ppmpoisson%planbz,__ZEROSI,&
      !!& ppmpoisson%fldzc1, &
      !!& info)

      !write(*,*) ppmpoisson%fldzc1
      #ifdef __NOPE
      if (ppm_rank .EQ. trank)   THEN
            write(*,*)
            write(*,*)
      !ppmpoisson%fldzc1 = 0.0_MK
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
      !AOK

      !!IF (isubl .EQ. 1) THEN
        !!ppmpoisson%fldzc1 = CMPLX(-1.0_MK,-1.0_MK,MK)
      !!ELSE
        !!ppmpoisson%fldzc1 = CMPLX(8.0_MK,8.0_MK,MK)
        !!DO isub=1,ppmpoisson%nsublistz
          !!isubl=ppmpoisson%isublistz(isub)
            !!DO k=1,ppmpoisson%ndataz(3,isubl)
              !!DO j=1,ppmpoisson%ndataz(2,isubl)
                !!DO i=1,ppmpoisson%ndataz(1,isubl)
                  !!ppmpoisson%fldzc1(1,i,j,k,isub) = REAL(k,MK)
                !!END DO
              !!END DO
            !!END DO
        !!END DO
      !!ENDIF
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

        !!ppmpoisson%fldxyc = CMPLX(-1.0_MK,-1.0_MK,MK)
      !Retrieve
      CALL ppm_map_field_pop(&
      & ppmpoisson%topoidxy, &
      & ppmpoisson%meshidxyc,ppmpoisson%fldxyc, &
      & __NCOM,__ZEROSI,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',info2)
        GOTO 9999
      ENDIF

      !!write(*,*) 'ifftz slab'
      !!DO isub=1,ppmpoisson%nsublistxy
        !!isubl=ppmpoisson%isublistxy(isub)
        !!IF (isubl.EQ.2) THEN
        !!DO k=1,ppmpoisson%ndataxy(3,isubl)-1
          !!write(*,*) 'z',k
            !!DO i=1,ppmpoisson%ndataxy(1,isubl)-1
          !!DO j=1,ppmpoisson%ndataxy(2,isubl)-1
              !!write(*,'(A,E12.4,E12.4,A,$)') '(',ppmpoisson%fldxyc(1,i,j,k,isub),')'
            !!END DO
            !!write(*,*)
          !!END DO
        !!END DO
      !!ENDIF
      !!END DO
      !-----------------------------------------------------------------------
      ! IFFT (XY) use the non-reduced topology
      !-----------------------------------------------------------------------
      CALL ppm_fft_execute_2d(ppmpoisson%topoidxy,&
      & ppmpoisson%meshidxy, ppmpoisson%planbxy, &
      & ppmpoisson%fldxyc, ppmpoisson%fldxyr, &
      & info)

      !!CALL ppm_fft_normalize(ppmpoisson%topoidxy,&
      !!& ppmpoisson%meshidxy, ppmpoisson%planbxy, __ZEROSI,&
      !!& ppmpoisson%fldxyr, &
      !!& info)

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

      !Push the data
      CALL ppm_map_field_push(&
      & ppmpoisson%topoidxy, &
      & ppmpoisson%meshidxy,ppmpoisson%fldxyr,__NCOM,info)
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

      !-------------------------------------------------------------------------
      ! FINAL RETRIEVE - Here we do different things depending on the task
      ! i.e. the receiver varies
      !-------------------------------------------------------------------------
      IF (ppmpoisson%derivatives .EQ. ppm_poisson_drv_curl_fd2 .AND.&
          presentcase .EQ. ppm_poisson_grn_pois_per) THEN
        CALL ppm_map_field_pop(&
        & topoid, &
        & meshid,ppmpoisson%drv_vr, &
        & __NCOM,gstw,info) !!! gst
        !-------------------------------------------------------------------------
        ! Ghost the temporary array for derivatives (drv_vr)
        !-------------------------------------------------------------------------
        CALL ppm_map_field_ghost_get(topoid,meshid,gstw,info)
        CALL ppm_map_field_push(topoid,meshid,ppmpoisson%drv_vr,__NCOM,info)
        CALL ppm_map_field_send(info)
        CALL ppm_map_field_pop(topoid,meshid,ppmpoisson%drv_vr,__NCOM,gstw,info)

      ELSE
        CALL ppm_map_field_pop(&
        & topoid, &
        & meshid,fieldout, &
        & __NCOM,gstw,info) !!! gst
      ENDIF
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',info2)
        GOTO 9999
      ENDIF


      !-------------------------------------------------------------------------
      ! Optionally do derivatives
      ! Perhaps make ppm_poisson_fd take _none as argument. Then maybe no
      ! if-statement is required
      !-------------------------------------------------------------------------
      IF (ppmpoisson%derivatives .EQ. ppm_poisson_drv_curl_fd2 .AND.&
          presentcase  .EQ. ppm_poisson_grn_pois_per) THEN
        CALL ppm_poisson_fd(topoid,meshid,ppmpoisson%drv_vr,fieldout,&
                          & ppm_poisson_drv_curl_fd2,info)
      ENDIF

      !-------------------------------------------------------------------------
      ! Finally ghost the velocity/stream function field before returning it
      !-------------------------------------------------------------------------
      CALL ppm_map_field_ghost_get(topoid,meshid,gstw,info)
      CALL ppm_map_field_push(topoid,meshid,fieldout,__NCOM,info)
      CALL ppm_map_field_send(info)
      CALL ppm_map_field_pop(topoid,meshid,fieldout,__NCOM,gstw,info)


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

