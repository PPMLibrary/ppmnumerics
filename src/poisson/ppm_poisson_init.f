      !-------------------------------------------------------------------------
      !  Subroutine   : ppm_poisson_init.f90
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      !-------------------------------------------------------------------------
      SUBROUTINE __ROUTINE(topoid,meshid,ppmpoisson,fieldin,fieldout,green,info&
                                &,bc,derive)
      !!! Routine to initialise Greens function solution of the Poisson
      !!! equation. green is the flag defining which Greens function to use:
      !!! * ppm_poisson_grn_pois_per - Poisson equation, periodic boundaries
      !!! * ppm_poisson_grn_pois_fre - Poisson equation, freespace boundaries
      !!! * ppm_poisson_grn_reprojec - Do vorticity reprojection to kill divergence
      !!!
      !!! [NOTE]
      !!! fieldin is not preserved by this routine!
      !!! fieldin and fieldout must NOT be the same fields. In-place FFTs have
      !!! not been implemented.
      !!!
      !!! Eventually the routine should be overloaded to accept custom Greens
      !!! functions such that more general convolutions can be performed.
      !!! green should be expanded to include more buildin Greens functions.
      !!!
      !!! The routine should eventually accept an optional flag to toggle 
      !!! deallocation of work arrays between calls to ppm_poisson_solve.


      USE ppm_module_mktopo
      USE ppm_module_topo_get
      USE ppm_module_mesh_define
      USE ppm_module_map_field
      USE ppm_module_map_field_global
      USE ppm_module_map


      IMPLICIT NONE
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN)                                         :: topoid
      !!! Topology ID
      INTEGER, INTENT(IN)                                         :: meshid
      !!! Mesh ID
      TYPE(ppm_poisson_plan),INTENT(INOUT)                        :: ppmpoisson
      !!! The PPM Poisson plan type (inspired by the FFTW plan)
      REAL(__PREC),DIMENSION(:,:,:,:,:),POINTER                   :: fieldin
      !!! Input data field. RHS to the Poisson equation/field to be convolved
      !@ strictly speaking fieldin is not being used in the init routine
      REAL(__PREC),DIMENSION(:,:,:,:,:),POINTER                   :: fieldout
      !!! Output data field
      INTEGER, INTENT(IN)                                         :: green
      !!!flag to select build-in Greens functions:
      !!!ppm_poisson_grn_pois_per - Poisson equation, periodic boundaries
      !!!ppm_poisson_grn_pois_fre - Poisson equation, freespace boundaries
      !!!ppm_poisson_grn_reprojec - Do vorticity reprojection to kill divergence
      !!!
      !!!Eventually this should also accept custom Greens function
      INTEGER, INTENT(OUT)                                        :: info
      INTEGER,INTENT(IN),OPTIONAL                                 :: bc
      !!!boundary condition for the convolution. Can be on of the following:
      !!!ppm_poisson_grn_pois_per, ppm_poisson_grn_pois_fre.
      !!!One could argue that this is redundant in the build-in case
      !@ Isnt this redundant because of green?
      INTEGER,INTENT(IN),OPTIONAL                                 :: derive
      !!!flag to toggle various derivatives of the solution (not to be used with
      !!!green=ppm_poisson_grn_reprojec):
      !!! * ppm_poisson_drv_none
      !!! * ppm_poisson_drv_curl_sp  (only for periodic BC)
      !!! * ppm_poisson_drv_grad_sp  (not implemented)
      !!! * ppm_poisson_drv_lapl_sp  (not implemented)
      !!! * ppm_poisson_drv_div_sp   (not implemented)
      !!! * ppm_poisson_drv_curl_fd2
      !!! * ppm_poisson_drv_grad_fd2 (not implemented)
      !!! * ppm_poisson_drv_lapl_fd2 (not implemented)
      !!! * ppm_poisson_drv_div_fd2  (not implemented)
      !!! * ppm_poisson_drv_curl_fd4
      !!! * ppm_poisson_drv_grad_fd4 (not implemented)
      !!! * ppm_poisson_drv_lapl_fd4 (not implemented)
      !!! * ppm_poisson_drv_div_fd4  (not implemented)
      !!!
      !!! curl=curl, grad=gradient,lapl=laplace operator,div=divergence
      !!! sp=spectrally, fd2=2nd order finite differences, fd4=4th order FD
      !!!
      !!!The spectral derivatives can only be computed in this routine. Since 
      !!!the flag exists finite difference derivatives have also been included.


      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(__PREC)                            :: t0
      REAL(__PREC),DIMENSION(:,:),POINTER     :: xp=>NULL()      !particle positions
      TYPE(ppm_t_topo),POINTER                :: topology=>NULL()
      TYPE(ppm_t_equi_mesh)                   :: mesh
      INTEGER ,DIMENSION(__DIM)               :: indl,indu
      INTEGER,PARAMETER                       :: MK = __PREC
      REAL(__PREC),PARAMETER                  :: PI=ACOS(-1.0_MK) !@ use ppm pi
      REAL(__PREC)                            :: normfac
      !!!factor for the Greens function, including FFT normalization
      INTEGER                                 :: i,j,k
      INTEGER                                 :: kx,ky,kz
      INTEGER                                 :: isubl,isub
      INTEGER,DIMENSION(__DIM*2)              :: bcdef
      INTEGER                                 :: assigning
      INTEGER                                 :: decomposition
      INTEGER,SAVE                            :: ttopoid
      INTEGER                                 :: tmeshid
      REAL(__PREC)                            :: dx,dy,dz
      REAL(__PREC)                            :: Lx2,Ly2,Lz2
      REAL(__PREC)                            :: rad !@
      REAL(__PREC)                            :: gzero

      REAL(__PREC),DIMENSION(__DIM)           :: tmpmin,tmpmax
      INTEGER, DIMENSION(__DIM)               :: maxndataxy,maxndataz
      INTEGER, DIMENSION(:  ), POINTER        :: dummynmxy,dummynmz



      !-------------------------------------------------------------------------
      ! Initialise routine
      !-------------------------------------------------------------------------
      CALL substart('ppm_poisson_init',t0,info)


      !-------------------------------------------------------------------------
      ! Investigate optional arguments, setup routine accordingly
      ! !@TODO: Also check if the input/output and derivatives match
      !-------------------------------------------------------------------------
      IF (green .EQ. ppm_poisson_grn_pois_per) THEN
        ppmpoisson%case  = ppm_poisson_grn_pois_per
      ELSE IF (green .EQ. ppm_poisson_grn_pois_fre) THEN
        ppmpoisson%case  = ppm_poisson_grn_pois_fre
      ELSE IF (green .EQ. ppm_poisson_grn_reprojec) THEN
        ppmpoisson%case  = ppm_poisson_grn_reprojec
      ENDIF


      !-------------------------------------------------------------------------
      ! Nullify pointers from the ppmpoisson plans and the fftplans
      !-------------------------------------------------------------------------
      NULLIFY(xp)
      NULLIFY(ppmpoisson%costxy)
      NULLIFY(ppmpoisson%istartxy)
      NULLIFY(ppmpoisson%ndataxy)
      NULLIFY(ppmpoisson%istartxyc)
      NULLIFY(ppmpoisson%ndataxyc)
      NULLIFY(ppmpoisson%costz)
      NULLIFY(ppmpoisson%istartz)
      NULLIFY(ppmpoisson%ndataz)
      NULLIFY(ppmpoisson%planfxy%plan)
      NULLIFY(ppmpoisson%planbxy%plan)
      NULLIFY(ppmpoisson%planfz%plan)
      NULLIFY(ppmpoisson%planbz%plan)
      NULLIFY(maxndataxy)
      NULLIFY(maxndataz)
      NULLIFY(dummynmxy)
      NULLIFY(dummynmz)

      !-------------------------------------------------------------------------
      ! Get topology and mesh values of input/output
      !-------------------------------------------------------------------------
      CALL ppm_topo_get(topoid,topology,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_poisson_init','Failed to get topology.',isub)
        GOTO 9999
      ENDIF
      mesh  = topology%mesh(meshid)


      !-------------------------------------------------------------------------
      ! Setup mesh sizes for intermediate meshes/topologies
      !-------------------------------------------------------------------------
      IF (green .EQ. ppm_poisson_grn_pois_per) THEN
        !size of real slabs
        ppmpoisson%nmxy (1) =  mesh%nm(1)
        ppmpoisson%nmxy (2) =  mesh%nm(2)
        ppmpoisson%nmxy (3) =  mesh%nm(3)
        !size of complex slabs
        ppmpoisson%nmxyc(1) = (mesh%nm(1)-1)/2+1
        !!ppmpoisson%nmxyc(1) =  mesh%nm(1)
        ppmpoisson%nmxyc(2) =  mesh%nm(2)
        ppmpoisson%nmxyc(3) =  mesh%nm(3)
        !size of complex pencils
        ppmpoisson%nmz  (1) = (ppmpoisson%nmxyc(1))
        ppmpoisson%nmz  (2) = (ppmpoisson%nmxyc(2))
        ppmpoisson%nmz  (3) = (ppmpoisson%nmxyc(3))
        !size of the fft
        ppmpoisson%nmfft(1) =  mesh%nm(1)-1
        ppmpoisson%nmfft(2) =  mesh%nm(2)-1
        ppmpoisson%nmfft(3) =  mesh%nm(3)-1
        !Inverse of the size of the domain squared
        Lx2 = 1.0_MK/(topology%max_physd(1)-topology%min_physd(1))**2
        Ly2 = 1.0_MK/(topology%max_physd(2)-topology%min_physd(2))**2
        Lz2 = 1.0_MK/(topology%max_physd(3)-topology%min_physd(3))**2
      ELSE IF (green .EQ. ppm_poisson_grn_pois_fre) THEN !vertex
        !size of real slabs
        ppmpoisson%nmxy (1) =  mesh%nm(1)*2
        ppmpoisson%nmxy (2) =  mesh%nm(2)*2
        ppmpoisson%nmxy (3) =  mesh%nm(3)*2
        !size of complex slabs
        ppmpoisson%nmxyc(1) = (mesh%nm(1)*2)/2+1
        !!ppmpoisson%nmxyc(1) =  mesh%nm(1)*2
        ppmpoisson%nmxyc(2) =  mesh%nm(2)*2
        ppmpoisson%nmxyc(3) =  mesh%nm(3)*2
        !size of complex pencils
        ppmpoisson%nmz  (1) = (ppmpoisson%nmxyc(1))
        ppmpoisson%nmz  (2) = (ppmpoisson%nmxyc(2))
        ppmpoisson%nmz  (3) = (ppmpoisson%nmxyc(3))
        !size of the fft
        ppmpoisson%nmfft(1) =  mesh%nm(1)*2
        ppmpoisson%nmfft(2) =  mesh%nm(2)*2
        ppmpoisson%nmfft(3) =  mesh%nm(3)*2
        !Determine the grid spacing !vertex
        dx = (topology%max_physd(1)-topology%min_physd(1))/(mesh%nm(1)-1)
        dy = (topology%max_physd(2)-topology%min_physd(2))/(mesh%nm(2)-1)
        dz = (topology%max_physd(3)-topology%min_physd(3))/(mesh%nm(3)-1)
      ENDIF


      !-------------------------------------------------------------------------
      ! Create temporary derivation arrays if necessary
      !-------------------------------------------------------------------------
      !!write(*,*) 'fieldout',ppm_rank, LBOUND(fieldout),UBOUND(fieldout)
      IF (PRESENT(derive)) THEN
        IF ((      derive .EQ. ppm_poisson_drv_curl_fd2  &
          &  .OR.  derive .EQ. ppm_poisson_drv_curl_fd4)) THEN
          ppmpoisson%derivatives = derive
          ALLOCATE(ppmpoisson%drv_vr(LBOUND(fieldout,1):UBOUND(fieldout,1),&
                                   & LBOUND(fieldout,2):UBOUND(fieldout,2),&
                                   & LBOUND(fieldout,3):UBOUND(fieldout,3),&
                                   & LBOUND(fieldout,4):UBOUND(fieldout,4),&
                                   & LBOUND(fieldout,5):UBOUND(fieldout,5)))
        ELSE IF (derive .EQ. ppm_poisson_drv_curl_sp) THEN
          ppmpoisson%derivatives = ppm_poisson_drv_curl_sp
          !IF (ppmpoisson%case .EQ. ppm_poisson_grn_pois_fre) THEN
            !CALL ppm_write(ppm_rank,'ppm_poisson_init', &
            !& 'WARNING: Spectral curl is not fully implemented in Freespace.',isub)
          !ENDIF
        ELSE IF (derive .EQ. ppm_poisson_drv_none) THEN
          ppmpoisson%derivatives = ppm_poisson_drv_none
        ELSE
          CALL ppm_write(ppm_rank,'ppm_poisson_init','Undefined derivation input.',isub)
          GOTO 9999
        ENDIF
      ELSE
        ppmpoisson%derivatives = ppm_poisson_drv_none
      ENDIF


      !-------------------------------------------------------------------------
      ! Create spectral scaling components always. Just in case some 
      ! reprojection comes up
      ! The conditionals need to be for not just the Poisson equation
      !-------------------------------------------------------------------------
      IF (ppmpoisson%case .EQ. ppm_poisson_grn_pois_per) THEN
        ppmpoisson%normkx = &
        & 2.0_MK*PI/(topology%max_physd(1)-topology%min_physd(1))
        ppmpoisson%normky = &
        & 2.0_MK*PI/(topology%max_physd(2)-topology%min_physd(2))
        ppmpoisson%normkz = &
        & 2.0_MK*PI/(topology%max_physd(3)-topology%min_physd(3))
      ELSE IF (ppmpoisson%case .EQ. ppm_poisson_grn_pois_fre) THEN
        ppmpoisson%normkx = &
        & 2.0_MK*PI/((topology%max_physd(1)-topology%min_physd(1))*2.0_MK)
        ppmpoisson%normky = &
        & 2.0_MK*PI/((topology%max_physd(2)-topology%min_physd(2))*2.0_MK)
        ppmpoisson%normkz = &
        & 2.0_MK*PI/((topology%max_physd(3)-topology%min_physd(3))*2.0_MK)
      ENDIF


      !-------------------------------------------------------------------------
      ! Create new slab topology
      !-------------------------------------------------------------------------
      ttopoid = 0
      tmeshid = -1
      decomposition       = ppm_param_decomp_xy_slab
      assigning           = ppm_param_assign_internal
      IF (ppmpoisson%case .EQ. ppm_poisson_grn_pois_per) THEN
        bcdef               = ppm_param_bcdef_periodic
      ELSE IF (ppmpoisson%case .EQ. ppm_poisson_grn_pois_fre) THEN
        bcdef               = ppm_param_bcdef_freespace
      ENDIF
      tmpmin              = topology%min_physd
      tmpmax              = topology%max_physd


      CALL ppm_mktopo(ttopoid,tmeshid,xp,0,&
      & decomposition,assigning,&
      & tmpmin,tmpmax,bcdef,&
      & __ZEROSI,ppmpoisson%costxy,&
      & ppmpoisson%nmxy,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_poisson_init','Failed to create xy-topology.',isub)
        GOTO 9999
      ENDIF
      ppmpoisson%topoidxy = ttopoid
      ppmpoisson%meshidxy = tmeshid
      !-------------------------------------------------------------------------
      ! Get additional xy-mesh information
      !-------------------------------------------------------------------------
      CALL ppm_topo_get_meshinfo(ppmpoisson%topoidxy,ppmpoisson%meshidxy, &
         & dummynmxy,ppmpoisson%istartxy,ppmpoisson%ndataxy,maxndataxy, &
         & ppmpoisson%isublistxy,ppmpoisson%nsublistxy,info)


      !-------------------------------------------------------------------------
      ! Create complex slab mesh
      !-------------------------------------------------------------------------
      ttopoid = ppmpoisson%topoidxy
      tmeshid = -1
      CALL ppm_mesh_define(ttopoid,tmeshid,&
      & ppmpoisson%nmxyc,ppmpoisson%istartxyc,ppmpoisson%ndataxyc,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_poisson_init','Failed to create complex xy mesh definition.',isub)
        GOTO 9999
      ENDIF
      ppmpoisson%meshidxyc = tmeshid


      !-------------------------------------------------------------------------
      ! Create new pencil topology
      !-------------------------------------------------------------------------
      ttopoid = 0
      tmeshid = -1
      IF (ppmpoisson%case .EQ. ppm_poisson_grn_pois_per) THEN
        bcdef               = ppm_param_bcdef_periodic
      ELSE IF (ppmpoisson%case .EQ. ppm_poisson_grn_pois_fre) THEN
        bcdef               = ppm_param_bcdef_freespace
      ENDIF
      assigning       = ppm_param_assign_internal
      decomposition   = ppm_param_decomp_zpencil

      CALL ppm_mktopo(ttopoid,tmeshid,xp,0,&
      & decomposition,assigning,&
      & tmpmin,tmpmax,bcdef,&
      & __ZEROSI,ppmpoisson%costz,&
      & ppmpoisson%nmz,info)
      IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_poisson_init','Failed to create z-topology.',isub)
        GOTO 9999
      ENDIF
      ppmpoisson%topoidz = ttopoid
      ppmpoisson%meshidz = tmeshid
      !-------------------------------------------------------------------------
      ! Get additional z-mesh information
      !-------------------------------------------------------------------------
      CALL ppm_topo_get_meshinfo(ppmpoisson%topoidz,ppmpoisson%meshidz, &
         & dummynmz,ppmpoisson%istartz,ppmpoisson%ndataz,maxndataz, &
         & ppmpoisson%isublistz,ppmpoisson%nsublistz,info)


      !-------------------------------------------------------------------------
      ! Set and get minimum and maximum indicies
      !-------------------------------------------------------------------------
      indl(1) = 1
      indl(2) = 1
      indl(3) = 1

      !-------------------------------------------------------------------------
      ! Allocate real xy slabs
      !-------------------------------------------------------------------------
      ALLOCATE(ppmpoisson%fldxyr(__DIM,&
      & indl(1):maxndataxy(1),indl(2):maxndataxy(2),indl(3):maxndataxy(3),&
      & 1:ppmpoisson%nsublistxy),stat=info)


      !-------------------------------------------------------------------------
      ! Set and get minimum and maximum indicies of COMPLEX xy slabs
      !-------------------------------------------------------------------------
      indl(1) = 1
      indl(2) = 1
      indl(3) = 1
      indu(1) = 0
      indu(2) = 0
      indu(3) = 0
      DO isub=1,ppmpoisson%nsublistxy
        isubl = ppmpoisson%isublistxy(isub)
        indu(1) = MAX(indu(1),ppmpoisson%ndataxyc(1,isubl))
        indu(2) = MAX(indu(2),ppmpoisson%ndataxyc(2,isubl))
        indu(3) = MAX(indu(3),ppmpoisson%ndataxyc(3,isubl))
      ENDDO


      !-------------------------------------------------------------------------
      ! Allocate complex xy slabs
      !-------------------------------------------------------------------------
      ALLOCATE(ppmpoisson%fldxyc(__DIM,&
      & indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
      & 1:ppmpoisson%nsublistxy),stat=info)


      !-------------------------------------------------------------------------
      ! Allocate two complex z pencils + Greens fcn array !@check return vars.
      !-------------------------------------------------------------------------
      ALLOCATE(ppmpoisson%fldzc1(__DIM,&
      & indl(1):maxndataz(1),indl(2):maxndataz(2),indl(3):maxndataz(3),&
      & 1:ppmpoisson%nsublistz),stat=info)

      ALLOCATE(ppmpoisson%fldzc2(__DIM,&
      & indl(1):maxndataz(1),indl(2):maxndataz(2),indl(3):maxndataz(3),&
      & 1:ppmpoisson%nsublistz),stat=info)


      !-------------------------------------------------------------------------
      ! The complex Greens function is always kept on the z-pencil topology
      !-------------------------------------------------------------------------
      IF (green .EQ. ppm_poisson_grn_pois_per) THEN
        !!ALLOCATE(ppmpoisson%fldgrnr(&
        !!& indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
        !!& 1:ppmpoisson%nsublistz),stat=info)
        ALLOCATE(ppmpoisson%fldgrnr(&
        & indl(1):maxndataz(1),indl(2):maxndataz(2),indl(3):maxndataz(3),&
        & 1:ppmpoisson%nsublistz),stat=info)
      ELSE IF (green .EQ. ppm_poisson_grn_pois_fre) THEN
        !!ALLOCATE(ppmpoisson%fldgrnc(&
        !!& indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
        !!& 1:ppmpoisson%nsublistz),stat=info)
        ALLOCATE(ppmpoisson%fldgrnc(&
        & indl(1):maxndataz(1),indl(2):maxndataz(2),indl(3):maxndataz(3),&
        & 1:ppmpoisson%nsublistz),stat=info)
      ENDIF


      !-------------------------------------------------------------------------
      ! Set up xy FFT plans
      ! The inverse plan takes the returning topology since it has the full size
      !-------------------------------------------------------------------------
      CALL ppm_fft_forward_2d(ppmpoisson%topoidxy,ppmpoisson%meshidxy,&
      & ppmpoisson%planfxy,ppmpoisson%fldxyr,&
      & ppmpoisson%fldxyc,info)

      CALL ppm_fft_backward_2d(ppmpoisson%topoidxy,ppmpoisson%meshidxy,&
      & ppmpoisson%planbxy,ppmpoisson%fldxyc,&
      & ppmpoisson%fldxyr,info)


      !-------------------------------------------------------------------------
      ! Set up z FFT plans
      !-------------------------------------------------------------------------
      CALL ppm_fft_forward_1d(ppmpoisson%topoidz,ppmpoisson%meshidz,&
      & ppmpoisson%planfz,ppmpoisson%fldzc1,&
      & ppmpoisson%fldzc2,info)

      CALL ppm_fft_backward_1d(ppmpoisson%topoidz,ppmpoisson%meshidz,&
      & ppmpoisson%planbz,ppmpoisson%fldzc2,&
      & ppmpoisson%fldzc1,info)


      !-------------------------------------------------------------------------
      ! Compute Greens function. Analytic, periodic
      !
      ! (d2_/dx2 + d2_/dy2 + d2_/dz2)psi = -omega     =>
      ! -4*pi2(kx2 + ky2 + kz2)PSI       = -OMEGA     =>
      ! PSI                              = 1/(4*pi2)*1/(kx2 + ky2 + kz2)OMEGA
      !-------------------------------------------------------------------------
      IF (green .EQ. ppm_poisson_grn_pois_per) THEN
        ! Scaling the spectral coefficients... 
        ! one minus due to (i*k)^2 and another due to the Poisson equation
        normfac = 1.0_MK/(4.0_MK*PI*PI * &
                !and normalisation of FFTs (full domain) !vertex
                & REAL((ppmpoisson%nmfft(1))* &
                &      (ppmpoisson%nmfft(2))* &
                &      (ppmpoisson%nmfft(3)),MK))
        DO isub=1,ppmpoisson%nsublistz
          isubl=ppmpoisson%isublistz(isub)
          DO k=1,ppmpoisson%ndataz(3,isubl)
            DO j=1,ppmpoisson%ndataz(2,isubl)
              DO i=1,ppmpoisson%ndataz(1,isubl)
                kx = i-1 + (ppmpoisson%istartz(1,isubl)-1)
                ky = j-1 + (ppmpoisson%istartz(2,isubl)-1)
                kz = k-1 + (ppmpoisson%istartz(3,isubl)-1)
                !This is a nasty way to do this but it is only done once so...:
                IF (kx .GT. (ppmpoisson%nmfft(1)/2)) kx = kx-(ppmpoisson%nmfft(1))
                IF (ky .GT. (ppmpoisson%nmfft(2)/2)) ky = ky-(ppmpoisson%nmfft(2))
                IF (kz .GT. (ppmpoisson%nmfft(3)/2)) kz = kz-(ppmpoisson%nmfft(3))
                ppmpoisson%fldgrnr(i,j,k,isub) = &
                      & normfac/(REAL(kx*kx,__PREC)*Lx2 &
                      &        + REAL(ky*ky,__PREC)*Ly2 &
                      &        + REAL(kz*kz,__PREC)*Lz2)
                !Take care of singularity
                !This is nasty as well
                IF ((kx*kx+ky*ky+kz*kz) .EQ. 0) THEN
                  ppmpoisson%fldgrnr(i,j,k,isub) = 0.0_MK
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO


      !-------------------------------------------------------------------------
      ! Compute real free space Greens function. Analytic
      ! The Greens function is initialised temporarily in the real xy slabs,
      ! then FFTed in those directions, mapped to z-pencils, FFTed in z and
      ! finally copied to ppmpoisson%fldgrnc. The real xy slabs have already
      ! been setup for FFTs etc so they offer a convenient container for the
      ! FFTing the Greens function instead of setting up the whole apparatus
      ! for this one-time affair.
      ! These loops must run over the padded(extended) domain thus %ndataxy
      ! \nabla \Psi = -\omega
      ! The Greens function is 1/4piR (R is distance to origo) and includes the
      ! minus of the RHS
      !-------------------------------------------------------------------------
      ELSE IF (green .EQ. ppm_poisson_grn_pois_fre) THEN
#define __GFREESPACE
!#define __GBLOB4

!#define __GSPEC
        !-----------------------------------------------------------------------
        ! First initialise the real Greens function
        !@alternatively this could come from as input
        !-----------------------------------------------------------------------
        !there should NOT be a minus here since THIS Greens function takes
        !the minus of the Poisson equation into account
        normfac = 1.0_MK/(4.0_MK*PI* &
        !remembering FFT normalization of ALL points: !vertex
                & REAL((ppmpoisson%nmfft(1))* &
                &      (ppmpoisson%nmfft(2))* &
                &      (ppmpoisson%nmfft(3)),MK))*dx*dy*dz
        DO isub=1,ppmpoisson%nsublistxy
          isubl=ppmpoisson%isublistxy(isub)
          DO k=1,ppmpoisson%ndataxy(3,isubl)
            DO j=1,ppmpoisson%ndataxy(2,isubl)
              DO i=1,ppmpoisson%ndataxy(1,isubl)
                !kx,ky,kz implies that they are spectral coordinates. they are not
                kx = i-1 + (ppmpoisson%istartxy(1,isubl)-1)
                ky = j-1 + (ppmpoisson%istartxy(2,isubl)-1)
                kz = k-1 + (ppmpoisson%istartxy(3,isubl)-1)
                !This is a nasty way to do this but it is only done once so...:
                IF (kx .GT. (ppmpoisson%nmfft(1))/2) kx = kx-(ppmpoisson%nmfft(1))
                IF (ky .GT. (ppmpoisson%nmfft(2))/2) ky = ky-(ppmpoisson%nmfft(2))
                IF (kz .GT. (ppmpoisson%nmfft(3))/2) kz = kz-(ppmpoisson%nmfft(3))
                rad = SQRT( REAL(kx*kx,__PREC)*dx*dx+ &
                          & REAL(ky*ky,__PREC)*dy*dy+ &
                          & REAL(kz*kz,__PREC)*dz*dz)
                !----FREESPACE GREENS FUNCTION----
#ifdef __GFREESPACE
                ppmpoisson%fldxyr(1,i,j,k,isub) = &
                        & normfac/(SQRT( REAL(kx*kx,__PREC)*dx*dx+ &
                        &                REAL(ky*ky,__PREC)*dy*dy+ &
                        &                REAL(kz*kz,__PREC)*dz*dz))
#endif
                !----BLOB----2nd ORDER---:
#ifdef __GBLOB2
                ppmpoisson%fldxyr(1,i,j,k,isub) = &
                        & normfac/(rad)*ERF(rad/(SQRT(2.0E0_MK)*dx))
                        !& normfac/(rad)*ERF(rad/(1.0_MK*dx))
                        !& normfac/(rad)*ERF(rad/1.0E-2_MK)
#endif
                !----BLOB----4th ORDER---cell-normalised:
#ifdef __GBLOB4
                ppmpoisson%fldxyr(1,i,j,k,isub) = &
                        & normfac/(rad)* &
                        & (ERF(rad/(dx*1.0_MK*SQRT(2.0_MK)))+ &
                        & 1.0_MK/SQRT(PI)*rad/(dx*1.0_MK*SQRT(2.0_MK))* &
                        & EXP(-rad**2/(dx*1.0_MK*SQRT(2.0_MK))**2))

                        !---4th ORDER---unnormalised:
                        !& (ERF(rad/SQRT(2.0_MK))+1.0_MK/SQRT(2.0_MK*PI)*rad* &
                        !& EXP(-rad**2/2))
#endif
                ppmpoisson%fldxyr(2,i,j,k,isub) = 0.0_MK
                ppmpoisson%fldxyr(3,i,j,k,isub) = 0.0_MK
                !Take care of singularity (This is nasty as well)
                !Simply zero
                !gzero = 0.0_MK
                !Simply one (H&E style)
                !gzero = 1.0_MK*normfac
                !Equal to the neighbour point (first order)
                !gzero = normfac/(dx)
                !Linear extrapolation into zero (2nd order)
                !gzero =   2.0_MK*normfac/(       dx) &
                        !& - 1.0_MK*normfac/(2.0_MK*dx)
                !extrapolation 3rd order
                !gzero =   3.0_MK*normfac/(       dx) &
                        !& - 3.0_MK*normfac/(2.0_MK*dx) &
                        !& + 1.0_MK*normfac/(3.0_MK*dx)
                !extrapolation 4th order
                !gzero =   4.0_MK*normfac/(       dx) &
                        !& - 6.0_MK*normfac/(2.0_MK*dx) &
                        !& + 4.0_MK*normfac/(3.0_MK*dx) &
                        !& - 1.0_MK*normfac/(4.0_MK*dx)
                !extrapolation 6th order
                !gzero =   6.0_MK*normfac/(       dx) &
                        !& -15.0_MK*normfac/(2.0_MK*dx) &
                        !& +20.0_MK*normfac/(3.0_MK*dx) &
                        !& -15.0_MK*normfac/(4.0_MK*dx) &
                        !& + 6.0_MK*normfac/(5.0_MK*dx) &
                        !& - 1.0_MK*normfac/(6.0_MK*dx)
                !one thousand!!!
                !gzero = normfac*1000.0_MK
                !Chatelain:2010
                gzero = &
                      & 0.5_MK*( 3.0_MK*dx*dy*dz/(4.0_MK*PI) ) ** (2.0_MK/3.0_MK) * &
                      & normfac *4.0_MK*PI
                !bi-section galore!!!!
                !gzero = normfac * &
                !Chatelain for dx=2/64:
                  !& 0.00236131633889131295_MK
                  !& 0.001_MK
                  !& 0.01_MK
                  !& 0.1_MK
                  !& 1.0_MK
                  !& 10.0_MK
                  !& 100.0_MK
                  !& 1.0E-4_MK
                  !& 1.0E-5_MK
                  !& 0.0_MK
                  !& -0.001_MK
                  !& -0.1_MK
                  !& -10.0_MK
                  !& -1000.0_MK
                  !& -100.0_MK
                  !& -20.0_MK
                  !& -30.0_MK
                !Matching the blob (both 2nd and 4th order)
#if defined __GBLOB2 || defined __GBLOB4
                gzero = 4.0_MK*PI*normfac/ &
                        & (2.0_MK*PI**(3.0_MK/2.0_MK)* &
                        !& SQRT(2.0_MK))
                        & 1.0_MK*SQRT(2.0_MK)*dx)
                        !& 0.01_MK)
                !!write(*,*) ERF(1.0_MK)
#endif
                IF ((kx*kx+ky*ky+kz*kz) .EQ. 0) THEN
                  !Professor Nutbutter is out.
                  ppmpoisson%fldxyr(1,i,j,k,isub) = gzero
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        !-------------------------------------------------------------------------
        ! FOURIER TRANSFORM AND MAP GALORE
        ! This part should be used both for freespace and a custom Greens function
        !-------------------------------------------------------------------------
        !-----------------------------------------------------------------------
        ! Do slab FFT (XY)
        !-----------------------------------------------------------------------
        CALL ppm_fft_execute_2d(ppmpoisson%topoidxy,&
        & ppmpoisson%meshidxy, ppmpoisson%planfxy, &
        & ppmpoisson%fldxyr, ppmpoisson%fldxyc, &
        & info)
        !-----------------------------------------------------------------------
        ! Map to the pencils (Z)
        !-----------------------------------------------------------------------
        !Initialise
        CALL ppm_map_field_global(&
        & ppmpoisson%topoidxy, &
        & ppmpoisson%topoidz, &
        & ppmpoisson%meshidxyc, &
        & ppmpoisson%meshidz,info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to initialise field mapping.',isub)
          GOTO 9999
        ENDIF

        !Push the data
        CALL ppm_map_field_push(&
        & ppmpoisson%topoidxy, &
        & ppmpoisson%meshidxyc,ppmpoisson%fldxyc,__NCOM,info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to push vector field.',isub)
          GOTO 9999
        ENDIF

        !Send
        CALL ppm_map_field_send(info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to send field.',isub)
          GOTO 9999
        ENDIF

        !Retrieve
        CALL ppm_map_field_pop(&
        & ppmpoisson%topoidz, &
        & ppmpoisson%meshidz,ppmpoisson%fldzc1, &
        & __NCOM,__ZEROSI,info)
        IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank, 'ppm_poisson_solve','Failed to pop vector field.',isub)
          GOTO 9999
        ENDIF


        !-----------------------------------------------------------------------
        ! Do pencil FFT (Z)
        !-----------------------------------------------------------------------
        CALL ppm_fft_execute_1d(ppmpoisson%topoidz,&
        & ppmpoisson%meshidz, ppmpoisson%planfz, &
        & ppmpoisson%fldzc1, ppmpoisson%fldzc2, &
        & info)


        !-----------------------------------------------------------------------
        ! Copy first component of the Fourier transformed vector to fldgrnc
        !-----------------------------------------------------------------------
        DO isub=1,ppmpoisson%nsublistz
          isubl=ppmpoisson%isublistz(isub)
          DO k=1,ppmpoisson%ndataz(3,isubl)
            DO j=1,ppmpoisson%ndataz(2,isubl)
              DO i=1,ppmpoisson%ndataz(1,isubl)
#ifdef __GSPEC
                ppmpoisson%fldgrnc(i,j,k,isub) = ppmpoisson%fldzc2(1,i,j,k,isub) &
                & + (1.0_MK * normfac*4.0_MK*PI - gzero) !dirac minus actual value
#else
                ppmpoisson%fldgrnc(i,j,k,isub) = ppmpoisson%fldzc2(1,i,j,k,isub)
#endif
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      END IF
      !-------------------------------------------------------------------------
      ! Or alternatively FFT real Green from input
      !-------------------------------------------------------------------------


      !-------------------------------------------------------------------------
      ! Deallocate fields? !@
      !-------------------------------------------------------------------------


      !-------------------------------------------------------------------------
      ! Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_poisson_init',t0,info)
      RETURN

      END SUBROUTINE __ROUTINE

