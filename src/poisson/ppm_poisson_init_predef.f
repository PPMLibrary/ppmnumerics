  !-------------------------------------------------------------------------
  ! ppm_poisson_init_predef.f90
  !-------------------------------------------------------------------------
  ! Routine to initialise the convolution using build in Greens function
  !-------------------------------------------------------------------------
 SUBROUTINE __ROUTINE(topoid,meshid,ppmpoisson,fieldin,fieldout,green,info&
       &,bc,derive)

    USE ppm_module_mktopo
    USE ppm_module_topo_get
    USE ppm_module_mesh_define

    IMPLICIT NONE
    !-------------------------------------------------------------------------
    ! Arguments
    !-------------------------------------------------------------------------
    INTEGER, INTENT(IN)                                         :: topoid
    INTEGER, INTENT(IN)                                         :: meshid
    TYPE(ppm_poisson_plan),INTENT(INOUT)                        :: ppmpoisson
    !@ strictly speaking fieldin is not being used in the init routine
    REAL(__PREC),DIMENSION(:,:,:,:,:),POINTER                   :: fieldin
    REAL(__PREC),DIMENSION(:,:,:,:,:),POINTER                   :: fieldout
    !REAL(__PREC),DIMENSION(:,:,:,:,:),POINTER                   :: green
    INTEGER, INTENT(IN)                                         :: green
!!!flag to select build-in Greens functions:
!!!ppm_poisson_grn_pois_per - Poisson equation, periodic boundaries
!!!ppm_poisson_grn_pois_fre - Poisson equation, freespace boundaries (not implemented)
!!!ppm_poisson_grn_reprojec - Do vorticity reprojection to kill divergence
    INTEGER, INTENT(OUT)                                        :: info
    INTEGER,INTENT(IN),OPTIONAL                                 :: bc
!!!boundary condition for the convolution. Can be on of the following:
!!!ppm_poisson_grn_pois_per, ppm_poisson_grn_pois_fre.
!!!One could argue that this is redundant in the build-in case
!!!@ Isnt this redundant because of green?
    INTEGER,INTENT(IN),OPTIONAL                                 :: derive
!!!flag to toggle various derivatives of the solution (not to be used with
!!!green=ppm_poisson_grn_reprojec):
!!!ppm_poisson_drv_none
!!!ppm_poisson_drv_curl_fd (not implemented)
!!!ppm_poisson_drv_grad_fd (not implemented)
!!!ppm_poisson_drv_lapl_fd (not implemented)
!!!ppm_poisson_drv_curl_sp (not implemented)
!!!ppm_poisson_drv_grad_sp (not implemented)
!!!ppm_poisson_drv_lapl_sp (not implemented)
!!!The spectral derivatives can only be computed in this routine. Since 
!!!the flag exists finite difference derivatives have also been included.

    !-------------------------------------------------------------------------
    ! Local variables
    !-------------------------------------------------------------------------
    REAL(__PREC)                            :: t0
    REAL(__PREC),DIMENSION(:,:),POINTER     :: xp      !particle positions
    TYPE(ppm_t_topo),POINTER                :: topology,topologyxy,topologyxyc,topologyz
    TYPE(ppm_t_equi_mesh)                   :: mesh!!,meshtmp1,meshtmp2 !@
    INTEGER ,DIMENSION(__DIM)               :: indl,indu
    INTEGER,PARAMETER                       :: MK = __PREC
    REAL(__PREC),PARAMETER                  :: PI=ACOS(-1.0_MK) !@ use ppm pi
    !factor for the Green's function, including FFT normalization
    REAL(__PREC)                            :: normfac
    INTEGER                                 :: i,j,k
    INTEGER                                 :: kx,ky,kz
    INTEGER                                 :: isubl,isub
    INTEGER,DIMENSION(__DIM*2)              :: bcdef
    INTEGER                                 :: assigning
    INTEGER                                 :: decomposition
    INTEGER,SAVE                            :: ttopoid
    INTEGER                                 :: tmeshid

    REAL(__PREC),DIMENSION(__DIM)           :: tmpmin,tmpmax

    !-------------------------------------------------------------------------
    ! Initialise routine
    !-------------------------------------------------------------------------
    CALL substart('ppm_poisson_init_predef',t0,info)

    !-------------------------------------------------------------------------
    ! Investigate optional arguments, setup routine accordingly
    ! !@ Also check if the input/output and derivatives match
    !-------------------------------------------------------------------------
    IF (green .EQ. ppm_poisson_grn_pois_per) THEN
       ppmpoisson%case  = ppm_poisson_grn_pois_per
    ELSE IF (green .EQ. ppm_poisson_grn_pois_fre) THEN
       ppmpoisson%case  = ppm_poisson_grn_pois_fre
    ELSE IF (green .EQ. ppm_poisson_grn_reprojec) THEN
       ppmpoisson%case  = ppm_poisson_grn_reprojec
    ENDIF

    IF (PRESENT(derive)) THEN
       !Create temporary arrays if necessary
       IF (derive .EQ. ppm_poisson_drv_curl_fd2 .AND.&
            & ppmpoisson%case  .EQ. ppm_poisson_grn_pois_per) THEN
          ppmpoisson%derivatives = derive
          ALLOCATE(ppmpoisson%drv_vr(LBOUND(fieldout,1):UBOUND(fieldout,1),&
               & LBOUND(fieldout,2):UBOUND(fieldout,2),&
               & LBOUND(fieldout,3):UBOUND(fieldout,3),&
               & LBOUND(fieldout,4):UBOUND(fieldout,4),&
               & LBOUND(fieldout,5):UBOUND(fieldout,5)))
       ELSE IF (derive .EQ. ppm_poisson_drv_none) THEN
          ppmpoisson%derivatives = ppm_poisson_drv_none
       ELSE
          CALL ppm_write(ppm_rank,'ppm_poisson_init_predef','Undefined derivation input.',isub)
          GOTO 9999
       ENDIF
    ELSE
       ppmpoisson%derivatives = ppm_poisson_drv_none
    ENDIF

    !-------------------------------------------------------------------------
    ! Get topology and mesh values
    !-------------------------------------------------------------------------
    CALL ppm_topo_get(topoid,topology,info)
    IF (info .NE. 0) THEN
       CALL ppm_write(ppm_rank,'ppm_poisson_init_predef','Failed to get topology.',isub)
       GOTO 9999
    ENDIF
    !nsubs = topology%nsublist
    mesh  = topology%mesh(meshid)
    !Copy size of global mesh
    ppmpoisson%nmxy (1) =  mesh%nm(1)
    ppmpoisson%nmxy (2) =  mesh%nm(2)
    ppmpoisson%nmxy (3) =  mesh%nm(3)
    !ppmpoisson%nmxyc(1) = (mesh%nm(1)-1)/2+1
    ppmpoisson%nmxyc(1) = (mesh%nm(1)) !tmp until meshid>1 works
    ppmpoisson%nmxyc(2) =  mesh%nm(2)
    ppmpoisson%nmxyc(3) =  mesh%nm(3)
    !ppmpoisson%nmz  (1) = (mesh%nm(1)-1)/2+1!/2+1 !trying to reduce the number of points in matrix+1
    ppmpoisson%nmz  (1) = (mesh%nm(1)) !tmp until meshid>1 works
    ppmpoisson%nmz  (2) =  mesh%nm(2)
    ppmpoisson%nmz  (3) =  mesh%nm(3)


    !-------------------------------------------------------------------------
    ! Create new slab topology !@sofar periodic
    !-------------------------------------------------------------------------
    ttopoid = 0
    tmeshid = -1
    decomposition       = ppm_param_decomp_xy_slab
    assigning           = ppm_param_assign_internal
    bcdef               = ppm_param_bcdef_periodic
    tmpmin              = topology%min_physd
    tmpmax              = topology%max_physd
    !!tmpmin2             = topology%min_physd
    !!tmpmax2             = topology%max_physd

    CALL ppm_mktopo(ttopoid,tmeshid,xp,0,&
         & decomposition,assigning,&
                                !& topology%min_physd,topology%max_physd,bcdef,&
         & tmpmin,tmpmax,bcdef,&
         & __ZEROSI,ppmpoisson%costxy,ppmpoisson%istartxy,ppmpoisson%ndataxy,&
         & ppmpoisson%nmxy,info)
    IF (info .NE. 0) THEN
       CALL ppm_write(ppm_rank,'ppm_poisson_init_predef','Failed to create xy-topology.',isub)
       GOTO 9999
    ENDIF
    ppmpoisson%topoidxy = ttopoid
    ppmpoisson%meshidxy = tmeshid

    !-------------------------------------------------------------------------
    ! Create complex slab mesh
    ! !@ not used until meshid>1 works
    !-------------------------------------------------------------------------
    ttopoid = ppmpoisson%topoidxy
    tmeshid = -1
    CALL ppm_mesh_define(ttopoid,tmeshid,&
         & ppmpoisson%nmxy,ppmpoisson%istartxyc,ppmpoisson%ndataxyc,info) !@nmxyC
    IF (info .NE. 0) THEN
       CALL ppm_write(ppm_rank,'ppm_poisson_init_predef','Failed to create complex xy mesh definition.',isub)
       GOTO 9999
    ENDIF
    ppmpoisson%meshidxyc = tmeshid
    !write(*,*) ppmpoisson%meshidxy,tmeshid,ttopoid,ppmpoisson%meshidxyc
    !!write(*,*) ppmpoisson%istartxyc,ppmpoisson%ndataxyc


    !-------------------------------------------------------------------------
    ! Create new pencil topology !@sofar periodic
    !-------------------------------------------------------------------------
    ttopoid = 0
    tmeshid = -1
    bcdef           = ppm_param_bcdef_periodic
    assigning       = ppm_param_assign_internal
    decomposition   = ppm_param_decomp_zpencil
    !tmpmin              = topology%min_physd
    !tmpmax              = topology%max_physd

    CALL ppm_mktopo(ttopoid,tmeshid,xp,0,&
         & decomposition,assigning,&
                                !& topology%min_physd,topology%max_physd,bcdef,&
         & tmpmin,tmpmax,bcdef,&
         & __ZEROSI,ppmpoisson%costz,ppmpoisson%istartz,ppmpoisson%ndataz,&
         & ppmpoisson%nmz,info)
    IF (info .NE. 0) THEN
       CALL ppm_write(ppm_rank,'ppm_poisson_init_predef','Failed to create z-topology.',isub)
       GOTO 9999
    ENDIF
    ppmpoisson%topoidz = ttopoid
    ppmpoisson%meshidz = tmeshid

    !-------------------------------------------------------------------------
    ! Get the new topologies (slabs and pencils)
    !-------------------------------------------------------------------------
    CALL ppm_topo_get(ppmpoisson%topoidxy,topologyxy,info)
    IF (info .NE. 0) THEN
       CALL ppm_write(ppm_rank,'ppm_poisson_init_predef','Failed to get xy topology.',isub)
       GOTO 9999
    ENDIF

    !!CALL ppm_topo_get(ppmpoisson%topoidxyc,topologyxyc,info)
    !!IF (info .NE. 0) THEN
    !!CALL ppm_write(ppm_rank,'ppm_poisson_init_predef','Failed to get complex xy topology.',isub)
    !!GOTO 9999
    !!ENDIF

    CALL ppm_topo_get(ppmpoisson%topoidz,topologyz,info)
    IF (info .NE. 0) THEN
       CALL ppm_write(ppm_rank,'ppm_poisson_init_predef','Failed to get z topology.',isub)
       GOTO 9999
    ENDIF


    !-------------------------------------------------------------------------
    ! Copy data to ppm_poisson type
    !-------------------------------------------------------------------------
    ppmpoisson%nsublistxy  = topologyxy %nsublist
    !!ppmpoisson%nsublistxyc = topologyxyc%nsublist
    ppmpoisson%nsublistz   = topologyz  %nsublist

    ALLOCATE(ppmpoisson%isublistxy (ppmpoisson%nsublistxy))
    !!ALLOCATE(ppmpoisson%isublistxyc(ppmpoisson%nsublistxyc))
    ALLOCATE(ppmpoisson%isublistz  (ppmpoisson%nsublistz))

    DO isub=1,ppmpoisson%nsublistxy
       ppmpoisson%isublistxy(isub)  = topologyxy%isublist(isub)
    ENDDO
    !!DO isub=1,ppmpoisson%nsublistxyc
    !!ppmpoisson%isublistxyc(isub) = topologyxyc%isublist(isub)
    !!ENDDO
    DO isub=1,ppmpoisson%nsublistz
       ppmpoisson%isublistz(isub)   = topologyz%isublist(isub)
    ENDDO

    !-------------------------------------------------------------------------
    ! Set and get minimum and maximum indicies of xy slabs
    !-------------------------------------------------------------------------
    indl(1) = 1
    indl(2) = 1
    indl(3) = 1
    indu(1) = 0
    indu(2) = 0
    indu(3) = 0
    DO isub=1,ppmpoisson%nsublistxy
       isubl = ppmpoisson%isublistxy(isub)
       indu(1) = MAX(indu(1),ppmpoisson%ndataxy(1,isubl))
       indu(2) = MAX(indu(2),ppmpoisson%ndataxy(2,isubl))
       indu(3) = MAX(indu(3),ppmpoisson%ndataxy(3,isubl))
    ENDDO

    !-------------------------------------------------------------------------
    ! Allocate real xy slabs
    !-------------------------------------------------------------------------
    ALLOCATE(ppmpoisson%fldxyr(__DIM,&
         & indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
         & 1:ppmpoisson%nsublistxy),stat=info)

    !!ALLOCATE(ppmpoisson%fldxyc(__DIM,&
    !!& indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
    !!& 1:ppmpoisson%nsublistxy),stat=info)

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
    ! Allocate real and complex xy slabs
    !-------------------------------------------------------------------------
    !!ALLOCATE(ppmpoisson%fldxyr(__DIM,&
    !!& indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
    !!& 1:ppmpoisson%nsublistxy),stat=info)
    !!
    ALLOCATE(ppmpoisson%fldxyc(__DIM,&
         & indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
         & 1:ppmpoisson%nsublistxy),stat=info)

    !-------------------------------------------------------------------------
    ! Set and get minimum and maximum indicies of z pencils
    !-------------------------------------------------------------------------
    indl(1) = 1
    indl(2) = 1
    indl(3) = 1
    indu(1) = 0
    indu(2) = 0
    indu(3) = 0
    DO isub=1,ppmpoisson%nsublistz
       isubl = ppmpoisson%isublistz(isub)
       indu(1) = MAX(indu(1),ppmpoisson%ndataz(1,isubl))
       indu(2) = MAX(indu(2),ppmpoisson%ndataz(2,isubl))
       indu(3) = MAX(indu(3),ppmpoisson%ndataz(3,isubl))
    ENDDO

    !-------------------------------------------------------------------------
    ! Allocate two complex z pencils + Greens fcn array !@check return vars.
    !-------------------------------------------------------------------------
    ALLOCATE(ppmpoisson%fldzc1(__DIM,&
         & indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
         & 1:ppmpoisson%nsublistz),stat=info)

    ALLOCATE(ppmpoisson%fldzc2(__DIM,&
         & indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
         & 1:ppmpoisson%nsublistz),stat=info)

    IF (green .EQ. ppm_poisson_grn_pois_per) THEN
       ALLOCATE(ppmpoisson%fldgrnr(&
            & indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
            & 1:ppmpoisson%nsublistz),stat=info)
    ELSE IF (green .EQ. ppm_poisson_grn_pois_fre) THEN
       !@ freespace solutions are not even remotely implemented
       ALLOCATE(ppmpoisson%fldgrnc(3,&
            & indl(1):indu(1),indl(2):indu(2),indl(3):indu(3),&
            & 1:ppmpoisson%nsublistz),stat=info)
    ENDIF

    !-------------------------------------------------------------------------
    ! Set up xy FFT plans
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
    ! Scaling the spectral koefficients... and normalisation of FFTs
    normfac = -1.0_MK/(4.0_MK*PI*PI * &
         & REAL((mesh%nm(1)-1)*(mesh%nm(2)-1)*(mesh%nm(3)-1),MK))
    IF (green .EQ. ppm_poisson_grn_pois_per) THEN
       DO isub=1,ppmpoisson%nsublistz
          isubl=ppmpoisson%isublistz(isub)
          DO k=1,ppmpoisson%ndataz(3,isubl)
             DO j=1,ppmpoisson%ndataz(2,isubl)
                DO i=1,ppmpoisson%ndataz(1,isubl)
                   kx = i-1 + (ppmpoisson%istartz(1,isubl)-1)
                   ky = j-1 + (ppmpoisson%istartz(2,isubl)-1)
                   kz = k-1 + (ppmpoisson%istartz(3,isubl)-1)
                   !The integer division depends on mesh%nm to include N+1 points:
                   !This is a nasty way to do this but it is only done once so...:
                   IF (kx .GT. (mesh%nm(1)-1)/2) kx = kx-(mesh%nm(1)-1)
                   IF (ky .GT. (mesh%nm(2)-1)/2) ky = ky-(mesh%nm(2)-1)
                   IF (kz .GT. (mesh%nm(3)-1)/2) kz = kz-(mesh%nm(3)-1)
                   ppmpoisson%fldgrnr(i,j,k,isub) = &
                        & normfac/(REAL(kx*kx+ky*ky+kz*kz,__PREC))
                   !Take care of singularity
                   !This is nasty as well
                   IF ((kx*kx+ky*ky+kz*kz) .EQ. 0) THEN
                      ppmpoisson%fldgrnr(i,j,k,isub) = 0.0_MK
                   ENDIF
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
    CALL substop('ppm_poisson_init_predef',t0,info)
    RETURN

  END SUBROUTINE __ROUTINE

