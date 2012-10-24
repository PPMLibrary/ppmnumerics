test_suite ppm_module_ode_typedef

use ppm_module_core
use ppm_module_user_numerics
use ppm_module_ode_typedef

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = ppm_kind_double 
!integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
real(mk),parameter              :: tol=epsilon(1._mk)*1000._mk
real(mk),parameter              :: pi = ACOS(-1._mk)
real(mk),parameter              :: skin = 0._mk
integer,parameter               :: ndim=2
integer                         :: decomp,assig,tolexp
integer                         :: info,comm,rank,nproc,topoid
integer , parameter             :: np_global = 100
integer                         :: ode_scheme
real(mk)                        :: cutoff
real(mk)                        :: dtime
real(mk),dimension(:,:),pointer :: xp=>NULL()
real(mk),dimension(:  ),pointer :: min_phys=>NULL(),max_phys=>NULL()
real(mk),dimension(:  ),pointer :: len_phys=>NULL()
integer                         :: i,j,k,ip,wp_id
integer                         :: nstep
integer,dimension(3)            :: ldc
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost=>NULL()
integer, dimension(:  ),pointer :: nm => NULL()
integer, dimension(:  ),pointer :: ighostsize => NULL()
integer                         :: isymm = 0
real(mk)                        :: t0,t1,t2,t3
integer                         :: seedsize
integer,  dimension(:),allocatable :: seed
integer, dimension(:),pointer   :: nvlist=>NULL()
integer, dimension(:,:),pointer :: vlist=>NULL()
real(mk)                        :: err

integer, dimension(:), pointer                 :: wp_1i => NULL()
integer, dimension(:,:), pointer               :: wp_2i => NULL()
integer(ppm_kind_int64),dimension(:),  pointer :: wp_1li => NULL()
integer(ppm_kind_int64),dimension(:,:),pointer :: wp_2li => NULL()
real(mk), dimension(:),   pointer              :: wp_1r => NULL()
real(mk), dimension(:,:), pointer              :: wp_2r => NULL()
complex(mk), dimension(:),   pointer           :: wp_1c => NULL()
complex(mk), dimension(:,:), pointer           :: wp_2c => NULL()
logical, dimension(:),   pointer               :: wp_1l => NULL()

integer, dimension(:),allocatable              :: degree,order
real(ppm_kind_double),dimension(:),allocatable :: coeffs
integer                                        :: nterms
real(mk)                      :: last_err

    init

        use ppm_module_init
        use ppm_module_mktopo
        
        allocate(ighostsize(ndim),nm(ndim),min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        max_phys(ndim) = 1.4_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic
        ighostsize(1:ndim) = 2 
#ifdef __MPI
        comm = mpi_comm_world
        call mpi_comm_rank(comm,rank,info)
        call mpi_comm_size(comm,nproc,info)
#else
        rank = 0
        nproc = 1
#endif
        tolexp = int(log10(epsilon(1._mk)))+10
        call ppm_init(ndim,mk,tolexp,0,debug,info,99)

        call random_seed(size=seedsize)
        allocate(seed(seedsize))
        do i=1,seedsize
            seed(i)=10+i*i*(rank+1)
        enddo
        call random_seed(put=seed)

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0
        if (ndim.eq.2) then
            cutoff = 0.3_mk
        else
            cutoff = 0.35_mk
        endif
        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,cutoff,cost,info)
    end init


    finalize
        use ppm_module_finalize

        call ppm_finalize(info)

        deallocate(min_phys,max_phys,len_phys)

    end finalize


    setup


    end setup
        

    teardown
        
    end teardown

    test ode_wp_step
        use ppm_module_io_vtk
        type(ppm_t_particles_d), target :: Part1
        class(ppm_t_field_), pointer    :: Field1
        class(ppm_t_main_abstr), pointer :: el
        type(ppm_t_ode)                 :: ode 
        integer                         :: np
        CLASS(ppm_t_operator_discr),POINTER   :: DCop => NULL()
        CLASS(ppm_t_operator_discr),POINTER   :: PSEop => NULL()
        class(ppm_t_neighlist_d_),POINTER :: Nlist => NULL()
        class(ppm_t_discr_data),POINTER :: prop => NULL()
        class(ppm_v_main_abstr),pointer :: fields
        class(ppm_v_var_discr_pair),pointer :: rhs_fields
        real(mk)                        :: t,dt
        procedure(ppm_p_rhsfunc_d),pointer :: rhsptr
        class(ppm_t_var_discr_pair), pointer :: pair

        !--------------------------
        !Define Fields
        !--------------------------
        allocate(ppm_t_field::Field1,stat=info)
        Assert_Equal(info,0)
        call Field1%create(1,info,name="Concentration") !vector field

        np = np_global
        call Part1%initialize(np,info,topoid=topoid,name="Part1")
        Assert_Equal(info,0)

!  print particles to a VTK file
!        CALL ppm_vtk_particles("part_test",Part1,info)
!        Assert_Equal(info,0)

        !call Part1%set_cutoff(1.0_mk * Part1%h_avg,info)
        !Assert_Equal(info,0)

        call Part1%map(info,global=.true.,topoid=topoid)
        Assert_Equal(info,0)

        call Part1%map_ghosts(info)
        Assert_Equal(info,0)

        call Field1%discretize_on(Part1,info)
        Assert_Equal(info,0)

        call Part1%map_ghosts(info)
        Assert_Equal(info,0)

        call Part1%get_field(Field1,wp_1r,info)
        Assert_Equal(info,0)
        wp_1r(:) = 1.0_mk

        allocate(fields,stat=info)
        Assert_Equal(info,0)
        allocate(rhs_fields,stat=info)
        Assert_Equal(info,0)
        allocate(pair,stat=info)
        Assert_Equal(info,0)
        el => Field1
        call fields%push(el,info)
        pair%var => Field1
        pair%discr => Part1
        call rhs_fields%push(pair,info)

        rhsptr => rhs_test1
        call ode%create(ppm_param_ode_scheme_eulerf,fields,rhsptr,rhs_fields,info)
        Assert_Equal(info,0)
        t = 0.0_mk
        dt = 0.1_mk
        call ode%step(t,dt,1,info)
        Assert_Equal(info,0)
        
        call Part1%get_field(Field1,wp_1r,info)
        Assert_Equal(info,0)
        Assert_True(((minval(wp_1r).eq.1.2_mk).and.(maxval(wp_1r).eq.1.2_mk)))
        

        call ode%destroy(info)
        Assert_Equal(info,0)
        call Part1%destroy(info)
        Assert_Equal(info,0)
        call Field1%destroy(info)
        Assert_Equal(info,0)
        deallocate(Field1,fields,rhs_fields,pair)
    end test
    
    test ode_xp_step
        use ppm_module_io_vtk
        type(ppm_t_particles_d), target :: Part1
        class(ppm_t_field_), pointer    :: Field1
        class(ppm_t_main_abstr), pointer :: el
        type(ppm_t_ode)                 :: ode 
        integer                         :: np
        CLASS(ppm_t_operator_discr),POINTER   :: DCop => NULL()
        CLASS(ppm_t_operator_discr),POINTER   :: PSEop => NULL()
        class(ppm_t_neighlist_d_),POINTER :: Nlist => NULL()
        class(ppm_t_discr_data),POINTER :: prop => NULL()
        class(ppm_v_main_abstr),pointer :: fields
        class(ppm_v_var_discr_pair),pointer :: rhs_fields
        real(mk)                        :: t,dt
        procedure(ppm_p_rhsfunc_d),pointer :: rhsptr
        class(ppm_t_var_discr_pair), pointer :: pair
        real(mk),dimension(:,:),pointer :: moved_xp=>NULL()

        start_subroutine("ode_xp_step")
        !--------------------------
        ! Define Fields
        !--------------------------
        allocate(ppm_t_field::Field1,stat=info)
        Assert_Equal(info,0)
        call Field1%create(1,info,name="Concentration") !vector field
        np = np_global
        call Part1%initialize(np,info,topoid=topoid,name="Part1")
        Assert_Equal(info,0)

!  print particles to a VTK file
!        CALL ppm_vtk_particles("part_test",Part1,info)
!        Assert_Equal(info,0)

        !call Part1%set_cutoff(1.0_mk * Part1%h_avg,info)
        !Assert_Equal(info,0)

        call Part1%map(info,global=.true.,topoid=topoid)
        Assert_Equal(info,0)

        call Part1%map_ghosts(info)
        Assert_Equal(info,0)

        call Field1%discretize_on(Part1,info)
        Assert_Equal(info,0)

        call Part1%map_ghosts(info)
        Assert_Equal(info,0)

        call Part1%get_field(Field1,wp_1r,info)
        Assert_Equal(info,0)
        wp_1r(:) = 1.0_mk

        call Part1%get_xp(moved_xp,info)
        allocate(xp(ppm_dim,Part1%Npart))
        xp(:,:) = moved_xp(:,:)
        
        allocate(fields,stat=info)
        Assert_Equal(info,0)
        allocate(rhs_fields,stat=info)
        Assert_Equal(info,0)
        allocate(pair,stat=info)
        Assert_Equal(info,0)
        el => Part1
        call fields%push(el,info)
        pair%var => Field1
        pair%discr => Part1
        call rhs_fields%push(pair,info)

        rhsptr => rhs_test2
        call ode%create(ppm_param_ode_scheme_eulerf,fields,rhsptr,rhs_fields,info)
        Assert_Equal(info,0)
        t = 0.0_mk
        dt = 0.1_mk
        call ode%step(t,dt,1,info)
        Assert_Equal(info,0)
       
        call Part1%get_xp(moved_xp,info)
        do i=1,Part1%Npart
          xp(:,i) = moved_xp(:,i) - xp(:,i)
        enddo
        Assert_True(((abs(minval(xp(:,1:Part1%Npart))-0.2_mk).le.tol).and.(abs(maxval(xp(:,1:Part1%Npart))-0.2_mk).le.tol)))

        call ode%destroy(info)
        Assert_Equal(info,0)
        call Part1%destroy(info)
        Assert_Equal(info,0)
        call Field1%destroy(info)
        Assert_Equal(info,0)
        deallocate(Field1,fields,rhs_fields,pair,xp)

        end_subroutine()
    end test
    
    test ode_wp_xp_step
        use ppm_module_io_vtk
        type(ppm_t_particles_d), target :: Part1
        class(ppm_t_field_), pointer    :: Field1
        class(ppm_t_main_abstr), pointer :: el
        type(ppm_t_ode)                 :: ode 
        integer                         :: np
        CLASS(ppm_t_operator_discr),POINTER   :: DCop => NULL()
        CLASS(ppm_t_operator_discr),POINTER   :: PSEop => NULL()
        class(ppm_t_neighlist_d_),POINTER :: Nlist => NULL()
        class(ppm_t_discr_data),POINTER :: prop => NULL()
        class(ppm_v_main_abstr),pointer :: fields
        class(ppm_v_var_discr_pair),pointer :: rhs_fields
        real(mk)                        :: t,dt
        procedure(ppm_p_rhsfunc_d),pointer :: rhsptr
        class(ppm_t_var_discr_pair), pointer :: pair
        real(mk),dimension(:,:),pointer :: moved_xp=>NULL()

        start_subroutine("ode_xp_wp,step")
        !--------------------------
        ! Define Fields
        !--------------------------
        allocate(ppm_t_field::Field1,stat=info)
        Assert_Equal(info,0)
        call Field1%create(1,info,name="Concentration") !vector field
        np = np_global
        call Part1%initialize(np,info,topoid=topoid,name="Part1")
        Assert_Equal(info,0)

!  print particles to a VTK file
!        CALL ppm_vtk_particles("part_test",Part1,info)
!        Assert_Equal(info,0)

        !call Part1%set_cutoff(1.0_mk * Part1%h_avg,info)
        !Assert_Equal(info,0)

        call Part1%map(info,global=.true.,topoid=topoid)
        Assert_Equal(info,0)

        call Part1%map_ghosts(info)
        Assert_Equal(info,0)

        call Field1%discretize_on(Part1,info)
        Assert_Equal(info,0)

        call Part1%map_ghosts(info)
        Assert_Equal(info,0)

        call Part1%get_field(Field1,wp_1r,info)
        Assert_Equal(info,0)
        wp_1r(:) = 1.0_mk

        call Part1%get_xp(moved_xp,info)
        allocate(xp(ppm_dim,Part1%Npart))
        xp(:,:) = moved_xp(:,:)
        
        allocate(fields,stat=info)
        Assert_Equal(info,0)
        allocate(rhs_fields,stat=info)
        Assert_Equal(info,0)
        allocate(pair,stat=info)
        Assert_Equal(info,0)
        el => Part1
        call fields%push(el,info)
        Assert_Equal(info,0)
        el => Field1
        call fields%push(el,info)
        Assert_Equal(info,0)
        pair%var => Field1
        pair%discr => Part1
        call rhs_fields%push(pair,info)
        Assert_Equal(info,0)

        rhsptr => rhs_test_xp_wp
        call ode%create(ppm_param_ode_scheme_eulerf,fields,rhsptr,rhs_fields,info)
        Assert_Equal(info,0)
        t = 0.0_mk
        dt = 0.1_mk
        call ode%step(t,dt,1,info)
        Assert_Equal(info,0)
       
        call Part1%get_xp(moved_xp,info)
        do i=1,Part1%Npart
          xp(:,i) = moved_xp(:,i) - xp(:,i)
        enddo
        Assert_True(((abs(minval(xp(:,1:Part1%Npart))-0.2_mk).le.tol).and.(abs(maxval(xp(:,1:Part1%Npart))-0.2_mk).le.tol)))

        call ode%destroy(info)
        Assert_Equal(info,0)
        call Part1%destroy(info)
        Assert_Equal(info,0)
        call Field1%destroy(info)
        Assert_Equal(info,0)
        deallocate(Field1,fields,rhs_fields,pair,xp)

        end_subroutine()
    end test
    
    test ode_mesh_step
        class(ppm_t_field_), POINTER        :: Field1,Field2
        class(ppm_t_equi_mesh),POINTER       :: Mesh1
        type(ppm_t_ode)                     :: ode 
        integer                             :: p_idx, nb_errors
        CLASS(ppm_t_discr_info_),POINTER    :: dinfo => NULL()
        logical                             :: assoc
        real(mk),dimension(2*ndim)          :: my_patch
        real(mk),dimension(ndim)            :: offset
        real(mk)                            :: t,dt
        procedure(ppm_p_rhsfunc_d),pointer    :: rhsptr
        class(ppm_v_main_abstr),pointer     :: fields
        class(ppm_v_var_discr_pair),pointer :: rhs_fields
        class(ppm_t_var_discr_pair),pointer :: fpair1
        class(ppm_t_var_discr_pair),pointer :: fpair2
        class(ppm_t_main_abstr), pointer :: el

        start_subroutine("ode_mesh_step")

        Nm = 25
        Nm(ndim) = 45

        allocate(Mesh1,stat=info)
        call Mesh1%create(topoid,offset,info,Nm=Nm,&
            ghostsize=ighostsize,name='Test_Mesh_1')
        Assert_Equal(info,0)

        if (ndim.eq.2) then
            my_patch(1:4) = (/0.15_mk,0.10_mk,0.99_mk,0.7_mk/)
        else
            my_patch(1:6) = (/0.15_mk,0.10_mk,0.51_mk,0.99_mk,0.7_mk,0.78_mk/)
        endif

        call Mesh1%def_patch(my_patch,info) 
        Assert_Equal(info,0)

        allocate(ppm_t_field::Field1,stat=info)
        Assert_Equal(info,0)
        allocate(ppm_t_field::Field2,stat=info)
        Assert_Equal(info,0)
        call Field1%create(2,info,name='vecField') 
            Assert_Equal(info,0)
        call Field1%discretize_on(Mesh1,info)
            Assert_Equal(info,0)
        call Field2%create(1,info,name='scaField') 
            Assert_Equal(info,0)
        call Field2%discretize_on(Mesh1,info)
            Assert_Equal(info,0)

        IF (ndim.EQ.2) THEN
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j)
            for real
                Field1_n(:) = 1._mk
                Field2_n    = 2._mk
        end foreach
        ELSE
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j,k)
            for real
                Field1_n(:) = 1._mk
                Field2_n    = 2._mk
        end foreach
        ENDIF

        !Do a ghost mapping
        call Mesh1%map_ghost_get(info)
            Assert_Equal(info,0)

        call Field1%map_ghost_push(Mesh1,info)
            Assert_Equal(info,0)
        call Field2%map_ghost_push(Mesh1,info)
            Assert_Equal(info,0)

        call Mesh1%map_send(info)
            Assert_Equal(info,0)

        call Field2%map_ghost_pop(Mesh1,info)
            Assert_Equal(info,0)
        call Field1%map_ghost_pop(Mesh1,info)
            Assert_Equal(info,0)
        
        allocate(fields,stat=info)
        Assert_Equal(info,0)
        allocate(rhs_fields,stat=info)
        Assert_Equal(info,0)
        el => Field2
        call fields%push(el,info)
        allocate(fpair1,stat=info)
        Assert_Equal(info,0)
        fpair1%var => Field1
        fpair1%discr => Mesh1
        call rhs_fields%push(fpair1,info)
        allocate(fpair2,stat=info)
        Assert_Equal(info,0)
        fpair2%var => Field2
        fpair2%discr => Mesh1
        call rhs_fields%push(fpair2,info)

        rhsptr => rhs_test3
        call ode%create(ppm_param_ode_scheme_eulerf,fields,rhsptr,rhs_fields,info)
        Assert_Equal(info,0)
        t = 0.0_mk
        dt = 0.1_mk
        call ode%step(t,dt,1,info)
        Assert_Equal(info,0)
        
        IF (ndim.EQ.2) THEN
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j)
            for real
                Assert_Equal(Field2_n,2.6_mk)
        end foreach
        ELSE
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2) vec_fields(Field1) indices(i,j,k)
            for real
                Assert_Equal(Field2_n,2.6_mk)
        end foreach
        ENDIF

        call Mesh1%destroy(info)
            Assert_Equal(info,0)
        call Field1%destroy(info)
            Assert_Equal(info,0)
        call Field2%destroy(info)
            Assert_Equal(info,0)

        end_subroutine()
    end test

!-------------------------------------------------------------
! right hand sides
!-------------------------------------------------------------

real(mk) function rhs_test1(fields_discr,t,changes)
  class(ppm_v_var_discr_pair),   pointer    :: fields_discr
  real(ppm_kind_double)                     :: t
  class(ppm_v_main_abstr),       pointer    :: changes
  class(ppm_t_main_abstr),       pointer    :: c
  class(ppm_t_var_discr_pair),   pointer    :: pair
  real(mk), dimension(:),        pointer    :: wp => null()
  real(mk), dimension(:),        pointer    :: dwp => null()
  class(ppm_t_field_),           pointer    :: field
  class(ppm_t_field_),           pointer    :: df
  class(ppm_t_particles_d),      pointer    :: pset => null()
  start_subroutine("rhs_test1")
  
  pair => fields_discr%at(1)
  select type(d => pair%discr)
  class is (ppm_t_particles_d)
    pset => d
  end select
  select type(f => pair%var)
  class is (ppm_t_field_)
    field => f
  end select

  check_associated(pset,"type mismatch")
  c => changes%at(1)
  select type(c)
  class is (ppm_t_field_)
    df => c
  end select
  
  foreach p in particles(pset) with sca_fields(w=field,dw=df)
    dw_p = 2.0_mk*w_p
  end foreach

  rhs_test1 = 0 
  end_subroutine()
end function rhs_test1

real(mk) function rhs_test2(fields_discr,t,changes)
  class(ppm_v_var_discr_pair), pointer    :: fields_discr
  real(ppm_kind_double)                     :: t
  class(ppm_v_main_abstr),      pointer          :: changes
  class(ppm_t_main_abstr), pointer          :: c
  class(ppm_t_var_discr_pair), pointer      :: pair
  real(mk), dimension(:), pointer           :: wp => null()
  real(mk), dimension(:,:), pointer         :: dxp => null()
  class(ppm_t_field_),     pointer          :: field
  class(ppm_t_part_prop_d_),     pointer    :: df
  class(ppm_t_discr_info_), pointer         :: di => null()
  class(ppm_t_particles_d), pointer         :: pset => null()
  start_subroutine("rhs_test2")
  
  pair => fields_discr%at(1)
  select type(d => pair%discr)
  class is (ppm_t_particles_d)
    pset => d
  end select
  select type(f => pair%var)
  class is (ppm_t_field_)
    field => f
  end select

  check_associated(pset,"type mismatch")
  c => changes%at(1)
  select type(c)
  class is (ppm_t_part_prop_d_)
    df => c
  end select

  foreach p in particles(pset) with sca_fields(w=field) vec_props(dx=df)
    dx_p(:) = 2.0_mk*w_p
  end foreach

  rhs_test2 = 0 
  end_subroutine()
end function rhs_test2

real(mk) function rhs_test_xp_wp(fields_discr,t,changes)
  class(ppm_v_var_discr_pair), pointer    :: fields_discr
  real(ppm_kind_double)                     :: t
  class(ppm_v_main_abstr),      pointer          :: changes
  class(ppm_t_main_abstr), pointer          :: c
  class(ppm_t_var_discr_pair), pointer      :: pair
  real(mk), dimension(:), pointer           :: wp => null()
  real(mk), dimension(:,:), pointer         :: dxp => null()
  class(ppm_t_field_),     pointer          :: field
  class(ppm_t_part_prop_d_),     pointer    :: dxi
  class(ppm_t_field_),     pointer          :: df
  class(ppm_t_discr_info_), pointer         :: di => null()
  class(ppm_t_particles_d), pointer         :: pset => null()
  start_subroutine("rhs_test_xp_wp")
  
  pair => fields_discr%at(1)
  select type(d => pair%discr)
  class is (ppm_t_particles_d)
    pset => d
  end select
  select type(f => pair%var)
  class is (ppm_t_field_)
    field => f
  end select
  check_associated(pset,"type mismatch")

  c => changes%at(1)
  select type(c)
  class is (ppm_t_part_prop_d_)
    dxi => c
  end select
  c => changes%at(2)
  select type(c)
  class is (ppm_t_field_)
    df => c
  end select
  
  foreach p in particles(pset) with sca_fields(w=field) vec_props(dx=dxi)
    dx_p(:) = 2.0_mk*w_p
  end foreach
  
  foreach p in particles(pset) with sca_fields(dw=df)
    dw_p = 1.2_mk
  end foreach

  rhs_test_xp_wp = 0 
  end_subroutine()
end function rhs_test_xp_wp

real(mk) function rhs_test3(fields_discr,t,changes)
  class(ppm_v_var_discr_pair), pointer    :: fields_discr
  real(ppm_kind_double)                     :: t
  class(ppm_v_main_abstr),      pointer          :: changes
  class(ppm_t_main_abstr), pointer          :: c
  class(ppm_t_var_discr_pair), pointer    :: pair
  class(ppm_t_field_),     pointer          :: field1,field2,dfield
  class(ppm_t_equi_mesh), pointer           :: mesh => null()
  start_subroutine("rhs_test3")
  
  pair => fields_discr%at(1)
  select type(f => pair%var)
  class is (ppm_t_field_)
    field1 => f
  end select
  select type(d => pair%discr)
  class is (ppm_t_equi_mesh)
    mesh => d
  end select
  check_associated(mesh,"type mismatch")
  pair => fields_discr%at(2)
  select type(f => pair%var)
  class is (ppm_t_field_)
    field2 => f
  end select
  c => changes%at(1)
  select type(c)
  class is (ppm_t_field_)
    dfield => c
  end select
        
  IF (ndim.EQ.2) THEN
    foreach n in equi_mesh(mesh) with sca_fields(field2,dfield) vec_fields(field1) indices(i,j)
      for real
        dfield_n = 3._mk*field2_n
    end foreach
  ELSE
    foreach n in equi_mesh(mesh) with sca_fields(field2,dfield) vec_fields(field1) indices(i,j,k)
      for real
        dfield_n = 3._mk*field2_n
    end foreach
  ENDIF
  

  rhs_test3 = 0 
  end_subroutine()
end function rhs_test3

!-------------------------------------------------------------
! ODE convergence tests
!-------------------------------------------------------------
    
    test test_ode_const({ode_scheme: [ppm_param_ode_scheme_eulerf,ppm_param_ode_scheme_tvdrk2,ppm_param_ode_scheme_rk4]})
        use ppm_module_io_vtk
        type(ppm_t_particles_d), target :: Part1
        class(ppm_t_field_), pointer    :: Field1
        class(ppm_t_main_abstr), pointer :: el
        type(ppm_t_ode)                 :: ode 
        class(ppm_v_main_abstr),pointer :: fields
        class(ppm_v_var_discr_pair),pointer :: rhs_fields
        real(mk)                        :: t,dt
        real(mk)         , parameter    :: ts = 3.0_mk
        real(mk)         , parameter    :: te = 3.5_mk
        integer                         :: istage
        integer                         :: np
        procedure(ppm_p_rhsfunc_d),pointer :: rhsptr
        class(ppm_t_field_discr_pair), pointer :: pair

        !--------------------------
        !Define Fields
        !--------------------------
        allocate(ppm_t_field::Field1,stat=info)
        Assert_Equal(info,0)
        call Field1%create(1,info,name="Concentration") !vector field

        np = np_global
        call Part1%initialize(np,info,topoid=topoid,name="Part1")
        Assert_Equal(info,0)

!  print particles to a VTK file
!        CALL ppm_vtk_particles("part_test",Part1,info)
!        Assert_Equal(info,0)

        !call Part1%set_cutoff(1.0_mk * Part1%h_avg,info)
        !Assert_Equal(info,0)

        call Part1%map(info,global=.true.,topoid=topoid)
        Assert_Equal(info,0)

        call Part1%map_ghosts(info)
        Assert_Equal(info,0)

        call Field1%discretize_on(Part1,info)
        Assert_Equal(info,0)

        call Part1%map_ghosts(info)
        Assert_Equal(info,0)

        call Part1%get_field(Field1,wp_1r,info)
        Assert_Equal(info,0)
        wp_1r(:) = 1.1_mk

        allocate(fields,stat=info)
        Assert_Equal(info,0)
        allocate(rhs_fields,stat=info)
        Assert_Equal(info,0)
        allocate(pair,stat=info)
        Assert_Equal(info,0)
        el => Field1
        call fields%push(el,info)

        rhsptr => const_ode
        call ode%create(ode_scheme,fields,rhsptr,rhs_fields,info)
        Assert_Equal(info,0)
        t = ts
        dt = 0.1_mk
        do while (t.lt.te-dt/2._mk)
          do istage=1,ode%integrator%scheme_nstages
            call ode%step(t,dt,istage,info)
          end do
        end do
        Assert_Equal(info,0)
        
        call Part1%get_field(Field1,wp_1r,info)
        Assert_Equal(info,0)
        Assert_True(((abs(minval(wp_1r)-(1.7_mk*0.5_mk + 1.1_mk)).lt.tol).and.(abs(maxval(wp_1r)-(1.7_mk*0.5_mk + 1.1_mk)).lt.tol)))
        

        call ode%destroy(info)
        Assert_Equal(info,0)
        call Part1%destroy(info)
        Assert_Equal(info,0)
        call Field1%destroy(info)
        Assert_Equal(info,0)
        deallocate(Field1,fields,rhs_fields,pair)
    end test

real(mk) function const_ode(fields_discr,t,changes)
  class(ppm_v_var_discr_pair), pointer    :: fields_discr
  real(ppm_kind_double)                   :: t
  class(ppm_v_main_abstr),     pointer    :: changes
  class(ppm_t_field_),         pointer    :: df
  class(ppm_t_main_abstr), pointer          :: c
  class(ppm_t_discr_info_),    pointer    :: di => null()
  class(ppm_t_particles_d),    pointer    :: pset => null()
  start_subroutine("const_ode")

  c => changes%at(1)
  select type(c)
  class is (ppm_t_field_)
    df => c
  end select
  di => df%discr_info%begin()
  select type(disc => di%discr_ptr)
  class is (ppm_t_particles_d)
    pset => disc
  end select
  
  foreach p in particles(pset) with sca_fields(dw=df)
    dw_p = 1.7_mk
  end foreach

  const_ode = 0 
  end_subroutine()
end function const_ode
    
     test test_ode_linear({dtime: [0.1_mk,0.01_mk,0.001_mk,0.0001_mk], ode_scheme: [ppm_param_ode_scheme_eulerf,ppm_param_ode_scheme_tvdrk2,ppm_param_ode_scheme_midrk2,ppm_param_ode_scheme_rk4]})
        use ppm_module_io_vtk
        type(ppm_t_particles_d), target :: Part1
        class(ppm_t_field_), pointer    :: Field1
        class(ppm_t_main_abstr), pointer :: el
        type(ppm_t_ode)                 :: ode 
        class(ppm_v_main_abstr),pointer :: fields
        class(ppm_v_var_discr_pair),pointer :: rhs_fields
        real(mk)                        :: t,dt
        real(mk)         , parameter    :: ts = 0.0_mk
        real(mk)         , parameter    :: te = 1.0_mk
        integer                         :: istage
        integer                         :: np
        real(mk)        , parameter     :: symres = 0.5_mk
        procedure(ppm_p_rhsfunc_d),pointer :: rhsptr

        !--------------------------
        !Define Fields
        !--------------------------
        allocate(ppm_t_field::Field1,stat=info)
        Assert_Equal(info,0)
        call Field1%create(1,info,name="Concentration") !vector field

        np = np_global
        call Part1%initialize(np,info,topoid=topoid,name="Part1")
        Assert_Equal(info,0)

!  print particles to a VTK file
!        CALL ppm_vtk_particles("part_test",Part1,info)
!        Assert_Equal(info,0)

        !call Part1%set_cutoff(1.0_mk * Part1%h_avg,info)
        !Assert_Equal(info,0)

        call Part1%map(info,global=.true.,topoid=topoid)
        Assert_Equal(info,0)

        call Part1%map_ghosts(info)
        Assert_Equal(info,0)

        call Field1%discretize_on(Part1,info)
        Assert_Equal(info,0)

        call Part1%map_ghosts(info)
        Assert_Equal(info,0)

        call Part1%get_field(Field1,wp_1r,info)
        Assert_Equal(info,0)
        wp_1r(:) = 0.0_mk

        allocate(fields,stat=info)
        Assert_Equal(info,0)
        allocate(rhs_fields,stat=info)
        Assert_Equal(info,0)
        el => Field1
        call fields%push(el,info)

        rhsptr => linear_ode
        call ode%create(ode_scheme,fields,rhsptr,rhs_fields,info)
        Assert_Equal(info,0)
        t = ts
        dt = dtime
        do while (t.lt.te-dt/2._mk)
          do istage=1,ode%integrator%scheme_nstages
            call ode%step(t,dt,istage,info)
          end do
          !call Part1%get_field(Field1,wp_1r,info)
          !Assert_Equal(info,0)
          !print *,t,wp_1r(1)
        end do
        Assert_Equal(info,0)
        
        call Part1%get_field(Field1,wp_1r,info)
        Assert_Equal(info,0)
        if (ode%integrator%scheme_order.gt.1) then
          Assert_True(abs(symres - maxval(wp_1r))/symres.lt.tol)
        else
          Assert_True((dtime.eq.0.1_mk).or.(log10(last_err / abs(symres - maxval(wp_1r))/symres).ge.1.0_mk))
          last_err = abs(symres - maxval(wp_1r))/symres
        end if
        !Assert_True(((abs(minval(wp_1r)-(1.7_mk*0.5_mk + 1.1_mk)).lt.tol).and.(abs(maxval(wp_1r)-(1.7_mk*0.5_mk + 1.1_mk)).lt.tol)))
        

        call ode%destroy(info)
        Assert_Equal(info,0)
        call Part1%destroy(info)
        Assert_Equal(info,0)
        call Field1%destroy(info)
        Assert_Equal(info,0)
        deallocate(Field1,fields,rhs_fields)
    end test

real(mk) function linear_ode(fields_discr,t,changes)
  class(ppm_v_var_discr_pair),   pointer    :: fields_discr
  real(ppm_kind_double)                     :: t
  class(ppm_v_main_abstr),       pointer    :: changes
  class(ppm_t_field_),           pointer    :: df
  class(ppm_t_main_abstr), pointer          :: c
  class(ppm_t_discr_info_),      pointer    :: di => null()
  class(ppm_t_particles_d),      pointer    :: pset => null()
  start_subroutine("linear_ode")

  c => changes%at(1)
  select type(c)
  class is (ppm_t_field_)
    df => c
  end select
  di => df%discr_info%begin()
  select type(disc => di%discr_ptr)
  class is (ppm_t_particles_d)
    pset => disc
  end select
  
  foreach p in particles(pset) with sca_fields(dw=df)
    dw_p = t
  end foreach

  linear_ode = 0 
  end_subroutine()
end function linear_ode
     
     test test_ode_exp({dtime: [0.1_mk,0.01_mk,0.001_mk,0.0001_mk], ode_scheme: [ppm_param_ode_scheme_eulerf,ppm_param_ode_scheme_tvdrk2,ppm_param_ode_scheme_midrk2,ppm_param_ode_scheme_rk4]})
        use ppm_module_io_vtk
        type(ppm_t_particles_d), target :: Part1
        class(ppm_t_field_), pointer    :: Field1
        class(ppm_t_main_abstr), pointer :: el
        type(ppm_t_ode)                 :: ode 
        class(ppm_v_main_abstr),pointer :: fields
        class(ppm_v_var_discr_pair),pointer :: rhs_fields
        real(mk)                        :: t,dt
        real(mk)         , parameter    :: ts = 0.0_mk
        real(mk)         , parameter    :: te = 1.0_mk
        integer                         :: istage
        integer                         :: np
        procedure(ppm_p_rhsfunc_d),pointer :: rhsptr
        real(mk)        , parameter     :: symres = exp(1.0_mk)-1.0_mk
        real(mk)                        :: rel_max_err 

        !--------------------------
        !Define Fields
        !--------------------------
        allocate(ppm_t_field::Field1,stat=info)
        Assert_Equal(info,0)
        call Field1%create(1,info,name="Concentration") !vector field

        np = np_global
        call Part1%initialize(np,info,topoid=topoid,name="Part1")
        Assert_Equal(info,0)

!  print particles to a VTK file
!        CALL ppm_vtk_particles("part_test",Part1,info)
!        Assert_Equal(info,0)

        !call Part1%set_cutoff(1.0_mk * Part1%h_avg,info)
        !Assert_Equal(info,0)

        call Part1%map(info,global=.true.,topoid=topoid)
        Assert_Equal(info,0)

        call Part1%map_ghosts(info)
        Assert_Equal(info,0)

        call Field1%discretize_on(Part1,info)
        Assert_Equal(info,0)

        call Part1%map_ghosts(info)
        Assert_Equal(info,0)

        call Part1%get_field(Field1,wp_1r,info)
        Assert_Equal(info,0)
        wp_1r(:) = 0.0_mk

        allocate(fields,stat=info)
        Assert_Equal(info,0)
        allocate(rhs_fields,stat=info)
        Assert_Equal(info,0)
        el => Field1
        call fields%push(el,info)

        rhsptr => exp_ode
        call ode%create(ode_scheme,fields,rhsptr,rhs_fields,info)
        Assert_Equal(info,0)
        t = ts
        dt = dtime
        do while (t.lt.te-dt/2._mk)
          do istage=1,ode%integrator%scheme_nstages
            call ode%step(t,dt,istage,info)
          end do
          !call Part1%get_field(Field1,wp_1r,info)
          !Assert_Equal(info,0)
          !print *,t,wp_1r(1)
        end do
        Assert_Equal(info,0)
        
        call Part1%get_field(Field1,wp_1r,info)
        Assert_Equal(info,0)
        rel_max_err = abs(symres - maxval(wp_1r))/symres
        ! the 0.01 is an additonal tolerance used to account for the asymptotic
        ! convergence towards the order and numerical inaccuracies.
        Assert_True((dtime.eq.0.1_mk).or.(rel_max_err.lt.tol).or.(abs(log10(last_err / rel_max_err)-real(ode%integrator%scheme_order,mk)).lt.0.01_mk))
        last_err =rel_max_err 
        

        call ode%destroy(info)
        Assert_Equal(info,0)
        call Part1%destroy(info)
        Assert_Equal(info,0)
        call Field1%destroy(info)
        Assert_Equal(info,0)
        deallocate(Field1,fields,rhs_fields)
    end test

real(mk) function exp_ode(fields_discr,t,changes)
  class(ppm_v_var_discr_pair),   pointer    :: fields_discr
  real(ppm_kind_double)                     :: t
  class(ppm_v_main_abstr),      pointer     :: changes
  class(ppm_t_field_),          pointer     :: df
  class(ppm_t_main_abstr), pointer          :: c
  class(ppm_t_discr_info_),     pointer     :: di => null()
  class(ppm_t_particles_d),     pointer     :: pset => null()
  start_subroutine("exp_ode")

  c => changes%at(1)
  select type(c)
  class is (ppm_t_field_)
    df => c
  end select
  di => df%discr_info%begin()
  select type(disc => di%discr_ptr)
  class is (ppm_t_particles_d)
    pset => disc
  end select
  
  foreach p in particles(pset) with sca_fields(dw=df)
    dw_p = exp(t)
  end foreach

  exp_ode = 0 
  end_subroutine()
end function exp_ode

!-------------------------------------------------------------
! ODE convergence test for meshes
!-------------------------------------------------------------
     
      test test_ode_mesh_exp({dtime: [0.1_mk,0.01_mk,0.001_mk,0.0001_mk], ode_scheme: [ppm_param_ode_scheme_eulerf,ppm_param_ode_scheme_tvdrk2,ppm_param_ode_scheme_midrk2,ppm_param_ode_scheme_rk4]})
        use ppm_module_io_vtk
        class(ppm_t_equi_mesh),pointer  :: Mesh1
        class(ppm_t_field_), pointer    :: Field1
        class(ppm_t_main_abstr), pointer:: el
        type(ppm_t_ode)                 :: ode 
        class(ppm_v_main_abstr),pointer :: fields
        class(ppm_v_var_discr_pair),pointer :: rhs_fields
        real(mk)                        :: t,dt
        real(mk)         , parameter    :: ts = 0.0_mk
        real(mk)         , parameter    :: te = 1.0_mk
        real(mk),dimension(2*ndim)          :: my_patch
        real(mk),dimension(ndim)            :: offset
        integer                         :: istage
        integer                         :: np
        procedure(ppm_p_rhsfunc_d),pointer :: rhsptr
        real(mk)        , parameter     :: symres = exp(1.0_mk)-1.0_mk
        real(mk)                        :: rel_max_err 
  
        start_subroutine("test_ode_mesh_exp")
        Nm = 25
        Nm(ndim) = 45

        allocate(Mesh1,stat=info)
        call Mesh1%create(topoid,offset,info,Nm=Nm,&
            ghostsize=ighostsize,name='Test_Mesh_1')
        Assert_Equal(info,0)

        if (ndim.eq.2) then
            my_patch(1:4) = (/0.15_mk,0.10_mk,0.99_mk,0.7_mk/)
        else
            my_patch(1:6) = (/0.15_mk,0.10_mk,0.51_mk,0.99_mk,0.7_mk,0.78_mk/)
        endif

        call Mesh1%def_patch(my_patch,info) 
        Assert_Equal(info,0)
        !--------------------------
        !Define Fields
        !--------------------------
        allocate(ppm_t_field::Field1,stat=info)
        Assert_Equal(info,0)
        call Field1%create(1,info,name="Concentration") !vector field


        call Field1%discretize_on(Mesh1,info)
        Assert_Equal(info,0)
        
        IF (ndim.EQ.2) THEN
        foreach n in equi_mesh(Mesh1) with sca_fields(Field1) indices(i,j)
            for real
                Field1_n    = 0._mk
        end foreach
        ELSE
        foreach n in equi_mesh(Mesh1) with sca_fields(Field1) indices(i,j,k)
            for real
                Field1_n    = 0._mk
        end foreach
        ENDIF

        !Do a ghost mapping
        call Mesh1%map_ghost_get(info)
            Assert_Equal(info,0)

        call Field1%map_ghost_push(Mesh1,info)
            Assert_Equal(info,0)

        call Mesh1%map_send(info)
            Assert_Equal(info,0)

        call Field1%map_ghost_pop(Mesh1,info)
            Assert_Equal(info,0)

        allocate(fields,stat=info)
        Assert_Equal(info,0)
        allocate(rhs_fields,stat=info)
        Assert_Equal(info,0)
        el => Field1
        call fields%push(el,info)

        rhsptr => exp_mesh_ode
        call ode%create(ode_scheme,fields,rhsptr,rhs_fields,info)
        Assert_Equal(info,0)
        t = ts
        dt = dtime
        do while (t.lt.te-dt/2._mk)
          do istage=1,ode%integrator%scheme_nstages
            call ode%step(t,dt,istage,info)
          end do
          !call Part1%get_field(Field1,wp_1r,info)
          !Assert_Equal(info,0)
          !print *,t,wp_1r(1)
        end do
        Assert_Equal(info,0)
        
        rel_max_err = 0.0_mk
        IF (ndim.EQ.2) THEN
        foreach n in equi_mesh(Mesh1) with sca_fields(Field1) indices(i,j)
            for real
              if (rel_max_err.lt.Field1_n) then
                rel_max_err = Field1_n
              end if
        end foreach
        ELSE
        foreach n in equi_mesh(Mesh1) with sca_fields(Field1) indices(i,j,k)
            for real
              if (rel_max_err.lt.Field1_n) then
                rel_max_err = Field1_n
              end if
        end foreach
        ENDIF
        rel_max_err = abs(symres - rel_max_err)/symres

        ! the 0.01 is an additonal tolerance used to account for the asymptotic
        ! convergence towards the order and numerical inaccuracies.
        Assert_True((dtime.eq.0.1_mk).or.(rel_max_err.lt.tol).or.(abs(log10(last_err / rel_max_err)-real(ode%integrator%scheme_order,mk)).lt.0.01_mk))
        last_err =rel_max_err 
        

        call ode%destroy(info)
        Assert_Equal(info,0)
        call Mesh1%destroy(info)
        Assert_Equal(info,0)
        call Field1%destroy(info)
        Assert_Equal(info,0)
        deallocate(Field1,fields,rhs_fields)

        end_subroutine()
    end test

real(mk) function exp_mesh_ode(fields_discr,t,changes)
  class(ppm_v_var_discr_pair), pointer    :: fields_discr
  real(ppm_kind_double)                   :: t
  class(ppm_v_main_abstr),     pointer    :: changes
  class(ppm_t_field_),         pointer    :: df
  class(ppm_t_main_abstr), pointer          :: c
  class(ppm_t_discr_info_),    pointer    :: di => null()
  class(ppm_t_equi_mesh),      pointer    :: mesh => null()
  start_subroutine("exp_mesh_ode")

  c => changes%at(1)
  select type(c)
  class is (ppm_t_field_)
    df => c
  end select
  di => df%discr_info%begin()
  select type(disc => di%discr_ptr)
  class is (ppm_t_equi_mesh)
    mesh => disc
  end select
  
  IF (ndim.EQ.2) THEN
  foreach n in equi_mesh(mesh) with sca_fields(df) indices(i,j)
      for real
        df_n = exp(t)
  end foreach
  ELSE
  foreach n in equi_mesh(mesh) with sca_fields(df) indices(i,j,k)
      for real
        df_n = exp(t)
  end foreach
  ENDIF

  exp_mesh_ode = 0 
  end_subroutine()
end function exp_mesh_ode
    

end test_suite
