test_suite ppm_module_ode_typedef

use ppm_module_core
use ppm_module_user_numerics

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
real(mk),parameter              :: tol=epsilon(1._mk)*100
real(mk),parameter              :: pi = ACOS(-1._mk)
real(mk),parameter              :: skin = 0._mk
integer,parameter               :: ndim=2
integer                         :: decomp,assig,tolexp
integer                         :: info,comm,rank,nproc,topoid
integer                         :: np_global = 3000
real(mk)                        :: cutoff
real(mk),dimension(:,:),pointer :: xp=>NULL()
real(mk),dimension(:  ),pointer :: min_phys=>NULL(),max_phys=>NULL()
real(mk),dimension(:  ),pointer :: len_phys=>NULL()
integer                         :: i,j,k,ip,wp_id
integer                         :: nstep
integer,dimension(3)            :: ldc
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost=>NULL()
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

    init

        use ppm_module_init
        use ppm_module_mktopo
        
        allocate(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        max_phys(ndim) = 1.4_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic
        
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
            cutoff = 0.15_mk
        else
            cutoff = 0.25_mk
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

    test PSE_client
        use ppm_module_io_vtk
        type(ppm_t_particles_d)         :: Part1
        type(ppm_t_field)               :: Field1
        type(ppm_t_field)               :: Field2
        type(ppm_t_operator)            :: Laplacian
        CLASS(ppm_t_operator_discr),POINTER   :: DCop => NULL()
        CLASS(ppm_t_operator_discr),POINTER   :: PSEop => NULL()
        class(ppm_t_neighlist_d_),POINTER :: Nlist => NULL()
        class(ppm_t_discr_data),POINTER :: prop => NULL()

        !--------------------------
        !Define Fields
        !--------------------------
        call Field1%create(5,info,name="Concentration") !vector field
        Assert_Equal(info,0)

        call Part1%initialize(np_global,info,topoid=topoid,name="Part1")
        Assert_Equal(info,0)

!  print particles to a VTK file
!        CALL ppm_vtk_particles("part_test",Part1,info)
!        Assert_Equal(info,0)

        call Part1%set_cutoff(3._mk * Part1%h_avg,info)
        Assert_Equal(info,0)

        allocate(wp_2r(ndim,Part1%Npart))
        call random_number(wp_2r)
        wp_2r = (wp_2r - 0.5_mk) * Part1%h_avg * 0.15_mk
        call Part1%move(wp_2r,info)
        Assert_Equal(info,0)
        deallocate(wp_2r)

        Assert_true(Part1%has_neighlist(Part1))
        call Part1%apply_bc(info)
        Assert_Equal(info,0)

        call Part1%map(info,global=.true.,topoid=topoid)
        Assert_Equal(info,0)

        call Part1%map_ghosts(info)
        Assert_Equal(info,0)

        call Field1%discretize_on(Part1,info)
        Assert_Equal(info,0)


        call Part1%set_xp(xp,info,read_only=.true.)
        Assert_Equal(info,0)

        call Part1%map_ghosts(info)
        Assert_Equal(info,0)





        call Part1%destroy(info)
        Assert_Equal(info,0)
        call Field1%destroy(info)
        Assert_Equal(info,0)
        deallocate(degree,coeffs,order)
    end test
!-------------------------------------------------------------
! test function
!-------------------------------------------------------------
pure function f0_test(pos,ndim)

    real(mk)                              :: f0_test
    integer                 ,  intent(in) :: ndim
    real(mk), dimension(ndim), intent(in) :: pos

    f0_test =  sin(2._mk*pi*pos(1)) * cos(2._mk*pi*pos(2)) * &
        & sin(2._mk*pi*pos(ndim))

    return

end function f0_test

    

end test_suite
