
module ppm_module_integrator_typedef

!----------------------------------------------------------------------
!  Modules
!----------------------------------------------------------------------
use ppm_module_core
use ppm_module_numerics_interfaces

implicit none

!-----------------------------------------------------------------------
! Internal Parameters
!-----------------------------------------------------------------------

! TODO: replace this temporary hack by ppmcg templates
  integer, parameter :: mk = ppm_kind_double

  private mk
!----------------------------------------------------------------------
!  Types
!----------------------------------------------------------------------

type,abstract,extends(ppm_t_integrator_) :: ppm_t_integrator
  contains
  procedure :: create => integrator_create
  procedure :: destroy => integrator_destroy
end type ppm_t_integrator

type,extends(ppm_t_integrator) :: ppm_t_eulerf
  contains
  procedure :: create => eulerf_create
  procedure :: step => eulerf_step
end type ppm_t_eulerf

type,extends(ppm_t_integrator) :: ppm_t_sts
  real(mk), dimension(20)        :: stsnu  = 0.0_mk
  integer                        :: nsts   = 1
  contains
  procedure :: create => sts_create
  procedure :: step => sts_step
end type ppm_t_sts

type,extends(ppm_t_options) :: ppm_t_sts_options
  integer       :: nsts = 1
  real(mk)      :: nu   = 0.0_mk
end type ppm_t_sts_options

type,extends(ppm_t_integrator) :: ppm_t_tvdrk2
  contains
  procedure :: create => tvdrk2_create
  procedure :: step => tvdrk2_step
end type ppm_t_tvdrk2

type,extends(ppm_t_integrator) :: ppm_t_rk4
  contains
  procedure :: create => rk4_create
  procedure :: step => rk4_step
end type ppm_t_rk4

type,extends(ppm_t_integrator) :: ppm_t_midrk2
  contains
  procedure :: create => midrk2_create
  procedure :: step => midrk2_step
end type ppm_t_midrk2

!----------------------------------------------------------------------
!  Type-bound procedures
!----------------------------------------------------------------------
contains


subroutine integrator_create(this,fields,rhsfunc,rhs_fields_discr,info,options)
  implicit none
  !----------------------------------------------------------------------
  !  Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_integrator)                :: this
  class(ppm_v_main_abstr)                :: fields
  !!! This vector holds all entities to be updated by the integrator the
  !!! elements may be fields or particle discretizations.
  !!!
  !!! fields must either have one discretization or wrapped together with the
  !!! intended discretization in a ppm_t_pair object
  procedure(ppm_p_rhsfunc)               :: rhsfunc
  !!! The right hand side function to be executed by the eulerf::step function
  class(ppm_v_field_discr_pair)          :: rhs_fields_discr
  !!! The fields to be passed to the right hand side function.
  !!!
  !!! The elements of the container can be fields, or pairs continaing a field
  !!! and its intended discretization.
  integer,                 intent(  out) :: info
  class(ppm_t_options),target,optional,intent(in   ) :: options
  !----------------------------------------------------------------------
  !  Variables
  !----------------------------------------------------------------------
  class(ppm_t_main_abstr), pointer       :: el => null()
  class(ppm_t_main_abstr), pointer       :: temp => null()
  class(ppm_t_field_), pointer           :: c => null()
  class(ppm_t_field_), pointer           :: buf => null()
  class(ppm_t_field_), pointer           :: el_f => null()
  class(ppm_t_discr_info_), pointer      :: di => null()
  class(ppm_t_discr_kind), pointer       :: d => null()
  class(ppm_t_field_discr_pair), pointer :: el_p => null()
  integer                                :: ibuf
  integer                                :: ifield
  character(len=16)                      :: bname
  character(len=16)                      :: cname

  start_subroutine("integrator_create")

  allocate(this%fields,stat=info)
  or_fail_alloc("this%fields")
  allocate(this%buffers,stat=info)
  or_fail_alloc("this%buffers")
  allocate(this%changes,stat=info)
  or_fail_alloc("this%changes")
  allocate(this%discretizations,stat=info)
  or_fail_alloc("this%discretizations")
  el => fields%begin()
  ifield = 0
  do while (associated(el))
    ifield = ifield + 1
    allocate(ppm_t_field::c,stat=info)
    allocate(ppm_t_field::buf,stat=info)
    write(cname,'(A,I0)') 'ode_change_',ifield
    write(bname,'(A,I0)') 'ode_buffer_',ifield
    select type(el)
    class is (ppm_t_field_)
      el_f => el
      di => el_f%discr_info%begin()
      call c%create(el_f%lda,info,name=cname)
      call c%discretize_on(di%discr_ptr,info)
      or_fail("Discretizing change failed")
      call this%discretizations%push(di%discr_ptr,info)
      or_fail("Pushing change discretization failed")
      temp => el
      call this%fields%push(temp,info)
      or_fail("Pushing field failed")
      if (el_f%lda*this%scheme_memsize.ne.0) then
        call buf%create(el_f%lda*this%scheme_memsize,info,name=bname)
        call buf%discretize_on(di%discr_ptr,info,with_ghosts=.false.)
      endif
    class is (ppm_t_discr_kind)
      d => el
      call c%create(ppm_dim,info,name=cname)
      call c%discretize_on(d,info)
      or_fail("Discretizing change failed")
      call this%discretizations%push(d,info)
      or_fail("Pushing change discretization failed")
      temp => d
      call this%fields%push(temp,info)
      or_fail("Pushing field failed")
      if (ppm_dim*this%scheme_memsize.ne.0) then
        call buf%create(ppm_dim*this%scheme_memsize,info,name=bname)
        call buf%discretize_on(d,info,with_ghosts=.false.)
      endif
    class is (ppm_t_field_discr_pair)
      el_p => el
      call c%create(el_p%field%lda,info,name=cname)
      call c%discretize_on(el_p%discretization,info)
      or_fail("Discretizing change failed")
      call this%discretizations%push(el_p%discretization,info)
      or_fail("Pushing change discretization failed")
      temp => el_p%field
      call this%fields%push(temp,info)
      or_fail("Pushing field failed")
      if (ppm_dim*this%scheme_memsize.ne.0) then
        call buf%create(el_p%field%lda*this%scheme_memsize,info,name=bname)
        call buf%discretize_on(el_p%discretization,info,with_ghosts=.false.)
      endif
    class default
      fail("fields should only contain fields discrs and pairs",ppm_err_argument)
  end select
  call this%changes%push(c,info)
  call this%buffers%push(buf,info)
  el => fields%next()
  end do

  ! copy the rhs_fields_discr vector
  allocate(this%rhs_fields_discr,stat=info)
  or_fail_alloc("this%rhs_fields_discr")
  el_p => rhs_fields_discr%begin()
  do while (associated(el_p))
    call this%rhs_fields_discr%push(el_p,info)
    el_p => rhs_fields_discr%next()
  end do
 
  ! set the right hand side function
  this%rhsfunc => rhsfunc
 
  end_subroutine()
end subroutine integrator_create

subroutine integrator_destroy(this,info)
  implicit none
  !----------------------------------------------------------------------
  ! Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_integrator)      :: this
  integer    , intent(  out)   :: info
  !----------------------------------------------------------------------
  ! Variables
  !----------------------------------------------------------------------
  class(ppm_t_field_), pointer               :: change => null()
  class(ppm_t_field_), pointer               :: buffer => null()
  start_subroutine("eulerf_destroy")
  
  change => this%changes%begin()
  buffer => this%buffers%begin()
  do while (associated(change))
    check_associated(buffer,"Buffers and changes should have same length")
    call change%destroy(info) 
     or_fail("Destroying change failed")
    change => this%changes%next()
    call buffer%destroy(info) 
     or_fail("Destroying buffer failed")
    buffer => this%buffers%next()
  end do
  deallocate(this%changes,this%buffers,stat=info)
 
  end_subroutine()
end subroutine integrator_destroy

subroutine eulerf_create(this,fields,rhsfunc,rhs_fields_discr,info,options)
  implicit none
  !----------------------------------------------------------------------
  !  Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_eulerf)                    :: this
  class(ppm_v_main_abstr)                :: fields
  !!! This vector holds all entities to be updated by the integrator the
  !!! elements may be fields or particle discretizations.
  !!!
  !!! fields must either have one discretization or wrapped together with the
  !!! intended discretization in a ppm_t_pair object
  procedure(ppm_p_rhsfunc)               :: rhsfunc
  !!! The right hand side function to be executed by the eulerf::step function
  class(ppm_v_field_discr_pair)          :: rhs_fields_discr
  !!! The fields to be passed to the right hand side function.
  !!!
  !!! The elements of the container can be fields, or pairs continaing a field
  !!! and its intended discretization.
  integer,                 intent(  out) :: info
  class(ppm_t_options),target,optional,intent(in   ) :: options
  
  start_subroutine("eulerf_create")
  
  this%scheme_order   = 1
  this%scheme_memsize = 0
  this%scheme_nstages = 1
  this%scheme_kickoff = ppm_param_ode_scheme_eulerf
  
  
  call integrator_create(this,fields,rhsfunc,rhs_fields_discr,info)


  end_subroutine()

end subroutine eulerf_create

subroutine eulerf_step(this,t,dt,istage,info)
  implicit none
  ! FIXME: must be templated
  integer, parameter  :: mk = ppm_kind_double
  !----------------------------------------------------------------------
  ! Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_eulerf)                        :: this
  real(ppm_kind_double)            ,intent(inout)    :: t
  real(ppm_kind_double)            ,intent(in   )    :: dt
  integer                          ,intent(in   )    :: istage
  integer                          ,intent(  out)    :: info
  !----------------------------------------------------------------------
  ! Variables
  !----------------------------------------------------------------------
  integer                                    :: rhs_info
  class(ppm_t_main_abstr), pointer           :: field => null()
  class(ppm_t_field_), pointer               :: change => null()
  class(ppm_t_discr_kind), pointer           :: discr  => null()
  class(ppm_t_particles_d), pointer          :: pset_d => null()
  class(ppm_t_particles_s), pointer          :: pset_s => null()
  class(ppm_t_equi_mesh), pointer            :: mesh => null()
  start_subroutine("eulerf_step")

  rhs_info = this%rhsfunc(this%rhs_fields_discr,t,this%changes)

  field => this%fields%begin()
  change => this%changes%begin()
  discr => this%discretizations%begin()
  do while (associated(field))
    check_associated(change,"Fields and changes should have same length")
    check_associated(discr,"Fields and discretizations should have same length")

    select type (field)
    class is (ppm_t_field)
      select type (discr)
      class is (ppm_t_particles_s)
        !pset_s => discr
        !foreach p in particles(pset_s) with fields(f=field,df=change) types(f=scalar,df=scalar)
        !  f_p = f_p + dt*df_p
        !end foreach
      class is (ppm_t_particles_d)
        pset_d => discr
        if (field%lda.eq.1) then
          foreach p in particles(pset_d) with sca_fields(f=field,df=change)
            f_p = f_p + dt*df_p
          end foreach
        else
          foreach p in particles(pset_d) with vec_fields(f=field,df=change)
            f_p(:) = f_p(:) + dt*df_p(:)
          end foreach
        endif
      class is (ppm_t_equi_mesh)
        mesh => discr
        if (ppm_dim.eq.2) then
          if (field%lda.eq.1) then
            foreach n in equi_mesh(mesh) with sca_fields(field,change) indices(i,j)
              for real
                field_n = field_n + dt*change_n
            end foreach
          else
            foreach n in equi_mesh(mesh) with vec_fields(field,change) indices(i,j)
              for real
                field_n(:) = field_n(:) + dt*change_n(:)
            end foreach
          endif
        else
          if (field%lda.eq.1) then
            foreach n in equi_mesh(mesh) with sca_fields(field,change) indices(i,j,k)
              for real
                field_n = field_n + dt*change_n
            end foreach
          else
            foreach n in equi_mesh(mesh) with vec_fields(field,change) indices(i,j,k)
              for real
                field_n(:) = field_n(:) + dt*change_n(:)
            end foreach
          endif
        endif
      end select
    class is (ppm_t_particles_s)
      !pset_s => field
      !foreach p in particles(pset_s) with positions(x) fields(dx=change) types(dx=vector)
      !  x_p(:) = x_p(:) + dt*dx_p(:)
      !end foreach
    class is (ppm_t_particles_d)
      pset_d => field
      foreach p in particles(pset_d) with positions(x) vec_fields(dx=change)
        x_p(:) = x_p(:) + dt*dx_p(:)
      end foreach
    end select
    field => this%fields%next()
    change => this%changes%next()
    discr => this%discretizations%next()
  end do
  check_false("associated(change)","Fields and changes should have same length")
  check_false("associated(discr)","Fields and discretizations should have same length")

  t  = t  + dt
  ! how much to save of this
  ! TODO FIXME
  ! ppm_ode_sent(mid) = 0

  end_subroutine()
end subroutine eulerf_step


subroutine sts_create(this,fields,rhsfunc,rhs_fields_discr,info,options)
  implicit none
  !----------------------------------------------------------------------
  !  Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_sts)                       :: this
  class(ppm_v_main_abstr)                :: fields
  !!! This vector holds all entities to be updated by the integrator the
  !!! elements may be fields or particle discretizations.
  !!!
  !!! fields must either have one discretization or wrapped together with the
  !!! intended discretization in a ppm_t_pair object
  procedure(ppm_p_rhsfunc)               :: rhsfunc
  !!! The right hand side function to be executed by the eulerf::step function
  class(ppm_v_field_discr_pair)          :: rhs_fields_discr
  !!! The fields to be passed to the right hand side function.
  !!!
  !!! The elements of the container can be fields, or pairs continaing a field
  !!! and its intended discretization.
  integer,                 intent(  out) :: info
  class(ppm_t_options),target,optional,intent(in   ) :: options
  !----------------------------------------------------------------------
  ! Variables
  !----------------------------------------------------------------------
  class(ppm_t_sts_options), pointer           :: sts_options
  start_subroutine("sts_create")
  
  this%scheme_order   = 1
  this%scheme_memsize = 0
  this%scheme_nstages = 1000
  this%scheme_kickoff = ppm_param_ode_scheme_sts
 
  this%stsnu(:)  = 0.0_mk
  this%stsnu(1)  = 0.0_mk
  this%stsnu(5)  = 0.04_mk
  this%stsnu(7)  = 0.0015_mk
  this%stsnu(9)  = 0.04_mk
  this%stsnu(20) = 0.006_mk
 
  select type(options)
  class is (ppm_t_sts_options)
  sts_options => options
  this%nsts      = sts_options%nsts
  this%stsnu(this%nsts) = sts_options%nu
  end select
  
  call integrator_create(this,fields,rhsfunc,rhs_fields_discr,info)


  end_subroutine()

end subroutine sts_create


subroutine sts_step(this,t,dt,istage,info)
  implicit none
  ! FIXME: must be templated
  integer, parameter  :: mk = ppm_kind_double
  !----------------------------------------------------------------------
  ! Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_sts)                                   :: this
  real(ppm_kind_double)            ,intent(inout)    :: t
  real(ppm_kind_double)            ,intent(in   )    :: dt
  integer                          ,intent(in   )    :: istage
  integer                          ,intent(  out)    :: info
  !----------------------------------------------------------------------
  ! Variables
  !----------------------------------------------------------------------
  integer                                    :: rhs_info
  class(ppm_t_main_abstr), pointer           :: field => null()
  class(ppm_t_field_), pointer               :: change => null()
  class(ppm_t_discr_kind), pointer           :: discr  => null()
  class(ppm_t_particles_d), pointer          :: pset_d => null()
  class(ppm_t_particles_s), pointer          :: pset_s => null()
  class(ppm_t_equi_mesh), pointer            :: mesh => null()
  real(mk)                                   :: tau
  real(mk)                                   :: m_pi
  start_subroutine("sts_step")
           
  m_pi = ACOS(-1.0_MK)
  !-----------------------------------------------------
  !  compute the new dt
  !-----------------------------------------------------
  tau = dt/((this%stsnu(this%nsts) - 1.0_mk) * &
  &     cos((2.0_mk * real(istage,mk) - 1.0_mk)/real(this%nsts,mk) * M_PI * 0.5_mk) + &
  &     1.0_mk + this%stsnu(this%nsts))

  rhs_info = this%rhsfunc(this%rhs_fields_discr,t,this%changes)

  field => this%fields%begin()
  change => this%changes%begin()
  discr => this%discretizations%begin()
  do while (associated(field))
    check_associated(change,"Fields and changes should have same length")
    check_associated(discr,"Fields and discretizations should have same length")

    select type (field)
    class is (ppm_t_field)
      select type (discr)
      class is (ppm_t_particles_s)
        !pset_s => discr
        !foreach p in particles(pset_s) with fields(f=field,df=change) types(f=scalar,df=scalar)
        !  f_p = f_p + tau*df_p
        !end foreach
      class is (ppm_t_particles_d)
        pset_d => discr
        if (field%lda.eq.1) then
          foreach p in particles(pset_d) with sca_fields(f=field,df=change)
            f_p = f_p + tau*df_p
          end foreach
        else
          foreach p in particles(pset_d) with vec_fields(f=field,df=change)
            f_p(:) = f_p(:) + tau*df_p(:)
          end foreach
        endif
      class is (ppm_t_equi_mesh)
        mesh => discr
        if (ppm_dim.eq.2) then
          if (field%lda.eq.1) then
            foreach n in equi_mesh(mesh) with sca_fields(field,change) indices(i,j)
              for real
                field_n = field_n + tau*change_n
            end foreach
          else
            foreach n in equi_mesh(mesh) with vec_fields(field,change) indices(i,j)
              for real
                field_n(:) = field_n(:) + tau*change_n(:)
            end foreach
          endif
        else
          if (field%lda.eq.1) then
            foreach n in equi_mesh(mesh) with sca_fields(field,change) indices(i,j,k)
              for real
                field_n = field_n + tau*change_n
            end foreach
          else
            foreach n in equi_mesh(mesh) with vec_fields(field,change) indices(i,j,k)
              for real
                field_n(:) = field_n(:) + tau*change_n(:)
            end foreach
          endif
        endif
      end select
    class is (ppm_t_particles_s)
      !pset_s => field
      !foreach p in particles(pset_s) with positions(x) fields(dx=change) types(dx=vector)
      !  x_p(:) = x_p(:) + tau*dx_p(:)
      !end foreach
    class is (ppm_t_particles_d)
      pset_d => field
      foreach p in particles(pset_d) with positions(x) vec_fields(dx=change)
        x_p(:) = x_p(:) + tau*dx_p(:)
      end foreach
    end select
    field => this%fields%next()
    change => this%changes%next()
    discr => this%discretizations%next()
  end do
  check_false("associated(change)","Fields and changes should have same length")
  check_false("associated(discr)","Fields and discretizations should have same length")

  t  = t  + tau

  end_subroutine()
end subroutine sts_step

subroutine tvdrk2_create(this,fields,rhsfunc,rhs_fields_discr,info,options)
  implicit none
  !----------------------------------------------------------------------
  !  Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_tvdrk2)                    :: this
  class(ppm_v_main_abstr)                :: fields
  !!! This vector holds all entities to be updated by the integrator the
  !!! elements may be fields or particle discretizations.
  !!!
  !!! fields must either have one discretization or wrapped together with the
  !!! intended discretization in a ppm_t_pair object
  procedure(ppm_p_rhsfunc)               :: rhsfunc
  !!! The right hand side function to be executed by the eulerf::step function
  class(ppm_v_field_discr_pair)          :: rhs_fields_discr
  !!! The fields to be passed to the right hand side function.
  !!!
  !!! The elements of the container can be fields, or pairs continaing a field
  !!! and its intended discretization.
  integer,                 intent(  out) :: info
  class(ppm_t_options),target,optional,intent(in   ) :: options
  start_subroutine("tvdrk2_create")
  
  this%scheme_order   = 2
  this%scheme_memsize = 1
  this%scheme_nstages = 2
  this%scheme_kickoff = ppm_param_ode_scheme_tvdrk2
 
  ! TODO allocate buffer   
  call integrator_create(this,fields,rhsfunc,rhs_fields_discr,info)


  end_subroutine()

end subroutine tvdrk2_create


subroutine tvdrk2_step(this,t,dt,istage,info)
  implicit none
  ! FIXME: must be templated
  integer, parameter  :: mk = ppm_kind_double
  !----------------------------------------------------------------------
  ! Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_tvdrk2)                                :: this
  real(ppm_kind_double)            ,intent(inout)    :: t
  real(ppm_kind_double)            ,intent(in   )    :: dt
  integer                          ,intent(in   )    :: istage
  integer                          ,intent(  out)    :: info
  !----------------------------------------------------------------------
  ! Variables
  !----------------------------------------------------------------------
  integer                                    :: rhs_info
  class(ppm_t_main_abstr), pointer           :: field => null()
  class(ppm_t_field_), pointer               :: change => null()
  class(ppm_t_field_), pointer               :: buffer => null()
  class(ppm_t_discr_kind), pointer           :: discr  => null()
  class(ppm_t_particles_d), pointer          :: pset_d => null()
  class(ppm_t_particles_s), pointer          :: pset_s => null()
  class(ppm_t_equi_mesh), pointer            :: mesh => null()
  start_subroutine("tvdrk2_step")
           


  ! TODO check if buffers where mapped
  select case (istage)
  case (1)
    !----------------------------------------------------------------------
    ! STAGE 1: Store the current field in a secondary buffer and update it
    ! with dt*df
    !----------------------------------------------------------------------
    rhs_info = this%rhsfunc(this%rhs_fields_discr,t,this%changes)

    field => this%fields%begin()
    change => this%changes%begin()
    discr => this%discretizations%begin()
    buffer => this%buffers%begin()
    do while (associated(field))
      check_associated(change,"Fields and changes should have same length")
      check_associated(discr,"Fields and discretizations should have same length")
      check_associated(buffer,"Fields and buffers should have same length")

      select type (field)
      class is (ppm_t_field)
        select type (discr)
        class is (ppm_t_particles_s)
          !pset_s => discr
          !foreach p in particles(pset_s) with fields(f=field,df=change) types(f=scalar,df=scalar)
          !  f_p = f_p + tau*df_p
          !end foreach
        class is (ppm_t_particles_d)
          pset_d => discr
          if (field%lda.eq.1) then
            foreach p in particles(pset_d) with sca_fields(f=field,df=change,bfr=buffer)
              bfr_p = f_p
              f_p = f_p + dt*df_p
            end foreach
          else
            foreach p in particles(pset_d) with vec_fields(f=field,df=change,bfr=buffer)
              bfr_p(:) = f_p(:)
              f_p(:) = f_p(:) + dt*df_p(:)
            end foreach
          endif
        class is (ppm_t_equi_mesh)
          mesh => discr
          if (ppm_dim.eq.2) then
            if (field%lda.eq.1) then
              foreach n in equi_mesh(mesh) with sca_fields(field,change,buffer) indices(i,j)
                for real
                  buffer_n = field_n
                  field_n = field_n + dt*change_n
              end foreach
            else
              foreach n in equi_mesh(mesh) with vec_fields(field,change,buffer) indices(i,j)
                for real
                  buffer_n(:) = field_n(:)
                  field_n(:) = field_n(:) + dt*change_n(:)
              end foreach
            endif
          else
            if (field%lda.eq.1) then
              foreach n in equi_mesh(mesh) with sca_fields(field,change,buffer) indices(i,j,k)
                for real
                  buffer_n = field_n
                  field_n = field_n + dt*change_n
              end foreach
            else
              foreach n in equi_mesh(mesh) with vec_fields(field,change,buffer) indices(i,j,k)
                for real
                  buffer_n(:) = field_n(:)
                  field_n(:) = field_n(:) + dt*change_n(:)
              end foreach
            endif
          endif
        end select
      class is (ppm_t_particles_s)
        !pset_s => field
        !foreach p in particles(pset_s) with positions(x) fields(dx=change) types(dx=vector)
        !  buffer_p(:) = x_p(:)
        !  x_p(:) = x_p(:) + dt*dx_p(:)
        !end foreach
      class is (ppm_t_particles_d)
        pset_d => field
        foreach p in particles(pset_d) with positions(x) vec_fields(dx=change,bfr=buffer)
          bfr_p(:) = x_p(:)
          x_p(:) = x_p(:) + dt*dx_p(:)
        end foreach
      end select
      field => this%fields%next()
      change => this%changes%next()
      buffer => this%buffers%next()
      discr => this%discretizations%next()
    end do
    check_false("associated(change)","Fields and changes should have same length")
    check_false("associated(discr)","Fields and discretizations should have same length")
    check_false("associated(buffer)","Fields and buffers should have same length")
  
  case (2)
    !----------------------------------------------------------------------
    ! STAGE 2: Do another update and interpolate with previously saved data
    !----------------------------------------------------------------------
    rhs_info = this%rhsfunc(this%rhs_fields_discr,t+dt,this%changes)
    
    field => this%fields%begin()
    change => this%changes%begin()
    buffer => this%buffers%begin()
    discr => this%discretizations%begin()
    do while (associated(field))
      check_associated(change,"Fields and changes should have same length")
      check_associated(discr,"Fields and discretizations should have same length")
      check_associated(buffer,"Fields and buffers should have same length")

      select type (field)
      class is (ppm_t_field)
        select type (discr)
        class is (ppm_t_particles_s)
          !pset_s => discr
          !foreach p in particles(pset_s) with fields(f=field,df=change) types(f=scalar,df=scalar)
          !  f_p = f_p + dt*df_p
          !  f_p  = 0.5_mk * (f_p + buffer_p)
          !end foreach
        class is (ppm_t_particles_d)
          pset_d => discr
          if (field%lda.eq.1) then
            foreach p in particles(pset_d) with sca_fields(f=field,df=change,bfr=buffer)
              f_p = f_p + dt*df_p
              f_p  = 0.5_mk * (f_p + bfr_p)
            end foreach
          else
            foreach p in particles(pset_d) with vec_fields(f=field,df=change,bfr=buffer)
              f_p(:) = f_p(:) + dt*df_p(:)
              f_p(:)  = 0.5_mk * (f_p(:) + bfr_p(:))
            end foreach
          endif
        class is (ppm_t_equi_mesh)
          mesh => discr
          if (ppm_dim.eq.2) then
            if (field%lda.eq.1) then
              foreach n in equi_mesh(mesh) with sca_fields(field,change,buffer) indices(i,j)
                for real
                  field_n = field_n + dt*change_n
                  field_n  = 0.5_mk * (field_n + buffer_n)
              end foreach
            else
              foreach n in equi_mesh(mesh) with vec_fields(field,change,buffer) indices(i,j)
                for real
                  field_n(:) = field_n(:) + dt*change_n(:)
                  field_n(:)  = 0.5_mk * (field_n(:) + buffer_n(:))
              end foreach
            endif
          else
            if (field%lda.eq.1) then
              foreach n in equi_mesh(mesh) with sca_fields(field,change,buffer) indices(i,j,k)
                for real
                  field_n = field_n + dt*change_n
                  field_n  = 0.5_mk * (field_n + buffer_n)
              end foreach
            else
              foreach n in equi_mesh(mesh) with vec_fields(field,change,buffer) indices(i,j,k)
                for real
                  field_n(:) = field_n(:) + dt*change_n(:)
                  field_n(:)  = 0.5_mk * (field_n(:) + buffer_n(:))
              end foreach
            endif
          endif
        end select
      class is (ppm_t_particles_s)
        !pset_s => field
        !foreach p in particles(pset_s) with positions(x) fields(dx=change) types(dx=vector)
        !  x_p(:) = x_p(:) + dt*dx_p(:)
        !  x_p(:)  = 0.5_mk * (x_p(:) + buffer_p(:))
        !end foreach
      class is (ppm_t_particles_d)
        pset_d => field
        foreach p in particles(pset_d) with positions(x) vec_fields(dx=change,bfr=buffer)
          x_p(:) = x_p(:) + dt*dx_p(:)
          x_p(:)  = 0.5_mk * (x_p(:) + bfr_p(:))
        end foreach
      end select
      field => this%fields%next()
      change => this%changes%next()
      buffer => this%buffers%next()
      discr => this%discretizations%next()
    end do
    check_false("associated(change)","Fields and changes should have same length")
    check_false("associated(discr)","Fields and discretizations should have same length")
    check_false("associated(buffer)","Fields and buffers should have same length")
  
    t  = t  + dt
  end select


  end_subroutine()
end subroutine tvdrk2_step


subroutine midrk2_create(this,fields,rhsfunc,rhs_fields_discr,info,options)
  implicit none
  !----------------------------------------------------------------------
  !  Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_midrk2)                    :: this
  class(ppm_v_main_abstr)                :: fields
  !!! This vector holds all entities to be updated by the integrator the
  !!! elements may be fields or particle discretizations.
  !!!
  !!! fields must either have one discretization or wrapped together with the
  !!! intended discretization in a ppm_t_pair object
  procedure(ppm_p_rhsfunc)               :: rhsfunc
  !!! The right hand side function to be executed by the eulerf::step function
  class(ppm_v_field_discr_pair)          :: rhs_fields_discr
  !!! The fields to be passed to the right hand side function.
  !!!
  !!! The elements of the container can be fields, or pairs continaing a field
  !!! and its intended discretization.
  integer,                 intent(  out) :: info
  class(ppm_t_options),target,optional,intent(in   ) :: options
  start_subroutine("midrk2_create")
  
  this%scheme_order   = 2
  this%scheme_memsize = 1
  this%scheme_nstages = 2
  this%scheme_kickoff = ppm_param_ode_scheme_midrk2
 
  ! TODO allocate buffer   
  call integrator_create(this,fields,rhsfunc,rhs_fields_discr,info)


  end_subroutine()

end subroutine midrk2_create


subroutine midrk2_step(this,t,dt,istage,info)
  implicit none
  ! FIXME: must be templated
  integer, parameter  :: mk = ppm_kind_double
  !----------------------------------------------------------------------
  ! Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_midrk2)                                :: this
  real(ppm_kind_double)            ,intent(inout)    :: t
  real(ppm_kind_double)            ,intent(in   )    :: dt
  integer                          ,intent(in   )    :: istage
  integer                          ,intent(  out)    :: info
  !----------------------------------------------------------------------
  ! Variables
  !----------------------------------------------------------------------
  integer                                    :: rhs_info
  class(ppm_t_main_abstr), pointer           :: field => null()
  class(ppm_t_field_), pointer               :: change => null()
  class(ppm_t_field_), pointer               :: buffer => null()
  class(ppm_t_discr_kind), pointer           :: discr  => null()
  class(ppm_t_particles_d), pointer          :: pset_d => null()
  class(ppm_t_particles_s), pointer          :: pset_s => null()
  class(ppm_t_equi_mesh), pointer            :: mesh => null()
  start_subroutine("midrk2_step")
           


  ! TODO check if buffers where mapped
  select case (istage)
  case (1)
    !----------------------------------------------------------------------
    ! STAGE 1: Compute midpoint
    ! Store the current field in a secondary buffer and update it
    ! with 0.5*dt*df
    !----------------------------------------------------------------------
    rhs_info = this%rhsfunc(this%rhs_fields_discr,t,this%changes)

    field => this%fields%begin()
    change => this%changes%begin()
    discr => this%discretizations%begin()
    buffer => this%buffers%begin()
    do while (associated(field))
      check_associated(change,"Fields and changes should have same length")
      check_associated(discr,"Fields and discretizations should have same length")
      check_associated(buffer,"Fields and buffers should have same length")

      select type (field)
      class is (ppm_t_field)
        select type (discr)
        class is (ppm_t_particles_s)
          !pset_s => discr
          !foreach p in particles(pset_s) with fields(f=field,df=change) types(f=scalar,df=scalar)
          !  f_p = f_p + tau*df_p
          !end foreach
        class is (ppm_t_particles_d)
          pset_d => discr
          if (field%lda.eq.1) then
            foreach p in particles(pset_d) with sca_fields(f=field,df=change,bfr=buffer)
              bfr_p = f_p
              f_p = f_p + 0.5_mk*dt*df_p
            end foreach
          else
            foreach p in particles(pset_d) with vec_fields(f=field,df=change,bfr=buffer)
              bfr_p(:) = f_p(:)
              f_p(:) = f_p(:) + 0.5_mk*dt*df_p(:)
            end foreach
          endif
        class is (ppm_t_equi_mesh)
          mesh => discr
          if (ppm_dim.eq.2) then
            if (field%lda.eq.1) then
              foreach n in equi_mesh(mesh) with sca_fields(field,change,buffer) indices(i,j)
                for real
                  buffer_n = field_n
                  field_n = field_n + 0.5_mk*dt*change_n
              end foreach
            else
              foreach n in equi_mesh(mesh) with vec_fields(field,change,buffer) indices(i,j)
                for real
                  buffer_n(:) = field_n(:)
                  field_n(:) = field_n(:) + 0.5_mk*dt*change_n(:)
              end foreach
            endif
          else
            if (field%lda.eq.1) then
              foreach n in equi_mesh(mesh) with sca_fields(field,change,buffer) indices(i,j,k)
                for real
                  buffer_n = field_n
                  field_n = field_n + 0.5_mk*dt*change_n
              end foreach
            else
              foreach n in equi_mesh(mesh) with vec_fields(field,change,buffer) indices(i,j,k)
                for real
                  buffer_n(:) = field_n(:)
                  field_n(:) = field_n(:) + 0.5_mk*dt*change_n(:)
              end foreach
            endif
          endif
        end select
      class is (ppm_t_particles_s)
        !pset_s => field
        !foreach p in particles(pset_s) with positions(x) fields(dx=change) types(dx=vector)
        !  bfr_p(:) = x_p(:)
        !  x_p(:) = x_p(:) + 0.5_mk*dt*dx_p(:)
        !end foreach
      class is (ppm_t_particles_d)
        pset_d => field
        foreach p in particles(pset_d) with positions(x) vec_fields(dx=change,bfr=buffer)
          bfr_p(:) = x_p(:)
          x_p(:) = x_p(:) + 0.5_mk*dt*dx_p(:)
        end foreach
      end select
      field => this%fields%next()
      change => this%changes%next()
      buffer => this%buffers%next()
      discr => this%discretizations%next()
    end do
    check_false("associated(change)","Fields and changes should have same length")
    check_false("associated(discr)","Fields and discretizations should have same length")
    check_false("associated(buffer)","Fields and buffers should have same length")
  
  case (2)
    !----------------------------------------------------------------------
    ! STAGE 2: Do another update and interpolate with previously saved data
    !----------------------------------------------------------------------
    rhs_info = this%rhsfunc(this%rhs_fields_discr,t+0.5*dt,this%changes)
    
    field => this%fields%begin()
    change => this%changes%begin()
    buffer => this%buffers%begin()
    discr => this%discretizations%begin()
    do while (associated(field))
      check_associated(change,"Fields and changes should have same length")
      check_associated(discr,"Fields and discretizations should have same length")
      check_associated(buffer,"Fields and buffers should have same length")

      select type (field)
      class is (ppm_t_field)
        select type (discr)
        class is (ppm_t_particles_s)
          !pset_s => discr
          !foreach p in particles(pset_s) with fields(f=field,df=change) types(f=scalar,df=scalar)
          !  f_p = bfr_p + dt*df_p
          !end foreach
        class is (ppm_t_particles_d)
          pset_d => discr
          if (field%lda.eq.1) then
            foreach p in particles(pset_d) with sca_fields(f=field,df=change,bfr=buffer)
              f_p = bfr_p + dt*df_p
            end foreach
          else
            foreach p in particles(pset_d) with vec_fields(f=field,df=change,bfr=buffer)
              f_p(:) = bfr_p(:) + dt*df_p(:)
            end foreach
          endif
        class is (ppm_t_equi_mesh)
          mesh => discr
          if (ppm_dim.eq.2) then
            if (field%lda.eq.1) then
              foreach n in equi_mesh(mesh) with sca_fields(field,change,buffer) indices(i,j)
                for real
                  field_n = buffer_n + dt*change_n
              end foreach
            else
              foreach n in equi_mesh(mesh) with vec_fields(field,change,buffer) indices(i,j)
                for real
                  field_n(:) = buffer_n(:) + dt*change_n(:)
              end foreach
            endif
          else
            if (field%lda.eq.1) then
              foreach n in equi_mesh(mesh) with sca_fields(field,change,buffer) indices(i,j,k)
                for real
                  field_n = buffer_n + dt*change_n
              end foreach
            else
              foreach n in equi_mesh(mesh) with vec_fields(field,change,buffer) indices(i,j,k)
                for real
                  field_n(:) = buffer_n(:) + dt*change_n(:)
              end foreach
            endif
          endif
        end select
      class is (ppm_t_particles_s)
        !pset_s => field
        !foreach p in particles(pset_s) with positions(x) fields(dx=change) types(dx=vector)
        !  x_p(:) = bfr_p(:) + dt*dx_p(:)
        !end foreach
      class is (ppm_t_particles_d)
        pset_d => field
        foreach p in particles(pset_d) with positions(x) vec_fields(dx=change,bfr=buffer)
          x_p(:) = bfr_p(:) + dt*dx_p(:)
        end foreach
      end select
      field => this%fields%next()
      change => this%changes%next()
      buffer => this%buffers%next()
      discr => this%discretizations%next()
    end do
    check_false("associated(change)","Fields and changes should have same length")
    check_false("associated(discr)","Fields and discretizations should have same length")
    check_false("associated(buffer)","Fields and buffers should have same length")
  
    t  = t  + dt

  end select


  end_subroutine()
end subroutine midrk2_step


subroutine rk4_create(this,fields,rhsfunc,rhs_fields_discr,info,options)
  implicit none
  !----------------------------------------------------------------------
  !  Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_rk4)                    :: this
  class(ppm_v_main_abstr)                :: fields
  !!! This vector holds all entities to be updated by the integrator the
  !!! elements may be fields or particle discretizations.
  !!!
  !!! fields must either have one discretization or wrapped together with the
  !!! intended discretization in a ppm_t_pair object
  procedure(ppm_p_rhsfunc)               :: rhsfunc
  !!! The right hand side function to be executed by the eulerf::step function
  class(ppm_v_field_discr_pair)          :: rhs_fields_discr
  !!! The fields to be passed to the right hand side function.
  !!!
  !!! The elements of the container can be fields, or pairs continaing a field
  !!! and its intended discretization.
  integer,                 intent(  out) :: info
  class(ppm_t_options),target,optional,intent(in   ) :: options
  start_subroutine("rk4_create")
  
  this%scheme_order   = 4
  this%scheme_memsize = 4
  this%scheme_nstages = 4
  this%scheme_kickoff = ppm_param_ode_scheme_rk4
 
  ! TODO allocate buffer   
  call integrator_create(this,fields,rhsfunc,rhs_fields_discr,info)


  end_subroutine()

end subroutine rk4_create


subroutine rk4_step(this,t,dt,istage,info)
  implicit none
  ! FIXME: must be templated
  integer, parameter  :: mk = ppm_kind_double
  !----------------------------------------------------------------------
  ! Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_rk4)                                :: this
  real(ppm_kind_double)            ,intent(inout)    :: t
  real(ppm_kind_double)            ,intent(in   )    :: dt
  integer                          ,intent(in   )    :: istage
  integer                          ,intent(  out)    :: info
  !----------------------------------------------------------------------
  ! Variables
  !----------------------------------------------------------------------
  integer                                    :: rhs_info
  class(ppm_t_main_abstr), pointer           :: field => null()
  class(ppm_t_field_), pointer               :: change => null()
  class(ppm_t_field_), pointer               :: buffer => null()
  class(ppm_t_discr_kind), pointer           :: discr  => null()
  class(ppm_t_particles_d), pointer          :: pset_d => null()
  class(ppm_t_particles_s), pointer          :: pset_s => null()
  class(ppm_t_equi_mesh), pointer            :: mesh => null()
  integer                                    :: bs
  start_subroutine("rk4_step")
           


  ! TODO check if buffers where mapped
  select case (istage)
  case (1)
    !----------------------------------------------------------------------
    ! STAGE 1: x_n + 1/2 dt k1
    !----------------------------------------------------------------------
    rhs_info = this%rhsfunc(this%rhs_fields_discr,t,this%changes)
    
    field => this%fields%begin()
    change => this%changes%begin()
    discr => this%discretizations%begin()
    buffer => this%buffers%begin()
    do while (associated(field))
      check_associated(change,"Fields and changes should have same length")
      check_associated(discr,"Fields and discretizations should have same length")
      check_associated(buffer,"Fields and buffers should have same length")

      select type (field)
      class is (ppm_t_field)
        bs = field%lda
        select type (discr)
        class is (ppm_t_particles_s)
          !pset_s => discr
          !foreach p in particles(pset_s) with fields(f=field,df=change) types(f=scalar,df=scalar)
          !  f_p = f_p + tau*df_p
          !end foreach
        class is (ppm_t_particles_d)
          pset_d => discr
          if (field%lda.eq.1) then
            foreach p in particles(pset_d) with sca_fields(f=field,df=change) vec_fields(bfr=buffer)
              bfr_p(1) = f_p
              bfr_p(2) = df_p
              f_p = f_p + 0.5_mk*dt*df_p
            end foreach
          else
            foreach p in particles(pset_d) with vec_fields(f=field,df=change,bfr=buffer)
              bfr_p(1:bs) = f_p(:)
              bfr_p(bs+1:2*bs) = df_p(:)
              f_p(:) = f_p(:) + 0.5_mk*dt*df_p(:)
            end foreach
          endif
        class is (ppm_t_equi_mesh)
          mesh => discr
          if (ppm_dim.eq.2) then
            if (field%lda.eq.1) then
              foreach n in equi_mesh(mesh) with sca_fields(field,change) vec_fields(buffer) indices(i,j)
                for real
                  buffer_n(1) = field_n
                  buffer_n(2) = change_n
                  field_n = field_n + 0.5_mk*dt*change_n
              end foreach
            else
              foreach n in equi_mesh(mesh) with vec_fields(field,change,buffer) indices(i,j)
                for real
                  buffer_n(1:bs) = field_n(:)
                  buffer_n(bs+1:2*bs) = change_n(:)
                  field_n(:) = field_n(:) + 0.5_mk*dt*change_n(:)
              end foreach
            endif
          else
            if (field%lda.eq.1) then
              foreach n in equi_mesh(mesh) with sca_fields(field,change) vec_fields(buffer) indices(i,j,k)
                for real
                  buffer_n(1) = field_n
                  buffer_n(2) = change_n
                  field_n = field_n + 0.5_mk*dt*change_n
              end foreach
            else
              foreach n in equi_mesh(mesh) with vec_fields(field,change,buffer) indices(i,j,k)
                for real
                  buffer_n(1:bs) = field_n(:)
                  buffer_n(bs+1:2*bs) = change_n(:)
                  field_n(:) = field_n(:) + 0.5_mk*dt*change_n(:)
              end foreach
            endif
          endif
        end select
      class is (ppm_t_particles_s)
        bs = ppm_dim
        !pset_s => field
        !foreach p in particles(pset_s) with positions(x) fields(dx=change) types(dx=vector)
        !  bfr_p(:) = x_p(:)
        !  x_p(:) = x_p(:) + 0.5_mk*dt*dx_p(:)
        !end foreach
      class is (ppm_t_particles_d)
        bs = ppm_dim
        pset_d => field
        foreach p in particles(pset_d) with positions(x) vec_fields(dx=change,bfr=buffer)
          bfr_p(1:bs) = x_p(:)
          bfr_p(bs+1:2*bs) = dx_p(:)
          x_p(:) = x_p(:) + 0.5_mk*dt*dx_p(:)
        end foreach
      end select
      field => this%fields%next()
      change => this%changes%next()
      buffer => this%buffers%next()
      discr => this%discretizations%next()
    end do
    check_false("associated(change)","Fields and changes should have same length")
    check_false("associated(discr)","Fields and discretizations should have same length")
    check_false("associated(buffer)","Fields and buffers should have same length")
  
  case (2)
    !----------------------------------------------------------------------
    ! STAGE 2: x_n + 1/2 dt k2
    !----------------------------------------------------------------------
    rhs_info = this%rhsfunc(this%rhs_fields_discr,t+0.5_mk*dt,this%changes)
    
    field => this%fields%begin()
    change => this%changes%begin()
    buffer => this%buffers%begin()
    discr => this%discretizations%begin()
    do while (associated(field))
      check_associated(change,"Fields and changes should have same length")
      check_associated(discr,"Fields and discretizations should have same length")
      check_associated(buffer,"Fields and buffers should have same length")

      select type (field)
      class is (ppm_t_field)
        bs = field%lda
        select type (discr)
        class is (ppm_t_particles_s)
          !pset_s => discr
          !foreach p in particles(pset_s) with fields(f=field,df=change) types(f=scalar,df=scalar)
          !   bfr_p(2) = df_p
          !  f_p = bfr_p + dt*df_p
          !end foreach
        class is (ppm_t_particles_d)
          pset_d => discr
          if (field%lda.eq.1) then
            foreach p in particles(pset_d) with sca_fields(f=field,df=change) vec_fields(bfr=buffer)
              bfr_p(3) = df_p
              f_p = bfr_p(1) + 0.5_mk*dt*df_p
            end foreach
          else
            foreach p in particles(pset_d) with vec_fields(f=field,df=change,bfr=buffer)
              bfr_p(2*bs+1:3*bs) = df_p(:)
              f_p(:) = bfr_p(1:bs) + 0.5_mk*dt*df_p(:)
            end foreach
          endif
        class is (ppm_t_equi_mesh)
          mesh => discr
          if (ppm_dim.eq.2) then
            if (field%lda.eq.1) then
              foreach n in equi_mesh(mesh) with sca_fields(field,change) vec_fields(buffer) indices(i,j)
                for real
                  buffer_n(3) = change_n
                  field_n = buffer_n(1) + 0.5_mk*dt*change_n
              end foreach
            else
              foreach n in equi_mesh(mesh) with vec_fields(field,change,buffer) indices(i,j)
                for real
                  buffer_n(2*bs+1:3*bs) = change_n(:)
                  field_n(:) = buffer_n(1:bs) + 0.5_mk*dt*change_n(:)
              end foreach
            endif
          else
            if (field%lda.eq.1) then
              foreach n in equi_mesh(mesh) with sca_fields(field,change) vec_fields(buffer) indices(i,j,k)
                for real
                  buffer_n(3) = change_n
                  field_n = buffer_n(1) + 0.5_mk*dt*change_n
              end foreach
            else
              foreach n in equi_mesh(mesh) with vec_fields(field,change,buffer) indices(i,j,k)
                for real
                  buffer_n(2*bs+1:3*bs) = change_n(:)
                  field_n(:) = buffer_n(1:bs) + 0.5_mk*dt*change_n(:)
              end foreach
            endif
          endif
        end select
      class is (ppm_t_particles_s)
        bs = ppm_dim
        !pset_s => field
        !foreach p in particles(pset_s) with positions(x) fields(dx=change) types(dx=vector)
        !  x_p(:) = bfr_p(:) + dt*dx_p(:)
        !end foreach
      class is (ppm_t_particles_d)
        bs = ppm_dim
        pset_d => field
        foreach p in particles(pset_d) with positions(x) vec_fields(dx=change,bfr=buffer)
          bfr_p(2*bs+1:3*bs) = dx_p(:)
          x_p(:) = bfr_p(1:bs) + 0.5_mk*dt*dx_p(:)
        end foreach
      end select
      field => this%fields%next()
      change => this%changes%next()
      buffer => this%buffers%next()
      discr => this%discretizations%next()
    end do
    check_false("associated(change)","Fields and changes should have same length")
    check_false("associated(discr)","Fields and discretizations should have same length")
    check_false("associated(buffer)","Fields and buffers should have same length")
  
  case (3)
    !----------------------------------------------------------------------
    ! STAGE 3: x_n + dt k3
    !----------------------------------------------------------------------
    rhs_info = this%rhsfunc(this%rhs_fields_discr,t+0.5_mk*dt,this%changes)
    
    field => this%fields%begin()
    change => this%changes%begin()
    buffer => this%buffers%begin()
    discr => this%discretizations%begin()
    do while (associated(field))
      check_associated(change,"Fields and changes should have same length")
      check_associated(discr,"Fields and discretizations should have same length")
      check_associated(buffer,"Fields and buffers should have same length")

      select type (field)
      class is (ppm_t_field)
        bs = field%lda
        select type (discr)
        class is (ppm_t_particles_s)
          !pset_s => discr
          !foreach p in particles(pset_s) with fields(f=field,df=change) types(f=scalar,df=scalar)
          !   bfr_p(2) = df_p
          !  f_p = bfr_p + dt*df_p
          !end foreach
        class is (ppm_t_particles_d)
          pset_d => discr
          if (field%lda.eq.1) then
            foreach p in particles(pset_d) with sca_fields(f=field,df=change) vec_fields(bfr=buffer)
              bfr_p(4) = df_p
              f_p = bfr_p(1) + dt*df_p
            end foreach
          else
            foreach p in particles(pset_d) with vec_fields(f=field,df=change,bfr=buffer)
              bfr_p(3*bs+1:4*bs) = df_p(:)
              f_p(:) = bfr_p(1:bs) + dt*df_p(:)
            end foreach
          endif
        class is (ppm_t_equi_mesh)
          mesh => discr
          if (ppm_dim.eq.2) then
            if (field%lda.eq.1) then
              foreach n in equi_mesh(mesh) with sca_fields(field,change) vec_fields(buffer) indices(i,j)
                for real
                  buffer_n(4) = change_n
                  field_n = buffer_n(1) + dt*change_n
              end foreach
            else
              foreach n in equi_mesh(mesh) with vec_fields(field,change,buffer) indices(i,j)
                for real
                  buffer_n(3*bs+1:4*bs) = change_n(:)
                  field_n(:) = buffer_n(1:bs) + dt*change_n(:)
              end foreach
            endif
          else
            if (field%lda.eq.1) then
              foreach n in equi_mesh(mesh) with sca_fields(field,change) vec_fields(buffer) indices(i,j,k)
                for real
                  buffer_n(4) = change_n
                  field_n = buffer_n(1) + dt*change_n
              end foreach
            else
              foreach n in equi_mesh(mesh) with vec_fields(field,change,buffer) indices(i,j,k)
                for real
                  buffer_n(3*bs+1:4*bs) = change_n(:)
                  field_n(:) = buffer_n(1:bs) + dt*change_n(:)
              end foreach
            endif
          endif
        end select
      class is (ppm_t_particles_s)
        bs = ppm_dim
        !pset_s => field
        !foreach p in particles(pset_s) with positions(x) fields(dx=change) types(dx=vector)
        !  x_p(:) = bfr_p(:) + dt*dx_p(:)
        !end foreach
      class is (ppm_t_particles_d)
        bs = ppm_dim
        pset_d => field
        foreach p in particles(pset_d) with positions(x) vec_fields(dx=change,bfr=buffer)
          bfr_p(3*bs+1:4*bs) = dx_p(:)
          x_p(:) = bfr_p(1:bs) + dt*dx_p(:)
        end foreach
      end select
      field => this%fields%next()
      change => this%changes%next()
      buffer => this%buffers%next()
      discr => this%discretizations%next()
    end do
    check_false("associated(change)","Fields and changes should have same length")
    check_false("associated(discr)","Fields and discretizations should have same length")
    check_false("associated(buffer)","Fields and buffers should have same length")
  
  case (4)
    !----------------------------------------------------------------------
    ! STAGE 4: Final step, interpolating the previous results
    ! x_n + 1/6 dt (k1 + 2 k2 + 2k3 +k4)
    !----------------------------------------------------------------------
    rhs_info = this%rhsfunc(this%rhs_fields_discr,t+dt,this%changes)
    
    field => this%fields%begin()
    change => this%changes%begin()
    buffer => this%buffers%begin()
    discr => this%discretizations%begin()
    do while (associated(field))
      check_associated(change,"Fields and changes should have same length")
      check_associated(discr,"Fields and discretizations should have same length")
      check_associated(buffer,"Fields and buffers should have same length")

      select type (field)
      class is (ppm_t_field)
        bs = field%lda
        select type (discr)
        class is (ppm_t_particles_s)
          !pset_s => discr
          !foreach p in particles(pset_s) with fields(f=field,df=change) types(f=scalar,df=scalar)
          !   bfr_p(2) = df_p
          !  f_p = bfr_p + dt*df_p
          !end foreach
        class is (ppm_t_particles_d)
          pset_d => discr
          if (field%lda.eq.1) then
            foreach p in particles(pset_d) with sca_fields(f=field,df=change) vec_fields(bfr=buffer)
              f_p = bfr_p(1) + 1.0_mk/6.0_mk * dt * &
              &     (bfr_p(2) + 2.0_mk*bfr_p(3) + 2.0_mk*bfr_p(4) + df_p)
            end foreach
          else
            foreach p in particles(pset_d) with vec_fields(f=field,df=change,bfr=buffer)
              f_p(:) = bfr_p(1:bs) + 1.0_mk/6.0_mk * dt * &
              &     (bfr_p(bs+1:2*bs) + 2.0_mk*bfr_p(2*bs+1:3*bs) + &
              &      2.0_mk*bfr_p(3*bs+1:4*bs) + df_p(:))
            end foreach
          endif
        class is (ppm_t_equi_mesh)
          mesh => discr
          if (ppm_dim.eq.2) then
            if (field%lda.eq.1) then
              foreach n in equi_mesh(mesh) with sca_fields(field,change) vec_fields(buffer) indices(i,j)
                for real
                  field_n = buffer_n(1) + 1.0_mk/6.0_mk * dt * &
                  &     (buffer_n(2) + 2.0_mk*buffer_n(3) + 2.0_mk*buffer_n(4) + change_n)
              end foreach
            else
              foreach n in equi_mesh(mesh) with vec_fields(field,change,buffer) indices(i,j)
                for real
                  field_n(:) = buffer_n(1:bs) + 1.0_mk/6.0_mk * dt * &
                  &     (buffer_n(bs+1:2*bs) + 2.0_mk*buffer_n(2*bs+1:3*bs) + &
                  &      2.0_mk*buffer_n(3*bs+1:4*bs) + change_n(:))
              end foreach
            endif
          else
            if (field%lda.eq.1) then
              foreach n in equi_mesh(mesh) with sca_fields(field,change) vec_fields(buffer) indices(i,j,k)
                for real
                  field_n = buffer_n(1) + 1.0_mk/6.0_mk * dt * &
                  &     (buffer_n(2) + 2.0_mk*buffer_n(3) + 2.0_mk*buffer_n(4) + change_n)
              end foreach
            else
              foreach n in equi_mesh(mesh) with vec_fields(field,change,buffer) indices(i,j,k)
                for real
                  field_n(:) = buffer_n(1:bs) + 1.0_mk/6.0_mk * dt * &
                  &     (buffer_n(bs+1:2*bs) + 2.0_mk*buffer_n(2*bs+1:3*bs) + &
                  &      2.0_mk*buffer_n(3*bs+1:4*bs) + change_n(:))
              end foreach
            endif
          endif
        end select
      class is (ppm_t_particles_s)
        !pset_s => field
        !foreach p in particles(pset_s) with positions(x) fields(dx=change) types(dx=vector)
        !  x_p(:) = bfr_p(:) + dt*dx_p(:)
        !end foreach
      class is (ppm_t_particles_d)
        pset_d => field
        foreach p in particles(pset_d) with positions(x) vec_fields(dx=change,bfr=buffer)
          x_p(:) = bfr_p(1:bs) + 1.0_mk/6.0_mk * dt * &
          &     (bfr_p(bs+1:2*bs) + 2.0_mk*bfr_p(2*bs+1:3*bs) + &
          &      2.0_mk*bfr_p(3*bs+1:4*bs) + dx_p(:))
        end foreach
      end select
      field => this%fields%next()
      change => this%changes%next()
      buffer => this%buffers%next()
      discr => this%discretizations%next()
    end do
    check_false("associated(change)","Fields and changes should have same length")
    check_false("associated(discr)","Fields and discretizations should have same length")
    check_false("associated(buffer)","Fields and buffers should have same length")
  
    t  = t  + dt

  end select


  end_subroutine()
end subroutine rk4_step



end module ppm_module_integrator_typedef
