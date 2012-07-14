
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

!----------------------------------------------------------------------
!  Types
!----------------------------------------------------------------------

type,extends(ppm_t_integrator_) :: ppm_t_integrator
  contains
  procedure :: create_s => integrator_create_s
  procedure :: create_d => integrator_create_d
  procedure :: destroy  => integrator_destroy
  procedure :: step_s => integrator_step_s
  procedure :: step_d => integrator_step_d
end type ppm_t_integrator

type,extends(ppm_t_integrator) :: ppm_t_eulerf
  contains
  procedure :: create_s => eulerf_create_s
  procedure :: create_d => eulerf_create_d
  procedure :: step_s => eulerf_step_s
  procedure :: step_d => eulerf_step_d
end type ppm_t_eulerf

type,extends(ppm_t_integrator) :: ppm_t_sts
  real(ppm_kind_single), dimension(20) :: stsnu_s
  real(ppm_kind_double), dimension(20) :: stsnu_d
  integer                              :: nsts = 1
  contains
  procedure :: create_s => sts_create_s
  procedure :: create_d => sts_create_d
  procedure :: step_s => sts_step_s
  procedure :: step_d => sts_step_d
end type ppm_t_sts

type,extends(ppm_t_options) :: ppm_t_sts_options_s
  integer               :: nsts = 1
  real(ppm_kind_single) :: nu   = 0.0_4
end type ppm_t_sts_options_s

type,extends(ppm_t_options) :: ppm_t_sts_options_d
  integer               :: nsts = 1
  real(ppm_kind_double) :: nu   = 0.0_8
end type ppm_t_sts_options_d

type,extends(ppm_t_integrator) :: ppm_t_tvdrk2
  contains
  procedure :: create_s => tvdrk2_create_s
  procedure :: create_d => tvdrk2_create_d
  procedure :: step_s => tvdrk2_step_s
  procedure :: step_d => tvdrk2_step_d
end type ppm_t_tvdrk2

type,extends(ppm_t_integrator) :: ppm_t_midrk2
  contains
  procedure :: create_s => midrk2_create_s
  procedure :: create_d => midrk2_create_d
  procedure :: step_s => midrk2_step_s
  procedure :: step_d => midrk2_step_d
end type ppm_t_midrk2

type,extends(ppm_t_integrator) :: ppm_t_rk4
  contains
  procedure :: create_s => rk4_create_s
  procedure :: create_d => rk4_create_d
  procedure :: step_s => rk4_step_s
  procedure :: step_d => rk4_step_d
end type ppm_t_rk4

!----------------------------------------------------------------------
!  Type-bound procedures
!----------------------------------------------------------------------
contains


template <RHST:[ppm_p_rhsfunc_s,ppm_p_rhsfunc_d],RHSV:[rhsfunc_s,rhsfunc_d]> nointerface suffixes [s,d]
subroutine integrator_create(this,variables,rhsfunc,rhs_variables,info,options)
  implicit none
  !----------------------------------------------------------------------
  !  Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_integrator)                :: this
  class(ppm_v_main_abstr)                :: variables
  !!! This vector holds all entities to be updated by the integrator the
  !!! elements may be fields or particle discretizations.
  !!!
  !!! fields must either have one discretization or wrapped together with the
  !!! intended discretization in a ppm_t_pair object
  procedure(RHST)                      :: rhsfunc
  !!! The right hand side function to be executed by the eulerf::step function
  class(ppm_v_var_discr_pair)          :: rhs_variables
  !!! The fields or particle properties (discr_data) to be passed to the right 
  !!! hand side function.
  !!!
  !!! The elements of the container are (var,discr) pairs where var can be
  !!! either a field or a discr_data and discr is the corresponding
  !!! discretization.
  integer,                 intent(  out) :: info
  class(ppm_t_options),target,optional,intent(in   ) :: options
  !----------------------------------------------------------------------
  !  Variables
  !----------------------------------------------------------------------
  class(ppm_t_main_abstr), pointer       :: el => null()
  class(ppm_t_main_abstr), pointer       :: temp => null()
  class(ppm_t_field_), pointer           :: cfield => null()
  class(ppm_t_field_), pointer           :: pfield => null()
  class(ppm_t_part_prop_s_), pointer     :: cprop_s => null()
  class(ppm_t_part_prop_d_), pointer     :: cprop_d => null()
  class(ppm_t_field_), pointer           :: buf => null()
  class(ppm_t_field_), pointer           :: el_f => null()
  class(ppm_t_discr_data), pointer       :: prop => null()
  class(ppm_t_discr_info_), pointer      :: di => null()
  class(ppm_t_discr_kind), pointer       :: d => null()
  class(ppm_t_field_discr_pair), pointer :: el_p => null()
  class(ppm_t_var_discr_pair), pointer   :: el_vp => null()
  integer                                :: ibuf
  integer                                :: ifield
  character(len=16)                      :: bname
  character(len=16)                      :: cname
  logical                                :: mkbuf = .true.
  start_subroutine("integrator_create")

  nullify(this%buffers)
  if (this%scheme_memsize.eq.0) then
    mkbuf = .false.
  else
    mkbuf = .true.
  end if

  allocate(this%variables,stat=info)
  or_fail_alloc("this%fields")
  if (mkbuf) then
    allocate(this%buffers,stat=info)
    or_fail_alloc("this%buffers")
  end if
  allocate(this%changes,stat=info)
  or_fail_alloc("this%changes")
  allocate(this%discretizations,stat=info)
  or_fail_alloc("this%discretizations")
  el => variables%begin()
  ifield = 0
  do while (associated(el))
    ifield = ifield + 1
    write(cname,'(A,I0)') 'ode_change_',ifield
    if (mkbuf) then
      allocate(ppm_t_field::buf,stat=info)
      write(bname,'(A,I0)') 'ode_buffer_',ifield
    end if
    select type(el)
    class is (ppm_t_field_)
      el_f => el
      di => el_f%discr_info%begin()
      allocate(ppm_t_field::cfield,stat=info)
      call cfield%create(el_f%lda,info,name=cname)
      call cfield%discretize_on(di%discr_ptr,info)
      or_fail("Discretizing change failed")
      temp => cfield
      call this%changes%push(temp,info)
      or_fail("Pushing change failed")
      call this%discretizations%push(di%discr_ptr,info)
      or_fail("Pushing change discretization failed")
      temp => el
      call this%variables%push(temp,info)
      or_fail("Pushing field failed")
      if (mkbuf) then
        call buf%create(el_f%lda*this%scheme_memsize,info,name=bname)
        call buf%discretize_on(di%discr_ptr,info,with_ghosts=.false.)
      endif
    class is (ppm_t_discr_data)
      prop => el
      select type (parts => prop%discr)
      class is (ppm_t_particles_s)
        call parts%create_prop(info,name=cname,lda=prop%lda,part_prop=cprop_s)
        or_fail("creating change property")
        temp => cprop_s
        call this%changes%push(temp,info)
      class is (ppm_t_particles_d)
        call parts%create_prop(info,name=cname,lda=prop%lda,part_prop=cprop_d)
        or_fail("creating change property")
        temp => cprop_d
        call this%changes%push(temp,info)
      class default
        fail("Only particle properties allowed",ppm_err_argument)
      end select
      or_fail("Pushing change failed")
      call this%discretizations%push(prop%discr,info)
      or_fail("Pushing change discretization failed")
      temp => el
      call this%variables%push(temp,info)
      or_fail("Pushing field failed")
      !if (mkbuf) then
      !  call buf%create(el_f%lda*this%scheme_memsize,info,name=bname)
      !  call buf%discretize_on(di%discr_ptr,info,with_ghosts=.false.)
      !endif
    class is (ppm_t_discr_kind)
      d => el
      select type (parts => d)
      class is (ppm_t_particles_s)
        call parts%create_prop(info,name=cname,lda=ppm_dim,part_prop=cprop_s)
        or_fail("creating change property")
        temp => cprop_s
        call this%changes%push(temp,info)
      class is (ppm_t_particles_d)
        call parts%create_prop(info,name=cname,lda=ppm_dim,part_prop=cprop_d)
        or_fail("creating change property")
        temp => cprop_d
        call this%changes%push(temp,info)
      class default
        fail("Only particles allowed",ppm_err_argument)
      end select
      or_fail("Pushing change failed")
      call this%discretizations%push(d,info)
      or_fail("Pushing change discretization failed")
      temp => d
      call this%variables%push(temp,info)
      or_fail("Pushing positions failed")
      if (mkbuf) then
        call buf%create(ppm_dim*this%scheme_memsize,info,name=bname)
        call buf%discretize_on(d,info,with_ghosts=.false.)
      endif
    class is (ppm_t_field_discr_pair)
      el_p => el
      allocate(ppm_t_field::cfield,stat=info)
      call cfield%create(el_p%field%lda,info,name=cname)
      call cfield%discretize_on(el_p%discr,info)
      or_fail("Discretizing change failed")
      temp => cfield
      call this%changes%push(temp,info)
      or_fail("Pushing change failed")
      call this%discretizations%push(el_p%discr,info)
      or_fail("Pushing change discretization failed")
      temp => el_p%field
      call this%variables%push(temp,info)
      or_fail("Pushing field failed")
      if (mkbuf) then
        call buf%create(el_p%field%lda*this%scheme_memsize,info,name=bname)
        call buf%discretize_on(el_p%discr,info,with_ghosts=.false.)
      endif
    class default
      fail("variables should only contain fields, props, discrs and pairs",ppm_err_argument)
  end select
  if (mkbuf) then
    call this%buffers%push(buf,info)
  endif
  el => variables%next()
  end do

  ! copy the rhs_variables vector
  allocate(this%rhs_variables,stat=info)
  or_fail_alloc("this%rhs_variables")
  el_vp => rhs_variables%begin()
  do while (associated(el_vp))
    call this%rhs_variables%push(el_vp,info)
    el_vp => rhs_variables%next()
  end do
 
  ! set the right hand side function
  this%RHSV => rhsfunc
 
  end_subroutine()
end subroutine

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
  class(ppm_t_main_abstr), pointer               :: change => null()
  class(ppm_t_main_abstr), pointer               :: buffer => null()
  start_subroutine("eulerf_destroy")
  
  change => this%changes%begin()
  do while (associated(change))
    !call change%destroy(info)  !FIXME
    ! or_fail("Destroying change failed")
    change => this%changes%next()
  end do
  deallocate(this%changes,stat=info)
  if (associated(this%buffers)) then
    buffer => this%buffers%begin()
    do while (associated(buffer))
      !call buffer%destroy(info) !FIXME
      ! or_fail("Destroying buffer failed")
      buffer => this%buffers%next()
    end do
  end if
 
  end_subroutine()
end subroutine integrator_destroy

template <MK:[ppm_kind_single,ppm_kind_double],part_type:[ppm_t_particles_s,ppm_t_particles_d],prop_type:[ppm_t_part_prop_s,ppm_t_part_prop_d],RHSV:[rhsfunc_s,rhsfunc_d]> nointerface suffixes [s,d]
subroutine integrator_step(this,t,dt,istage,info)
  implicit none
  !----------------------------------------------------------------------
  ! Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_integrator)                            :: this
  real(MK)                         ,intent(inout)    :: t
  real(MK)                         ,intent(in   )    :: dt
  integer                          ,intent(in   )    :: istage
  integer                          ,intent(  out)    :: info
  !----------------------------------------------------------------------
  ! Variables
  !----------------------------------------------------------------------
  real(MK)                                   :: rhs_info
  class(ppm_t_main_abstr), pointer           :: var => null()
  class(ppm_t_main_abstr), pointer           :: change => null()
  class(ppm_t_discr_kind), pointer           :: discr  => null()
  start_subroutine("integrator_step")

  ! nothing happening here

  end_subroutine()
end subroutine integrator_step

template <RHST:[ppm_p_rhsfunc_s,ppm_p_rhsfunc_d]> nointerface suffixes [s,d]
subroutine eulerf_create(this,variables,rhsfunc,rhs_variables,info,options)
  implicit none
  !----------------------------------------------------------------------
  !  Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_eulerf)                    :: this
  class(ppm_v_main_abstr)                :: variables
  !!! This vector holds all entities to be updated by the integrator the
  !!! elements may be fields or particle discretizations.
  !!!
  !!! fields must either have one discretization or wrapped together with the
  !!! intended discretization in a ppm_t_pair object
  procedure(RHST)                      :: rhsfunc
  !!! The right hand side function to be executed by the eulerf::step function
  class(ppm_v_var_discr_pair)          :: rhs_variables
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
  
  call this%ppm_t_integrator%create(variables,rhsfunc,rhs_variables,info)

  end_subroutine()

end subroutine eulerf_create

template <MK:[ppm_kind_single,ppm_kind_double],part_type:[ppm_t_particles_s,ppm_t_particles_d],prop_type:[ppm_t_part_prop_s,ppm_t_part_prop_d],RHSV:[rhsfunc_s,rhsfunc_d]> nointerface suffixes [s,d]
subroutine eulerf_step(this,t,dt,istage,info)
  implicit none
  !----------------------------------------------------------------------
  ! Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_eulerf)                                :: this
  real(MK)                         ,intent(inout)    :: t
  real(MK)                         ,intent(in   )    :: dt
  integer                          ,intent(in   )    :: istage
  integer                          ,intent(  out)    :: info
  !----------------------------------------------------------------------
  ! Variables
  !----------------------------------------------------------------------
  real(MK)                                   :: rhs_info
  class(ppm_t_main_abstr), pointer           :: var => null()
  class(ppm_t_main_abstr), pointer           :: change => null()
  class(ppm_t_discr_kind), pointer           :: discr  => null()
  start_subroutine("eulerf_step")

  rhs_info = this%RHSV(this%rhs_variables,t,this%changes)

  var => this%variables%begin()
  change => this%changes%begin()
  discr => this%discretizations%begin()
  do while (associated(var))
    check_associated(change,"Fields and changes should have same length")
    check_associated(discr,"Fields and discretizations should have same length")

    foreach e in discrp(discr) with vars([u=>var],[du=>change]) prec(MK) part_type(part_type) prop_type(prop_type)
      for sca
        u_e = u_e + dt*du_e
      for vec
        u_e(:) = u_e(:) + dt*du_e(:)
      for pos
        u_e(:) = u_e(:) + dt*du_e(:)
    end foreach

    var => this%variables%next()
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


template <stsnu:[stsnu_s,stsnu_d],MK:[ppm_kind_single,ppm_kind_double],RHST:[ppm_p_rhsfunc_s,ppm_p_rhsfunc_d]> nointerface suffixes [s,d]
subroutine sts_create(this,variables,rhsfunc,rhs_variables,info,options)
  implicit none
  parameter(mk,integer,MK)
  !----------------------------------------------------------------------
  !  Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_sts)                       :: this
  class(ppm_v_main_abstr)                :: variables
  !!! This vector holds all entities to be updated by the integrator the
  !!! elements may be fields or particle discretizations.
  !!!
  !!! fields must either have one discretization or wrapped together with the
  !!! intended discretization in a ppm_t_pair object
  procedure(RHST)                        :: rhsfunc
  !!! The right hand side function to be executed by the eulerf::step function
  class(ppm_v_var_discr_pair)          :: rhs_variables
  !!! The fields to be passed to the right hand side function.
  !!!
  !!! The elements of the container can be fields, or pairs continaing a field
  !!! and its intended discretization.
  integer,                 intent(  out) :: info
  class(ppm_t_options),target,optional,intent(in   ) :: options
  !----------------------------------------------------------------------
  ! Variables
  !----------------------------------------------------------------------
  class(ppm_t_sts_options_s), pointer           :: sts_options_s
  class(ppm_t_sts_options_d), pointer           :: sts_options_d
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
  class is (ppm_t_sts_options_s)
  sts_options_s => options
  this%nsts      = sts_options_s%nsts
  this%stsnu(this%nsts) = sts_options_s%nu
  class is (ppm_t_sts_options_d)
  sts_options_d => options
  this%nsts      = sts_options_d%nsts
  this%stsnu(this%nsts) = sts_options_d%nu
  end select
  
  call this%ppm_t_integrator%create(variables,rhsfunc,rhs_variables,info)


  end_subroutine()

end subroutine sts_create


template <MK:[ppm_kind_single,ppm_kind_double],part_type:[ppm_t_particles_s,ppm_t_particles_d],prop_type:[ppm_t_part_prop_s,ppm_t_part_prop_d],RHSV:[rhsfunc_s,rhsfunc_d]> nointerface suffixes [s,d]
subroutine sts_step(this,t,dt,istage,info)
  implicit none
  parameter(mk,integer,MK)
  !----------------------------------------------------------------------
  ! Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_sts)                                   :: this
  real(mk)                         ,intent(inout)    :: t
  real(mk)                         ,intent(in   )    :: dt
  integer                          ,intent(in   )    :: istage
  integer                          ,intent(  out)    :: info
  !----------------------------------------------------------------------
  ! Variables
  !----------------------------------------------------------------------
  real(mk)                                   :: rhs_info
  class(ppm_t_main_abstr), pointer           :: var => null()
  class(ppm_t_main_abstr), pointer           :: change => null()
  class(ppm_t_discr_kind), pointer           :: discr  => null()
  real(mk)                                   :: tau
  real(mk)                                   :: m_pi
  start_subroutine("sts_step")
           
  m_pi = ACOS(-1.0_mk)
  !-----------------------------------------------------
  !  compute the new dt
  !-----------------------------------------------------
  if (mk.eq.ppm_kind_single) then
    tau = dt/((this%stsnu_s(this%nsts) - 1.0_mk) * &
    &     cos((2.0_mk * real(istage,mk) - 1.0_mk)/real(this%nsts,mk) * M_PI * 0.5_mk) + &
    &     1.0_mk + this%stsnu_s(this%nsts))
  else
    tau = dt/((this%stsnu_d(this%nsts) - 1.0_mk) * &
    &     cos((2.0_mk * real(istage,mk) - 1.0_mk)/real(this%nsts,mk) * M_PI * 0.5_mk) + &
    &     1.0_mk + this%stsnu_d(this%nsts))
  end if
  
  rhs_info = this%RHSV(this%rhs_variables,t,this%changes)
  
  var => this%variables%begin()
  change => this%changes%begin()
  discr => this%discretizations%begin()
  do while (associated(var))
    check_associated(change,"Fields and changes should have same length")
    check_associated(discr,"Fields and discretizations should have same length")

    foreach e in discrp(discr) with vars([u=>var],[du=>change]) prec(MK) part_type(part_type) prop_type(prop_type)
      for sca
        u_e = u_e + tau*du_e
      for vec
        u_e(:) = u_e(:) + tau*du_e(:)
      for pos
        u_e(:) = u_e(:) + tau*du_e(:)
    end foreach

    var => this%variables%next()
    change => this%changes%next()
    discr => this%discretizations%next()
  end do
  check_false("associated(change)","Fields and changes should have same length")
  check_false("associated(discr)","Fields and discretizations should have same length")

  t  = t  + tau

  end_subroutine()
end subroutine sts_step

template <RHST:[ppm_p_rhsfunc_s,ppm_p_rhsfunc_d]> nointerface suffixes [s,d]
subroutine tvdrk2_create(this,variables,rhsfunc,rhs_variables,info,options)
  implicit none
  !----------------------------------------------------------------------
  !  Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_tvdrk2)                    :: this
  class(ppm_v_main_abstr)                :: variables
  !!! This vector holds all entities to be updated by the integrator the
  !!! elements may be fields or particle discretizations.
  !!!
  !!! fields must either have one discretization or wrapped together with the
  !!! intended discretization in a ppm_t_pair object
  procedure(RHST)                        :: rhsfunc
  !!! The right hand side function to be executed by the eulerf::step function
  class(ppm_v_var_discr_pair)          :: rhs_variables
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
  call this%ppm_t_integrator%create(variables,rhsfunc,rhs_variables,info)


  end_subroutine()

end subroutine tvdrk2_create


template <MK:[ppm_kind_single,ppm_kind_double],part_type:[ppm_t_particles_s,ppm_t_particles_d],prop_type:[ppm_t_part_prop_s,ppm_t_part_prop_d],RHSV:[rhsfunc_s,rhsfunc_d]> nointerface suffixes [s,d]
subroutine tvdrk2_step(this,t,dt,istage,info)
  implicit none
  parameter(mk,integer,MK)
  !----------------------------------------------------------------------
  ! Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_tvdrk2)                                :: this
  real(mk)                         ,intent(inout)    :: t
  real(mk)                         ,intent(in   )    :: dt
  integer                          ,intent(in   )    :: istage
  integer                          ,intent(  out)    :: info
  !----------------------------------------------------------------------
  ! Variables
  !----------------------------------------------------------------------
  real(MK)                                   :: rhs_info
  class(ppm_t_main_abstr), pointer           :: var => null()
  class(ppm_t_main_abstr), pointer           :: change => null()
  class(ppm_t_field_), pointer               :: buffer => null()
  class(ppm_t_discr_kind), pointer           :: discr  => null()
  start_subroutine("tvdrk2_step")
           


  ! TODO check if buffers where mapped
  select case (istage)
  case (1)
    !----------------------------------------------------------------------
    ! STAGE 1: Store the current field in a secondary buffer and update it
    ! with dt*du
    !----------------------------------------------------------------------
    rhs_info = this%RHSV(this%rhs_variables,t,this%changes)
    
    var => this%variables%begin()
    change => this%changes%begin()
    discr => this%discretizations%begin()
    buffer => this%buffers%begin()
    do while (associated(var))
      check_associated(change,"Fields and changes should have same length")
      check_associated(discr,"Fields and discretizations should have same length")

      foreach e in discrp(discr) with vars([u=>var],[du=>change]) buffer([bfr=>buffer]) prec(MK) part_type(part_type) prop_type(prop_type)
        for sca
          bfr_e = u_e
          u_e = u_e + dt*du_e
        for vec
          bfr_e(:) = u_e(:)
          u_e(:) = u_e(:) + dt*du_e(:)
        for pos
          bfr_e(:) = u_e(:)
          u_e(:) = u_e(:) + dt*du_e(:)
      end foreach

      var => this%variables%next()
      change => this%changes%next()
      discr => this%discretizations%next()
      buffer => this%buffers%next()
    end do
    check_false("associated(change)","Fields and changes should have same length")
    check_false("associated(discr)","Fields and discretizations should have same length")
    check_false("associated(buffer)","Fields and buffers should have same length")
  
  case (2)
    !----------------------------------------------------------------------
    ! STAGE 2: Do another update and interpolate with previously saved data
    !----------------------------------------------------------------------
    rhs_info = this%RHSV(this%rhs_variables,t+dt,this%changes)
    
    var => this%variables%begin()
    change => this%changes%begin()
    discr => this%discretizations%begin()
    buffer => this%buffers%begin()
    do while (associated(var))
      check_associated(change,"Fields and changes should have same length")
      check_associated(discr,"Fields and discretizations should have same length")

      foreach e in discrp(discr) with vars([u=>var],[du=>change]) buffer([bfr=>buffer]) prec(MK) part_type(part_type) prop_type(prop_type)
        for sca
          u_e = u_e + dt*du_e
          u_e  = 0.5_mk * (u_e + bfr_e)
        for vec
          u_e(:) = u_e(:) + dt*du_e(:)
          u_e(:)  = 0.5_mk * (u_e(:) + bfr_e(:))
        for pos
          u_e(:) = u_e(:) + dt*du_e(:)
          u_e(:)  = 0.5_mk * (u_e(:) + bfr_e(:))
      end foreach

      var => this%variables%next()
      change => this%changes%next()
      discr => this%discretizations%next()
      buffer => this%buffers%next()
    end do
    check_false("associated(change)","Fields and changes should have same length")
    check_false("associated(discr)","Fields and discretizations should have same length")
    check_false("associated(buffer)","Fields and buffers should have same length")
  
    t  = t  + dt
  end select


  end_subroutine()
end subroutine tvdrk2_step


template <RHST:[ppm_p_rhsfunc_s,ppm_p_rhsfunc_d]> nointerface suffixes [s,d]
subroutine midrk2_create(this,variables,rhsfunc,rhs_variables,info,options)
  implicit none
  !----------------------------------------------------------------------
  !  Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_midrk2)                    :: this
  class(ppm_v_main_abstr)                :: variables
  !!! This vector holds all entities to be updated by the integrator the
  !!! elements may be fields or particle discretizations.
  !!!
  !!! fields must either have one discretization or wrapped together with the
  !!! intended discretization in a ppm_t_pair object
  procedure(RHST)                        :: rhsfunc
  !!! The right hand side function to be executed by the eulerf::step function
  class(ppm_v_var_discr_pair)          :: rhs_variables
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
  call this%ppm_t_integrator%create(variables,rhsfunc,rhs_variables,info)


  end_subroutine()

end subroutine midrk2_create


template <MK:[ppm_kind_single,ppm_kind_double],part_type:[ppm_t_particles_s,ppm_t_particles_d],prop_type:[ppm_t_part_prop_s,ppm_t_part_prop_d],RHSV:[rhsfunc_s,rhsfunc_d]> nointerface suffixes [s,d]
subroutine midrk2_step(this,t,dt,istage,info)
  implicit none
  parameter(mk,integer,MK)
  !----------------------------------------------------------------------
  ! Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_midrk2)                                :: this
  real(mk)                         ,intent(inout)    :: t
  real(mk)                         ,intent(in   )    :: dt
  integer                          ,intent(in   )    :: istage
  integer                          ,intent(  out)    :: info
  !----------------------------------------------------------------------
  ! Variables
  !----------------------------------------------------------------------
  real(MK)                                   :: rhs_info
  class(ppm_t_main_abstr), pointer           :: var => null()
  class(ppm_t_main_abstr), pointer           :: change => null()
  class(ppm_t_field_), pointer               :: buffer => null()
  class(ppm_t_discr_kind), pointer           :: discr  => null()
  start_subroutine("midrk2_step")
           


  ! TODO check if buffers where mapped
  select case (istage)
  case (1)
    !----------------------------------------------------------------------
    ! STAGE 1: Compute midpoint
    ! Store the current field in a secondary buffer and update it
    ! with 0.5*dt*du
    !----------------------------------------------------------------------
    rhs_info = this%RHSV(this%rhs_variables,t,this%changes)
    
    var => this%variables%begin()
    change => this%changes%begin()
    discr => this%discretizations%begin()
    buffer => this%buffers%begin()
    do while (associated(var))
      check_associated(change,"Fields and changes should have same length")
      check_associated(discr,"Fields and discretizations should have same length")

      foreach e in discrp(discr) with vars([u=>var],[du=>change]) buffer([bfr=>buffer]) prec(MK) part_type(part_type) prop_type(prop_type)
        for sca
          bfr_e = u_e
          u_e = u_e + 0.5_mk*dt*du_e
        for vec
          bfr_e(:) = u_e(:)
          u_e(:) = u_e(:) + 0.5_mk*dt*du_e(:)
        for pos
          bfr_e(:) = u_e(:)
          u_e(:) = u_e(:) + 0.5_mk*dt*du_e(:)
      end foreach

      var => this%variables%next()
      change => this%changes%next()
      discr => this%discretizations%next()
      buffer => this%buffers%next()
    end do
    check_false("associated(change)","Fields and changes should have same length")
    check_false("associated(discr)","Fields and discretizations should have same length")
    check_false("associated(buffer)","Fields and buffers should have same length")
  
  case (2)
    !----------------------------------------------------------------------
    ! STAGE 2: Do another update and interpolate with previously saved data
    !----------------------------------------------------------------------
    rhs_info = this%RHSV(this%rhs_variables,t+0.5_mk*dt,this%changes)
    
    var => this%variables%begin()
    change => this%changes%begin()
    discr => this%discretizations%begin()
    buffer => this%buffers%begin()
    do while (associated(var))
      check_associated(change,"Fields and changes should have same length")
      check_associated(discr,"Fields and discretizations should have same length")

      foreach e in discrp(discr) with vars([u=>var],[du=>change]) buffer([bfr=>buffer]) prec(MK) part_type(part_type) prop_type(prop_type)
        for sca
          u_e = bfr_e + dt*du_e
        for vec
          u_e(:) = bfr_e(:) + dt*du_e(:)
        for pos
          u_e(:) = bfr_e(:) + dt*du_e(:)
      end foreach

      var => this%variables%next()
      change => this%changes%next()
      discr => this%discretizations%next()
      buffer => this%buffers%next()
    end do
    check_false("associated(change)","Fields and changes should have same length")
    check_false("associated(discr)","Fields and discretizations should have same length")
    check_false("associated(buffer)","Fields and buffers should have same length")
  
    t  = t  + dt

  end select


  end_subroutine()
end subroutine midrk2_step


template <RHST:[ppm_p_rhsfunc_s,ppm_p_rhsfunc_d]> nointerface suffixes [s,d]
subroutine rk4_create(this,variables,rhsfunc,rhs_variables,info,options)
  implicit none
  !----------------------------------------------------------------------
  !  Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_rk4)                    :: this
  class(ppm_v_main_abstr)                :: variables
  !!! This vector holds all entities to be updated by the integrator the
  !!! elements may be fields or particle discretizations.
  !!!
  !!! fields must either have one discretization or wrapped together with the
  !!! intended discretization in a ppm_t_pair object
  procedure(RHST)                        :: rhsfunc
  !!! The right hand side function to be executed by the eulerf::step function
  class(ppm_v_var_discr_pair)          :: rhs_variables
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
  call this%ppm_t_integrator%create(variables,rhsfunc,rhs_variables,info)


  end_subroutine()

end subroutine rk4_create


template <MK:[ppm_kind_single,ppm_kind_double],part_type:[ppm_t_particles_s,ppm_t_particles_d],prop_type:[ppm_t_part_prop_s,ppm_t_part_prop_d],RHSV:[rhsfunc_s,rhsfunc_d]> nointerface suffixes [s,d]
subroutine rk4_step(this,t,dt,istage,info)
  implicit none
  parameter(mk,integer,MK)
  !----------------------------------------------------------------------
  ! Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_rk4)                                :: this
  real(mk)                         ,intent(inout)    :: t
  real(mk)                         ,intent(in   )    :: dt
  integer                          ,intent(in   )    :: istage
  integer                          ,intent(  out)    :: info
  !----------------------------------------------------------------------
  ! Variables
  !----------------------------------------------------------------------
  real(MK)                                   :: rhs_info
  class(ppm_t_main_abstr), pointer           :: var => null()
  class(ppm_t_main_abstr), pointer           :: change => null()
  class(ppm_t_field_), pointer               :: buffer => null()
  class(ppm_t_discr_kind), pointer           :: discr  => null()
  integer                                    :: bs
  start_subroutine("rk4_step")
           


  ! TODO check if buffers where mapped
  select case (istage)
  case (1)
    !----------------------------------------------------------------------
    ! STAGE 1: x_n + 1/2 dt k1
    !----------------------------------------------------------------------
    rhs_info = this%RHSV(this%rhs_variables,t,this%changes)
    
    var => this%variables%begin()
    change => this%changes%begin()
    discr => this%discretizations%begin()
    buffer => this%buffers%begin()
    do while (associated(var))
      check_associated(change,"Fields and changes should have same length")
      check_associated(discr,"Fields and discretizations should have same length")

      foreach e in discrp(discr) with vars([u=>var],[du=>change]) buffer([bfr=>buffer],vecbuf=true) prec(MK) part_type(part_type) prop_type(prop_type)
        for sca
          bfr_e(1) = u_e
          bfr_e(2) = du_e
          u_e = u_e + 0.5_mk*dt*du_e
        for vec
          bfr_e(1:ulda) = u_e(:)
          bfr_e(ulda+1:2*ulda) = du_e(:)
          u_e(:) = u_e(:) + 0.5_mk*dt*du_e(:)
        for pos
          bfr_e(1:ulda) = u_e(:)
          bfr_e(ulda+1:2*ulda) = du_e(:)
          u_e(:) = u_e(:) + 0.5_mk*dt*du_e(:)
      end foreach

      var => this%variables%next()
      change => this%changes%next()
      discr => this%discretizations%next()
      buffer => this%buffers%next()
    end do
    check_false("associated(change)","Fields and changes should have same length")
    check_false("associated(discr)","Fields and discretizations should have same length")
    check_false("associated(buffer)","Fields and buffers should have same length")
  
  case (2)
    !----------------------------------------------------------------------
    ! STAGE 2: x_n + 1/2 dt k2
    !----------------------------------------------------------------------
    rhs_info = this%RHSV(this%rhs_variables,t+0.5_mk*dt,this%changes)
    
    var => this%variables%begin()
    change => this%changes%begin()
    discr => this%discretizations%begin()
    buffer => this%buffers%begin()
    do while (associated(var))
      check_associated(change,"Fields and changes should have same length")
      check_associated(discr,"Fields and discretizations should have same length")

      foreach e in discrp(discr) with vars([u=>var],[du=>change]) buffer([bfr=>buffer],vecbuf=true) prec(MK) part_type(part_type) prop_type(prop_type)
        for sca
          bfr_e(3) = du_e
          u_e = bfr_e(1) + 0.5_mk*dt*du_e
        for vec
          bfr_e(2*ulda+1:3*ulda) = du_e(:)
          u_e(:) = bfr_e(1:ulda) + 0.5_mk*dt*du_e(:)
        for pos
          bfr_e(2*ulda+1:3*ulda) = du_e(:)
          u_e(:) = bfr_e(1:ulda) + 0.5_mk*dt*du_e(:)
      end foreach

      var => this%variables%next()
      change => this%changes%next()
      discr => this%discretizations%next()
      buffer => this%buffers%next()
    end do
    check_false("associated(change)","Fields and changes should have same length")
    check_false("associated(discr)","Fields and discretizations should have same length")
    check_false("associated(buffer)","Fields and buffers should have same length")
  
  case (3)
    !----------------------------------------------------------------------
    ! STAGE 3: x_n + dt k3
    !----------------------------------------------------------------------
    rhs_info = this%RHSV(this%rhs_variables,t+0.5_mk*dt,this%changes)
    
    var => this%variables%begin()
    change => this%changes%begin()
    discr => this%discretizations%begin()
    buffer => this%buffers%begin()
    do while (associated(var))
      check_associated(change,"Fields and changes should have same length")
      check_associated(discr,"Fields and discretizations should have same length")

      foreach e in discrp(discr) with vars([u=>var],[du=>change]) buffer([bfr=>buffer],vecbuf=true) prec(MK) part_type(part_type) prop_type(prop_type)
        for sca
          bfr_e(4) = du_e
          u_e = bfr_e(1) + dt*du_e
        for vec
          bfr_e(3*ulda+1:4*ulda) = du_e(:)
          u_e(:) = bfr_e(1:ulda) + dt*du_e(:)
        for pos
          bfr_e(3*ulda+1:4*ulda) = du_e(:)
          u_e(:) = bfr_e(1:ulda) + dt*du_e(:)
      end foreach

      var => this%variables%next()
      change => this%changes%next()
      discr => this%discretizations%next()
      buffer => this%buffers%next()
    end do
    check_false("associated(change)","Fields and changes should have same length")
    check_false("associated(discr)","Fields and discretizations should have same length")
    check_false("associated(buffer)","Fields and buffers should have same length")
  
  case (4)
    !----------------------------------------------------------------------
    ! STAGE 4: Final step, interpolating the previous results
    ! x_n + 1/6 dt (k1 + 2 k2 + 2k3 +k4)
    !----------------------------------------------------------------------
    rhs_info = this%RHSV(this%rhs_variables,t+dt,this%changes)
    
    var => this%variables%begin()
    change => this%changes%begin()
    discr => this%discretizations%begin()
    buffer => this%buffers%begin()
    do while (associated(var))
      check_associated(change,"Fields and changes should have same length")
      check_associated(discr,"Fields and discretizations should have same length")

      foreach e in discrp(discr) with vars([u=>var],[du=>change]) buffer([bfr=>buffer],vecbuf=true) prec(MK) part_type(part_type) prop_type(prop_type)
        for sca
          u_e = bfr_e(1) + 1.0_mk/6.0_mk * dt * &
          &     (bfr_e(2) + 2.0_mk*bfr_e(3) + 2.0_mk*bfr_e(4) + du_e)
        for vec
          u_e(:) = bfr_e(1:ulda) + 1.0_mk/6.0_mk * dt * &
          &     (bfr_e(ulda+1:2*ulda) + 2.0_mk*bfr_e(2*ulda+1:3*ulda) + &
          &      2.0_mk*bfr_e(3*ulda+1:4*ulda) + du_e(:))
        for pos
          u_e(:) = bfr_e(1:ulda) + 1.0_mk/6.0_mk * dt * &
          &     (bfr_e(ulda+1:2*ulda) + 2.0_mk*bfr_e(2*ulda+1:3*ulda) + &
          &      2.0_mk*bfr_e(3*ulda+1:4*ulda) + du_e(:))
      end foreach

      var => this%variables%next()
      change => this%changes%next()
      discr => this%discretizations%next()
      buffer => this%buffers%next()
    end do
    check_false("associated(change)","Fields and changes should have same length")
    check_false("associated(discr)","Fields and discretizations should have same length")
    check_false("associated(buffer)","Fields and buffers should have same length")
  
    t  = t  + dt

  end select


  end_subroutine()
end subroutine rk4_step



end module ppm_module_integrator_typedef
