
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

!----------------------------------------------------------------------
!  Types
!----------------------------------------------------------------------

type,extends(ppm_t_integrator_) :: ppm_t_eulerf
  contains
  procedure :: create => eulerf_create
  procedure :: destroy => eulerf_destroy
  procedure :: step => eulerf_step
end type ppm_t_eulerf

!----------------------------------------------------------------------
!  Type-bound procedures
!----------------------------------------------------------------------
contains


subroutine eulerf_create(this,fields,rhsfunc,rhs_fields,info)
  implicit none
  !----------------------------------------------------------------------
  !  Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_eulerf)                    :: this
  class(ppm_v_main_abstr), pointer       :: fields
  !!! This vector holds all entities to be updated by the integrator the
  !!! elements may be fields or particle discretizations.
  !!!
  !!! fields must either have one discretization or wrapped together with the
  !!! intended discretization in a ppm_t_pair object
  procedure(ppm_p_rhsfunc)               :: rhsfunc
  !!! The right hand side function to be executed by the eulerf::step function
  class(ppm_v_main_abstr), pointer       :: rhs_fields
  !!! The fields to be passed to the right hand side function.
  integer,                 intent(  out) :: info
  !----------------------------------------------------------------------
  !  Variables
  !----------------------------------------------------------------------
  class(ppm_t_main_abstr), pointer       :: el => null()
  class(ppm_t_main_abstr), pointer       :: temp => null()
  class(ppm_t_field_), pointer           :: c => null()
  class(ppm_t_field_), pointer           :: el_f => null()
  class(ppm_t_discr_info_), pointer      :: di => null()
  class(ppm_t_discr_kind), pointer       :: d => null()
  class(ppm_t_field_discr_pair), pointer :: el_p => null()
 
  start_subroutine("eulerf_create")
  this%scheme_order   = 1
  this%scheme_memsize = 1
  this%scheme_nstages = 1
  this%scheme_kickoff = 1
  

 
  allocate(this%fields,stat=info)
  allocate(this%changes,stat=info)
  allocate(this%discretizations,stat=info)
  el => fields%begin()
  do while (associated(el))
   allocate(ppm_t_field::c,stat=info)
   call c%create(ppm_dim,info)
   select type(el)
   class is (ppm_t_field_)
     el_f => el
     di => el_f%discr_info%begin()
     call c%discretize_on(di%discr_ptr,info)
     call this%discretizations%push(di%discr_ptr,info)
     temp => el
     call this%fields%push(temp,info)
   class is (ppm_t_discr_kind)
     d => el
     call c%discretize_on(d,info)
     call this%discretizations%push(d,info)
     temp => d
     call this%fields%push(temp,info)
   class is (ppm_t_field_discr_pair)
     el_p => el
     call c%discretize_on(el_p%discretization,info)
     call this%discretizations%push(el_p%discretization,info)
     temp => el_p%field
     call this%fields%push(temp,info)
   end select
   call this%changes%push(c,info)
    el => fields%next()
  end do

  allocate(this%rhs_fields,stat=info)
  allocate(this%rhs_discr,stat=info)
  el => rhs_fields%begin()
  do while (associated(el))
    select type(el)
    class is (ppm_t_field)
      el_f => el
      nullify(d)
      call this%rhs_fields%push(el_f,info)
      call this%rhs_discr%push(d,info)
    class is (ppm_t_field_discr_pair)
      el_p => el
      call this%rhs_fields%push(el_p%field,info)
      call this%rhs_discr%push(el_p%discretization,info)
    end select
    el => rhs_fields%next()
  end do
 
  end_subroutine()
end subroutine eulerf_create

subroutine eulerf_destroy(this,info)
  implicit none
  class(ppm_t_eulerf)          :: this
  integer    , intent(  out)   :: info
  start_subroutine("eulerf_destroy")
  
  end_subroutine()
end subroutine eulerf_destroy

subroutine eulerf_step(this,t,dt,info)
  implicit none
  !----------------------------------------------------------------------
  ! Arguments 
  !----------------------------------------------------------------------
  class(ppm_t_eulerf)                        :: this
  real(ppm_kind_double)            ,intent(inout)    :: t
  real(ppm_kind_double)            ,intent(in   )    :: dt
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


  rhs_info = this%rhsfunc(this%rhs_fields,this%rhs_discr,this%changes)

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
              for all
                field_n = field_n + dt*change_n
            end foreach
          else
            foreach n in equi_mesh(mesh) with vec_fields(field,change) indices(i,j)
              for all
                field_n(:) = field_n(:) + dt*change_n(:)
            end foreach
          endif
        else
          if (field%lda.eq.1) then
            foreach n in equi_mesh(mesh) with sca_fields(field,change) indices(i,j,k)
              for all
                field_n = field_n + dt*change_n
            end foreach
          else
            foreach n in equi_mesh(mesh) with vec_fields(field,change) indices(i,j,k)
              for all
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



end module ppm_module_integrator_typedef
