module post_storage_module
use parameters_module
use storage_module
use mesh_module
use FileModule
implicit none

integer       :: f         !< file index
integer       :: file_ini  !< first file  
integer       :: file_end  !< last file
character(dl) :: ifile

! file timing
real(rp) :: time_restart        !< restarting time
real(rp) :: time_new            !< time of file f
real(rp) :: time_old            !< time of file f-1
real(rp) :: file_dt             !< time step between file f and f-1

! periodic directions
character(1), dimension(3) :: periodic_dir
character(1), dimension(3) :: symmtric_dir

integer, dimension(2), parameter :: gnu_unit = (/15,16/)

type exact_field_type
  real(rp), dimension(:,:,:), allocatable :: rho
  real(rp), dimension(:,:,:), allocatable :: vel
  real(rp), dimension(:,:,:), allocatable :: tmp
end type exact_field_type

type error_field_type
  real(rp), dimension(:,:,:), allocatable :: rho
  real(rp), dimension(:,:,:), allocatable :: vel
  real(rp), dimension(:,:,:), allocatable :: tmp
end type error_field_type

type Rey_average_type
  real(rp), dimension(:,:,:), allocatable :: rho
  real(rp), dimension(:,:,:), allocatable :: rhu
  real(rp), dimension(:,:,:), allocatable :: rhv
  real(rp), dimension(:,:,:), allocatable :: rhw
  real(rp), dimension(:,:,:), allocatable :: ruu
  real(rp), dimension(:,:,:), allocatable :: rvv
  real(rp), dimension(:,:,:), allocatable :: rww
  real(rp), dimension(:,:,:), allocatable :: ruv
  real(rp), dimension(:,:,:), allocatable :: prs
  real(rp), dimension(:,:,:), allocatable :: Tmp
  real(rp), dimension(:,:,:), allocatable :: Vis
  real(rp), dimension(:,:,:), allocatable :: Mac
end type Rey_average_type

type bulk_type
  real(rp) :: rho
  real(rp) :: rhou
  real(rp) :: rhov
  real(rp) :: rhow
  real(rp) :: u
  real(rp) :: v
  real(rp) :: w
  real(rp) :: Ma
  real(rp) :: Re
end type bulk_type

real(rp), dimension(:,:,:), allocatable :: corr_v

! commit derived data 
type(DirType)          :: main_dir
type(exact_field_type) :: exact
type(error_field_type) :: error
type(Rey_average_type) :: Re_av

! variable that can be plotted
real(rp), dimension(:,:,:), allocatable :: Mach_x, Mach_y, Mach_z       !< mach number
real(rp), dimension(:,:,:), allocatable :: vor_x, vor_y, vor_z          !< vorticity
real(rp), dimension(:,:,:), allocatable :: L1, L2, L3, L4, L5           !< charactheristic waves
real(rp), dimension(:,:,:), allocatable :: s_density                    !< shielereen density
real(rp), dimension(:,:,:), allocatable :: abs_vor                      !< vorticity magnitude
real(rp), dimension(:,:,:), allocatable :: q_criterion                  !< Q criterion for vortices
real(rp), dimension(:,:,:), allocatable :: UV, UW                       !< mixed velocities
real(rp), dimension(:,:,:), allocatable :: wrk                          !< work array


contains

subroutine look_for_symmetry
        implicit none
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        
        ! PERIODICAL DIRECTIONS ==============================
        periodic_dir(:) = 'F'

        ! look for periodicity in x-direction
        if(bc(E) == bc(W) .and. bc(E) == 'periodic') then
                periodic_dir(1) = 'T'
        endif

        ! look for periodicity in y-direction
        if(bc(S) == bc(N) .and. bc(S) == 'periodic') then
                periodic_dir(2) = 'T'
        endif

        ! look for periodicity in z-direction
        if(bc(B) == bc(F) .and. bc(B) == 'periodic') then
                periodic_dir(3) = 'T'
        endif

        ! SYMMETRICAL DIRECTIONS =============================
        symmtric_dir(:) = 'F'

        ! look for symmetry in x-direction
        if(bc(E) == bc(W) .and. bc(E) /= 'periodic') then
                symmtric_dir(1) = 'T'
        endif

        ! look for symmetry in y-direction
        if(bc(S) == bc(N) .and. bc(S) /= 'periodic') then
                symmtric_dir(2) = 'T'
        endif

        ! look for symmetry in z-direction
        if(bc(B) == bc(F) .and. bc(B) /= 'periodic') then
                symmtric_dir(3) = 'T'
        endif

        return
end subroutine look_for_symmetry


subroutine init_cons_fields
! ----------------------------------------------------------------------------------
!       This subroutine allocates the conservative variable and set the problem size
! ----------------------------------------------------------------------------------
        implicit none

        sx = 1 ; ex = nx
        sy = 1 ; ey = ny

        if(dims == 2) nz = 1   !< improving speed
        sz = 1 ; ez = nz

        lbx = sx - GN ; ubx = ex + GN
        lby = sy - GN ; uby = ey + GN
        lbz = sz - GN ; ubz = ez + GN

        ! conservative variables 
        allocate(phi(lbx:ubx, lby:uby, lbz:ubz,eqs), STAT=err) ; if(err.ne.0) stop 'allocation error' ; phi = 0._rp

        return
end subroutine init_cons_fields


subroutine init_inst_fields
! ----------------------------------------------------------------------------------------------------------
!       This subroutine allocates the fields we want to post-treat
! ----------------------------------------------------------------------------------------------------------
        use allocate_module
        implicit none

        allocate(c_rk(3))
        c_rk = (/0.0_rp, 0.50_rp, 0.50_rp/)

        if(velocity) then
          call AllocateReal(U,lbx,ubx,lby,uby,lbz,ubz)
          call AllocateReal(V,lbx,ubx,lby,uby,lbz,ubz)
          call AllocateReal(W,lbx,ubx,lby,uby,lbz,ubz)
        endif
        
        if(pressure)    call AllocateReal(P,lbx,ubx,lby,uby,lbz,ubz)
        if(temperature) call AllocateReal(T,lbx,ubx,lby,uby,lbz,ubz)
        
        if(mach_) then
          call AllocateReal(Mach_x,lbx,ubx,lby,uby,lbz,ubz)
          call AllocateReal(Mach_y,lbx,ubx,lby,uby,lbz,ubz)
          call AllocateReal(Mach_z,lbx,ubx,lby,uby,lbz,ubz)
        endif

        if(vorticity) then
          call AllocateReal(vor_x,lbx,ubx,lby,uby,lbz,ubz)
          call AllocateReal(vor_y,lbx,ubx,lby,uby,lbz,ubz)
          call AllocateReal(vor_z,lbx,ubx,lby,uby,lbz,ubz)
        endif

        if(vorticity_magnitude) call AllocateReal(abs_vor,lbx,ubx,lby,uby,lbz,ubz)
        if(speed_div)           call AllocateReal(div,lbx,ubx,lby,uby,lbz,ubz)
        
        if(charactheristic) then
          call AllocateReal(L1,lbx,ubx,lby,uby,lbz,ubz)
          call AllocateReal(L2,lbx,ubx,lby,uby,lbz,ubz)
          call AllocateReal(L3,lbx,ubx,lby,uby,lbz,ubz)
          call AllocateReal(L4,lbx,ubx,lby,uby,lbz,ubz)
          call AllocateReal(L5,lbx,ubx,lby,uby,lbz,ubz)
        endif

        if(hybrid_weno) then
          call AllocateReal(SSENSOR,lbx,ubx,lby,uby,lbz,ubz)
          call AllocateInteger(weno%flag,lbx,ubx,lby,uby,lbz,ubz); weno%flag = weno%smooth
        endif

        if(sdensity)   call AllocateReal(s_density,lbx,ubx,lby,uby,lbz,ubz)
        if(qcriterion) call AllocateReal(q_criterion,lbx,ubx,lby,uby,lbz,ubz)
        
        if(viscous) then
          call AllocateReal(VIS,lbx,ubx,lby,uby,lbz,ubz)
          call AllocateReal(LMD,lbx,ubx,lby,uby,lbz,ubz)
        endif

        if(MixedVelocity) then
          call AllocateReal(UV,lbx,ubx,lby,uby,lbz,ubz)
          call AllocateReal(UW,lbx,ubx,lby,uby,lbz,ubz)
        endif

        if(ic == 'turbulent_channel') then
          allocate(corr_v(lby:uby,lbz:ubz,3)) ; if(err.ne.0) stop 'allocation error'; corr_v = 0.0_rp
        endif

        call AllocateReal(WMLES_DATA_LW,lbx,ubx,lbz,ubz,1,nvWmlesData)

        return
end subroutine init_inst_fields


subroutine init_exct_fields
! -------------------------------------------------------------------------------------
! allocatation of exact and error field for problems that admits analitycal solutions
! -------------------------------------------------------------------------------------
        implicit none

        select case(ic)

          case('periodic_euler_x'       , &
               'periodic_euler_y'       , &
               'periodic_euler_z'       , &
               'linear_ode'             , &
               'linear_advection'       , &
               'isentropic_vortex_x'    , &
               'shock_tube_x'           , & 
               'nscbc_perturbation')

          
                allocate(exact%rho(sx:ex,sy:ey,sz:ez), &
                         error%rho(sx:ex,sy:ey,sz:ez), stat=err)
                
                if(err.ne.0) stop 'allocation error' 

                exact%rho = 0.0_rp
                error%rho = 0.0_rp

          case('couette_x'          , 'couette_y'          , 'couette_z'          , &
               'I_stokes_problem_x' , 'I_stokes_problem_y' , 'I_stokes_problem_z' , &
               'II_stokes_problem_x', 'II_stokes_problem_y', 'II_stokes_problem_z', &
               'steady_couette_x'   , 'steady_couette_y'   , 'steady_couette_z'   , &
               'KolmogorovFlow')

                allocate(exact%vel(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1), &
                         error%vel(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1), stat=err) 

                if(err.ne.0) stop 'allocation error'

                exact%vel = 0.0_rp
                error%vel = 0.0_rp

        case('poiseuille_x', 'poiseuille_y', 'poiseuille_z', 'inflow_poiseuille')

                allocate(exact%vel(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1), &
                         exact%tmp(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1), &

                         error%vel(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1), &
                         error%tmp(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1), stat = err)

                if(err.ne.0) stop 'allocation error' 
                
                exact%vel = 0.0_rp
                exact%tmp = 0.0_rp
                error%vel = 0.0_rp
                error%tmp = 0.0_rp

        endselect

        
        return
end subroutine init_exct_fields


subroutine init_stat_fields
   implicit none
   integer :: err = 0

   if(Rey_average) then
   
     allocate(Re_av%rho(lbx:ubx,lby:uby,lbz:ubz), stat = err)
     if(err.ne.0) stop ' Allocation error'
     allocate(Re_av%rhu(lbx:ubx,lby:uby,lbz:ubz), stat = err)
     if(err.ne.0) stop ' Allocation error'
     allocate(Re_av%rhv(lbx:ubx,lby:uby,lbz:ubz), stat = err)
     if(err.ne.0) stop ' Allocation error'
     allocate(Re_av%rhw(lbx:ubx,lby:uby,lbz:ubz), stat = err)
     if(err.ne.0) stop ' Allocation error'
     allocate(Re_av%ruu(lbx:ubx,lby:uby,lbz:ubz), stat = err)
     if(err.ne.0) stop ' Allocation error'
     allocate(Re_av%rvv(lbx:ubx,lby:uby,lbz:ubz), stat = err)
     if(err.ne.0) stop ' Allocation error'
     allocate(Re_av%rww(lbx:ubx,lby:uby,lbz:ubz), stat = err)
     if(err.ne.0) stop ' Allocation error'
     allocate(Re_av%ruv(lbx:ubx,lby:uby,lbz:ubz), stat = err)
     if(err.ne.0) stop ' Allocation error'
     allocate(Re_av%prs(lbx:ubx,lby:uby,lbz:ubz), stat = err)
     if(err.ne.0) stop ' Allocation error'
     allocate(Re_av%Tmp(lbx:ubx,lby:uby,lbz:ubz), stat = err)
     if(err.ne.0) stop ' Allocation error'
     allocate(Re_av%Mac(lbx:ubx,lby:uby,lbz:ubz), stat = err)
     if(err.ne.0) stop ' Allocation error'

     Re_av%rho = 0.0_rp
     Re_av%rhu = 0.0_rp
     Re_av%rhv = 0.0_rp
     Re_av%rhw = 0.0_rp
     Re_av%ruu = 0.0_rp
     Re_av%rvv = 0.0_rp
     Re_av%rww = 0.0_rp
     Re_av%ruv = 0.0_rp
     Re_av%prs = 0.0_rp
     Re_av%Tmp = 0.0_rp
     Re_av%Mac = 0.0_rp

     if(viscous) then
       allocate(Re_av%vis(lbx:ubx,lby:uby,lbz:ubz), stat = err) 
       if (err.ne.0) stop ' Allocation err in statistics'
       Re_av%vis = 0.0_rp
     endif

   
   endif

   return
end subroutine init_stat_fields


subroutine destroy_variables
        
        if(allocated(phi)) deallocate(phi)

        if(allocated(U)) deallocate(U)
        if(allocated(V)) deallocate(V)
        if(allocated(W)) deallocate(W)

        if(allocated(P)) deallocate(P)
        if(allocated(T)) deallocate(T)

        if(allocated(Mach_x)) deallocate(Mach_x)
        if(allocated(Mach_y)) deallocate(Mach_y)
        if(allocated(Mach_z)) deallocate(Mach_z)

        if(allocated(vor_x)) deallocate(vor_x)
        if(allocated(vor_y)) deallocate(vor_y)
        if(allocated(vor_z)) deallocate(vor_z)

        if(allocated(abs_vor)) deallocate(abs_vor)

        if(allocated(div)) deallocate(div)

        call DeallocateReal(SSENSOR)

        if(allocated(weno%flag)) deallocate(weno%flag)
        if(allocated(weno%temp)) deallocate(weno%temp)

        if(allocated(s_density)) deallocate(s_density)
        if(allocated(q_criterion)) deallocate(q_criterion)
        
        ! deallocate Rey av
        if(allocated(Re_av%rho)) deallocate(Re_av%rho)
        if(allocated(Re_av%rhu)) deallocate(Re_av%rhu)
        if(allocated(Re_av%rhv)) deallocate(Re_av%rhv)
        if(allocated(Re_av%rhw)) deallocate(Re_av%rhw)
        if(allocated(Re_av%ruu)) deallocate(Re_av%ruu)
        if(allocated(Re_av%rvv)) deallocate(Re_av%rvv)
        if(allocated(Re_av%rww)) deallocate(Re_av%rww)
        if(allocated(Re_av%ruv)) deallocate(Re_av%ruv)
        if(allocated(Re_av%prs)) deallocate(Re_av%prs)
        if(allocated(Re_av%Tmp)) deallocate(Re_av%Tmp)
        if(allocated(Re_av%vis)) deallocate(Re_av%vis)

end subroutine destroy_variables





end module post_storage_module

