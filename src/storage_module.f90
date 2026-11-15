module storage_module
use parameters_module
use mpi_module
use mesh_module
use allocate_module
use iso_c_binding
implicit none

integer :: i, j, k, iv

integer :: err   = 0
integer :: istop = 0

integer, parameter :: nvAve1D      = 35
integer, parameter :: nvAve2D      = 100
character(3)       :: stat_kind    = 'FFF'

real(rp) :: time
real(rp) :: time_n
real(rp) :: t_restart

real(rp), allocatable, dimension(:,:,:,:) :: phi
real(rp), allocatable, dimension(:,:,:,:) :: phi_n
real(rp), allocatable, dimension(:,:,:,:) :: RHS

real(rp), allocatable, dimension(:,:,:) :: U
real(rp), allocatable, dimension(:,:,:) :: V
real(rp), allocatable, dimension(:,:,:) :: W
real(rp), allocatable, dimension(:,:,:) :: P
real(rp), allocatable, dimension(:,:,:) :: T
real(rp), allocatable, dimension(:,:,:) :: VIS
real(rp), allocatable, dimension(:,:,:) :: LMD
real(rp), allocatable, dimension(:,:,:) :: DIV
real(rp), allocatable, dimension(:,:,:) :: isoSGS

real(rp), allocatable, dimension(:,:,:) :: vfclo_x, vfchi_x, vfcsl_x, vfcsh_x
real(rp), allocatable, dimension(:,:,:) :: vfclo_y, vfchi_y, vfcsl_y, vfcsh_y
real(rp), allocatable, dimension(:,:,:) :: vfclo_z, vfchi_z, vfcsl_z, vfcsh_z

real(rp), allocatable, dimension(:,:,:) :: SSensor
real(rp), allocatable, dimension(:,:,:) :: GAMMA_T
real(rp), allocatable, dimension(:,:,:) :: ENTHALPY

real(rp), allocatable, dimension(:,:) :: pri_1D_x,  pri_1D_y,  pri_1D_z
real(rp), allocatable, dimension(:,:) :: phi_arr_x, phi_arr_y, phi_arr_z
real(rp), allocatable, dimension(:,:) :: flx_arr_x, flx_arr_y, flx_arr_z
real(rp), allocatable, dimension(:,:) :: flx_x,     flx_y,     flx_z

integer(1), allocatable, dimension(:)     :: ishock_x, ishock_y, ishock_z
integer(1), allocatable, dimension(:,:,:) :: weno_flag
integer(1), allocatable, dimension(:,:,:) :: weno_flag_xyz

integer(1), parameter :: weno_smooth = 1
integer(1), parameter :: weno_shock  = 0

real(rp), allocatable, dimension(:,:)   :: vmean1D
real(rp), allocatable, dimension(:,:,:) :: vmean2D

real(rp), allocatable, dimension(:,:,:) :: stat_lcl2D, stat_gbl2D
real(rp), allocatable, dimension(:,:,:) :: stat_pln1D
real(rp), allocatable, dimension(:,:)   :: stat_lcl1D, stat_gbl1D

contains

subroutine init_GL_variables
        implicit none

        time      = 0.0_rp
        t_restart = 0.0_rp

        call AllocateReal(phi  , lbx,ubx, lby,uby, lbz,ubz, 1,5)
        call AllocateReal(phi_n, lbx,ubx, lby,uby, lbz,ubz, 1,5)
        call AllocateReal(RHS  , lbx,ubx, lby,uby, lbz,ubz, 1,5)

        call AllocateReal(U       , lbx,ubx, lby,uby, lbz,ubz)
        call AllocateReal(V       , lbx,ubx, lby,uby, lbz,ubz)
        call AllocateReal(W       , lbx,ubx, lby,uby, lbz,ubz)
        call AllocateReal(P       , lbx,ubx, lby,uby, lbz,ubz)
        call AllocateReal(T       , lbx,ubx, lby,uby, lbz,ubz)
        call AllocateReal(DIV     , lbx,ubx, lby,uby, lbz,ubz)
        call AllocateReal(VIS     , lbx,ubx, lby,uby, lbz,ubz)
        call AllocateReal(LMD     , lbx,ubx, lby,uby, lbz,ubz)
        call AllocateReal(isoSGS  , lbx,ubx, lby,uby, lbz,ubz)
        call AllocateReal(ENTHALPY, lbx,ubx, lby,uby, lbz,ubz)
        call AllocateReal(GAMMA_T , lbx,ubx, lby,uby, lbz,ubz)
        GAMMA_T = gamma0

        call AllocateReal(pri_1D_x , lbx,ubx, 1,6)
        call AllocateReal(pri_1D_y , lby,uby, 1,6)
        call AllocateReal(pri_1D_z , lbz,ubz, 1,6)
        call AllocateReal(phi_arr_x, lbx,ubx, 1,5)
        call AllocateReal(phi_arr_y, lby,uby, 1,5)
        call AllocateReal(phi_arr_z, lbz,ubz, 1,5)
        call AllocateReal(flx_arr_x, lbx,ubx, 1,5)
        call AllocateReal(flx_arr_y, lby,uby, 1,5)
        call AllocateReal(flx_arr_z, lbz,ubz, 1,5)
        call AllocateReal(flx_x    , lbx,ubx, 1,5)
        call AllocateReal(flx_y    , lby,uby, 1,5)
        call AllocateReal(flx_z    , lbz,ubz, 1,5)

        call AllocateInteger(ishock_x, sx-1,ex)
        call AllocateInteger(ishock_y, sy-1,ey)
        call AllocateInteger(ishock_z, sz-1,ez)

        if(hybrid_weno) then
          call AllocateInteger(weno_flag    , lbx,ubx, lby,uby, lbz,ubz)
          call AllocateInteger(weno_flag_xyz, lbx,ubx, lby,uby, lbz,ubz)
          call AllocateReal   (SSensor      , lbx,ubx, lby,uby, lbz,ubz)
          weno_flag = weno_smooth
        endif

        if(trim(diffusion_scheme) == 'conservative') then
          call AllocateReal(vfclo_x, 1,8, sy,ey, sz,ez)
          call AllocateReal(vfchi_x, 1,8, sy,ey, sz,ez)
          call AllocateReal(vfcsl_x, 1,8, sy,ey, sz,ez)
          call AllocateReal(vfcsh_x, 1,8, sy,ey, sz,ez)
          call AllocateReal(vfclo_y, 1,8, sx,ex, sz,ez)
          call AllocateReal(vfchi_y, 1,8, sx,ex, sz,ez)
          call AllocateReal(vfcsl_y, 1,8, sx,ex, sz,ez)
          call AllocateReal(vfcsh_y, 1,8, sx,ex, sz,ez)
          call AllocateReal(vfclo_z, 1,8, sx,ex, sy,ey)
          call AllocateReal(vfchi_z, 1,8, sx,ex, sy,ey)
          call AllocateReal(vfcsl_z, 1,8, sx,ex, sy,ey)
          call AllocateReal(vfcsh_z, 1,8, sx,ex, sy,ey)
        endif

        return
end subroutine init_GL_variables

subroutine init_statistics_fields
        implicit none
        character(1), dimension(3) :: pchar

        if(StOut <= 0) return

        pchar = 'F'
        if(periods(1) .and. grid_invariant_along(yy(sx:ex,sy:ey)           , comm1Dx)) pchar(1) = 'T'
        if(periods(2) .and. grid_invariant_along(transpose(xx(sx:ex,sy:ey)), comm1Dy)) pchar(2) = 'T'
        if(periods(3)) pchar(3) = 'T'
        stat_kind = pchar(1)//pchar(2)//pchar(3)

        if(rank == root) write(*,'(1x,2A)') ' Statistics averaged along (x,y,z) = ', stat_kind

        selectcase(ic)
        case('turbulent_channel')
          if(stat_kind /= 'TFT') call stat_layout_stop('TFT (x and z homogeneous)')
        endselect

        selectcase(stat_kind)

        case('TFF')

        case('FTF')

        case('FFT')
          call AllocateReal(vmean2D   , lbx,ubx, lby,uby, 1,nvAve2D)
          call AllocateReal(stat_lcl2D, lbx,ubx, lby,uby, 1,nvAve2D)
          call AllocateReal(stat_gbl2D, lbx,ubx, lby,uby, 1,nvAve2D)

        case('TTF')

        case('FTT')
          call AllocateReal(vmean2D   , lbx,ubx, lby,uby, 1,nvAve2D)
          call AllocateReal(stat_lcl2D, lbx,ubx, lby,uby, 1,nvAve2D)
          call AllocateReal(stat_gbl2D, lbx,ubx, lby,uby, 1,nvAve2D)

        case('TFT')
          call AllocateReal(vmean1D   , lby,uby, 1,nvAve1D)
          call AllocateReal(stat_pln1D, lbx,ubx, lby,uby, 1,nvAve1D)
          call AllocateReal(stat_lcl1D, lby,uby, 1,nvAve1D)
          call AllocateReal(stat_gbl1D, lby,uby, 1,nvAve1D)

        case default

        endselect

        return
end subroutine init_statistics_fields

logical function grid_invariant_along(c, comm)
        implicit none
        real(rp), dimension(:,:), intent(in) :: c
        integer                 , intent(in) :: comm

        real(rp), dimension(size(c,2)) :: cmin, cmax
        real(rp)                       :: tol, gap
        real(rp), dimension(2)         :: ext

        cmin = minval(c, dim=1)
        cmax = maxval(c, dim=1)
        call MPI_allreduce(MPI_IN_PLACE, cmin, size(cmin), MPI_RP, MPI_MIN, comm, err)
        call MPI_allreduce(MPI_IN_PLACE, cmax, size(cmax), MPI_RP, MPI_MAX, comm, err)

        ext = (/maxval(cmax), -minval(cmin)/)
        call MPI_allreduce(MPI_IN_PLACE, ext, 2, MPI_RP, MPI_MAX, mpi_comm_cart, err)
        tol = 100.0_rp*epsilon(1.0_rp)*(ext(1) + ext(2))
        gap = maxval(cmax - cmin)
        call MPI_allreduce(MPI_IN_PLACE, gap, 1, MPI_RP, MPI_MAX, mpi_comm_cart, err)

        grid_invariant_along = (gap <= tol)

        return
end function grid_invariant_along

subroutine stat_layout_stop(needed)
        implicit none
        character(*), intent(in) :: needed

        if(rank == root) write(*,'(1x,6A)') ' ERROR: ic = ', trim(ic), ' needs statistics ', needed, &
                                             ', this grid gives ', stat_kind
        call secure_stop

end subroutine stat_layout_stop

subroutine init_FD_coefficients
        implicit none

        central_fd_order = fd_order
        bward_fd_order   = fd_order
        fward_fd_order   = fd_order

        selectcase(fd_order)
        case(2)
          fd_L = 0
          fd_R = 1
        case(4)
          fd_L = -1
          fd_R =  2
        case(6)
          fd_L = -2
          fd_R =  3
        case default
          if(rank == root) print*, ' Finite-difference order ', fd_order, ' is not implemented.'
          call secure_stop
        endselect

        call AllocateReal(bward_1       , -6               , 0)
        call AllocateReal(bward_2       , -bward_fd_order-1, 0)
        call AllocateReal(fward_1       , 0                , 6               )
        call AllocateReal(fward_2       , 0                , fward_fd_order+1)
        call AllocateReal(mid_point_lele, fd_L             , fd_R)

        central_1          = 0.0_rp
        central_2          = 0.0_rp
        central_1_one_half = 0.0_rp
        central_2_one_half = 0.0_rp

        visc_stag   = 0.0_rp
        visc_intp   = 0.0_rp
        visc_stag_b = 0.0_rp
        visc_intp_b = 0.0_rp
        if(fd_order >= 6) then

          visc_stag(1) =   75.0_rp/ 64.0_rp
          visc_stag(2) = - 25.0_rp/384.0_rp
          visc_stag(3) =    3.0_rp/640.0_rp

          visc_intp(1) =   75.0_rp/128.0_rp
          visc_intp(2) = - 25.0_rp/256.0_rp
          visc_intp(3) =    3.0_rp/256.0_rp

        elseif(fd_order >= 4) then

          visc_stag(1) =    9.0_rp/  8.0_rp
          visc_stag(2) = -  1.0_rp/ 24.0_rp

          visc_intp(1) =    9.0_rp/ 16.0_rp
          visc_intp(2) = -  1.0_rp/ 16.0_rp

        else

          visc_stag(1) =    1.0_rp
          visc_intp(1) =    0.5_rp

        endif
        visc_divs = visc_stag

        if(fd_order >= 4) then

          visc_stag_b(1) =    9.0_rp/  8.0_rp
          visc_stag_b(2) = -  1.0_rp/ 24.0_rp

          visc_intp_b(1) =    9.0_rp/ 16.0_rp
          visc_intp_b(2) = -  1.0_rp/ 16.0_rp

        else

          visc_stag_b(1) =    1.0_rp
          visc_intp_b(1) =    0.5_rp

        endif

        if(trim(diffusion_scheme) == 'conservative') then
          if(visc_face_reach > GN) then
            if(rank == root) print*, ' ERROR: diffusion_scheme conservative needs GN >=', &
                                     visc_face_reach, ' but GN =', GN
            call secure_stop
          endif
          if(visc_div_reach - (GN - visc_face_reach + 1) > 1) then
            if(rank == root) print*, ' ERROR: viscous_conservative_module trades ONE face plane', &
                                     ' per side, but the node sum overshoots by', &
                                     visc_div_reach - (GN - visc_face_reach + 1)
            call secure_stop
          endif
          if(visc_div_reach + visc_bnd_reach - 1 > GN) then
            if(rank == root) print*, ' ERROR: the closing face of a physical border needs', &
                                     ' Ld+Rb-1 <= GN, but it is', visc_div_reach + visc_bnd_reach - 1
            call secure_stop
          endif
          if(fd_order > 6 .and. rank == root) then
            print*, ' NOTE: diffusion_scheme conservative runs at order 6 (fd_order =', &
                    fd_order, '); the halo GN =', GN, ' caps it.'
          endif
        endif
        fward_1            = 0.0_rp
        bward_1            = 0.0_rp

        selectcase(fd_order)

        case(2)

          central_1(-1) = - 0.5_rp
          central_1( 0) =   0.0_rp
          central_1( 1) = - central_1(-1)

          central_2(-1) =   1.0_rp
          central_2( 0) = - 2.0_rp
          central_2( 1) =   central_2(-1)

          central_1_one_half(0) = 0.5_rp
          central_1_one_half(1) = central_1_one_half(0)

          central_2_one_half(0) = - 1.0_rp
          central_2_one_half(1) = - central_2_one_half(0)

          mid_point_lele(0) = 0.5_rp
          mid_point_lele(1) = mid_point_lele(0)

          bward_1(-2) =   0.5_rp
          bward_1(-1) = - 2.0_rp
          bward_1( 0) =   1.5_rp

          bward_2(-3) = - 1.0_rp
          bward_2(-2) =   4.0_rp
          bward_2(-1) = - 5.0_rp
          bward_2( 0) =   2.0_rp

          fward_1( 0) = - bward_1( 0)
          fward_1( 1) = - bward_1(-1)
          fward_1( 2) = - bward_1(-2)

          fward_2( 0) =   bward_2( 0)
          fward_2( 1) =   bward_2(-1)
          fward_2( 2) =   bward_2(-2)
          fward_2( 3) =   bward_2(-3)

        case(4)

          central_1(-2) =   1.0_rp/12.0_rp
          central_1(-1) = - 2.0_rp/ 3.0_rp
          central_1( 0) =   0.0_rp
          central_1( 1) = - central_1(-1)
          central_1( 2) = - central_1(-2)

          central_2(-2) = - 1.0_rp/12.0_rp
          central_2(-1) =   4.0_rp/ 3.0_rp
          central_2( 0) = - 5.0_rp/ 2.0_rp
          central_2( 1) =   central_2(-1)
          central_2( 2) =   central_2(-2)

          central_1_one_half(-1) = - 1.0_rp/12.0_rp
          central_1_one_half( 0) =   7.0_rp/12.0_rp
          central_1_one_half( 1) =   central_1_one_half( 0)
          central_1_one_half( 2) =   central_1_one_half(-1)

          central_2_one_half(-1) =   1.0_rp/12.0_rp
          central_2_one_half( 0) = - 5.0_rp/ 4.0_rp
          central_2_one_half( 1) = - central_2_one_half( 0)
          central_2_one_half( 2) = - central_2_one_half(-1)

          mid_point_lele(-1) = - 1.0_rp/16.0_rp
          mid_point_lele( 0) =   9.0_rp/16.0_rp
          mid_point_lele( 1) =   mid_point_lele( 0)
          mid_point_lele( 2) =   mid_point_lele(-1)

          bward_1(-4) =    3.0_rp/12.0_rp
          bward_1(-3) = - 16.0_rp/12.0_rp
          bward_1(-2) =   36.0_rp/12.0_rp
          bward_1(-1) = - 48.0_rp/12.0_rp
          bward_1( 0) =   25.0_rp/12.0_rp

          bward_2(-5) = - 10.0_rp/12.0_rp
          bward_2(-4) =   61.0_rp/12.0_rp
          bward_2(-3) = -156.0_rp/12.0_rp
          bward_2(-2) =  214.0_rp/12.0_rp
          bward_2(-1) = -154.0_rp/12.0_rp
          bward_2( 0) =   45.0_rp/12.0_rp

          fward_1( 0) = - bward_1( 0)
          fward_1( 1) = - bward_1(-1)
          fward_1( 2) = - bward_1(-2)
          fward_1( 3) = - bward_1(-3)
          fward_1( 4) = - bward_1(-4)

          fward_2( 0) =   bward_2( 0)
          fward_2( 1) =   bward_2(-1)
          fward_2( 2) =   bward_2(-2)
          fward_2( 3) =   bward_2(-3)
          fward_2( 4) =   bward_2(-4)
          fward_2( 5) =   bward_2(-5)

        case(6)

          central_1(-3) = - 1.0_rp/60.0_rp
          central_1(-2) =   3.0_rp/20.0_rp
          central_1(-1) = - 3.0_rp/ 4.0_rp
          central_1( 0) =   0.0_rp
          central_1( 1) = - central_1(-1)
          central_1( 2) = - central_1(-2)
          central_1( 3) = - central_1(-3)

          central_2(-3) =    1.0_rp/90.0_rp
          central_2(-2) = -  3.0_rp/20.0_rp
          central_2(-1) =    3.0_rp/ 2.0_rp
          central_2( 0) = - 49.0_rp/18.0_rp
          central_2( 1) =   central_2(-1)
          central_2( 2) =   central_2(-2)
          central_2( 3) =   central_2(-3)

          central_1_one_half(-2) =   1.0_rp/60.0_rp
          central_1_one_half(-1) = - 2.0_rp/15.0_rp
          central_1_one_half( 0) =  37.0_rp/60.0_rp
          central_1_one_half( 1) =   central_1_one_half( 0)
          central_1_one_half( 2) =   central_1_one_half(-1)
          central_1_one_half( 3) =   central_1_one_half(-2)

          central_2_one_half(-2) = -  1.0_rp/90.0_rp
          central_2_one_half(-1) =    5.0_rp/36.0_rp
          central_2_one_half( 0) = - 49.0_rp/36.0_rp
          central_2_one_half( 1) = - central_2_one_half( 0)
          central_2_one_half( 2) = - central_2_one_half(-1)
          central_2_one_half( 3) = - central_2_one_half(-2)

          mid_point_lele(-2) =   3.0_rp/256.0_rp
          mid_point_lele(-1) = -25.0_rp/256.0_rp
          mid_point_lele( 0) =  75.0_rp/128.0_rp
          mid_point_lele( 1) =  mid_point_lele( 0)
          mid_point_lele( 2) =  mid_point_lele(-1)
          mid_point_lele( 3) =  mid_point_lele(-2)

          bward_1(-6) =    1.0_rp/  6.0_rp
          bward_1(-5) = -  6.0_rp/  5.0_rp
          bward_1(-4) =   15.0_rp/  4.0_rp
          bward_1(-3) = - 20.0_rp/  3.0_rp
          bward_1(-2) =   15.0_rp/  2.0_rp
          bward_1(-1) = -  6.0_rp
          bward_1( 0) =   49.0_rp/ 20.0_rp

          bward_2(-7) = -  7.0_rp/ 10.0_rp
          bward_2(-6) = 1019.0_rp/180.0_rp
          bward_2(-5) = -201.0_rp/ 10.0_rp
          bward_2(-4) =   41.0_rp
          bward_2(-3) = -949.0_rp/ 18.0_rp
          bward_2(-2) =  879.0_rp/ 20.0_rp
          bward_2(-1) = -223.0_rp/ 10.0_rp
          bward_2( 0) =  469.0_rp/ 90.0_rp

          fward_1( 0) = - bward_1( 0)
          fward_1( 1) = - bward_1(-1)
          fward_1( 2) = - bward_1(-2)
          fward_1( 3) = - bward_1(-3)
          fward_1( 4) = - bward_1(-4)
          fward_1( 5) = - bward_1(-5)
          fward_1( 6) = - bward_1(-6)

          fward_2( 0) =   bward_2( 0)
          fward_2( 1) =   bward_2(-1)
          fward_2( 2) =   bward_2(-2)
          fward_2( 3) =   bward_2(-3)
          fward_2( 4) =   bward_2(-4)
          fward_2( 5) =   bward_2(-5)
          fward_2( 6) =   bward_2(-6)
          fward_2( 7) =   bward_2(-7)

        endselect

        c1_wbias = 0.0_rp
        c1_wbias(0,0:6) = fward_1(0:6)

        selectcase(fd_order)
        case(2)
          c1_wbias(1,0) = central_1(-1); c1_wbias(1,1) = central_1(0); c1_wbias(1,2) = central_1(1)
          c1_wbias(2,1) = central_1(-1); c1_wbias(2,2) = central_1(0); c1_wbias(2,3) = central_1(1)

        case(4)
          c1_wbias(1,0) = - 1.0_rp/ 4.0_rp
          c1_wbias(1,1) = - 5.0_rp/ 6.0_rp
          c1_wbias(1,2) =   3.0_rp/ 2.0_rp
          c1_wbias(1,3) = - 1.0_rp/ 2.0_rp
          c1_wbias(1,4) =   1.0_rp/12.0_rp
          c1_wbias(2,0) = central_1(-2); c1_wbias(2,1) = central_1(-1); c1_wbias(2,2) = central_1(0)
          c1_wbias(2,3) = central_1( 1); c1_wbias(2,4) = central_1( 2)

        case(6)
          c1_wbias(1,0) = - 1.0_rp/ 6.0_rp
          c1_wbias(1,1) = -77.0_rp/60.0_rp
          c1_wbias(1,2) =   5.0_rp/ 2.0_rp
          c1_wbias(1,3) = - 5.0_rp/ 3.0_rp
          c1_wbias(1,4) =   5.0_rp/ 6.0_rp
          c1_wbias(1,5) = - 1.0_rp/ 4.0_rp
          c1_wbias(1,6) =   1.0_rp/30.0_rp

          c1_wbias(2,0) =   1.0_rp/30.0_rp
          c1_wbias(2,1) = - 2.0_rp/ 5.0_rp
          c1_wbias(2,2) = - 7.0_rp/12.0_rp
          c1_wbias(2,3) =   4.0_rp/ 3.0_rp
          c1_wbias(2,4) = - 1.0_rp/ 2.0_rp
          c1_wbias(2,5) =   2.0_rp/15.0_rp
          c1_wbias(2,6) = - 1.0_rp/60.0_rp
        endselect

        c1_ebias = - c1_wbias(:,6:0:-1)

        if(any(abs(sum(c1_wbias,dim=2)) > toll_equality)) stop ' c1_wbias gradient stencils are unbalanced'
        if(any(abs(sum(c1_ebias,dim=2)) > toll_equality)) stop ' c1_ebias gradient stencils are unbalanced'

        if(abs(sum(central_1))                 > toll_equality) stop ' central_1 FD coefficients are unbalanced'
        if(abs(sum(central_2))                 > toll_equality) stop ' central_2 FD coefficients are unbalanced'
        if(abs(sum(central_1_one_half)-1.0_rp) > toll_equality) stop ' central_1_one_half are unbalanced'
        if(abs(sum(central_2_one_half))        > toll_equality) stop ' central_2_one_half are unbalanced'
        if(abs(sum(bward_1))                   > toll_equality) stop ' bward_1 FD coefficients are unbalanced'
        if(abs(sum(bward_2))                   > toll_equality) stop ' bward_2 FD coefficients are unbalanced'
        if(abs(sum(fward_1))                   > toll_equality) stop ' fward_1 FD coefficients are unbalanced'
        if(abs(sum(fward_2))                   > toll_equality) stop ' fward_2 FD coefficients are unbalanced'
        if(abs(sum(mid_point_lele)-1.0_rp)     > toll_equality) stop ' mid_point_lele are unbalanced'

        return
end subroutine init_FD_coefficients

subroutine init_RK_coefficients
        implicit none

        call AllocateReal(a_rk, 1,n_step)
        call AllocateReal(b_rk, 1,n_step)
        call AllocateReal(c_rk, 1,n_step)

        selectcase(n_step)
        case(3)
          a_rk = (/ 0.0_rp, 0.75_rp, 1.0_rp/3.0_rp /)
          b_rk = (/ 1.0_rp, 0.25_rp, 2.0_rp/3.0_rp /)
          c_rk = (/ 1.0_rp, 0.50_rp, 1.0_rp        /)
        case default
          if(rank == root) print*, ' Runge-Kutta coefficients for ', n_step, ' steps are not implemented.'
          call secure_stop
        endselect

        return
end subroutine init_RK_coefficients

subroutine end_all_variables
        implicit none

        call DeallocateReal(bward_1)
        call DeallocateReal(bward_2)
        call DeallocateReal(fward_1)
        call DeallocateReal(fward_2)
        call DeallocateReal(mid_point_lele)
        call DeallocateReal(a_rk)
        call DeallocateReal(b_rk)
        call DeallocateReal(c_rk)

        call DeallocateReal(phi)
        call DeallocateReal(phi_n)
        call DeallocateReal(RHS)

        call DeallocateReal(U)
        call DeallocateReal(V)
        call DeallocateReal(W)
        call DeallocateReal(P)
        call DeallocateReal(T)
        call DeallocateReal(VIS)
        call DeallocateReal(LMD)
        call DeallocateReal(DIV)
        call DeallocateReal(isoSGS)
        call DeallocateReal(SSensor)

        call DeallocateReal(vfclo_x)
        call DeallocateReal(vfchi_x)
        call DeallocateReal(vfcsl_x)
        call DeallocateReal(vfcsh_x)
        call DeallocateReal(vfclo_y)
        call DeallocateReal(vfchi_y)
        call DeallocateReal(vfcsl_y)
        call DeallocateReal(vfcsh_y)
        call DeallocateReal(vfclo_z)
        call DeallocateReal(vfchi_z)
        call DeallocateReal(vfcsl_z)
        call DeallocateReal(vfcsh_z)
        call DeallocateReal(GAMMA_T)
        call DeallocateReal(ENTHALPY)

        call DeallocateReal(pri_1D_x);  call DeallocateReal(pri_1D_y);  call DeallocateReal(pri_1D_z)
        call DeallocateReal(phi_arr_x); call DeallocateReal(phi_arr_y); call DeallocateReal(phi_arr_z)
        call DeallocateReal(flx_arr_x); call DeallocateReal(flx_arr_y); call DeallocateReal(flx_arr_z)
        call DeallocateReal(flx_x);     call DeallocateReal(flx_y);     call DeallocateReal(flx_z)

        call DeallocateInteger(ishock_x)
        call DeallocateInteger(ishock_y)
        call DeallocateInteger(ishock_z)
        call DeallocateInteger(weno_flag)
        call DeallocateInteger(weno_flag_xyz)

        call DeallocateReal(vmean1D)
        call DeallocateReal(vmean2D)
        call DeallocateReal(stat_lcl2D)
        call DeallocateReal(stat_gbl2D)
        call DeallocateReal(stat_pln1D)
        call DeallocateReal(stat_lcl1D)
        call DeallocateReal(stat_gbl1D)

        call destroy_grid

        return
end subroutine end_all_variables

end module storage_module
