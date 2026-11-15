module parameters_module

        use mpi

        implicit none

        integer , parameter :: sp  = selected_real_kind(6,37)
        integer , parameter :: dp  = selected_real_kind(15,307)
#ifdef SINGLE_PRECISION
        integer , parameter :: rp     = sp
        integer , parameter :: MPI_RP = MPI_SINGLE_PRECISION
#else
        integer , parameter :: rp     = dp
        integer , parameter :: MPI_RP = MPI_DOUBLE_PRECISION
#endif

        integer , parameter :: msm = MPI_SUM
        integer , parameter :: ip  = 4
        integer , parameter :: dl  = 600
        real(rp), parameter :: pi  = 2.0_rp*asin(1.0_rp)

        real(rp), parameter :: gamma0 = 1.4_rp
        real(rp), parameter :: cv   = 1._rp /(gamma0-1.0_rp)
        real(rp), parameter :: cp   = gamma0/(gamma0-1.0_rp)
        real(rp), parameter :: cv_i = 1._rp/cv
        real(rp), parameter :: cp_i = 1._rp/cp
        real(rp), parameter :: toll_equality = 1.0E-13_rp

        integer             :: thermally_perfect_gas = 0

        integer :: mpi_opt_level = 2
        logical :: cuda_aware = .true.

        integer, dimension(3) :: cart_dims
        real(rp),dimension(3) :: gbl_min_step
        real(rp)              :: xmin, xmax, Lx
        real(rp)              :: ymin, ymax, Ly
        real(rp)              :: zmin, zmax, Lz
        integer               :: nx, ny, nz
        integer               :: dims
        character(dl)         :: gridpoint_x
        character(dl)         :: gridpoint_y
        character(dl)         :: gridpoint_z
        character(dl)         :: mesh_file = ''
        real(rp)              :: str_it_time, end_it_time

        real(rp) :: turbulent_intensity = 0.03_rp

        real(rp)   :: tmax, Dt, CFL, time_from_restart

        real(rp)   :: Fourier = 0.1_rp
        integer    :: it,itmax, itStat = 0
        integer    :: n_step, ik = 1, itout, rsOut = 10
        integer    :: StOut = 0
        logical    :: printres = .true.
        integer    :: itPrb = 0
        integer    :: itWallMap = 0
        logical    :: restart_flag  = .false.
        logical    :: inflow        = .true.
        logical    :: smooth_inflow = .false.
        logical    :: turb_inflow   = .false.
        logical    :: StFlg         = .true.
        logical    :: restartStat   = .false.
        real(rp)   :: s_cpu_time, e_cpu_time

        integer :: central_fd_order
        integer :: bward_fd_order
        integer :: fward_fd_order
        integer :: weno_num = 3
        integer :: weno_order = 5
        integer :: fd_L, fd_R

        real(rp), dimension(-3:3) :: central_1
        real(rp), dimension(-3:3) :: central_2
        real(rp), dimension(-2:3) :: central_1_one_half
        real(rp), dimension(-2:3) :: central_2_one_half
        integer , parameter       :: visc_div_reach  = 3
        integer , parameter       :: visc_face_reach = 3
        integer , parameter       :: visc_bnd_reach  = 2
        real(rp), dimension(1:visc_div_reach ) :: visc_divs
        real(rp), dimension(1:visc_face_reach) :: visc_stag
        real(rp), dimension(1:visc_face_reach) :: visc_intp
        real(rp), dimension(1:visc_face_reach) :: visc_stag_b
        real(rp), dimension(1:visc_face_reach) :: visc_intp_b
        real(rp), dimension(:), allocatable :: bward_1  , bward_2
        real(rp), dimension(:), allocatable :: fward_1  , fward_2

        real(rp), dimension(0:2,0:6) :: c1_wbias, c1_ebias

        real(rp), dimension(:),   allocatable :: mid_point_lele
        real(rp), dimension(:,:), allocatable :: mid_point_lele_x
        real(rp), dimension(:,:), allocatable :: mid_point_lele_y
        real(rp), dimension(:,:), allocatable :: mid_point_lele_z

        real(rp), dimension(0:1,0:1), parameter :: aweno3 = reshape( &
                        1.0_rp/2.0_rp * [ &
                                          - 1.0_rp ,   3.0_rp, &
                                            1.0_rp ,   1.0_rp  &
                                        ], shape(aweno3), order=[2,1])
        real(rp), dimension(0:1)    , parameter :: cweno3 = &
                        1.0_rp/3.0_rp * [1.0_rp, 2.0_rp]

        real(rp), dimension(0:2,0:2), parameter :: aweno5 = reshape( &
                        1.0_rp/6.0_rp * [ &
                                            2.0_rp , - 7.0_rp,  11.0_rp, &
                                          - 1.0_rp ,   5.0_rp,   2.0_rp, &
                                            2.0_rp ,   5.0_rp, - 1.0_rp  &
                                        ], shape(aweno5), order=[2,1])
        real(rp), dimension(0:2)    , parameter :: cweno5 = &
                        1.0_rp/10.0_rp * [1.0_rp, 6.0_rp, 3.0_rp]

        real(rp), dimension(0:3,0:3), parameter :: aweno7 = reshape( &
                        1.0_rp/12.0_rp * [ &
                                          - 3.0_rp ,  13.0_rp, -23.0_rp,  25.0_rp, &
                                            1.0_rp ,  -5.0_rp,  13.0_rp,   3.0_rp, &
                                          - 1.0_rp ,   7.0_rp,   7.0_rp, - 1.0_rp, &
                                            3.0_rp ,  13.0_rp, - 5.0_rp,   1.0_rp  &
                                        ], shape(aweno7), order=[2,1])
        real(rp), dimension(0:3)    , parameter :: cweno7 = &
                        1.0_rp/35.0_rp * [ 1.0_rp, 12.0_rp, 18.0_rp,  4.0_rp]

        type face_type
          integer      :: node
          integer      :: norm
          character(1) :: face
        endtype face_type

        real(rp), dimension(:), allocatable :: a_rk, b_rk, c_rk

        character(50)               :: ic
        character(50), dimension(6) :: bc
        logical      , dimension(6) :: extrapolation_bc = .false.
        real(rp)                    :: Trat  = 1.0_rp
        real(rp)                    :: Tref  = 273.15_rp
        character(10)               :: species_name  = "air"
        integer                     :: vis_flag = 1

        !$acc declare create(Tref)

        real(rp) :: Reynolds, Mach, Prandtl, mu_inf, u_inf, k_inf, iPrandtl, q_inf, ReTau
        real(rp) :: T_recovery
        real(rp) :: force_turbulent_channel

        character(dl) :: output_file_name, data_dir, inflow_profile
        character(dl) :: restart_file

        logical       :: LES   = .false.
        logical       :: yoshizawa = .true.
        character(40) :: sgs_model

        logical       :: logical_CFL = .false.
        logical       :: hybrid_weno = .false.
        character(50) :: diffusion_scheme = 'staggered'
        real(rp)      :: kolmogorov_forcing = 0.1_rp
        real(rp)      :: mu_jump_ratio = 1.0_rp
        real(rp)      :: mu_jump_loc   = 0.25_rp
        character(50) :: scheme
        character(20) :: sensor
        real(rp)      :: ducrosToll = 0.1_rp
        integer       :: fd_order

        logical       :: enable_limiter    = .false.
        real(rp)      :: limiter_rho_floor = 1.0e-3_rp
        real(rp)      :: limiter_e_floor   = 1.0e-4_rp
        real(rp)      :: limiter_abort_frac = 1.0e-4_rp
        character(10) :: ep_type = 'JP'
        integer(1)    :: ep_id   =   1

        logical :: viscous = .false.

        integer , parameter :: prb_max = 16
        integer  :: prb_nwall  = 0, prb_nfield = 0
        real(rp) :: prb_xc(prb_max)   = 0.0_rp
        integer  :: prb_side(prb_max) = 0
        real(rp) :: prb_x(prb_max)    = huge(1.0_rp)
        real(rp) :: prb_y(prb_max)    = huge(1.0_rp)

        namelist /input_list/              &
              &   xmin,                &
              &   xmax,                &
              &   ymin,                &
              &   ymax,                &
              &   zmin,                &
              &   zmax,                &
              &   gridpoint_x,         &
              &   gridpoint_y,         &
              &   gridpoint_z,         &
              &   mesh_file,           &
              &   dims,                &
              &   tmax,                &
              &   itmax,               &
              &   nx,                  &
              &   ny,                  &
              &   nz,                  &
              &   cart_dims,           &
              &   mpi_opt_level,       &
              &   cuda_aware,          &
              &   CFL,                 &
              &   Fourier,             &
              &   Dt,                  &
              &   bc,                  &
              &   Trat,                &
              &   Tref,                &
              &   species_name,        &
              &   vis_flag,            &
              &   thermally_perfect_gas,   &
              &   inflow_profile,      &
              &   smooth_inflow,       &
              &   turb_inflow,         &
              &   turbulent_intensity, &
              &   Reynolds,            &
              &   Mach,                &
              &   Prandtl,             &
              &   logical_CFL,         &
              &   scheme,              &
              &   ep_type,             &
              &   diffusion_scheme,    &
              &   kolmogorov_forcing,  &
              &   mu_jump_ratio,       &
              &   mu_jump_loc,         &
              &   sensor,              &
              &   ducrosToll,          &
              &   enable_limiter,      &
              &   limiter_rho_floor,   &
              &   limiter_e_floor,     &
              &   limiter_abort_frac,  &
              &   fd_order,            &
              &   weno_order,          &
              &   LES,                 &
              &   sgs_model,           &
              &   yoshizawa,           &
              &   viscous,             &
              &   n_step,              &
              &   printres,            &
              &   itout,               &
              &   StOut,               &
              &   rsOut,               &
              &   itPrb,               &
              &   itWallMap,           &
              &   StFlg,               &
              &   restartStat,         &
              &   ic,                  &
              &   output_file_name,    &
              &   data_dir,            &
              &   prb_nwall,           &
              &   prb_xc,              &
              &   prb_side,            &
              &   prb_nfield,          &
              &   prb_x,               &
              &   prb_y

end module parameters_module
