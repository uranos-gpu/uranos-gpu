module parameters_module

        use mpi

        implicit none

        integer , parameter :: sp  = selected_real_kind(6,37)  !< single
        integer , parameter :: dp  = selected_real_kind(15,307)  !< double
#ifdef SINGLE_PRECISION
        integer , parameter :: rp     = sp
        integer , parameter :: MPI_RP = MPI_SINGLE_PRECISION
#else
        integer , parameter :: rp     = dp
        integer , parameter :: MPI_RP = MPI_DOUBLE_PRECISION
#endif

        integer , parameter :: msm = MPI_SUM
        integer , parameter :: ip  = 4                   !< integer precision
        integer , parameter :: dl  = 600                 !< character default lenght
        integer , parameter :: eqs = 5                   !< number of equations
        real(rp), parameter :: pi  = 2.0_rp*asin(1.0_rp)       !< pi

        ! gamma and related properties
        real(rp), parameter :: gamma0 = 1.4_rp               !< gamma 
        real(rp), parameter :: cv   = 1._rp /(gamma0-1.0_rp) !< cv
        real(rp), parameter :: cp   = gamma0/(gamma0-1.0_rp) !< cv
        real(rp), parameter :: cv_i = 1._rp/cv               !< inverse of cv
        real(rp), parameter :: cp_i = 1._rp/cp               !< inverse of cp
        real(rp), parameter :: toll_equality = 1.0E-13_rp    !< check for equality
        real(rp)            :: suthc = 110.4_rp/273.15_rp    !< sutherland adimensianal coefficient

        ! shock bc parameter
        real(rp) :: wedge_angle = 8.0_rp
        real(rp) :: impin_point = 54.0_rp

        ! input MPI
        ! mpi_opt_level = 1 > sendrecv
        ! mpi_opt_level = 2 > isend/irecv + buffers
        ! mpi_opt_level = 3 > isend/irecv + derived data types
        integer :: mpi_opt_level = 3
        logical :: cuda_aware = .false.

        ! input per costruire la griglia di calcolo
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
        real(rp)              :: stretching_par
        
        ! input per la discretizzazione temporale
        real(rp)   :: tmax, Dt, CFL, time_from_restart
        integer    :: it,itmax, itStat = 0
        integer    :: n_step, ik = 1, itout, StOut
        integer    :: itPrb = 0
        logical    :: restart_flag  = .false.
        logical    :: inflow        = .true.
        logical    :: smooth_inflow = .false.
        logical    :: turb_inflow   = .false.
        logical    :: StFlg         = .false.
        logical    :: restartStat   = .false.
        real(rp)   :: s_cpu_time, e_cpu_time
        
        ! finite difference orders
        integer :: central_fd_order 
        integer :: bward_fd_order 
        integer :: fward_fd_order
        integer :: weno_num = 3
        integer :: weno_order = 5
        integer :: fd_L, fd_R

        ! finite difference coefficients
        real(rp), dimension(:), allocatable :: central_1, central_2
        real(rp), dimension(:), allocatable :: bward_1  , bward_2
        real(rp), dimension(:), allocatable :: fward_1  , fward_2
        real(rp), dimension(:), allocatable :: central_1_one_half
        real(rp), dimension(:), allocatable :: central_2_one_half

        real(rp), dimension(:), allocatable :: mid_point_lele
        real(rp), dimension(:,:), allocatable :: mid_point_lele_x
        real(rp), dimension(:,:), allocatable :: mid_point_lele_y
        real(rp), dimension(:,:), allocatable :: mid_point_lele_z

        real(rp), dimension(:,:) , allocatable :: aweno
        real(rp), dimension(:)   , allocatable :: cweno

        ! face type
        type face_type
          integer      :: node
          integer      :: norm
          character(1) :: face
        endtype face_type

        ! probes
        type prb
          integer , dimension(:,:), allocatable :: i
          real(rp), dimension(:,:), allocatable :: x
          integer                               :: n
        endtype
        type(prb) :: Probe

        ! Runge-Kutta coefficients
        real(rp), dimension(:), allocatable :: a_rk, b_rk, c_rk
        
        ! variabili per la scelta delle condizioni al contorno
        character(50)               :: ic     !< initial  condition
        character(50), dimension(6) :: bc     !< boundary condition
        logical      , dimension(6) :: sponge = .false.   
        real(rp)                    :: Trat  = 1.0_rp ! T_adi/T_rec ratio
        
        ! gruppi adimensionali
        real(rp) :: Reynolds, Mach, Prandtl, mu_inf, u_inf, k_inf, q_inf, ReTau
        real(rp) :: force_turbulent_channel
        
        ! nome dei file di output
        character(dl) :: output_file_name, data_dir, inflow_profile
        character(dl) :: restart_file

        ! large eddy simulation
        logical       :: LES   = .false. 
        logical       :: WMLES = .false.
        character(40) :: sgs_model
        
        ! condizione logica sul CFL e suo default
        logical       :: logical_CFL = .false.
        logical       :: hybrid_weno = .false.
        character(50) :: scheme
        character(20) :: sensor
        real(rp)      :: ducrosToll = 0.05_rp
        character(20) :: L1_wave
        integer       :: fd_order
        
        ! logical condition of viscous therms
        logical :: viscous = .false.
       
        ! skip file in post process
        integer :: skip_file = 1

        ! variabili che si vuole far stampare e loro default
        logical :: density             = .false.
        logical :: pressure            = .false.
        logical :: velocity            = .false.
        logical :: temperature         = .false.
        logical :: mach_               = .false.
        logical :: sdensity            = .false.
        logical :: vorticity           = .false.
        logical :: speed_div           = .false.
        logical :: charactheristic     = .false.
        logical :: hw_sensor           = .false.
        logical :: hweno_flag          = .false.
        logical :: bulk_quantities     = .false.
        logical :: qcriterion          = .false.
        logical :: vorticity_magnitude = .false.
        logical :: molecular_viscosity = .false.
        logical :: turbulent_viscosity = .false.
        logical :: hTurb_parameters    = .false.
        logical :: MixedVelocity       = .false.

        ! statistic
        logical :: Rey_average         = .false.
        logical :: Fav_average         = .false.
        
        ! formati di stampa e loro default
        logical :: BINARY     = .false. 
        logical :: GNUPLOT    = .false. 
        logical :: VTK        = .false. 
        logical :: VTK2D      = .false. 
        logical :: LIFT_DRAG  = .false. 
        logical :: NPY        = .false.
        
        ! lista dei parametri di input
        namelist /input_list/          &
        
              &   xmin,                & 
              &   xmax,                & 
              &   ymin,                & 
              &   ymax,                & 
              &   zmin,                & 
              &   zmax,                & 
              &   gridpoint_x,         &
              &   gridpoint_y,         &
              &   gridpoint_z,         &
              &   stretching_par,      &
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
              &   Dt,                  & 
              &   bc,                  &
              &   wedge_angle,         &
              &   impin_point,         &
              &   sponge,              &
              &   Trat,                &
              &   suthc,               &
              &   inflow_profile,      &
              &   smooth_inflow,       &
              &   turb_inflow,         &
              &   Reynolds,            &
              &   Mach,                &
              &   Prandtl,             &
              &   logical_CFL,         &         
              &   scheme,              &         
              &   sensor,              &         
              &   ducrosToll,          &
              &   L1_wave,             &
              &   fd_order,            &         
              &   weno_order,          &
              &   LES,                 &
              &   sgs_model,           &
              &   viscous,             &
              &   n_step,              & 
              &   itout,               &         
              &   StOut,               &         
              &   itPrb,               &
              &   StFlg,               &
              &   restartStat,         &
              &   ic,                  &
              &   output_file_name,    & 
              &   data_dir,            & 
              &   skip_file,           &
              &   density,             & 
              &   pressure,            & 
              &   velocity,            & 
              &   temperature,         & 
              &   mach_,               & 
              &   sdensity,            & 
              &   vorticity,           &
              &   speed_div,           &
              &   charactheristic,     &
              &   hw_sensor,           &
              &   hweno_flag,          &
              &   bulk_quantities,     &
              &   qcriterion,          &
              &   Rey_average,         &
              &   Fav_average,         &
              &   vorticity_magnitude, &
              &   molecular_viscosity, &
              &   turbulent_viscosity, &
              &   MixedVelocity,       &
              &   hTurb_parameters,    &
              &   BINARY,              &
              &   GNUPLOT,             &
              &   NPY,                 &
              &   VTK,                 &
              &   VTK2D,               &
              &   LIFT_DRAG


end module parameters_module
