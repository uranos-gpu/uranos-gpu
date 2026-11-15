!================================================================================
!  ramp_mesh_gen : body-fitted curvilinear mesh generator for the
!                  SUPERSONIC TURBULENT COMPRESSION RAMP
!                  (FLEW paper, Soldati et al. Sect. 3.4 ; Wu & Martin 2007)
!--------------------------------------------------------------------------------
!  theta-deg compression corner at x = x_corner (= 0).
!
!  DESIGN: ADAMS (1998) two-boundary LINEAR BLEND -- Fortran port of
!  martin_mesh_explorer.py / gen_martin_mesh.py (SAME equations, validated
!  value-identical to ~3e-14).  This is the body-fitted Martin/FLEW domain:
!  vertical inlet, perpendicular outlet, rounded corner, top follows the wall.
!
!    * lower wall  z_l(x) = lncosh ramp, SHARP radius Rlo = corner_length_bottom,
!                  corner centred at x_corner.
!    * upper bound z_u(x) = lncosh ramp, radius Rhi = corner_length_top, corner
!                  centred at X_TOP (the top fold/radius start -- independent of
!                  the wall corner), shifted so the inlet edge is VERTICAL and
!                  the outlet edge is PERPENDICULAR to the ramp.
!    * columns = straight lines lower->upper, linear blend r = eta/ymax:
!         xx = (1-r) x_l + r x_u ,   yy = (1-r) z_l + r z_u .
!      No folding, no hourglass; corner skew ~minimal.
!
!    * Streamwise grid: 'uniform' (FLEW Dx+ ~ const) or 'rampCornerSinh'.
!      Wall-normal grid: 'hyperbolicSin' (clustered, NyWr=0 -> auto-fit so the
!      geometric tail lands on ymax, like gen .py fit_nywr) / 'uniform'.
!      z is NOT written: URANOS extrudes it as a uniform 1D grid.
!
!  ----- TUNABLE PARAMETERS (same as martin_mesh_explorer.py / gen .py) -----
!     ramp_angle_deg        compression angle theta            [deg]
!     corner_length_bottom  Rlo  lower-wall radius (low=sharp)  [delta]
!     corner_length_top     Rhi  upper-boundary radius          [delta]
!     x_corner              wall (bottom) corner location       [delta]
!     x_top                 top fold/radius start (>0 keeps top flat longer)
!     ymax                  top offset H (wall-normal gap)      [delta]
!     LyWr,YPlus,NyWr       wall-normal hyperbolicSin clustering (NyWr=0 -> auto)
!     x_cluster             streamwise sinh factor (rampCornerSinh only)
!     corner_blend_height   UNUSED by the Adams blend (kept for namelist parsing)
!
!  Output layout (MUST match mesh_module::ReadCurvilinearGrid2D, ng_file = 7):
!     line 1 : nx  ny  1                           <- header: INTERIOR counts only
!     block X: (nx+2*GN)*(ny+2*GN) reals, i fastest then j   ; xx(i,j)
!     block Y: (nx+2*GN)*(ny+2*GN) reals, i fastest then j   ; yy(i,j)
!
!  Build : gfortran -O2 ramp_mesh_gen.f90 -o ramp_mesh_gen.exe
!  Run   : ./ramp_mesh_gen.exe input_ramp_curv.dat
!================================================================================
program ramp_mesh_gen
  implicit none
  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: dl = 1000
  integer, parameter :: GN = 7        !< MUST equal ng_file in ReadCurvilinearGrid2D
  real(dp),parameter :: pi = 3.14159265358979323846_dp

  ! --- namelist inputs ---
  real(dp)      :: xmin, xmax, ymin, ymax, zmin, zmax
  integer       :: nx, ny, nz
  character(dl) :: GridTypeX, GridTypeY, GridTypeZ
  real(dp)      :: ReTau, Mach
  ! wall-normal clustering (hyperbolicSin) controls
  real(dp)      :: LyWr, YPlus
  integer       :: NyWr
  ! streamwise (corner) controls
  real(dp)      :: x_corner          !< streamwise location of the ramp (wall) corner
  real(dp)      :: x_top             !< streamwise location where the UPPER boundary turns
                                     !< (top fold/radius start; = martin_mesh_explorer x_top)
  real(dp)      :: x_cluster         !< sinh clustering factor (GridTypeX='rampCornerSinh')
  ! ramp geometry (blend mapping -- same knobs as mesh_explorer.py)
  real(dp)      :: ramp_angle_deg       !< compression angle theta [deg]   (FLEW: 24)
  real(dp)      :: corner_length_bottom !< SHARP  wall turn length Lb [delta]
  real(dp)      :: corner_length_top    !< GENTLE top  turn length Lt [delta]  (Lt > Lb)
  real(dp)      :: corner_blend_height  !< sharp->gentle blend height Hb [delta] (< ramp height)
  ! --- legacy knobs kept only so older .dat files still parse (UNUSED) ---
  real(dp)      :: R_fillet, R_top, column_turn_length
  ! misc
  character(dl) :: ResDir
  character(dl) :: out_file

  namelist /input_list/ xmin,xmax,ymin,ymax,zmin,zmax, nx,ny,nz, &
       GridTypeX,GridTypeY,GridTypeZ, ReTau,Mach, &
       LyWr,YPlus,NyWr, x_corner,x_top,x_cluster, ramp_angle_deg, &
       corner_length_bottom, corner_length_top, corner_blend_height, &
       R_fillet, R_top, column_turn_length, ResDir, out_file

  real(dp), allocatable :: x(:), eta(:), z(:)
  real(dp), allocatable :: ywb(:), psib(:)    !< SHARP  wall height & deflection angle
  real(dp), allocatable :: ywt(:), psit(:)    !< GENTLE top  height & deflection angle
  real(dp), allocatable :: xx(:,:), yy(:,:)
  real(dp) :: theta                            !< ramp angle [rad]
  character(dl) :: infile
  integer :: ios, ufile
  logical :: hsg_quiet = .false.               !< suppress hyperbolic_sin_grid prints (auto-NyWr search)
  integer :: nw_try, nw_best                   !< auto-NyWr search (mirror gen .py fit_nywr)
  real(dp) :: em, err_best

  ! ---- defaults (mirror mesh_explorer.py slider initial values) ----
  GridTypeX = 'uniform'; GridTypeY = 'hyperbolicSin'; GridTypeZ = 'uniform'
  Mach = 2.9_dp; ResDir = 'temp'
  x_corner = 0.0_dp; x_top = 0.0_dp; x_cluster = 2.8_dp
  ramp_angle_deg       = 24.0_dp
  corner_length_bottom = 0.2_dp                ! Rlo (sharp lower wall)  -- gen .py default
  corner_length_top    = 0.2_dp                ! Rhi (upper boundary)    -- gen .py default
  corner_blend_height  = 2.5_dp
  R_fillet = -1.0_dp; R_top = -1.0_dp; column_turn_length = -1.0_dp   ! legacy, unused
  LyWr = 0.0_dp; YPlus = 0.5_dp; NyWr = 0
  out_file = 'mesh_ramp_curv.txt'

  ! ---- read input file ----
  call get_command_argument(1, infile)
  if (len_trim(infile) == 0) stop ' ERROR: missing input file (usage: ./ramp_mesh_gen.exe input.dat)'
  open(newunit=ufile, file=trim(infile), status='old', iostat=ios)
  if (ios /= 0) stop ' ERROR: could not open input file'
  read(ufile, nml=input_list, iostat=ios)
  if (ios /= 0) stop ' ERROR: could not parse &input_list namelist'
  close(ufile)

  theta = ramp_angle_deg*pi/180.0_dp

  write(*,*) ' --------------------------------------------------------- '
  write(*,*) '  URANOS curvilinear COMPRESSION-RAMP mesh generator'
  write(*,*) '  (Adams 1998 two-boundary blend; x_top = top fold start)'
  write(*,*) ' --------------------------------------------------------- '
  write(*,'(A,3I6)')    '  nx ny nz             : ', nx, ny, nz
  write(*,'(A,f8.2,A)') '  ramp angle           : ', ramp_angle_deg, ' deg'
  write(*,'(A,f8.3,A)') '  corner len bottom Lb : ', corner_length_bottom, ' delta (SHARP wall turn)'
  write(*,'(A,f8.3,A)') '  corner len top    Lt : ', corner_length_top,    ' delta (GENTLE top turn, starts upstream)'
  write(*,'(A,f8.3,A)') '  blend height      Hb : ', corner_blend_height,  ' delta (sharp below -> gentle above)'
  write(*,'(A,f10.4)')  '  corner location  x_c : ', x_corner
  write(*,'(A,f10.4)')  '  top fold start  x_top: ', x_top
  write(*,'(A,f10.3)')  '  ReTau                : ', ReTau

  allocate(x(1-GN:nx+GN), eta(1-GN:ny+GN), z(1-GN:nz+GN))
  allocate(ywb(1-GN:nx+GN), psib(1-GN:nx+GN))
  allocate(ywt(1-GN:nx+GN), psit(1-GN:nx+GN))
  x = 0.0_dp; eta = 0.0_dp; z = 0.0_dp

  ! ============================ streamwise grid x(i) =========================
  select case (trim(GridTypeX))
  case ('uniform');        call uniform_grid(x, nx, xmin, xmax, 'x')
  case ('rampCornerSinh'); call ramp_corner_grid(x)
  case default;            stop ' GridTypeX not supported (uniform / rampCornerSinh)'
  end select

  ! ============== (Adams two-boundary blend uses lncosh in build_mesh) ========
  ! corner_length_bottom -> Rlo (sharp lower-wall radius)
  ! corner_length_top    -> Rhi (gentle upper-boundary radius)

  ! ===================== wall-normal clustering eta(j) in [ymin,ymax] ========
  ! auto-size NyWr so the geometric tail lands exactly on ymax (mirror of the
  ! gen_martin_mesh.py fit_nywr): only when NyWr<=0 in the namelist.
  if (trim(GridTypeY) == 'hyperbolicSin' .and. NyWr <= 0) then
     err_best = huge(1.0_dp); nw_best = 60; hsg_quiet = .true.
     do nw_try = 60, ny-20
        NyWr = nw_try
        call hyperbolic_sin_grid(eta)
        em = maxval(eta(1:ny))
        if (abs(em - ymax) < err_best) then
           err_best = abs(em - ymax); nw_best = nw_try
        end if
     end do
     hsg_quiet = .false.; NyWr = nw_best
     write(*,'(A,I0)') '  auto NyWr            : ', NyWr
  end if

  select case (trim(GridTypeY))
  case ('hyperbolicSin');         call hyperbolic_sin_grid(eta)
  case ('hyperbolicSinConstant'); call hyperbolic_sin_constant_grid(eta)
  case ('uniform');               call uniform_grid(eta, ny, ymin, ymax, 'y')
  case default; stop ' GridTypeY not supported (hyperbolicSin / hyperbolicSinConstant / uniform)'
  end select

  ! ===================== spanwise (report only; URANOS extrudes z) ===========
  select case (trim(GridTypeZ))
  case ('uniform'); call uniform_grid(z, nz, zmin, zmax, 'z')
  case default;     stop ' GridTypeZ not supported (only "uniform")'
  end select

  ! ============================ build + check + write ========================
  allocate(xx(1-GN:nx+GN,1-GN:ny+GN), yy(1-GN:nx+GN,1-GN:ny+GN))
  call build_mesh(x, eta, xx, yy)
  call check_jacobian(xx, yy)
  call write_uranos_plot3d(xx, yy)

  write(*,*)
  write(*,'(A,A)')     '  Mesh written to      : ', trim(out_file)
  write(*,'(A,2I6,A)') '  (header dims         : ', nx, ny, '   1 )  interior counts (no ghosts)'
  write(*,'(A,2I6)')   '  (data rows per block : ', nx+2*GN, ny+2*GN
  write(*,'(A)')       '  -> set  mesh_file = '''//trim(out_file)//'''  in the URANOS .dat'
  write(*,*)

contains

  !------------------------------------------------------------------
  ! smootherstep : C2 Perlin ramp  S(t) = 6t^5 - 15t^4 + 10t^3,
  ! clamped to [0,1].  S(0)=0, S(1)=1, S'(0)=S'(1)=S''(0)=S''(1)=0.
  !------------------------------------------------------------------
  pure elemental function smootherstep(t) result(s)
    real(dp), intent(in) :: t
    real(dp) :: s, u
    u = max(0.0_dp, min(1.0_dp, t))
    s = u*u*u*(u*(u*6.0_dp - 15.0_dp) + 10.0_dp)
  end function smootherstep

  !------------------------------------------------------------------
  ! deflection_profile : wall deflection psi(x) ramping 0 -> theta with a
  ! smootherstep over a streamwise turn length L, centred on x_corner (so the
  ! turn spans [x_corner-L, x_corner]).  A SHORT L = sharp; a LONG L starts
  ! the turn further UPSTREAM.  The wall elevation y_w(x) = int tan(psi) dx'
  ! (trapezoidal, from the upstream ghost where psi=0 so y_w=0 on the flat).
  !------------------------------------------------------------------
  subroutine deflection_profile(L, psi, yw)
    real(dp), intent(in)  :: L
    real(dp), intent(out) :: psi(1-GN:), yw(1-GN:)
    integer  :: i
    real(dp) :: Le, t
    Le = max(L, 1.0e-9_dp)
    do i = 1-GN, nx+GN
       t      = (x(i) - (x_corner - 0.5_dp*Le)) / Le
       psi(i) = theta*smootherstep(t)
    enddo
    yw(1-GN) = 0.0_dp
    do i = 1-GN+1, nx+GN
       yw(i) = yw(i-1) + 0.5_dp*(tan(psi(i)) + tan(psi(i-1)))*(x(i) - x(i-1))
    enddo
  end subroutine deflection_profile

  !------------------------------------------------------------------
  ! Numerically-stable log(cosh(t)) = |t| + log1p(exp(-2|t|)) - log 2.
  !------------------------------------------------------------------
  pure function logcosh(t) result(r)
    real(dp), intent(in) :: t
    real(dp) :: r, a
    a = abs(t)
    r = a + log(1.0_dp + exp(-2.0_dp*a)) - log(2.0_dp)
  end function logcosh

  !------------------------------------------------------------------
  ! Adams (1998) eq.(2) smooth ramp: flat (z=0) upstream, slope tan(theta)
  ! downstream, finite corner radius Rc at (x_corner,0).  Bigger Rc = gentler.
  ! Anchored so the upstream flat (at xmin) is z=0.
  !------------------------------------------------------------------
  pure function lncosh_wall(xq, th, Rc, xc) result(z)
    real(dp), intent(in) :: xq, th, Rc, xc      !< xc = streamwise corner location
    real(dp) :: z, d1, d2, z0
    d2 = tan(th)/2.0_dp
    d1 = (1.0_dp + d2*d2)**1.5_dp / max(Rc, 1.0e-9_dp)
    z  = d2*(xq   + (1.0_dp/d1)*logcosh(d1*(xq   - xc)))
    z0 = d2*(xmin + (1.0_dp/d1)*logcosh(d1*(xmin - xc)))
    z  = z - z0
  end function lncosh_wall

  !------------------------------------------------------------------
  ! BUILD MESH : ADAMS (1998) two-boundary LINEAR BLEND.
  ! (Theor. Comput. Fluid Dyn. 12:109, Sect. 2.3.)  Linear interpolation between
  ! a SHARP lower wall and a GENTLE upper boundary; columns are STRAIGHT lines
  ! connecting corresponding lower/upper points:
  !     xx = (1-r) x_l + r x_u ,  yy = (1-r) z_l + r z_u ,   r = eta/ymax
  !   * lower z_l : ln(cosh) ramp, SHARP corner radius Rlo (=corner_length_bottom)
  !   * upper     : ln(cosh) ramp, GENTLE radius Rhi (=corner_length_top), placed so
  !     the inlet edge is VERTICAL and the outlet edge is PERPENDICULAR to the ramp
  !     (upper outlet offset by H=ymax along the ramp normal).
  ! Gentle upper -> straight columns rotate smoothly vertical->perp: minimal corner
  ! skew (~15 deg), NEVER folds, NO hourglass / dip.  (corner_blend_height unused.)
  !------------------------------------------------------------------
  subroutine build_mesh(x, eta, xx, yy)
    real(dp), intent(in)  :: x(1-GN:), eta(1-GN:)
    real(dp), intent(out) :: xx(1-GN:,1-GN:), yy(1-GN:,1-GN:)
    integer  :: i, j
    real(dp) :: H, th, sphi, cphi, s, xl, zl, xu, zu, zsh, r
    real(dp) :: zl_out, zsh_out, scale
    th   = theta
    H    = ymax
    sphi = sin(th);  cphi = cos(th)
    zl_out  = lncosh_wall(xmax,            th, corner_length_bottom, x_corner)
    zsh_out = lncosh_wall(xmax - H*sphi,   th, corner_length_top,    x_top)
    if (abs(zsh_out) < 1.0e-12_dp) zsh_out = 1.0e-12_dp
    scale   = (zl_out + H*cphi - H) / zsh_out
    do i = 1-GN, nx+GN
       s   = (x(i) - xmin) / (xmax - xmin)            ! streamwise parameter
       xl  = x(i)
       zl  = lncosh_wall(xl,            th, corner_length_bottom, x_corner)   ! SHARP lower (wall corner)
       xu  = xl - H*sphi*s                            ! perp-offset outlet
       zsh = lncosh_wall(xu,            th, corner_length_top,    x_top)      ! GENTLE upper (own corner = x_top)
       zu  = H + zsh*scale
       do j = 1-GN, ny+GN
          r = eta(j)/H
          xx(i,j) = (1.0_dp - r)*xl + r*xu
          yy(i,j) = (1.0_dp - r)*zl + r*zu
       enddo
    enddo
  end subroutine build_mesh

  !------------------------------------------------------------------
  ! Streamwise grid clustered at the corner (symmetric sinh on each side).
  ! OPTIONAL -- the FLEW reference uses a ~uniform streamwise grid, so this
  ! is only used when GridTypeX = 'rampCornerSinh'.  Cell-centered.
  !------------------------------------------------------------------
  subroutine ramp_corner_grid(g)
    real(dp), intent(inout) :: g(1-GN:)
    integer  :: i, nL, nR
    real(dp) :: LL, LR, s, beta, dxe
    LL = x_corner - xmin
    LR = xmax - x_corner
    nL = nint(real(nx,dp)*LL/(LL+LR))
    nR = nx - nL
    beta = max(x_cluster, 1.0e-6_dp)
    do i = 1, nL
       s    = (real(i,dp) - 0.5_dp)/real(nL,dp)
       g(i) = x_corner - LL*sinh(beta*(1.0_dp - s))/sinh(beta)
    end do
    do i = 1, nR
       s        = (real(i,dp) - 0.5_dp)/real(nR,dp)
       g(nL+i)  = x_corner + LR*sinh(beta*s)/sinh(beta)
    end do
    dxe = g(2) - g(1)
    do i = 0, 1-GN, -1
       g(i) = g(i+1) - dxe
    end do
    dxe = g(nx) - g(nx-1)
    do i = nx+1, nx+GN
       g(i) = g(i-1) + dxe
    end do
    call smooth_grid(g, 2000)
    write(*,'(A,2I6)')   '  Streamwise pts up/down: ', nL, nR
    write(*,'(A,2f12.5)')'  DeltaX+ corner/coarse : ', &
         min(g(nL+1)-g(nL), g(nL)-g(nL-1))*ReTau, (g(2)-g(1))*ReTau
  end subroutine ramp_corner_grid

  !------------------------------------------------------------------
  ! Cell-centered uniform grid (copied from BL_grid_module::uniform*Grid)
  !------------------------------------------------------------------
  subroutine uniform_grid(g, n, gmin, gmax, lbl)
    real(dp), intent(inout) :: g(1-GN:)
    integer,  intent(in)    :: n
    real(dp), intent(in)    :: gmin, gmax
    character(*), intent(in):: lbl
    integer  :: i
    real(dp) :: xi, L
    L = gmax - gmin
    do i = 1-GN, n+GN
       xi   = (i - 0.5_dp)/real(n,dp)
       g(i) = gmin + L*xi
    enddo
    write(*,'(A,A,A,f14.6)') '  Delta', lbl, '+ (uniform) : ', (g(2)-g(1))*ReTau
  end subroutine uniform_grid

  !------------------------------------------------------------------
  ! HyperbolicSin wall-normal grid (copied VERBATIM from
  ! BL_grid_module::HyperbolicSinYGrid via curv_mesh_gen): sinh-clustered
  ! well-resolved region [0,LyWr] with NyWr cells, geometric-progression
  ! tail to ymax, antisymmetric wall ghosts, then Laplacian smoothing.
  ! Returns eta in [ymin=0, ymax]  (the wall-normal offset, NOT y).
  !------------------------------------------------------------------
  subroutine hyperbolic_sin_grid(g)
    real(dp), intent(inout) :: g(1-GN:)
    real(dp), parameter :: toll  = 1.0e-15_dp
    integer , parameter :: itmax = 1000
    integer , parameter :: nsmooth = 100
    real(dp) :: err, alpha, alpha_old, xi, yp_target, yp_c
    real(dp) :: dyFinal, Ly, LyRest, r, r_old
    real(dp) :: mean_dyP, maxj_dyP
    integer  :: iter, j, jg, ji, NyRest

    yp_target = YPlus

    alpha = 1.0_dp; err = 2.0_dp*toll; iter = 0
    do while (err > toll .and. iter < itmax)
       alpha_old = alpha
       xi    = (1.0_dp - 0.5_dp)/real(NyWr,dp)
       alpha = asinh(ReTau*LyWr*sinh(alpha_old*xi)/yp_target)
       yp_c  = ReTau*LyWr*sinh(alpha*xi)/sinh(alpha)
       err   = abs(yp_c - yp_target)
       iter  = iter + 1
    enddo
    if (alpha < 1.0_dp) then
       if (hsg_quiet) then
          g = huge(1.0_dp); return        ! bad NyWr during auto-search: reject, don't stop
       else
          stop ' ERROR: hyperbolicSin generation fails (alpha < 1)'
       end if
    end if

    g = 0.0_dp
    do j = 1, NyWr+1
       xi   = (j - 0.5_dp)/real(NyWr,dp)
       g(j) = LyWr*sinh(alpha*xi)/sinh(alpha)
    enddo

    do j = 1, GN
       jg = 1 - j;  ji = j
       g(jg) = -g(ji)
    enddo
    dyFinal = g(NyWr+1) - g(NyWr)

    Ly     = ymax - ymin
    LyRest = Ly - LyWr
    NyRest = ny - NyWr
    err = 2.0_dp*toll; iter = 0; r = 1.01_dp
    do while (err > toll .and. iter < itmax)
       r_old = r
       r = 1.0_dp - (1.0_dp - r)/dyFinal*LyRest
       r = r**(1.0_dp/real(NyRest,dp))
       iter = iter + 1
       err  = abs(r - r_old)
    enddo
    do j = 1, NyRest+GN
       g(NyWr+j) = g(NyWr+j-1) + dyFinal*r**(j-1.0_dp)
    enddo

    call smooth_grid(g, nsmooth)

    mean_dyP = 0.0_dp; maxj_dyP = 0.0_dp
    do j = 1, ny
       mean_dyP = mean_dyP + (g(j)-g(j-1))*ReTau
       maxj_dyP = max((g(j)-g(j-1))*ReTau, maxj_dyP)
    enddo
    mean_dyP = mean_dyP/real(ny,dp)
    if (.not. hsg_quiet) then
       write(*,'(A,3f14.6)') '  DeltaY+ wall/max/mean : ', g(1)*ReTau, maxj_dyP, mean_dyP
       write(*,'(A,f12.6,A,e12.4)') '  alpha = ', alpha, '   r = ', r
    end if
  end subroutine hyperbolic_sin_grid

  !------------------------------------------------------------------
  ! HyperbolicSinConstant wall-normal grid (copied VERBATIM from
  ! BL_grid_module::HyperbolicSinConstantYGrid): same sinh-clustered
  ! well-resolved region [0,LyWr] with NyWr cells (first cell -> YPlus),
  ! but the OUTER region [LyWr,ymax] is filled with CONSTANT spacing
  ! dyFinal = (ymax-LyWr)/(ny-NyWr) instead of a geometric progression.
  ! -> uniform resolution in the free stream / outer region (good for a
  !    shock or shear layer that lives ABOVE the wall-resolved band), at
  !    the price of a spacing jump at the LyWr interface (softened by the
  !    100 Laplacian smooths) when (ymax-LyWr)/(ny-NyWr) >> last sinh cell.
  !------------------------------------------------------------------
  subroutine hyperbolic_sin_constant_grid(g)
    real(dp), intent(inout) :: g(1-GN:)
    real(dp), parameter :: toll  = 1.0e-15_dp
    integer , parameter :: itmax = 1000
    integer , parameter :: nsmooth = 100
    real(dp) :: err, alpha, alpha_old, xi, yp_target, yp_c
    real(dp) :: dyFinal, Ly, LyRest
    real(dp) :: mean_dyP, maxj_dyP
    integer  :: iter, j, jg, ji, NyRest

    yp_target = YPlus

    alpha = 1.0_dp; err = 2.0_dp*toll; iter = 0
    do while (err > toll .and. iter < itmax)
       alpha_old = alpha
       xi    = (1.0_dp - 0.5_dp)/real(NyWr,dp)
       alpha = asinh(ReTau*LyWr*sinh(alpha_old*xi)/yp_target)
       yp_c  = ReTau*LyWr*sinh(alpha*xi)/sinh(alpha)
       err   = abs(yp_c - yp_target)
       iter  = iter + 1
    enddo
    if (alpha < 1.0_dp) stop ' ERROR: hyperbolicSinConstant generation fails (alpha < 1)'

    g = 0.0_dp
    do j = 1, NyWr+1
       xi   = (j - 0.5_dp)/real(NyWr,dp)
       g(j) = LyWr*sinh(alpha*xi)/sinh(alpha)
    enddo

    do j = 1, GN
       jg = 1 - j;  ji = j
       g(jg) = -g(ji)
    enddo

    ! --- outer region: CONSTANT spacing ---
    Ly     = ymax - ymin
    LyRest = Ly - LyWr
    NyRest = ny - NyWr
    dyFinal = LyRest/real(NyRest,dp)
    do j = 1, NyRest+GN
       g(NyWr+j) = g(NyWr+j-1) + dyFinal
    enddo

    call smooth_grid(g, nsmooth)

    mean_dyP = 0.0_dp; maxj_dyP = 0.0_dp
    do j = 1, ny
       mean_dyP = mean_dyP + (g(j)-g(j-1))*ReTau
       maxj_dyP = max((g(j)-g(j-1))*ReTau, maxj_dyP)
    enddo
    mean_dyP = mean_dyP/real(ny,dp)
    write(*,'(A,3f14.6)') '  DeltaY+ wall/max/mean : ', g(1)*ReTau, maxj_dyP, mean_dyP
    write(*,'(A,f12.6,A,f12.6)') '  alpha = ', alpha, '   dy_outer = ', dyFinal
  end subroutine hyperbolic_sin_constant_grid

  !------------------------------------------------------------------
  ! Laplacian smoothing (copied from BL_grid_module::SmoothTheGrid)
  !------------------------------------------------------------------
  subroutine smooth_grid(g, nsmooth)
    real(dp), intent(inout) :: g(1-GN:)
    integer,  intent(in)    :: nsmooth
    real(dp), allocatable :: tmp(:)
    integer :: j, is, lb, ub
    lb = lbound(g,1); ub = ubound(g,1)
    allocate(tmp(lb:ub))
    do is = 1, nsmooth
       tmp = g
       do j = lb+1, ub-1
          g(j) = (1.0_dp/6.0_dp)*(tmp(j-1) + 4.0_dp*tmp(j) + tmp(j+1))
       enddo
    enddo
    deallocate(tmp)
  end subroutine smooth_grid

  !------------------------------------------------------------------
  ! Cell Jacobian scan : flags grid folding (J <= 0 anywhere in the
  ! interior).  J = x_xi*y_eta - x_eta*y_xi (central differences in i,j).
  !------------------------------------------------------------------
  subroutine check_jacobian(xx, yy)
    real(dp), intent(in) :: xx(1-GN:,1-GN:), yy(1-GN:,1-GN:)
    integer  :: i, j, nbad
    real(dp) :: xxi, yxi, xet, yet, jac, jmin
    jmin = huge(1.0_dp); nbad = 0
    do j = 1, ny
       do i = 1, nx
          xxi = 0.5_dp*(xx(i+1,j) - xx(i-1,j))
          yxi = 0.5_dp*(yy(i+1,j) - yy(i-1,j))
          xet = 0.5_dp*(xx(i,j+1) - xx(i,j-1))
          yet = 0.5_dp*(yy(i,j+1) - yy(i,j-1))
          jac = xxi*yet - xet*yxi
          jmin = min(jmin, jac)
          if (jac <= 0.0_dp) nbad = nbad + 1
       enddo
    enddo
    write(*,'(A,e13.5)') '  min cell Jacobian     : ', jmin
    if (nbad > 0) then
       write(*,'(A,I0,A)') '  *** GRID FOLDED in ', nbad, &
            ' cells -> raise corner_length_top (gentler) or lower corner_blend_height ***'
    else
       write(*,'(A)')      '  Jacobian OK (mesh is valid, no folding).'
    end if
  end subroutine check_jacobian

  !------------------------------------------------------------------
  ! Write the URANOS curvilinear PLOT3D block file (i fastest then j).
  !------------------------------------------------------------------
  subroutine write_uranos_plot3d(xx, yy)
    real(dp), intent(in) :: xx(1-GN:,1-GN:), yy(1-GN:,1-GN:)
    integer :: i, j, u
    open(newunit=u, file=trim(out_file), status='replace', action='write')
    ! URANOS (mesh_module::read_curvilinear_mesh) validates the header against the
    ! INTERIOR counts (nxh/=nx .or. nyh/=ny), so the first row carries nx,ny WITHOUT
    ! ghost nodes.  The GN ghost layers below are still written (the reader consumes
    ! ng = GN_solver + mesh_gh = 7 ghost layers per side from the data blocks).
    write(u,'(3I8)') nx, ny, 1
    do j = 1-GN, ny+GN
       do i = 1-GN, nx+GN
          write(u,'(ES25.16)') xx(i,j)
       enddo
    enddo
    do j = 1-GN, ny+GN
       do i = 1-GN, nx+GN
          write(u,'(ES25.16)') yy(i,j)
       enddo
    enddo
    close(u)
  end subroutine write_uranos_plot3d

end program ramp_mesh_gen
