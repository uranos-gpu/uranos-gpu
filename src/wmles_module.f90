module wmles_module
use parameters_module, only: rp
implicit none

real(rp), parameter , public :: wmles_toll = 1.0E-08_rp    !< WMLES tolerance to convergency (1.0E-14_rp)
real(rp), parameter , public :: wmles_intr = 40.0_rp       !< WMLES DYNAMIC INTERFACE
integer , parameter , public :: wmles_indx = 4          !< WMLES STATIC  INTERFACE
integer , parameter , public :: wmles_npts = 30         !< WMLES NUMBER OF POINTS (50)
integer , parameter , public :: wmles_imax = 100        !< WMLES NUMBER OF ITERATIONS
integer , parameter , public :: wmles_strI = 4          !< WMLES STARTING INDEX

real(rp), parameter , public :: xpt_wm = 40.0_rp           !< WR/WM interface along x
real(rp), parameter , public :: ypt_wm =  5.0_rp           !< WR/WM interface along y
real(rp), parameter , public :: zpt_wm = 20.0_rp           !< WR/WM interface along z

real(rp), parameter , public :: omg_wm = 0.05_rp           !< Vorticity magnitude threshold
integer , parameter , public :: wmles_max_id = 6
integer , parameter , public :: wmles_min_id = wmles_strI

contains

subroutine Get_WMLES_Grid(nw,hw,ReTau,yf,yc,dy)
        implicit none

        integer , intent(in) :: nw      !< grid point number
        real(rp), intent(in) :: hw      !< height of the wmles grid
        real(rp), intent(in) :: ReTau   !< local friction Reynolds number

        real(rp), dimension(0:nw), intent(out) :: yf   !> face grid
        real(rp), dimension(1:nw), intent(out) :: yc   !> cent grid
        real(rp), dimension(0:nw), intent(out) :: dy   !> step grid
        
        ! local declarations
        real(rp), parameter :: toll = wmles_toll
        integer , parameter :: itmx = wmles_imax
        real(rp)            :: dy_wall
        real(rp)            :: rsOld, rs, rerr
        integer             :: j, iter

        ! build the face grid
        dy_wall = 0.1_rp/Retau

        ! get the stretching parameter
        rsOld = 1.01_rp
        rs    = rsOld
        iter  = 0
        rerr  = 2*toll
        do while(rerr > toll .and. iter < itmx)
           rs = (1.0_rp+(rs-1.0_rp)/dy_wall*hw)**(1.0_rp/nw)

           iter = iter + 1
           rerr = abs(rs - rsOld)
           rsOld = rs
        enddo

        yf(0) = 0.0_rp
        do j = 1,nw
           yf(j) = dy_wall*(rs**j-1.0_rp)/(rs-1.0_rp)
        enddo

        ! build the centroids grid
        do j = 1,nw
           yc(j) = 0.5_rp*(yf(j) + yf(j-1))
        enddo
        ! build the between centroids
        dy(0) = yc(1)
        do j = 1,nw-1
           dy(j) = yc(j+1)-yc(j)
        enddo
        dy(nw) = yf(nw)-yc(nw)


        return
end subroutine Get_WMLES_Grid



subroutine OdeWMLES(gamma0,Prandtl,mu_inf,hw,ReTau,u_w,u_p,T_w,T_h,p_h,u_wm,T_wm,tauW,qauW)
        !$acc routine seq
        
        use fluid_functions_module , only: Sutherland, getmuT
        use matrix_inversion_module, only: tdma

        implicit none
        
        integer, parameter :: nw = wmles_npts

        real(rp)                 , intent(in) :: ReTau   !< local friction Reynolds number
        real(rp)                 , intent(in) :: gamma0, Prandtl, mu_inf
        real(rp)                 , intent(in) :: hw, u_p, T_h, p_h, u_w, T_w
        real(rp)                 , intent(out):: tauW, qauW
        real(rp), dimension(1:nw), intent(out):: u_wm, T_wm
        
        ! local declarations
        real(rp), dimension(1:nw) :: a_,b_,c_,d_
        real(rp), dimension(0:nw) :: k_u, k_T

        real(rp) :: r_w, mlW, klW
        real(rp) :: r_m, T_m, mLm, mTm
        real(rp) :: u_hm, u_hp ,u_dm, u_dp

        real(rp), parameter :: toll = wmles_toll
        integer , parameter :: itmx = wmles_imax
        real(rp)            :: PrL, PrT, c_p
        real(rp)            :: tauW_old, qauW_old, rerr
        integer             :: j, iter

        ! local declarations
        real(rp), dimension(0:nw) :: yf
        real(rp), dimension(1:nw) :: yc
        real(rp), dimension(0:nw) :: dy
        real(rp)                  :: dy_wall
        real(rp)                  :: rsOld, rs

        !
        ! GET THE GRID
        !
        ! build the face grid
        dy_wall = 0.1_rp/Retau

        ! get the stretching parameter
        rsOld = 1.01_rp
        rs    = rsOld
        iter  = 0
        rerr  = 2*toll
        do while(rerr > toll .and. iter < itmx)
           rs = (1.0_rp+(rs-1.0_rp)/dy_wall*hw)**(1.0_rp/nw)

           iter = iter + 1
           rerr = abs(rs - rsOld)
           rsOld = rs
        enddo

        yf(0) = 0.0_rp
        do j = 1,nw
           yf(j) = dy_wall*(rs**j-1.0_rp)/(rs-1.0_rp)
        enddo

        ! build the centroids grid
        do j = 1,nw
           yc(j) = 0.5_rp*(yf(j) + yf(j-1))
        enddo
        ! build the between centroids
        dy(0) = yc(1)
        do j = 1,nw-1
           dy(j) = yc(j+1)-yc(j)
        enddo
        dy(nw) = yf(nw)-yc(nw)



        !
        ! SOLVE EQUATIONS
        !
        c_p = gamma0/(gamma0-1.0_rp)
        PrL = Prandtl
        PrT = 0.9_rp
        
        ! get missing wall info
        r_w = p_h/T_w
        mlW = mu_inf*Sutherland(T_w)
        klW = c_p*mlW/PrL

        ! init velocities and temperature profiles
        do j = 1,nw
           u_wm(j) = u_w + yc(j)*(u_p-u_w)/hw
           T_wm(j) = T_w + yc(j)*(T_h-T_w)/hw
        enddo

        ! start iterating 
        iter     = 0
        rerr     = 2*toll
        tauW     = 0.0_rp
        qauW     = 0.0_rp
        tauW_old = 0.0_rp
        qauW_old = 0.0_rp
        do while(rerr > toll .and. iter < itmx)
        
           tauW_old = tauW
           qauW_old = qauW
        
           tauW = mlW*(u_wm(1) - u_w)/dy(0)
           qauW = klW*(T_wm(1) - T_w)/dy(0)
        
           ! get momentum coefficients
           j = 0
           T_m    = T_w
           r_m    = p_h/T_m
           mLm    = mu_inf*Sutherland(T_m)
           mTm    = getMuT(r_m,tauW,yf(j),mLm)
           k_u(j) = (mLm + mTm)/dy(j)
           k_T(j) = c_p*(mLm/PrL + mTm/PrT)/dy(j)

           do j = 1,nw-1
              T_m = 0.5_rp*(T_wm(j) + T_wm(j+1))
              r_m    = p_h/T_m
              mLm    = mu_inf*Sutherland(T_m)
              mTm    = getMuT(r_m,tauW,yf(j),mLm)
              k_u(j) = (mLm + mTm)/dy(j)
              k_T(j) = c_p*(mLm/PrL + mTm/PrT)/dy(j)
           enddo

           j = nw
           T_m    = T_h
           r_m    = p_h/T_m
           mLm    = mu_inf*Sutherland(T_m)
           mTm    = getMuT(r_m,tauW,yf(j),mLm)
           k_u(j) = (mLm + mTm)/dy(j)
           k_T(j) = c_p*(mLm/PrL + mTm/PrT)/dy(j)
        
           ! building tri-diagonal matrix for momentum equation
           do j = 1,nw 
              a_(j) =   k_u(j-1) 
              b_(j) = - k_u(j)-k_u(j-1)
              c_(j) =   k_u(j)
              d_(j) =   0.0_rp
           enddo
           d_(1 ) = - k_u(0 )*u_w
           d_(nw) = - k_u(nw)*u_p
        
           call tdma(a_,b_,c_,d_,u_wm,1,nw)
        
           ! building tri-diagonal matrix for temperature
           do j = 2,nw-1
        
              a_(j) =  k_T(j-1)
              b_(j) = -k_T(j) - k_T(j-1)
              c_(j) =  k_T(j)
        
              u_hm = (u_wm(j) + u_wm(j-1))*0.5_rp
              u_hp = (u_wm(j+1) + u_wm(j))*0.5_rp
              u_dm = (u_wm(j) - u_wm(j-1))
              u_dp = (u_wm(j+1) - u_wm(j))
              d_(j) = - k_u(j)*u_hp*u_dp + k_u(j-1)*u_hm*u_dm
        
           enddo
        
           j = 1
           a_(j) =  k_T(j-1)
           b_(j) = -k_T(j) - k_T(j-1)
           c_(j) =  k_T(j)
        
           u_hm = u_w!(u_wm(j) + u_wm(j-1))*0.5_rp
           u_hp = (u_wm(j+1) + u_wm(j))*0.5_rp
           u_dm = (u_wm(j) - u_w)
           u_dp = (u_wm(j+1) - u_wm(j))
           d_(j) = -k_T(0)*T_w - k_u(j)*u_hp*u_dp + k_u(j-1)*u_hm*u_dm
        
           j = nw
           a_(j) =  k_T(j-1)
           b_(j) = -k_T(j) - k_T(j-1)
           c_(j) =  k_T(j)
        
           u_hm = (u_wm(j) + u_wm(j-1))*0.5_rp
           u_hp = u_p!(u_wm(j+1) + u_wm(j))*0.5_rp
           u_dm = (u_wm(j) - u_wm(j-1))
           u_dp = (u_p - u_wm(j))
           d_(j) = -k_T(nw)*T_h - k_u(j)*u_hp*u_dp + k_u(j-1)*u_hm*u_dm
        
           call tdma(a_,b_,c_,d_,T_wm,1,nw)
        
           iter = iter + 1 
           rerr = max(abs(tauW - tauW_old), abs(qauW - qauW_old))

        enddo

        return
end subroutine OdeWMLES





subroutine LogLawWMLES(mu_inf,d_h,T_W,u_h,p_h,tauWall)
        !$acc routine seq

        use fluid_functions_module, only: sutherland, iReichardtLaw

        implicit none
        real(rp), intent(in)  :: mu_inf
        real(rp), intent(in)  :: T_W
        real(rp), intent(in)  :: d_h,u_h,p_h
        real(rp), intent(out) :: tauWall

        ! local declarations
        real(rp), parameter :: toll = wmles_toll
        integer , parameter :: imax = wmles_imax
        real(rp) :: r_w, m_w, n_w, t_W0, u_t0, u_T

        real(rp) :: x0   
        real(rp) :: xnew
        integer  :: iter                     !< iter variable
        real(rp) :: ctoll                    !< calculated tolerance
        real(rp) :: xold1, xold2, hk
        

        ! get wall quantities
        r_w = p_h/T_w
        m_w = mu_inf*Sutherland(T_W)
        n_w = m_w/r_w

        t_W0 = m_w*u_h/d_h
        u_t0 = sqrt(t_w0/r_w)

        ! initializin the loop
        x0    = u_t0
        xold1 = x0
        xold2 = x0 + 1.E-14_rp
        ctoll = 2._rp*toll
        iter  = 0

        do while ((ctoll.ge.toll).and.(iter.le.imax))
           iter = iter + 1
           
           ! calculting finite difference
           hk = (iReichardtLaw(xold1,d_h,u_h,n_w) - &
                 iReichardtLaw(xold2,d_h,u_h,n_w))/(xold1 - xold2)
           
           ! update x with Newton-Raphson formula
           xnew = xold1 - iReichardtLaw(xold1,d_h,u_h,n_w)/hk

           ! update tollerance
           ctoll = abs(xnew-xold1)

           ! update loop points
           xold2 = xold1
           xold1 = xnew
        enddo
        
        u_t = xnew

        tauWall = r_w*u_t**2

        return
end subroutine LogLawWMLES





subroutine compute_tau_wall_wmles1D(jW,phi_mean,tWall,tauW)
        
        use fluid_functions_module, only: Sutherland
        use parameters_module     , only: rp, gamma0, mu_inf, Lx, Lz, Prandtl
        use mpi_module            , only: sy, nx, nz
        use mesh_module           , only: y

        implicit none
        real(rp), dimension(:,:), allocatable, intent(in)  :: phi_mean
        real(rp)                             , intent(in)  :: tWall
        integer                              , intent(in)  :: jW
        real(rp)                             , intent(out) :: tauW
        real(rp)                                           :: qauW

        ! local declaration
        integer , parameter       :: nw = wmles_npts
        real(rp), dimension(1:nw) :: u_wm, T_wm
        
        real(rp) :: yw , hw , u_w, T_w, u_p, mW
        real(rp) :: r_h, irh, u_h, v_h, w_h, ekh, p_h, T_h, mlh
        real(rp) :: yPlWall0, yPlWalln, xPlWall0, zPlWall0, inuWall0
        real(rp) :: rhoWall0, utaWall0
        real(rp) :: dltWall0, tauWall0, h1, hl
        integer  :: j, jl, ji, j0,  j1, jInt
             
        jl = jW+4
        j1 = jW
        j0 = jW-1
        yW = 0.5_rp*(y(j1) + y(j0)) + 1.0_rp

        ! get wall quantities
        T_w = tWall
        u_w = 0.0_rp
        mW = Sutherland(t_w)

        ! get les quantities at the jl node
        r_h = phi_mean(jl,1)
        irh = 1.0_rp/r_h
        u_h = phi_mean(jl,2)*irh
        v_h = phi_mean(jl,3)*irh
        w_h = phi_mean(jl,4)*irh
        ekh = 0.5_rp*(u_h*u_h + v_h*v_h + w_h*w_h)
        p_h = (gamma0-1.0_rp) * (phi_mean(jl,5) - r_h*ekh)
        T_h = p_h*irh
        mlh = mu_inf*Sutherland(T_h)

        u_p = sqrt(u_h*u_h + w_h*w_h)

        !
        ! === compute yPlusW with Reichardt's law
        !
        h1 = abs(y(j1)-yW)
        hl = abs(y(jl)-yW)
        call LogLawWMLES(mu_inf,hl,T_w,u_p,p_h,tauWall0)
        rhoWall0 = p_h/T_w
        utaWall0 = sqrt(tauWall0/rhoWall0)

        inuWall0 = rhoWall0*utaWall0/(mW*mu_inf)

        xPlWall0 = Lx/real(nx,rp)*inuWall0
        yPlWall0 = h1         *inuWall0
        zPlWall0 = Lz/real(nz,rp)*inuWall0

        dltWall0 = yPlWall0/h1

        !if(yPlWall0 < 5.0_rp) then ! OLD VERSION
        if(xPlWall0 < xpt_wm .and. yPlWall0 < ypt_wm .and. zPlWall0 < zpt_wm) then
          tauW = rhoWall0*utaWall0**2

        else

          j = jW + wmles_strI
          do
            ji = jW + (j-1)
            yPlWalln = yPlWall0*(y(ji)-yW)/(y(sy)-yW)
            if(yPlWalln > wmles_intr) exit
            j = j+1
          enddo
          jInt = ji
          hW = abs(y(jInt) - yW)

          ! get les quantities
          r_h = phi_mean(jInt,1)
          irh = 1.0_rp/r_h
          u_h = phi_mean(jInt,2)*irh
          v_h = phi_mean(jInt,3)*irh
          w_h = phi_mean(jInt,4)*irh
          ekh = 0.5_rp*(u_h*u_h + v_h*v_h + w_h*w_h)
          p_h = (gamma0-1.0_rp) * (phi_mean(jInt,5) - r_h*ekh)
          T_h = p_h*irh

          u_p = sqrt(u_h*u_h + w_h*w_h)
          
          call OdeWMLES(gamma0,Prandtl,mu_inf,hw,dltWall0,u_w,u_p,T_w,T_h,p_h,u_wm,T_wm,tauW,qauW)

        endif
        
        return
end subroutine compute_tau_wall_wmles1D



subroutine compute_tau_wall_wmles2D(i,phi_mean,tWall,yc,u_wm,T_wm,tauW)
       
        use fluid_functions_module, only: Sutherland
        use parameters_module     , only: rp, gamma0, mu_inf, Lx, Lz, Prandtl
        use mpi_module            , only: sy, nx, nz
        use mesh_module           , only: y

        implicit none
        integer , parameter       :: nw = wmles_npts
        integer , parameter       :: jl = wmles_indx

        integer                                , intent(in)  :: i
        real(rp), dimension(:,:,:), allocatable, intent(in)  :: phi_mean
        real(rp)                               , intent(in)  :: tWall
        real(rp)                               , intent(out) :: tauW
        real(rp), dimension(1:nw)              , intent(out) :: yc
        real(rp), dimension(1:nw)              , intent(out) :: u_wm, T_wm
        real(rp)                                             :: qauW

        ! local declaration
        
        real(rp) :: yw , hw , u_w, T_w, mw, u_p
        real(rp) :: r_h, irh, u_h, v_h, w_h, ekh, p_h, T_h, mlh
        real(rp) :: h1, hl, rhoWall0, utaWall0, yPlWall0,dltWall0,tauWall0
        real(rp) :: yPlWalln, xPlWall0, zPlWall0, inuWall0
        real(rp) :: tauW_WR, tauW_WM, xpR, ypR, zpR, A, B
        integer  :: j, ji, j0, j1, jInt
             
        j1 = sy
        j0 = sy-1
        yW = 0.5_rp*(y(j1) + y(j0))

        ! get wall quantities
        T_w = tWall
        u_w = 0.0_rp
        mW = Sutherland(t_w)

        ! get les quantities at the jl node
        r_h = phi_mean(i,jl,1)
        irh = 1.0_rp/r_h
        u_h = phi_mean(i,jl,2)*irh
        v_h = phi_mean(i,jl,3)*irh
        w_h = phi_mean(i,jl,4)*irh
        ekh = 0.5_rp*(u_h*u_h + v_h*v_h + w_h*w_h)
        p_h = (gamma0-1.0_rp) * (phi_mean(i,jl,5) - r_h*ekh)
        T_h = p_h*irh
        mlh = mu_inf*Sutherland(T_h)

        u_p = sqrt(u_h*u_h + w_h*w_h)
        
        ! estimate near-wall resolution with Reichardt's law
        h1 = abs(y(j1)-yW)
        hl = abs(y(jl)-yW)
        call LogLawWMLES(mu_inf,hl,T_w,u_p,p_h,tauWall0)
        rhoWall0 = p_h/T_w
        utaWall0 = sqrt(tauWall0/rhoWall0)

        inuWall0 = rhoWall0*utaWall0/(mW*mu_inf)

        xPlWall0 = Lx/real(nx,rp)*inuWall0
        yPlWall0 = h1         *inuWall0
        zPlWall0 = Lz/real(nz,rp)*inuWall0

        dltWall0 = yPlWall0/h1
        
        ! get the interface location
        j = wmles_strI
        do
          ji = sy + (j-1)
          yPlWalln = yPlWall0*(y(ji)-yW)/(y(sy)-yW)
          if(yPlWalln > wmles_intr) exit
          j = j+1
        enddo
        jInt = ji
        hW = abs(y(jInt) - yW)

        ! get les quantities
        r_h = phi_mean(i,jInt,1)
        irh = 1.0_rp/r_h
        u_h = phi_mean(i,jInt,2)*irh
        v_h = phi_mean(i,jInt,3)*irh
        w_h = phi_mean(i,jInt,4)*irh
        ekh = 0.5_rp*(u_h*u_h + v_h*v_h + w_h*w_h)
        p_h = (gamma0-1.0_rp) * (phi_mean(i,jInt,5) - r_h*ekh)
        T_h = p_h*irh

        ! compute tauW_WM
        u_p = sqrt(u_h*u_h + w_h*w_h)
        call OdeWMLES(gamma0,Prandtl,mu_inf,hw,dltWall0,u_w,u_p,T_w,T_h,p_h,u_wm,T_wm,tauW_WM,qauW)

        ! compute tauW_WR
        tauW_WR = mW*mu_inf*phi_mean(i,sy,2)/phi_mean(i,sy,1)/h1

        ! blending tauW_WR and tauW_WR
        xpR = xPlWall0/xpt_wm
        ypR = yPlWall0/ypt_wm
        zpR = zPlWall0/zpt_wm

        A = max(1.0_rp-xpR,0.0_rp) + max(1.0_rp-ypR,0.0_rp) + max(1.0_rp-zpR,0.0_rp)
        B = min(xpR,1.0_rp) + min(ypR,1.0_rp) + min(zpR,1.0_rp)
        A = A/3.0_rp
        B = B/3.0_rp

        tauW = tauW_WR*A + tauW_WM*B
        
        return
end subroutine compute_tau_wall_wmles2D





subroutine compute_tau_wall_improved_wmles2D(i,phi_mean,tWall,yc,u_wm,T_wm,tauW)
       
        use fluid_functions_module, only: Sutherland
        use parameters_module     , only: rp, gamma0, mu_inf, Prandtl
        use mpi_module            , only: sy, nx, nz
        use mesh_module           , only: y

        implicit none
        integer , parameter       :: nw = wmles_npts

        integer                                , intent(in)  :: i
        real(rp), dimension(:,:,:), allocatable, intent(in)  :: phi_mean
        real(rp)                               , intent(in)  :: tWall
        real(rp)                               , intent(out) :: tauW
        real(rp), dimension(1:nw)              , intent(out) :: yc
        real(rp), dimension(1:nw)              , intent(out) :: u_wm, T_wm
        real(rp)                                             :: qauW

        ! local declaration
        real(rp), dimension(0:nw) :: yf   !> face grid
        real(rp), dimension(0:nw) :: dy   !> step grid
        
        real(rp) :: yw , hw , u_w, T_w, mw, u_p
        real(rp) :: r_h, irh, u_h, v_h, w_h, ekh, p_h, T_h
        real(rp) :: r_1, ir1, u_1, v_1, w_1, ek1, p_1, T_1
        real(rp) :: rhoWall0, utaWall0, yPlWall0,dltWall0,tauWall0, dyh
        real(rp) :: yPlWalln, inuWall0
        real(rp) :: tauW_WR
        integer  :: j, ji, j0, j1, jInt
             
        j1 = sy
        j0 = sy-1
        yW = 0.5_rp*(y(j1) + y(j0))
        dyh= 0.5_rp*(y(j1) - y(j0))

        ! get wall quantities
        T_w = tWall
        u_w = 0.0_rp
        mW = Sutherland(t_w)

        !
        ! get the numerical shear stress and heat flux
        !
        r_1 = phi_mean(i,j1,1)
        ir1 = 1.0_rp/r_1
        u_1 = phi_mean(i,j1,2)*ir1
        v_1 = phi_mean(i,j1,3)*ir1
        w_1 = phi_mean(i,j1,4)*ir1
        ek1 = 0.5_rp*(u_1*u_1 + v_1*v_1 + w_1*w_1)
        p_1 = (gamma0-1.0_rp) * (phi_mean(i,j1,5) - r_1*ek1)
        T_1 = p_1*ir1

        tauW_WR = mu_inf*mW*u_1/dyh
        !
        ! get the viscous length
        !
        tauWall0 = tauW_WR!mu_inf*mW*u_1/dyh
        rhoWall0 = p_1/T_w
        utaWall0 = sqrt(abs(tauWall0)/rhoWall0)
        inuWall0 = rhoWall0*utaWall0/(mW*mu_inf)
        !
        ! get internal resolutions estimation
        !
        yPlWall0 = dyh*inuWall0
        dltWall0 = yPlWall0/dyh
        !
        ! === looking at the grid point where yPlus>wmles_interface
        !
        do j = wmles_min_id,wmles_max_id
          ji = sy + j-1
          yPlWalln = yPlWall0*(y(ji)-yW)/(y(j1)-yW)
          if(yPlWalln > wmles_intr) exit
        enddo
        jInt = j
        hw = abs(y(jInt) - yW)

        call Get_WMLES_Grid(nw,hw,dltWall0,yf,yc,dy)
        !
        ! === get the LES field at the interface location
        !
        r_h = phi_mean(i,jInt,1)
        irh = 1.0_rp/r_h
        u_h = phi_mean(i,jInt,2)*irh
        v_h = phi_mean(i,jInt,3)*irh
        w_h = phi_mean(i,jInt,4)*irh
        ekh = 0.5_rp*(u_h*u_h + v_h*v_h + w_h*w_h)
        p_h = (gamma0-1.0_rp) * (phi_mean(i,jInt,5) - r_h*ekh)
        T_h = p_h*irh

        u_p = sqrt(u_h*u_h + w_h*w_h)

        !
        ! === get the correct tauWall and qWall from ODE model
        !
        call OdeWMLES(gamma0,Prandtl,mu_inf,hw,dltWall0,u_w,u_p,T_w,T_h,p_h,u_wm,T_wm,tauW,qauW)

        
        return
end subroutine compute_tau_wall_improved_wmles2D








subroutine compute_tau_wall_wmles2D_static(i,phi_mean,tWall,yc,u_wm,T_wm,tauW)
       
        use fluid_functions_module, only: Sutherland
        use parameters_module     , only: rp, gamma0, mu_inf, Prandtl
        use mpi_module            , only: sy
        use mesh_module           , only: y

        implicit none
        integer , parameter :: nw = wmles_npts

        integer                                , intent(in)  :: i
        real(rp), dimension(:,:,:), allocatable, intent(in)  :: phi_mean
        real(rp)                               , intent(in)  :: tWall
        real(rp)                               , intent(out) :: tauW
        real(rp), dimension(1:nw)              , intent(out) :: yc, u_wm, T_wm
        real(rp)                                             :: qauW

        ! local declaration
        
        real(rp) :: yw , hw , u_w, T_w, mw, u_p
        real(rp) :: r_h, irh, u_h, v_h, w_h, ekh, p_h, T_h, mlh
        real(rp) :: hl, rhoWall0, utaWall0,dltWall0,tauWall0
        integer  :: j0, j1, jInt
             
        j1 = sy
        j0 = sy-1
        yW = 0.5_rp*(y(j1) + y(j0))
        jInt = wmles_strI

        ! get wall quantities
        T_w = tWall
        u_w = 0.0_rp
        mW = Sutherland(t_w)

        ! get les quantities at the jInt node
        r_h = phi_mean(i,jInt,1)
        irh = 1.0_rp/r_h
        u_h = phi_mean(i,jInt,2)*irh
        v_h = phi_mean(i,jInt,3)*irh
        w_h = phi_mean(i,jInt,4)*irh
        ekh = 0.5_rp*(u_h*u_h + v_h*v_h + w_h*w_h)
        p_h = (gamma0-1.0_rp) * (phi_mean(i,jInt,5) - r_h*ekh)
        T_h = p_h*irh
        mlh = mu_inf*Sutherland(T_h)

        u_p = sqrt(u_h*u_h + w_h*w_h)
        
        !
        ! === compute yPlusW with Reichardt's law
        !
        hl = abs(y(jInt)-yW)
        call LogLawWMLES(mu_inf,hl,T_w,u_p,p_h,tauWall0)
        rhoWall0 = p_h/T_w
        utaWall0 = sqrt(tauWall0/rhoWall0)
        dltWall0 = rhoWall0*utaWall0/(mW*mu_inf)

        hW = abs(y(jInt) - yW)
        
        ! get tauWall and qWall
        if(u_h<0.0_rp) u_p = -u_p
        call OdeWMLES(gamma0,Prandtl,mu_inf,hw,dltWall0,u_w,u_p,T_w,T_h,p_h,u_wm,T_wm,tauW,qauW)

        return
end subroutine compute_tau_wall_wmles2D_static





subroutine compute_tau_wall_wmles2D_vorticity(i,phi_mean,phi_mean_aux,tWall,yc,u_wm,T_wm,tauW)
       
        use fluid_functions_module, only: Sutherland
        use parameters_module     , only: rp, gamma0, mu_inf, Prandtl
        use real_to_integer_module, only: nearest_integer_opt
        use mpi_module            , only: sy, ey
        use mesh_module           , only: y

        implicit none
        integer , parameter :: nw = wmles_npts

        integer                                , intent(in)  :: i
        real(rp), dimension(:,:,:), allocatable, intent(in)  :: phi_mean
        real(rp), dimension(:,:,:), allocatable, intent(in)  :: phi_mean_aux
        real(rp)                               , intent(in)  :: tWall
        real(rp)                               , intent(out) :: tauW
        real(rp), dimension(1:nw)              , intent(out) :: yc, u_wm, T_wm
        real(rp)                                             :: qauW

        ! local declaration
        
        real(rp) :: yw , hw , u_w, T_w, mw, u_p
        real(rp) :: r_h, irh, u_h, v_h, w_h, ekh, p_h, T_h, mlh
        real(rp) :: hl, rhoWall0, utaWall0,dltWall0,tauWall0
        real(rp) :: om, yTgT
        integer  :: j, j0, j1, jInt, jOme
             
        j1 = sy
        j0 = sy-1
        yW = 0.5_rp*(y(j1) + y(j0))
        
        ! get the interface location throught vorticity threshold
        jOme = sy
        do j = sy,ey
           om = phi_mean_aux(i,j,19)
           if(om < omg_wm) then
             jOme = j-1
             exit
           endif
        enddo
        yTgT = abs(y(jOme)-yW)/3.0_rp
        jInt = nearest_integer_opt(y,sy,ey,yTgT)
        if(jInt > wmles_max_id) jInt = wmles_max_id
        if(jInt < wmles_min_id) jInt = wmles_min_id

        ! get wall quantities
        T_w = tWall
        u_w = 0.0_rp
        mW = Sutherland(t_w)

        ! get les quantities at the jInt node
        r_h = phi_mean(i,jInt,1)
        irh = 1.0_rp/r_h
        u_h = phi_mean(i,jInt,2)*irh
        v_h = phi_mean(i,jInt,3)*irh
        w_h = phi_mean(i,jInt,4)*irh
        ekh = 0.5_rp*(u_h*u_h + v_h*v_h + w_h*w_h)
        p_h = (gamma0-1.0_rp) * (phi_mean(i,jInt,5) - r_h*ekh)
        T_h = p_h*irh
        mlh = mu_inf*Sutherland(T_h)

        u_p = sqrt(u_h*u_h + w_h*w_h)
        
        !
        ! === compute yPlusW with Reichardt's law
        !
        hl = abs(y(jInt)-yW)
        call LogLawWMLES(mu_inf,hl,T_w,u_p,p_h,tauWall0)
        rhoWall0 = p_h/T_w
        utaWall0 = sqrt(tauWall0/rhoWall0)
        dltWall0 = rhoWall0*utaWall0/(mW*mu_inf)

        hW = abs(y(jInt) - yW)
        
        ! get tauWall and qWall
        if(u_h<0.0_rp) u_p = -u_p
        call OdeWMLES(gamma0,Prandtl,mu_inf,hw,dltWall0,u_w,u_p,T_w,T_h,p_h,u_wm,T_wm,tauW,qauW)

        return
end subroutine compute_tau_wall_wmles2D_vorticity











end module wmles_module
