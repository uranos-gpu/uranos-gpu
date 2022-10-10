module fluid_functions_module
use parameters_module
use storage_module
implicit none
contains

elemental function Sutherland(t) result(mu_t)
! -------------------------------------------------------------
!       Sutherland's Law valid till T/T_ref = 3.5_rp
! -------------------------------------------------------------
        !$acc routine seq
        use parameters_module, only: rp

        implicit none

        real(rp), parameter :: suthc = 110.4_rp/273.15_rp

        real(rp), intent(in) :: T
        real(rp)             :: mu_t

        mu_t = T*sqrt(T) * (1.0_rp+suthc)/(T+suthc)

        return
end function Sutherland



elemental function getmuT(rho,tauw,yj,mul) result(mut)
        !$acc routine seq

        implicit none
        real(rp), intent(in) :: rho,tauw,yj,mul
        real(rp)             :: mut

        real(rp), parameter :: vkc = 0.41_rp !0.384_rp
        real(rp), parameter :: A_p = 17.0_rp !15.5_rp
        real(rp)            :: y_s
        real(rp)            :: DvD
        real(rp)            :: abstauw

        abstauw = abs(tauW)

        y_s = yj*sqrt(rho*abstauw)/mul
        DvD = (1.0_rp-exp(-y_s/A_p))**2
        mut = vkc*rho*sqrt(abstauw/rho)*yj*DvD

        return
end function getmuT


elemental function KarmanSchoenherr(ReThI) result(cfi)
        implicit none
        real(rp), intent(in) :: ReThI   !< inc. momentum thickness Reynolds number
        real(rp)             :: cfi     !> inc. friction coefficient

        cfi = 1.0_rp/(17.08_rp*(log10(ReThI))**2 + 25.11_rp*log10(ReThI) + 6.012_rp)

        return
end function KarmanSchoenherr

function ReThetaIncSpalding(uep) result(ReThInc)
        implicit none
        real(rp), intent(in) :: uep
        real(rp)             :: ReThInc

        real(rp), parameter :: k = 0.4_rp, E = 12.0_rp
        real(rp)            :: a1, a2, a3, a4, a5, a6, a7, a8

        !a1 = uep**2/6.0_rp
        !a2 = (1.0_rp-2.0_rp/(k*uep))*exp(k*uep)
        !a3 = 2.0_rp/(k*uep)
        !a4 = 1.0_rp
        !a5 = -(uep)**2/6.0_rp
        !a6 = -(uep)**2/12.0_rp
        !a7 = -(uep)**2/40.0_rp
        !a8 = -(uep)**2/180.0_rp

        a1 = uep**2/6
        a2 = (1.0_rp-2/(k*uep))*exp(k*uep)
        a3 = 2.0_rp/(k*uep)
        a4 = 1.0_rp
        a5 = -(k*uep)**2/6.0_rp
        a6 = -(k*uep)**3/12.0_rp
        a7 = -(k*uep)**4/40.0_rp
        a8 = -(k*uep)**5/180.0_rp

        ReThInc = a1 + (a2+a3+a4+a5+a6+a7+a8)/(k*E)


        return
end function ReThetaIncSpalding




elemental function divergency(i,j,k) result(div)
! -----------------------------------------------------------------
!       This function compute the speed divergency in a i,j,k point
! -----------------------------------------------------------------
        implicit none
        integer , intent(in) :: i,j,k
        real(rp)             :: div 

        ! local declarations
        real(rp), dimension(-1:1), parameter :: c1 = (/-0.5_rp,0.0_rp,0.5_rp/)
        real(rp)                             :: i_dx, i_dy, i_dz, cs
        real(rp)                             :: u_x, v_y, w_z
        real(rp)                             :: i_rho_is, i_rho_js, i_rho_ks
        integer                              :: s, is, js, ks

        i_dx = xstep_i(i)
        i_dy = ystep_i(j)
        i_dz = zstep_i(k)
        
        u_x = 0.0_rp
        v_y = 0.0_rp
        w_z = 0.0_rp
        do s = -1, 1

           cs = c1(s)

           is = i+s
           js = j+s
           ks = k+s

           i_rho_is = 1.0_rp/phi(is,j,k,1)
           i_rho_js = 1.0_rp/phi(i,js,k,1)
           i_rho_ks = 1.0_rp/phi(i,j,ks,1)

           u_x = u_x + i_dx * cs * phi(is,j,k,2) * i_rho_is
           v_y = v_y + i_dy * cs * phi(i,js,k,3) * i_rho_js
           w_z = w_z + i_dz * cs * phi(i,j,ks,4) * i_rho_ks

        enddo

        div = u_x + v_y + w_z
        
        return
end function divergency


function vort(i,j,k) result(vor)
! -----------------------------------------------------------------
!       This function compute the speed vorticity in a i,j,k point
! -----------------------------------------------------------------
        implicit none
        integer, intent(in)    :: i,j,k
        real(rp), dimension(3) :: vor

        ! local declaration
        real(rp), dimension(-1:1), parameter :: c1 = (/-0.5_rp,0.0_rp,0.5_rp/)
        real(rp)                             :: i_dx, i_dy, i_dz, cs
        real(rp)                             :: u_y, u_z, v_x, v_z, w_x, w_y
        real(rp)                             :: i_rho_is, i_rho_js, i_rho_ks
        integer                              :: s, is, js, ks

        i_dx = xstep_i(i)
        i_dy = ystep_i(j)
        i_dz = zstep_i(k)
        
        u_y = 0.0_rp ; u_z = 0.0_rp
        v_x = 0.0_rp ; v_z = 0.0_rp
        w_x = 0.0_rp ; w_y = 0.0_rp
        do s = -1, 1

           is = i+s
           js = j+s
           ks = k+s

           cs = c1(s)

           i_rho_is = 1.0_rp/phi(is,j,k,1)
           i_rho_js = 1.0_rp/phi(i,js,k,1)
           i_rho_ks = 1.0_rp/phi(i,j,ks,1)

           u_y = u_y + i_dy * cs * phi(i,js,k,2) * i_rho_js
           u_z = u_z + i_dz * cs * phi(i,j,ks,2) * i_rho_ks

           v_x = v_x + i_dx * cs * phi(is,j,k,3) * i_rho_is
           v_z = v_z + i_dz * cs * phi(i,j,ks,3) * i_rho_ks

           w_x = w_x + i_dx * cs * phi(is,j,k,4) * i_rho_is
           w_y = w_y + i_dy * cs * phi(i,js,k,4) * i_rho_js

        enddo

        vor(1) = w_y - v_z
        vor(2) = u_z - w_x
        vor(3) = v_x - u_y

        return
end function vort


subroutine Rankine_Hugoniot(rho1, p1, T1, M1, rho2, p2, T2 ,M2)
! --------------------------------------------------------------------------------
!
!       This subroutine applies Rankine-Hugoniot relations in order to compute
!       a post-shock termo-dynamical state
!
! --------------------------------------------------------------------------------
        implicit none
        real(rp), intent(in)  :: rho1, p1, T1, M1        !< pre -shock state
        real(rp)              :: phi_r, phi_p, phi_T     !< shock quotients
        real(rp), intent(out) :: rho2, p2, T2, M2        !< post-shock state
        
        ! shock quotients
        phi_p = (1.0_rp + (2*gamma0)/(gamma0+1.0_rp)*(M1**2 - 1.0_rp))

        phi_r = ((gamma0+1.0_rp)*M1**2)/(2.0_rp + (gamma0-1.0_rp)*M1**2)

        phi_T = phi_p / phi_r

        ! post shock state
        M2   = sqrt((1.0_rp+0.5_rp*(gamma0-1._rp)*M1**2)/(gamma0*M1**2 - 0.5_rp*(gamma0-1.0_rp)))
        rho2 = phi_r * rho1
        p2   = phi_p * p1
        T2   = phi_T * T1
return
end subroutine Rankine_Hugoniot


function inverse_rankine_hugoniot(M1) result(f_M1)
! ----------------------------------------------------------------------------------------
!
!       This function computes the shock-wave Mach number based on 
!       R-H equilibrium, fixing the post-shocked condition, instead of the pre-shocked one
!
! ----------------------------------------------------------------------------------------

        use parameters_module, only: rp, Mach

        implicit none
        real(rp), intent(in) :: M1
        real(rp)             :: f_M1

        ! local declarations
        real(rp) :: phi_p, phi_r, phi_T, M2

        phi_p = (1.0_rp + (2*gamma0)/(gamma0+1.0_rp)*(M1**2 - 1.0_rp))

        phi_r = ((gamma0+1.0_rp)*M1**2)/(2.0_rp + (gamma0-1.0_rp)*M1**2)

        phi_T = phi_p / phi_r

        M2   = sqrt((1.0_rp+0.5_rp*(gamma0-1._rp)*M1**2)/(gamma0*M1**2 - 0.5_rp*(gamma0-1.0_rp)))
        
        ! equation to be solved
        f_M1 = Mach - phi_T**(-0.5_rp)*M1 + M2

        return
end function inverse_rankine_hugoniot




function ThetaBetaMach(beta,theta,mach) result(fbeta)
! -------------------------------------------------------------------
!
!       This function implements the Theta-Beta-Mach function
!       for an oblique shock over a wedge
!
!       INPUT:  beta     !< shock angle
!               theta    !< wedge angle
!               Mach     !< free stream Mach number
!
!       OUTPUT: fbeta    !< Theta-Beta-Mach function
!
! -------------------------------------------------------------------
        implicit none
        real(rp), intent(in) :: theta
        real(rp), intent(in) :: beta
        real(rp), intent(in) :: Mach
        real(rp)             :: fbeta

        ! local declaration
        real(rp), parameter :: gamma = 1.4_rp
        real(rp)            :: cot_beta, nnum, dnum

        cot_beta = 1.0_rp/tan(beta)

        nnum = Mach**2 * sin(beta)**2 - 1.0_rp
        dnum = Mach**2 * (gamma + cos(2*beta)) + 2

        fbeta = tan(theta) - 2*cot_beta * nnum/dnum 

        return
end function ThetaBetaMach

function direct_PrandtlMeyer(M) result(vM)
! --------------------------------------------------------
!       This function provides a direct computation of the 
!       Prandtl-Meyer function
!
!       INPUT : M        !< Mach number
!       OUTPUT: vM       !> v(M) Prandlt-Meyer function
!
! ------------------------------------------------
        implicit none
        real(rp), intent(in) :: M
        real(rp)             :: vM

        !local declarations
        real(rp), parameter  :: g = 1.4_rp
        real(rp)             :: par, i_par, x
        
        par   = sqrt((g + 1.0_rp)/(g - 1.0_rp))
        i_par = 1.0_rp/par

        x     = sqrt(M**2 - 1.0_rp)
        vM    = par * atan(i_par * x) - atan(x)

        return
end function direct_PrandtlMeyer

function inverse_PrandtlMeyer(M2,M1,theta) result(fM2)
! -----------------------------------------------------
!       This function solves the Prandlt-Meyer equation
!       f(M2) = v(M2) - theta - v(M1) = 0
! -----------------------------------------------------
        implicit none
        real(rp), intent(in) :: M2
        real(rp), intent(in) :: M1
        real(rp), intent(in) :: theta 
        real(rp)             :: fM2

        ! local declaration
        real(rp)             :: vM1, vM2

        vM1 = direct_PrandtlMeyer(M1)
        
        vM2 = direct_PrandtlMeyer(M2)

        fM2 = vM2 - theta - vM1

        return
end function inverse_PrandtlMeyer



subroutine oblique_shock(r1,p1,T1,M1,theta,r2,p2,T2,M2,beta)
! ----------------------------------------------------------------
!       
!       This subroutine provides the compatibility relation for an
!       oblique shock wave over a wedge
!
!       INPUT : r1       !< density state 1
!               p1       !< pressure state 1
!               T1       !< temperature state 1
!               M1       !< Mach state 1
!               theta    !< wedge angle
!       OUTPUT: r2       !> density state 2
!               p2       !> pressure state 2
!               T2       !> Temperature state 2
!               M2       !> Mach state 2
!               beta     !> shock angle
!
! ----------------------------------------------------------------

        use math_tools_module, only: newton_raphson

        implicit none
        real(rp), intent(in ) :: r1
        real(rp), intent(in ) :: p1
        real(rp), intent(in ) :: T1
        real(rp), intent(in ) :: M1
        real(rp), intent(in ) :: theta           

        real(rp), intent(out) :: r2
        real(rp), intent(out) :: p2
        real(rp), intent(out) :: T2
        real(rp), intent(out) :: M2
        real(rp), intent(out) :: beta

        ! local declarations
        real(rp)            :: Mn1, Mn2
        real(rp), parameter :: toll   = 1.0E-14_rp
        integer , parameter :: itmax  = 1000
        
        ! compute the shock angle beta
        call newton_raphson(ThetaBetaMach,theta,M1,itmax,toll,theta,beta)

        ! normal component of the Mach number at state 1
        Mn1 = M1 * sin(beta)

        ! STATE 2
        call Rankine_Hugoniot(r1, p1, T1, Mn1, r2, p2, T2, Mn2)

        M2 = Mn2 / sin(beta - theta)

        return
end subroutine oblique_shock




subroutine ComputeStress2D(xstep_i,ystep_i,U,V,VIS,sigma)
        
        use parameters_module, only: rp, mu_inf
        use storage_module   , only: central_fd_order, central_1

        implicit none
        real(rp), dimension(:)      , allocatable, intent(in)    :: xstep_i, ystep_i
        real(rp), dimension(:,:,:)  , allocatable, intent(in)    :: U,V,VIS
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: sigma

        ! local declarations
        real(rp), dimension(2,2) :: grad_v, grad_vT, div_vI, stress_tensor
        real(rp), dimension(2)   :: du, dv
        real(rp), dimension(2)   :: i_st, cl1
        real(rp), parameter      :: two3 = 2.0_rp/3.0_rp
        real(rp)                 :: mu
        integer                  :: i,j,k,l, fR

        fR = central_fd_order/2

        k = 1
        !$omp parallel do collapse(2) default(private), shared(sx,ex,sy,ey,xstep_i,ystep_i), &
        !$omp shared(U,V,VIS,sigma,k,mu_inf,central_1,fR)
        do j    = sy,ey
           do i = sx,ex

              ! grid steps
              i_st(1) = xstep_i(i)
              i_st(2) = ystep_i(j)

              ! viscosity
              mu = VIS(i,j,k)

              ! === first derivatives
              du = 0.0_rp
              dv = 0.0_rp
              do l = 1,fR
                 cl1(:) = i_st(:) * central_1(l)

                 du(1) = du(1) + cl1(1) * (U(i+l,j,k) - U(i-l,j,k))
                 du(2) = du(2) + cl1(2) * (U(i,j+l,k) - U(i,j-l,k))

                 dv(1) = dv(1) + cl1(1) * (V(i+l,j,k) - V(i-l,j,k))
                 dv(2) = dv(2) + cl1(2) * (V(i,j+l,k) - V(i,j-l,k))

              enddo

              ! === calculating stress tensor
              grad_v(1,:)  = du
              grad_v(2,:)  = dv

              grad_vT(:,1) = du
              grad_vT(:,2) = dv
              
              div_vI = 0.0_rp
              div_vI(1,1) = du(1) + dv(2)
              div_vI(2,2) = du(1) + dv(2)

              stress_tensor = mu_inf*mu*(grad_v + grad_vT - two3 * div_vI)
        
              ! === storing sigma
              sigma(i,j,k,1) = stress_tensor(1,1)
              sigma(i,j,k,2) = stress_tensor(1,2)
              sigma(i,j,k,3) = stress_tensor(2,1)
              sigma(i,j,k,4) = stress_tensor(2,2)

          enddo
        enddo
        !$omp end parallel do


        return
end subroutine ComputeStress2D



subroutine hTurbEpsilon(xstep,ystep,zstep,phi,U,V,W,VIS,gbl_rh_mean,&
                gbl_eps_sol,gbl_eps_dil,gbl_eps_tot,gbl_eps_str,mpi_flag)
! --------------------------------------------------------------------------------------------------
!
!       This function computes the Turbulent dissipation components for an hTurb case
!
!       INPUT:  xstep           !< grid x step
!               ystep           !< grid x step
!               zstep           !< grid x step
!               U,V,W           !< velocity fields
!               gbl_rh_mean     !< gbobal mean density 
!               mpi_flag        !< flag handling mpi
!
!       OUTPUT: gbl_eps_dil     !> dilatation epsilon component
!               gbl_eps_sol     !> solenoidal epsilon component
!               gbl_eps_tot     !> total dissipation 
!               gbl_eps_str     !> total dissipation (strain calculus)
!
! --------------------------------------------------------------------------------------------------

        use parameters_module, only: rp
        use storage_module   , only: central_fd_order, central_1

        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(in)  :: phi
        real(rp), dimension(:,:,:)  , allocatable, intent(in)  :: U,V,W,VIS
        real(rp), dimension(:)      , allocatable, intent(in)  :: xstep, ystep, zstep
        real(rp)                                 , intent(in)  :: gbl_rh_mean
        logical                                  , intent(in)  :: mpi_flag
        real(rp)                                 , intent(out) :: gbl_eps_sol
        real(rp)                                 , intent(out) :: gbl_eps_dil
        real(rp)                                 , intent(out) :: gbl_eps_tot
        real(rp)                                 , intent(out) :: gbl_eps_str
        
        !local declarations
        real(rp)                 :: lcl_eps_sol, lcl_eps_dil
        real(rp), dimension(3,3) :: lcl_S2_mean, gbl_S2_mean
        real(rp)                 :: lcl_Vtot_sp
        real(rp)                 :: gbl_Vtot_sp

        real(rp), dimension(3,3) :: S, grad_v, grad_vT
        real(rp), dimension(3)   :: st, i_st, cl1
        real(rp), dimension(3)   :: du, dv, dw
        real(rp), dimension(3)   :: vor
        real(rp)                 :: vor2, div2, mu, nu
        real(rp)                 :: dVol, irhVtot, iVtot
        integer                  :: fR, l
        integer                  :: i,j,k
        integer, dimension(4)    :: err = 0

        fR   = central_fd_order/2


        lcl_Vtot_sp = 0.0_rp
        lcl_eps_sol = 0.0_rp
        lcl_eps_dil = 0.0_rp
        lcl_S2_mean = 0.0_rp
        
        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex
                
                 ! grid data
                 st(1) = xstep(i)
                 st(2) = ystep(j)
                 st(3) = zstep(k)

                 i_st(:) = 1.0_rp/st(:)

                 dVol = st(1)*st(2)*st(3)
                 lcl_Vtot_sp = lcl_Vtot_sp + dVol

                 ! viscosities
                 mu = mu_inf * VIS(i,j,k)
                 nu = mu/phi(i,j,k,1)

                 ! compute derivatives 
                 du(:) = 0.0_rp
                 dv(:) = 0.0_rp
                 dw(:) = 0.0_rp
                 do l = 1 , fR

                    cl1(:) = i_st(:) * central_1(l)

                    du(1) = du(1) + cl1(1) * (U(i+l,j,k) - U(i-l,j,k))
                    du(2) = du(2) + cl1(2) * (U(i,j+l,k) - U(i,j-l,k))
                    du(3) = du(3) + cl1(3) * (U(i,j,k+l) - U(i,j,k-l))

                    dv(1) = dv(1) + cl1(1) * (V(i+l,j,k) - V(i-l,j,k))
                    dv(2) = dv(2) + cl1(2) * (V(i,j+l,k) - V(i,j-l,k))
                    dv(3) = dv(3) + cl1(3) * (V(i,j,k+l) - V(i,j,k-l))

                    dw(1) = dw(1) + cl1(1) * (W(i+l,j,k) - W(i-l,j,k))
                    dw(2) = dw(2) + cl1(2) * (W(i,j+l,k) - W(i,j-l,k))
                    dw(3) = dw(3) + cl1(3) * (W(i,j,k+l) - W(i,j,k-l))

                 enddo

                 ! Strain tensor
                 grad_v (1,:) = du(:)
                 grad_v (2,:) = dv(:)
                 grad_v (3,:) = dw(:)

                 grad_vT(:,1) = du(:)
                 grad_vT(:,2) = dv(:)
                 grad_vT(:,3) = dw(:)

                 S = 0.5_rp*(grad_v + grad_vT)

                 ! mean S components
                 lcl_S2_mean(:,:) = lcl_S2_mean(:,:) + nu * S(:,:)*S(:,:) * dVol

                 ! vorticity module
                 vor(1) = dw(2) - dv(3)
                 vor(2) = du(3) - dw(1)
                 vor(3) = dv(1) - du(2)

                 vor2 = vor(1)*vor(1) + vor(2)*vor(2) + vor(3)*vor(3)

                 ! divergency module
                 div2 = (du(1) + dv(2) + dw(3))**2

                 ! compute locol dissapation
                 lcl_eps_sol = lcl_eps_sol + mu*vor2*dVol
                 lcl_eps_dil = lcl_eps_dil + mu*div2*dVol

              enddo
           enddo
        enddo
        

        if(mpi_flag) then
          call MPI_allreduce(lcl_Vtot_sp, gbl_Vtot_sp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_cart, err(1))
          call MPI_allreduce(lcl_eps_sol, gbl_eps_sol, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_cart, err(2))
          call MPI_allreduce(lcl_eps_dil, gbl_eps_dil, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_cart, err(3))
          call MPI_allreduce(lcl_S2_mean, gbl_S2_mean, 9, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_cart, err(4))
          if(sum(err) .ne. 0) stop ' MPI error in hTurbEpsilon'

        else
          gbl_Vtot_sp = lcl_Vtot_sp
          gbl_eps_sol = lcl_eps_sol
          gbl_eps_dil = lcl_eps_dil
          gbl_S2_mean = lcl_S2_mean

        endif
        
        iVtot       = 1.0_rp/gbl_Vtot_sp
        gbl_S2_mean = gbl_S2_mean * iVtot

        irhVtot     = 1.0_rp/(gbl_rh_mean * gbl_Vtot_sp)

        gbl_eps_sol = gbl_eps_sol*irhVtot
        gbl_eps_dil = gbl_eps_dil*irhVtot*4.0_rp/3.0_rp
        gbl_eps_tot = gbl_eps_sol + gbl_eps_dil
        gbl_eps_str = 2*sum(gbl_S2_mean)

        return
end subroutine hTurbEpsilon



subroutine BlasiusProfile(BlasiusVel, BlasiusTmp,  BlasiusRho)

        use interpolation_module

        implicit none
        real(rp), dimension(1-GN:ny+GN), intent(out) :: BlasiusVel
        real(rp), dimension(1-GN:ny+GN), intent(out) :: BlasiusTmp
        real(rp), dimension(1-GN:ny+GN), intent(out) :: BlasiusRho

        real(rp), parameter :: eta_max = 10.0_rp
        integer             :: j, err = 0, ntot

        real(rp), dimension(:), allocatable :: etaC
        real(rp), dimension(:), allocatable :: f0_lb, f1_lb, f2_lb
        real(rp), dimension(:), allocatable :: f0_ub, f1_ub, f2_ub
        real(rp), dimension(:), allocatable :: f0_mb, f1_mb, f2_mb

        real(rp), dimension(:), allocatable :: g0_mb, g1_mb
        real(rp), dimension(:), allocatable :: g0_ub, g1_ub
        real(rp), dimension(:), allocatable :: g0_lb, g1_lb
        
        real(rp), dimension(:), allocatable :: y_tmp, y_tmp1
        real(rp)                            :: erri, d_eta, yo, B
        real(rp)                            :: f1_int, g0_int
        real(rp), parameter                 :: toll = 1.0E-14_rp
        real(rp), parameter                 :: alp  = 1.0_rp
        integer                             :: iter, itermax = 1000, nhuge = 10000

        real(rp) :: dtheta, theta, theta1, theta2
        integer  :: n1, ntheta

        
        ! ==== allocate f blasius functions for momentum
        err = 0
        allocate(etaC(1:nhuge)                                 ,stat=err)
        if(err .ne. 0) stop ' Allocation error in BlasiusProfile'
        allocate(f0_lb(1:nhuge), f1_lb(1:nhuge), f2_lb(1:nhuge),stat=err)
        if(err .ne. 0) stop ' Allocation error in BlasiusProfile'
        allocate(f0_ub(1:nhuge), f1_ub(1:nhuge), f2_ub(1:nhuge),stat=err)
        if(err .ne. 0) stop ' Allocation error in BlasiusProfile'
        allocate(f0_mb(1:nhuge), f1_mb(1:nhuge), f2_mb(1:nhuge),stat=err)
        if(err .ne. 0) stop ' Allocation error in BlasiusProfile'

        f0_lb = 0.0_rp; f1_lb = 0.0_rp; f2_lb = 0.0_rp
        f0_ub = 0.0_rp; f1_ub = 0.0_rp; f2_ub = 0.0_rp
        f0_mb = 0.0_rp; f1_mb = 0.0_rp; f2_mb = 0.0_rp

        ! ==== allocate g blasius functions for temperature
        allocate(g0_mb(1:nhuge), g1_mb(1:nhuge),stat=err)
        if(err .ne. 0) stop ' Allocation error in BlasiusProfile'
        allocate(g0_lb(1:nhuge), g1_lb(1:nhuge),stat=err)
        if(err .ne. 0) stop ' Allocation error in BlasiusProfile'
        allocate(g0_ub(1:nhuge), g1_ub(1:nhuge),stat=err)
        if(err .ne. 0) stop ' Allocation error in BlasiusProfile'

        g0_lb = 0.0_rp; g1_lb = 0.0_rp
        g0_ub = 0.0_rp; g1_ub = 0.0_rp
        g0_mb = 0.0_rp; g1_mb = 0.0_rp

        ! === create a uniform grid of nhuge points
        do j = 1, nhuge
           etaC(j) = eta_max*(j-1)/real(nhuge,rp)
        enddo
        !
        ! === compute blasius solution
        !
        f0_lb(1) = 0.0_rp
        f0_ub(1) = 0.0_rp
        f0_mb(1) = 0.0_rp

        f1_lb(1) = 0.0_rp
        f1_ub(1) = 0.0_rp
        f1_mb(1) = 0.0_rp

        f2_lb(1) = 0.1_rp
        f2_ub(1) = 1.0_rp

        erri = 2*toll
        iter = 1
        
        d_eta =  etaC(2) - etaC(1)
        do j = 2, nhuge

           f0_lb(j) = f0_lb(j-1) +     f1_lb(j-1)            * d_eta
           f1_lb(j) = f1_lb(j-1) +     f2_lb(j-1)            * d_eta
           f2_lb(j) = f2_lb(j-1) - alp*f0_lb(j-1)*f2_lb(j-1) * d_eta

           f0_ub(j) = f0_ub(j-1) +     f1_ub(j-1)            * d_eta
           f1_ub(j) = f1_ub(j-1) +     f2_ub(j-1)            * d_eta
           f2_ub(j) = f2_ub(j-1) - alp*f0_ub(j-1)*f2_ub(j-1) * d_eta
        enddo

        ! === compute blasius solution for momentum
        do while (erri > toll .and. iter < itermax)
        
           f2_mb(1) = 0.5_rp * (f2_lb(1) + f2_ub(1))

           do j = 2, nhuge

              f0_mb(j) = f0_mb(j-1) +     f1_mb(j-1)            * d_eta
              f1_mb(j) = f1_mb(j-1) +     f2_mb(j-1)            * d_eta
              f2_mb(j) = f2_mb(j-1) - alp*f0_mb(j-1)*f2_mb(j-1) * d_eta
           enddo
           
           if    ((f1_lb(nhuge) - 1.0_rp) * (f1_mb(nhuge) - 1.0_rp) > 0.0_rp) then 
           ! the solution is in between of mb ub
                 f2_lb(1)     = f2_mb(1)
                 f1_lb(nhuge) = f1_mb(nhuge)

           elseif((f1_lb(nhuge) - 1.0_rp) * (f1_mb(nhuge) - 1.0_rp) < 0.0_rp) then 
           ! the solution is in between of lb mb
                 f2_ub(1)    = f2_mb(1)

           endif

           erri = abs(f1_mb(nhuge) - 1.0_rp)
           iter = iter + 1

        enddo
       
        ! === compute blasius solution for temperature
        g0_lb(1) = Trat*0.9_rp*(1.0_rp + sqrt(Prandtl)*0.5_rp*(gamma0-1.0_rp)*Mach**2)
        g0_ub(1) = Trat*1.1_rp*(1.0_rp + sqrt(Prandtl)*0.5_rp*(gamma0-1.0_rp)*Mach**2)

        g1_lb(1) = 0.0_rp
        g1_ub(1) = 0.0_rp
        g1_mb(1) = 0.0_rp

        B = Prandtl*(gamma0-1.0_rp)*Mach**2

        erri = 2*toll
        iter = 1

        do j = 2, nhuge

           g0_lb(j) = g0_lb(j-1) + g1_lb(j-1)                                       *d_eta
           g0_ub(j) = g0_ub(j-1) + g1_ub(j-1)                                       *d_eta

           g1_lb(j) = g1_lb(j-1) - (Prandtl*f0_mb(j-1)*g1_lb(j-1) + B*f2_mb(j-1)**2)*d_eta
           g1_ub(j) = g1_ub(j-1) - (Prandtl*f0_mb(j-1)*g1_ub(j-1) + B*f2_mb(j-1)**2)*d_eta

        enddo

        ! === compute blasius solution for momentum
        do while (erri > toll .and. iter < itermax)
        
           g0_mb(1) = 0.5_rp * (g0_lb(1) + g0_ub(1))

           do j = 2, nhuge

              g0_mb(j) = g0_mb(j-1) + g1_mb(j-1)                                       *d_eta
              g1_mb(j) = g1_mb(j-1) - (Prandtl*f0_mb(j-1)*g1_mb(j-1) + B*f2_mb(j-1)**2)*d_eta

           enddo

           if    ((g0_lb(nhuge) - 1.0_rp) * (g0_mb(nhuge) - 1.0_rp) > 0.0_rp) then 
           ! the solution is in between of mb ub
                 g0_lb(1)     = g0_mb(1)
                 g0_lb(nhuge) = g0_mb(nhuge)

           elseif((g0_lb(nhuge) - 1.0_rp) * (g0_mb(nhuge) - 1.0_rp) < 0.0_rp) then 
           ! the solution is in between of lb mb
                 g0_ub(1)    = g0_mb(1)
           endif

           erri = abs(g0_mb(nhuge) - 1.0_rp)
           iter = iter + 1

        enddo


        !
        ! cycle setting delta99 = 1
        !
        allocate(y_tmp(1-GN:ny+GN))
        call compute_grid_point(y_tmp,ny,lbound(y_tmp,1),ubound(y_tmp,1),ymin,ymax,stretching_par,gridpoint_y)
       
        ntot = 0
        linearSearch0:do j = 1,ny
           if(y_tmp(j) > eta_max) then
             ntot = j
             exit linearSearch0
           endif
        enddo linearSearch0

        allocate(y_tmp1(1:ntot))
        call compute_grid_point(y_tmp1,ny,lbound(y_tmp1,1),ubound(y_tmp1,1),ymin,ymax,stretching_par,gridpoint_y)


        erri = 2.0_rp * toll
        dtheta = 1.0_rp
        iter = 0

        theta1 = 1.0_rp
        theta2 = 5.0_rp
        theta  = 0.5_rp*(theta1+theta2)

        do while(erri > toll .and. dtheta > 1.0E-15_rp .and. iter < itermax)


           ntheta = 0
           linearSearch1: do j = 1, ny
            if(y_tmp(j) > theta) then
              ntheta = j
              exit linearSearch1
            endif
           enddo linearSearch1
           
           etaC = etaC * y_tmp1(ntheta)/etaC(nhuge)
       
           n1 = 0
           linearSearch2: do j = 1, ntheta
             f1_int = interpl_field_1D(y_tmp(j),etaC,f1_mb)
             if(abs(1.0_rp-f1_int)<toll) then
               n1 = j
               exit linearSearch2
             endif
           enddo linearSearch2
        
           yo = y_tmp(n1)

           if (yo > 1.0_rp) then
              theta2 = theta
           elseif (yo < 1.0_rp) then
              theta1 = theta
           else
              iter = itermax
           end if

           theta = (theta1 + theta2) * 0.5_rp

           dtheta = abs(theta2 - theta1)

           erri = abs(yo - 1.0_rp)
           iter = iter + 1

        enddo
        
        !
        ! === assign blasius profile
        !
        do j = 1-GN,ny+GN

           ! output grid point point
           yo = y_tmp(j)

           if(yo <= eta_max) then
             f1_int = interpl_field_1D(yo,etaC,f1_mb)
             g0_int = interpl_field_1D(yo,etaC,g0_mb)

             BlasiusVel(j) = f1_int * u_inf
             BlasiusTmp(j) = g0_int 
             BlasiusRho(j) = 1.0_rp/g0_int

           else
             BlasiusVel(j) = u_inf
             BlasiusTmp(j) = 1.0_rp
             BlasiusRho(j) = 1.0_rp

           endif

        enddo

        deallocate(etaC, f0_lb, f1_lb, f2_lb, f0_ub, f1_ub, f2_ub, f0_mb, f1_mb, f2_mb)
        deallocate(      g0_lb, g1_lb,        g0_ub, g1_ub       , g0_mb, g1_mb)
        deallocate(y_tmp)

        return
end subroutine BlasiusProfile




subroutine PohlhausenProfile(PohlhausenVel)
! ---------------------------------------------------------------
!
!       This subroutine computes the Pohlhausen velocity profile
!       according to a third order Blasius approximation. 
!
! ---------------------------------------------------------------
        use mesh_module      , only: compute_grid_point
        use parameters_module, only: u_inf

        implicit none
        real(rp), dimension(1-GN:ny+GN) :: PohlhausenVel

        ! local declarations
        real(rp), parameter :: c1 = 1.5_rp, c2 = - 0.5_rp
        real(rp), dimension(:), allocatable :: y_tmp
        real(rp) :: yj, etaj
        integer  :: j

        allocate(y_tmp(1-GN:ny+GN))
        call compute_grid_point(y_tmp,ny,lbound(y_tmp,1),ubound(y_tmp,1),ymin,ymax,stretching_par,gridpoint_y)

        do j = 1-GN, ny+GN

           yj = y_tmp(j)
           etaj = min(yj,1.0_rp)

           PohlhausenVel(j) = u_inf*(c1*etaj +c2*etaj**3)

        enddo

        deallocate(y_tmp)

        return
end subroutine PohlhausenProfile


subroutine MuskerProfile(ny,y_tmp,ReTau,MuskerUPlus)
! ------------------------------------------------------------------------------------------------
!
!       This function provides an explicit formulation the velocity near wall profile folling:
!
!       Musker, "Explicit Expression fortehSmooth Wall Velocity Distribution in a
!                Turbulent Boundary Layer", AIAA 1979
!
!       INPUT : ReTau    !< nominal shear Reynolds
!       OUTPUT: uPlus    !> uPlus velocity
!
! ------------------------------------------------------------------------------------------------
        implicit none
        integer                            , intent(in)  :: ny
        real(rp), dimension(:), allocatable, intent(in)  :: y_tmp
        real(rp)                           , intent(in)  :: ReTau
        real(rp), dimension(0:ny)          , intent(out) :: MuskerUPlus

        ! local declarations
        real(rp), parameter :: pi_wake = 0.434_rp
        real(rp)            :: yj, eta, yp, up, iReTau
        integer             :: j
        

        iRetau = 1.0_rp/ReTau
        MuskerUPlus(0) = 0.0_rp
        do j = 1,ny
        
           ! compute y+ coordinate
           yj  = abs(y_tmp(j) - ymin)
           yp  = yj * ReTau

           eta = yp*iReTau
           eta = min(1.0_rp,eta)
           yp  = eta*retau

           up = 5.424_rp*atan((2*yp-8.15_rp)/16.7_rp)+log10((yp+10.6_rp)**9.6_rp/(yp**2-8.15_rp*yp+86)**2)-3.51132976630723_rp+&
                2.44_rp*(pi_wake*(6*eta**2-4*eta**3)+(eta**2*(1-eta)))

           MuskerUPlus(j) = up
        
        enddo

        return
end subroutine MuskerProfile


















subroutine CompressibleCorrection(ny,turbFlag,Mach,Prandtl,Trat,Tw,Rw,u_inc,u_cmp,RhoY,TmpY)

        implicit none
        integer                    , intent(in)  :: ny
        logical                    , intent(in)  :: turbFlag
        real(rp)                   , intent(in)  :: Mach
        real(rp)                   , intent(in)  :: Prandtl
        real(rp)                   , intent(in)  :: Trat
        real(rp), dimension(0:ny), intent(in)  :: u_inc
        real(rp), dimension(0:ny), intent(out) :: u_cmp
        real(rp), dimension(0:ny), intent(out) :: RhoY
        real(rp), dimension(0:ny), intent(out) :: TmpY
        real(rp)                 , intent(out) :: Tw              !< Wall Temperature
        real(rp)                 , intent(out) :: Rw              !< Wall Density

        ! local declarations
        real(rp), dimension(0:ny) :: u_cmp_old
        real(rp)                  :: Tr                !< recovery temperature
        real(rp)                  :: uu
        real(rp)                  :: du, uci, err
        integer                   :: itr

        real(rp)                  :: RecFac            !< recovery factor
        real(rp)                  :: alf               !< See Zhang
        real(rp)                  :: fuu
        real(rp), parameter       :: s      = 1.1_rp      !< Reynolds analogy factor
        integer , parameter       :: itmx   = 100
        real(rp), parameter       :: toll   = 1.0E-14_rp
        real(rp), parameter       :: gamma0 = 1.4_rp      

        if(turbFlag) then
          RecFac = Prandtl**(1.0_rp/3.0_rp)
        else
          RecFac = sqrt(Prandtl)
        endif
       
        Tr  = 1.0_rp + 0.5_rp*RecFac*(gamma0-1.0_rp)*Mach**2
        Tw  = TRat*Tr
        Rw  = 1.0_rp/Tw
        alf = s*Prandtl

        ! init the loop
        itr   = 0
        err   = huge(1.0_rp)
        u_cmp = u_inc 

        do while(err > toll .and. itr < itmx)

           u_cmp_old = u_cmp
                
           ! compute Temperature and Density profiles
           do j = 0,ny
              uu      = u_cmp(j)/u_cmp(ny)
              fuu     = alf*uu+(1.0_rp-alf)*uu**2
              TmpY(j) = tw+(tr-tw)*fuu+(1.0_rp-tr)*uu**2           !< Zhang
              !TmpY(j) = Tw + (Tr-Tw)*uu + (1.0_rp-Tr)*uu**2       !< Crocco-Busemann
              RhoY(j) = 1.0_rp/TmpY(j)
           enddo
        
           ! apply the correction to uPlus_i according to VAN DRIEST
           do j = 1,ny
              du       = u_inc(j)-u_inc(j-1)
              uci      = 0.5_rp*(sqrt(Rw/RhoY(j))+sqrt(Rw/RhoY(j-1))) 
              u_cmp(j) = u_cmp(j-1)+uci*du
           enddo
           !do j = ny+1,ny
           !   TmpY(j) = TmpY(ny)
           !   RhoY(j) = RhoY(ny)
           !   u_cmp(j) = u_cmp(ny)
           !enddo
        
           err = sum(abs(u_cmp - u_cmp_old))
           itr = itr + 1

        enddo

        return
end subroutine CompressibleCorrection





subroutine BoundaryLayerQuantities(y_tmp,RhoY,UPlusCmp,rWall,theta,deltaStar,cf)
        implicit none
        real(rp), dimension(:), allocatable, intent(in)  :: y_tmp
        real(rp), dimension(0:ny)          , intent(in)  :: UPlusCmp
        real(rp), dimension(0:ny)          , intent(in)  :: RhoY
        real(rp)                           , intent(in)  :: rWall
        real(rp)                           , intent(out) :: theta
        real(rp)                           , intent(out) :: deltaStar
        real(rp)                           , intent(out) :: cf

        ! local declarations
        real(rp) :: i_ue, u00, um1, r00, rm1, dy
        integer  :: j

        i_ue = 1.0_rp/uPlusCmp(ny)
        
        ! first point
        deltaStar = 0.0_rp
        theta     = 0.0_rp
        do j = 1,ny
        
           dy  = y_tmp(j) - y_tmp(j-1)

           u00 = uPlusCmp(j)   * i_ue
           um1 = uPlusCmp(j-1) * i_ue

           r00 = RhoY(j)
           rm1 = RhoY(j-1)

           deltaStar = deltaStar + 0.5_rp*((1.0_rp-r00*u00)     + (1.0_rp-rm1*um1)    )*dy
           theta     = theta     + 0.5_rp*(r00*u00*(1.0_rp-u00) + rm1*um1*(1.0_rp-um1))*dy

        enddo

        cf = 2.0_rp*rWall*i_ue**2

        return
end subroutine BoundaryLayerQuantities


function logLaw(ut,y,nW) result(u)

        implicit none
        real(rp), intent(in) :: ut
        real(rp), intent(in) :: y
        real(rp), intent(in) :: nW
        real(rp)             :: u
        
        real(rp), parameter :: ivk = 1.0_rp/0.41_rp
        real(rp), parameter :: B   = 5.2_rp
        real(rp)            :: yPlus

        yPlus = y*ut/nW

        u = ut*(ivk*log(yPlus) + B)

        return
end function logLaw




function ilogLaw(ut,y,u,nW) result(fut)

        implicit none
        real(rp), intent(in) :: ut
        real(rp), intent(in) :: y
        real(rp), intent(in) :: u
        real(rp), intent(in) :: nW
        real(rp)             :: fut
        
        real(rp), parameter :: ivk = 1.0_rp/0.41_rp
        real(rp), parameter :: B   = 5.2_rp
        real(rp)            :: yPlus

        yPlus = y*ut/nW

        fut = u/ut - ivk*log(yPlus) - B

        return
end function ilogLaw



function ReichardtLaw(yPlus) result(uPlus)

        implicit none
        real(rp), intent(in) :: yPlus
        real(rp)             :: uPlus
        
        ! local declarations
        real(rp), parameter :: one11 = 1.0_rp/11.0_rp
        real(rp)            :: ePlus

        ePlus = exp(-yPlus*one11)
        
        uPlus = 2.5_rp*log(1.0_rp+0.4_rp*yPlus) &
              + 7.8_rp*(1.0_rp - ePlus - yPlus*one11*exp(-0.33_rp*yPlus))

        return
end function ReichardtLaw


function iReichardtLaw(ut,y,u,nW) result(fut)
        !$acc routine seq

        implicit none
        real(rp), intent(in) :: ut
        real(rp), intent(in) :: y
        real(rp), intent(in) :: u
        real(rp), intent(in) :: nW
        real(rp)             :: fut
        
        ! local declarations
        real(rp), parameter :: Cvk = 7.8_rp
        real(rp), parameter :: i11 = 1.0_rp/11.0_rp
        real(rp)            :: ePlus, yPlus

        yPlus = y*ut/nW
        ePlus = exp(-yPlus*i11)
        
        fut = u/(ut+1.0E-10_rp)           &
              - 2.5_rp*log(1.0_rp+0.4_rp*yPlus) &
              - Cvk*(1.0_rp - ePlus - yPlus*i11*exp(-0.33_rp*yPlus))

        return
end function iReichardtLaw




function hyperbolic_step(yj,y1,y2) result(fy)
        implicit none
        real(rp), intent(in) :: yj, y1, y2
        real(rp)             :: fy

        real(rp) :: htan1, htan2
        real(rp), parameter :: i_a = 1.0_rp/0.03_rp

        htan1 = tanh((yj-y1)*i_a)
        htan2 = tanh((yj-y2)*i_a)

        fy = htan1 - htan2 + 1.0_rp
        
        return
endfunction hyperbolic_step
























end module fluid_functions_module
