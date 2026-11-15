module fluid_functions_module
use parameters_module
use storage_module
use thermodynamics_module, only: suth_S
implicit none

real(rp), parameter :: zhang_s_rey = 1.1_rp

contains

elemental function laminar_viscosity(t,tref,vis_flag) result(mu_t)
        !$acc routine seq
        use parameters_module, only: rp

        implicit none
        real(rp), intent(in) :: T
        real(rp), intent(in) :: Tref
        integer , intent(in) :: vis_flag

        real(rp)             :: a,b,mu_t

        if    (vis_flag == 0) then

          a     = sqrt(T)
          mu_t  = a*sqrt(a)

        elseif(vis_flag == 1) then

          a = suth_S/Tref
          b = 1.0_rp + a

          mu_t = T*sqrt(T) * b/(T+a)

        elseif(vis_flag == 2) then

          mu_t = 1.0_rp

        else

          mu_t = 0.0_rp

        endif

        return
end function laminar_viscosity

elemental function KarmanSchoenherr(ReThI) result(cfi)
        implicit none
        real(rp), intent(in) :: ReThI
        real(rp)             :: cfi

        cfi = 1.0_rp/(17.08_rp*(log10(ReThI))**2 + 25.11_rp*log10(ReThI) + 6.012_rp)

        return
end function KarmanSchoenherr

function ReThetaIncSpalding(uep) result(ReThInc)
        implicit none
        real(rp), intent(in) :: uep
        real(rp)             :: ReThInc

        real(rp), parameter :: k = 0.4_rp, E = 12.0_rp
        real(rp)            :: a1, a2, a3, a4, a5, a6, a7, a8

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
        implicit none
        integer , intent(in) :: i,j,k
        real(rp)             :: div

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

function ThetaBetaMach(beta,theta,mach) result(fbeta)
        implicit none
        real(rp), intent(in) :: theta
        real(rp), intent(in) :: beta
        real(rp), intent(in) :: Mach
        real(rp)             :: fbeta

        real(rp), parameter :: gamma = 1.4_rp
        real(rp)            :: cot_beta, nnum, dnum

        cot_beta = 1.0_rp/tan(beta)

        nnum = Mach**2 * sin(beta)**2 - 1.0_rp
        dnum = Mach**2 * (gamma + cos(2*beta)) + 2

        fbeta = tan(theta) - 2*cot_beta * nnum/dnum

        return
end function ThetaBetaMach

function direct_PrandtlMeyer(M) result(vM)
        implicit none
        real(rp), intent(in) :: M
        real(rp)             :: vM

        real(rp), parameter  :: g = 1.4_rp
        real(rp)             :: par, i_par, x

        par   = sqrt((g + 1.0_rp)/(g - 1.0_rp))
        i_par = 1.0_rp/par

        x     = sqrt(M**2 - 1.0_rp)
        vM    = par * atan(i_par * x) - atan(x)

        return
end function direct_PrandtlMeyer

function inverse_PrandtlMeyer(M2,M1,theta) result(fM2)
        implicit none
        real(rp), intent(in) :: M2
        real(rp), intent(in) :: M1
        real(rp), intent(in) :: theta
        real(rp)             :: fM2

        real(rp)             :: vM1, vM2

        vM1 = direct_PrandtlMeyer(M1)

        vM2 = direct_PrandtlMeyer(M2)

        fM2 = vM2 - theta - vM1

        return
end function inverse_PrandtlMeyer

subroutine hTurbEpsilon(xstep,ystep,zstep,phi,U,V,W,VIS,gbl_rh_mean,&
                gbl_eps_sol,gbl_eps_dil,gbl_eps_tot,gbl_eps_str,mpi_flag)

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

                 st(1) = xstep(i)
                 st(2) = ystep(j)
                 st(3) = zstep(k)

                 i_st(:) = 1.0_rp/st(:)

                 dVol = st(1)*st(2)*st(3)
                 lcl_Vtot_sp = lcl_Vtot_sp + dVol

                 mu = mu_inf * VIS(i,j,k)
                 nu = mu/phi(i,j,k,1)

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

                 grad_v (1,:) = du(:)
                 grad_v (2,:) = dv(:)
                 grad_v (3,:) = dw(:)

                 grad_vT(:,1) = du(:)
                 grad_vT(:,2) = dv(:)
                 grad_vT(:,3) = dw(:)

                 S = 0.5_rp*(grad_v + grad_vT)

                 lcl_S2_mean(:,:) = lcl_S2_mean(:,:) + nu * S(:,:)*S(:,:) * dVol

                 vor(1) = dw(2) - dv(3)
                 vor(2) = du(3) - dw(1)
                 vor(3) = dv(1) - du(2)

                 vor2 = vor(1)*vor(1) + vor(2)*vor(2) + vor(3)*vor(3)

                 div2 = (du(1) + dv(2) + dw(3))**2

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

subroutine MuskerProfile(ny,y_tmp,ReTau,MuskerUPlus)
        implicit none
        integer                            , intent(in)  :: ny
        real(rp), dimension(:), allocatable, intent(in)  :: y_tmp
        real(rp)                           , intent(in)  :: ReTau
        real(rp), dimension(0:ny)          , intent(out) :: MuskerUPlus

        real(rp), parameter :: pi_wake = 0.434_rp
        real(rp)            :: yj, eta, yp, up, iReTau
        integer             :: j

        iRetau = 1.0_rp/ReTau
        MuskerUPlus(0) = 0.0_rp
        do j = 1,ny

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

subroutine CompressibleCorrection(ny,u_inc,u_cmp,RhoY,TmpY,Tw,Rw)
        use parameters_module    , only: Mach, Prandtl, Trat, T_recovery, &
                                         thermally_perfect_gas
        use thermodynamics_module, only: enthalpy_t, heat_capacity_t

        implicit none
        integer                  , intent(in)  :: ny
        real(rp), dimension(0:ny), intent(in)  :: u_inc
        real(rp), dimension(0:ny), intent(out) :: u_cmp
        real(rp), dimension(0:ny), intent(out) :: RhoY
        real(rp), dimension(0:ny), intent(out) :: TmpY
        real(rp)                 , intent(out) :: Tw
        real(rp)                 , intent(out) :: Rw

        real(rp), parameter :: tollVD = 1.0E-14_rp
        real(rp), parameter :: tollTP = 1.0E-12_rp
        integer , parameter :: itmxVD = 100
        integer , parameter :: itmxTP = 50

        real(rp), dimension(0:ny) :: u_cmp_old
        real(rp) :: Tr, hw, hr, h_inf, h_j, cp_j
        real(rp) :: alpha, uu, fuu, du, uci, err, Tnode, dT
        real(rp) :: inv_ue
        integer  :: j, itr, itrNR

        Tr    = T_recovery
        Tw    = Trat*Tr
        Rw    = 1.0_rp/Tw
        alpha = zhang_s_rey*Prandtl

        if(thermally_perfect_gas > 0) then
           h_inf = enthalpy_t(1.0_rp)
           hr    = enthalpy_t(Tr)
           hw    = enthalpy_t(Tw)
        endif

        u_cmp = u_inc
        err   = huge(1.0_rp)
        itr   = 0

        picard: do while(err > tollVD .and. itr < itmxVD)

           u_cmp_old = u_cmp

           inv_ue = 1.0_rp/u_cmp(ny)
           do j = 0, ny
              uu  = u_cmp(j)*inv_ue
              fuu = alpha*uu + (1.0_rp - alpha)*uu*uu

              if(thermally_perfect_gas > 0) then
                 h_j   = hw + (hr - hw)*fuu + (h_inf - hr)*uu*uu
                 Tnode = Tw + (Tr - Tw)*fuu + (1.0_rp - Tr)*uu*uu
                 do itrNR = 1, itmxTP
                    cp_j  = heat_capacity_t(Tnode)
                    dT    = (h_j - enthalpy_t(Tnode))/cp_j
                    Tnode = Tnode + dT
                    if(abs(dT) < tollTP) exit
                 enddo
                 TmpY(j) = Tnode
              else
                 TmpY(j) = Tw + (Tr - Tw)*fuu + (1.0_rp - Tr)*uu*uu
              endif

              RhoY(j) = 1.0_rp/TmpY(j)
           enddo

           do j = 1, ny
              du       = u_inc(j) - u_inc(j-1)
              uci      = 0.5_rp*(sqrt(Rw/RhoY(j)) + sqrt(Rw/RhoY(j-1)))
              u_cmp(j) = u_cmp(j-1) + uci*du
           enddo

           err = sum(abs(u_cmp - u_cmp_old))
           itr = itr + 1

        enddo picard

        return
end subroutine CompressibleCorrection

function GRATemperatureSlope(uu, Tm) result(dTduu)
        use parameters_module    , only: Prandtl, Trat, T_recovery, thermally_perfect_gas
        use thermodynamics_module, only: enthalpy_t, heat_capacity_t

        implicit none
        real(rp), intent(in) :: uu, Tm
        real(rp)             :: dTduu, Tr, Tw, alpha, dfuu

        Tr    = T_recovery
        Tw    = Trat*Tr
        alpha = zhang_s_rey*Prandtl
        dfuu  = alpha + 2.0_rp*(1.0_rp - alpha)*uu

        if(thermally_perfect_gas > 0) then
           dTduu = ((enthalpy_t(Tr) - enthalpy_t(Tw))*dfuu                 &
                 + 2.0_rp*(enthalpy_t(1.0_rp) - enthalpy_t(Tr))*uu)/heat_capacity_t(Tm)
        else
           dTduu = (Tr - Tw)*dfuu + 2.0_rp*(1.0_rp - Tr)*uu
        endif

        return
end function GRATemperatureSlope

subroutine BoundaryLayerQuantities(y_tmp,RhoY,UPlusCmp,rWall,theta,deltaStar,cf)
        implicit none
        real(rp), dimension(:), allocatable, intent(in)  :: y_tmp
        real(rp), dimension(0:ny)          , intent(in)  :: UPlusCmp
        real(rp), dimension(0:ny)          , intent(in)  :: RhoY
        real(rp)                           , intent(in)  :: rWall
        real(rp)                           , intent(out) :: theta
        real(rp)                           , intent(out) :: deltaStar
        real(rp)                           , intent(out) :: cf

        real(rp) :: i_ue, u00, um1, r00, rm1, dy
        integer  :: j

        i_ue = 1.0_rp/uPlusCmp(ny)

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

function ReichardtLaw(yPlus) result(uPlus)

        implicit none
        real(rp), intent(in) :: yPlus
        real(rp)             :: uPlus

        real(rp), parameter :: one11 = 1.0_rp/11.0_rp
        real(rp)            :: ePlus

        ePlus = exp(-yPlus*one11)

        uPlus = 2.5_rp*log(1.0_rp+0.4_rp*yPlus) &
              + 7.8_rp*(1.0_rp - ePlus - yPlus*one11*exp(-0.33_rp*yPlus))

        return
end function ReichardtLaw

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
