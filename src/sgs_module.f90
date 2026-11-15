module sgs_module

        use parameters_module, only: rp, pi, sgs_model, dims, &
                                     central_1, mu_inf
        use mpi_module
        use mpi_comm_module
        use profiling_module , only: StartProfRange, EndProfRange
        use storage_module   , only: phi, U, V, W, GAMMA_T, VIS, LMD
        use mesh_module      , only: zstep_i, ztastep_i, iz_zta, Jac2D, &
                                     dxi_dx, dxi_dy, deta_dx, deta_dy
        use sgs_bc_module    , only: TrbVisBC
        use wale_module      , only: compute_wale

        implicit none
        private

        real(rp), parameter :: Pr_turb   = 0.9_rp
        real(rp), parameter :: i_Pr_turb = 1.0_rp/Pr_turb
        real(rp), parameter :: one3      = 1.0_rp/3.0_rp
        real(rp), parameter :: two3      = 2.0_rp/3.0_rp
        real(rp), parameter :: one6      = 1.0_rp/6.0_rp

        public :: compute_subgrid_model

contains

subroutine compute_subgrid_model
        implicit none

        call StartProfRange("compute_subgrid_model")
        call StartProfRange("sgs_model")

        if (dims /= 3) stop ' LES requires 3D!'

        selectcase(trim(sgs_model))

          case('Smagorinsky')
            call compute_smagorinsky(VIS,LMD)

          case('WALE')
            call compute_wale(VIS,LMD)

          case('sigma_model')
            call compute_sigma_model(VIS,LMD)

          case('MixedTimeScale')
            call compute_MixedTimeScale(VIS,LMD)

          case default
            write(*,'(a)')   ' FATAL ERROR!'
            write(*,'(a,a)') trim(sgs_model), ' is not implemented. The program will be stopped'
            call secure_stop

        endselect
        call EndProfRange

        call StartProfRange("sgs_mpi_calls")
        if (mpi_flag) then
           call mpi_share(mpi_comm_cart,type_send_prim,type_recv_prim,my_neighbour,dims, &
                   bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                   bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                   bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                   VIS)

           call mpi_share(mpi_comm_cart,type_send_prim,type_recv_prim,my_neighbour,dims, &
                   bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                   bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                   bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                   LMD)
        endif
        call EndProfRange

        call TrbVisBC(VIS,LMD)

        call EndProfRange

        return
end subroutine compute_subgrid_model

subroutine compute_smagorinsky(VIS,LMD)

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: VIS, LMD

        real(rp), parameter :: Cs2 = 0.12_rp**2

        real(rp) :: csi_x, csi_y, eta_x, eta_y, zeta_z
        real(rp) :: duc, due, duz
        real(rp) :: dvc, dve, dvz
        real(rp) :: dwc, dwe, dwz
        real(rp) :: du_x, du_y, du_z
        real(rp) :: dv_x, dv_y, dv_z
        real(rp) :: dw_x, dw_y, dw_z
        real(rp) :: S_xy, S_xz, S_yz
        real(rp) :: sumSS, delta2, sqrtSS
        real(rp) :: mu_turb, cp_loc, i_mu_inf
        integer  :: i, j, k, l

        i_mu_inf = 1.0_rp/mu_inf

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex

                 csi_x = dxi_dx (i,j); csi_y = dxi_dy (i,j)
                 eta_x = deta_dx(i,j); eta_y = deta_dy(i,j)
                 zeta_z = ztastep_i(k)*iz_zta(k)

                 duc = 0.0_rp; due = 0.0_rp; duz = 0.0_rp
                 dvc = 0.0_rp; dve = 0.0_rp; dvz = 0.0_rp
                 dwc = 0.0_rp; dwe = 0.0_rp; dwz = 0.0_rp

                 !$acc loop seq
                 do l = 1,3
                    duc = duc + central_1(l) * (U(i+l,j,k) - U(i-l,j,k))
                    due = due + central_1(l) * (U(i,j+l,k) - U(i,j-l,k))
                    duz = duz + central_1(l) * (U(i,j,k+l) - U(i,j,k-l))

                    dvc = dvc + central_1(l) * (V(i+l,j,k) - V(i-l,j,k))
                    dve = dve + central_1(l) * (V(i,j+l,k) - V(i,j-l,k))
                    dvz = dvz + central_1(l) * (V(i,j,k+l) - V(i,j,k-l))

                    dwc = dwc + central_1(l) * (W(i+l,j,k) - W(i-l,j,k))
                    dwe = dwe + central_1(l) * (W(i,j+l,k) - W(i,j-l,k))
                    dwz = dwz + central_1(l) * (W(i,j,k+l) - W(i,j,k-l))
                 enddo

                 du_x = csi_x*duc + eta_x*due; du_y = csi_y*duc + eta_y*due; du_z = zeta_z*duz
                 dv_x = csi_x*dvc + eta_x*dve; dv_y = csi_y*dvc + eta_y*dve; dv_z = zeta_z*dvz
                 dw_x = csi_x*dwc + eta_x*dwe; dw_y = csi_y*dwc + eta_y*dwe; dw_z = zeta_z*dwz

                 S_xy = 0.5_rp*(du_y + dv_x)
                 S_xz = 0.5_rp*(du_z + dw_x)
                 S_yz = 0.5_rp*(dv_z + dw_y)

                 sumSS = du_x*du_x + dv_y*dv_y + dw_z*dw_z + &
                         2.0_rp*(S_xy*S_xy + S_xz*S_xz + S_yz*S_yz)
                 sqrtSS = sqrt(2.0_rp*sumSS)

                 delta2 = 1.0_rp/(Jac2D(i,j)*zstep_i(k))**two3

                 mu_turb = i_mu_inf * phi(i,j,k,1) * Cs2 * delta2 * sqrtSS
                 cp_loc  = GAMMA_T(i,j,k)/(GAMMA_T(i,j,k) - 1.0_rp)

                 VIS(i,j,k) = VIS(i,j,k) + mu_turb
                 LMD(i,j,k) = LMD(i,j,k) + cp_loc * i_Pr_turb * mu_turb

              enddo
           enddo
        enddo
        !$acc end parallel

        return
end subroutine compute_smagorinsky

subroutine compute_sigma_model(VIS,LMD)

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: VIS, LMD

        real(rp), parameter :: Cs2        = 1.0_rp**2
        real(rp), parameter :: toll       = 1.0E-14_rp
        real(rp), parameter :: twoPithird = 2.0_rp*pi*one3

        real(rp) :: csi_x, csi_y, eta_x, eta_y, zeta_z
        real(rp) :: duc, due, duz
        real(rp) :: dvc, dve, dvz
        real(rp) :: dwc, dwe, dwz
        real(rp) :: du_x, du_y, du_z
        real(rp) :: dv_x, dv_y, dv_z
        real(rp) :: dw_x, dw_y, dw_z
        real(rp) :: G_xx, G_yy, G_zz, G_xy, G_xz, G_yz
        real(rp) :: B_xx, B_yy, B_zz, B_xy, B_xz, B_yz
        real(rp) :: trac, q, p, p1, p2, ip, r, phi_
        real(rp) :: lambda_x, lambda_y, lambda_z
        real(rp) :: s_x, s_y, s_z
        real(rp) :: delta2, mu_turb, cp_loc, i_mu_inf
        integer  :: i, j, k, l

        i_mu_inf = 1.0_rp/mu_inf

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex

                 csi_x = dxi_dx (i,j); csi_y = dxi_dy (i,j)
                 eta_x = deta_dx(i,j); eta_y = deta_dy(i,j)
                 zeta_z = ztastep_i(k)*iz_zta(k)

                 duc = 0.0_rp; due = 0.0_rp; duz = 0.0_rp
                 dvc = 0.0_rp; dve = 0.0_rp; dvz = 0.0_rp
                 dwc = 0.0_rp; dwe = 0.0_rp; dwz = 0.0_rp

                 !$acc loop seq
                 do l = 1,3
                    duc = duc + central_1(l) * (U(i+l,j,k) - U(i-l,j,k))
                    due = due + central_1(l) * (U(i,j+l,k) - U(i,j-l,k))
                    duz = duz + central_1(l) * (U(i,j,k+l) - U(i,j,k-l))

                    dvc = dvc + central_1(l) * (V(i+l,j,k) - V(i-l,j,k))
                    dve = dve + central_1(l) * (V(i,j+l,k) - V(i,j-l,k))
                    dvz = dvz + central_1(l) * (V(i,j,k+l) - V(i,j,k-l))

                    dwc = dwc + central_1(l) * (W(i+l,j,k) - W(i-l,j,k))
                    dwe = dwe + central_1(l) * (W(i,j+l,k) - W(i,j-l,k))
                    dwz = dwz + central_1(l) * (W(i,j,k+l) - W(i,j,k-l))
                 enddo

                 du_x = csi_x*duc + eta_x*due; du_y = csi_y*duc + eta_y*due; du_z = zeta_z*duz
                 dv_x = csi_x*dvc + eta_x*dve; dv_y = csi_y*dvc + eta_y*dve; dv_z = zeta_z*dvz
                 dw_x = csi_x*dwc + eta_x*dwe; dw_y = csi_y*dwc + eta_y*dwe; dw_z = zeta_z*dwz

                 G_xx = du_x*du_x + dv_x*dv_x + dw_x*dw_x
                 G_yy = du_y*du_y + dv_y*dv_y + dw_y*dw_y
                 G_zz = du_z*du_z + dv_z*dv_z + dw_z*dw_z
                 G_xy = du_x*du_y + dv_x*dv_y + dw_x*dw_y
                 G_xz = du_x*du_z + dv_x*dv_z + dw_x*dw_z
                 G_yz = du_y*du_z + dv_y*dv_z + dw_y*dw_z

                 trac = G_xx + G_yy + G_zz
                 p1   = G_xy*G_xy + G_xz*G_xz + G_yz*G_yz

                 if (p1 < toll) then

                    lambda_x = max(G_xx, G_yy, G_zz)
                    lambda_z = min(G_xx, G_yy, G_zz)
                    lambda_y = trac - (lambda_x + lambda_z)

                 else

                    q  = trac*one3
                    p2 = G_xx*G_xx + G_yy*G_yy + G_zz*G_zz - 3.0_rp*q*q + 2.0_rp*p1
                    p  = sqrt(p2*one6)
                    ip = 1.0_rp/p

                    B_xx = ip*(G_xx - q)
                    B_yy = ip*(G_yy - q)
                    B_zz = ip*(G_zz - q)
                    B_xy = ip*G_xy
                    B_xz = ip*G_xz
                    B_yz = ip*G_yz

                    r = 0.5_rp*( B_xx*(B_yy*B_zz - B_yz*B_yz)   &
                               - B_xy*(B_xy*B_zz - B_yz*B_xz)   &
                               + B_xz*(B_xy*B_yz - B_yy*B_xz) )

                    phi_ = acos(max(-1.0_rp, min(1.0_rp, r)))*one3

                    lambda_x = q + 2.0_rp*p*cos(phi_)
                    lambda_z = q + 2.0_rp*p*cos(phi_ + twoPithird)
                    lambda_y = trac - (lambda_x + lambda_z)

                 endif

                 s_x = sqrt(abs(lambda_x))
                 s_y = sqrt(abs(lambda_y))
                 s_z = sqrt(abs(lambda_z))

                 delta2 = 1.0_rp/(Jac2D(i,j)*zstep_i(k))**two3

                 mu_turb = i_mu_inf * phi(i,j,k,1) * Cs2 * delta2 &
                                    * s_z * (s_x - s_y) * (s_y - s_z) / (s_x*s_x)
                 cp_loc  = GAMMA_T(i,j,k)/(GAMMA_T(i,j,k) - 1.0_rp)

                 VIS(i,j,k) = VIS(i,j,k) + mu_turb
                 LMD(i,j,k) = LMD(i,j,k) + cp_loc * i_Pr_turb * mu_turb

              enddo
           enddo
        enddo
        !$acc end parallel

        return
end subroutine compute_sigma_model

subroutine compute_MixedTimeScale(VIS,LMD)

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: VIS, LMD

        real(rp), parameter :: Cmts = 0.03_rp
        real(rp), parameter :: C_t  = 10.0_rp
        real(rp), parameter :: i_Ct = 1.0_rp/C_t

        real(rp) :: csi_x, csi_y, eta_x, eta_y, zeta_z
        real(rp) :: duc, due, duz
        real(rp) :: dvc, dve, dvz
        real(rp) :: dwc, dwe, dwz
        real(rp) :: du_x, du_y, du_z
        real(rp) :: dv_x, dv_y, dv_z
        real(rp) :: dw_x, dw_y, dw_z
        real(rp) :: S_xy, S_xz, S_yz
        real(rp) :: sum_r, sum_ru, sum_rv, sum_rw
        real(rp) :: imeanr, mean_u, mean_v, mean_w
        real(rp) :: sumSS, asS, kes, T_S, delta2
        real(rp) :: mu_turb, cp_loc, i_mu_inf
        integer  :: i, j, k, l

        i_mu_inf = 1.0_rp/mu_inf

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex

                 csi_x = dxi_dx (i,j); csi_y = dxi_dy (i,j)
                 eta_x = deta_dx(i,j); eta_y = deta_dy(i,j)
                 zeta_z = ztastep_i(k)*iz_zta(k)

                 duc = 0.0_rp; due = 0.0_rp; duz = 0.0_rp
                 dvc = 0.0_rp; dve = 0.0_rp; dvz = 0.0_rp
                 dwc = 0.0_rp; dwe = 0.0_rp; dwz = 0.0_rp

                 !$acc loop seq
                 do l = 1,3
                    duc = duc + central_1(l) * (U(i+l,j,k) - U(i-l,j,k))
                    due = due + central_1(l) * (U(i,j+l,k) - U(i,j-l,k))
                    duz = duz + central_1(l) * (U(i,j,k+l) - U(i,j,k-l))

                    dvc = dvc + central_1(l) * (V(i+l,j,k) - V(i-l,j,k))
                    dve = dve + central_1(l) * (V(i,j+l,k) - V(i,j-l,k))
                    dvz = dvz + central_1(l) * (V(i,j,k+l) - V(i,j,k-l))

                    dwc = dwc + central_1(l) * (W(i+l,j,k) - W(i-l,j,k))
                    dwe = dwe + central_1(l) * (W(i,j+l,k) - W(i,j-l,k))
                    dwz = dwz + central_1(l) * (W(i,j,k+l) - W(i,j,k-l))
                 enddo

                 du_x = csi_x*duc + eta_x*due; du_y = csi_y*duc + eta_y*due; du_z = zeta_z*duz
                 dv_x = csi_x*dvc + eta_x*dve; dv_y = csi_y*dvc + eta_y*dve; dv_z = zeta_z*dvz
                 dw_x = csi_x*dwc + eta_x*dwe; dw_y = csi_y*dwc + eta_y*dwe; dw_z = zeta_z*dwz

                 S_xy = 0.5_rp*(du_y + dv_x)
                 S_xz = 0.5_rp*(du_z + dw_x)
                 S_yz = 0.5_rp*(dv_z + dw_y)

                 sumSS = du_x*du_x + dv_y*dv_y + dw_z*dw_z + &
                         2.0_rp*(S_xy*S_xy + S_xz*S_xz + S_yz*S_yz)

                 sum_r  = phi(i  ,j  ,k  ,1) + &
                          phi(i+1,j  ,k  ,1) + phi(i-1,j  ,k  ,1) + &
                          phi(i  ,j+1,k  ,1) + phi(i  ,j-1,k  ,1) + &
                          phi(i  ,j  ,k+1,1) + phi(i  ,j  ,k-1,1)
                 sum_ru = phi(i  ,j  ,k  ,2) + &
                          phi(i+1,j  ,k  ,2) + phi(i-1,j  ,k  ,2) + &
                          phi(i  ,j+1,k  ,2) + phi(i  ,j-1,k  ,2) + &
                          phi(i  ,j  ,k+1,2) + phi(i  ,j  ,k-1,2)
                 sum_rv = phi(i  ,j  ,k  ,3) + &
                          phi(i+1,j  ,k  ,3) + phi(i-1,j  ,k  ,3) + &
                          phi(i  ,j+1,k  ,3) + phi(i  ,j-1,k  ,3) + &
                          phi(i  ,j  ,k+1,3) + phi(i  ,j  ,k-1,3)
                 sum_rw = phi(i  ,j  ,k  ,4) + &
                          phi(i+1,j  ,k  ,4) + phi(i-1,j  ,k  ,4) + &
                          phi(i  ,j+1,k  ,4) + phi(i  ,j-1,k  ,4) + &
                          phi(i  ,j  ,k+1,4) + phi(i  ,j  ,k-1,4)

                 imeanr = 1.0_rp/sum_r
                 mean_u = sum_ru*imeanr
                 mean_v = sum_rv*imeanr
                 mean_w = sum_rw*imeanr

                 kes = (U(i,j,k) - mean_u)**2 + &
                       (V(i,j,k) - mean_v)**2 + &
                       (W(i,j,k) - mean_w)**2

                 delta2 = 1.0_rp/(Jac2D(i,j)*zstep_i(k))**two3
                 asS    = sqrt(2.0_rp*sumSS)
                 T_S    = 1.0_rp/(sqrt(kes/delta2) + asS*i_Ct)

                 mu_turb = i_mu_inf * phi(i,j,k,1) * Cmts * kes * T_S
                 cp_loc  = GAMMA_T(i,j,k)/(GAMMA_T(i,j,k) - 1.0_rp)

                 VIS(i,j,k) = VIS(i,j,k) + mu_turb
                 LMD(i,j,k) = LMD(i,j,k) + cp_loc * i_Pr_turb * mu_turb

              enddo
           enddo
        enddo
        !$acc end parallel

        return
end subroutine compute_MixedTimeScale

end module sgs_module
