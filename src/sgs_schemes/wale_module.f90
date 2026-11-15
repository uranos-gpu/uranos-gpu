module wale_module

        use parameters_module, only: rp, central_1, mu_inf
        use mpi_module       , only: sx, ex, sy, ey, sz, ez
        use storage_module   , only: phi, U, V, W, GAMMA_T
        use mesh_module      , only: zstep_i, ztastep_i, iz_zta, Jac2D, &
                                     dxi_dx, dxi_dy, deta_dx, deta_dy

        implicit none
        private

        real(rp), parameter :: Pr_turb   = 0.9_rp
        real(rp), parameter :: i_Pr_turb = 1.0_rp/Pr_turb
        real(rp), parameter :: one3      = 1.0_rp/3.0_rp
        real(rp), parameter :: two3      = 2.0_rp/3.0_rp

        public :: compute_wale

contains

subroutine compute_wale(VIS,LMD)

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: VIS, LMD

        real(rp), parameter :: Cw2  = 0.325_rp**2
        real(rp), parameter :: toll = 1.0E-14_rp

        real(rp) :: csi_x, csi_y, eta_x, eta_y, zeta_z
        real(rp) :: duc, due, duz
        real(rp) :: dvc, dve, dvz
        real(rp) :: dwc, dwe, dwz
        real(rp) :: du_x, du_y, du_z
        real(rp) :: dv_x, dv_y, dv_z
        real(rp) :: dw_x, dw_y, dw_z
        real(rp) :: S_xy, S_xz, S_yz
        real(rp) :: DV_xx, DV_yy, DV_zz, SdTrace
        real(rp) :: Sd_xx, Sd_yy, Sd_zz, Sd_xy, Sd_xz, Sd_yz
        real(rp) :: sumSS, sumSdSd
        real(rp) :: SdSd_32, SdSd_54, S_S__52
        real(rp) :: mu_turb, cp_loc, delta2, Cw2_imu
        integer  :: i, j, k

        Cw2_imu = Cw2/mu_inf

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex

                 csi_x = dxi_dx (i,j); csi_y = dxi_dy (i,j)
                 eta_x = deta_dx(i,j); eta_y = deta_dy(i,j)
                 zeta_z = ztastep_i(k)*iz_zta(k)

                 duc = central_1(1) * (U(i+1,j,k) - U(i-1,j,k)) &
                     + central_1(2) * (U(i+2,j,k) - U(i-2,j,k)) &
                     + central_1(3) * (U(i+3,j,k) - U(i-3,j,k))
                 due = central_1(1) * (U(i,j+1,k) - U(i,j-1,k)) &
                     + central_1(2) * (U(i,j+2,k) - U(i,j-2,k)) &
                     + central_1(3) * (U(i,j+3,k) - U(i,j-3,k))
                 duz = central_1(1) * (U(i,j,k+1) - U(i,j,k-1)) &
                     + central_1(2) * (U(i,j,k+2) - U(i,j,k-2)) &
                     + central_1(3) * (U(i,j,k+3) - U(i,j,k-3))

                 dvc = central_1(1) * (V(i+1,j,k) - V(i-1,j,k)) &
                     + central_1(2) * (V(i+2,j,k) - V(i-2,j,k)) &
                     + central_1(3) * (V(i+3,j,k) - V(i-3,j,k))
                 dve = central_1(1) * (V(i,j+1,k) - V(i,j-1,k)) &
                     + central_1(2) * (V(i,j+2,k) - V(i,j-2,k)) &
                     + central_1(3) * (V(i,j+3,k) - V(i,j-3,k))
                 dvz = central_1(1) * (V(i,j,k+1) - V(i,j,k-1)) &
                     + central_1(2) * (V(i,j,k+2) - V(i,j,k-2)) &
                     + central_1(3) * (V(i,j,k+3) - V(i,j,k-3))

                 dwc = central_1(1) * (W(i+1,j,k) - W(i-1,j,k)) &
                     + central_1(2) * (W(i+2,j,k) - W(i-2,j,k)) &
                     + central_1(3) * (W(i+3,j,k) - W(i-3,j,k))
                 dwe = central_1(1) * (W(i,j+1,k) - W(i,j-1,k)) &
                     + central_1(2) * (W(i,j+2,k) - W(i,j-2,k)) &
                     + central_1(3) * (W(i,j+3,k) - W(i,j-3,k))
                 dwz = central_1(1) * (W(i,j,k+1) - W(i,j,k-1)) &
                     + central_1(2) * (W(i,j,k+2) - W(i,j,k-2)) &
                     + central_1(3) * (W(i,j,k+3) - W(i,j,k-3))

                 du_x = csi_x*duc + eta_x*due; du_y = csi_y*duc + eta_y*due; du_z = zeta_z*duz
                 dv_x = csi_x*dvc + eta_x*dve; dv_y = csi_y*dvc + eta_y*dve; dv_z = zeta_z*dvz
                 dw_x = csi_x*dwc + eta_x*dwe; dw_y = csi_y*dwc + eta_y*dwe; dw_z = zeta_z*dwz

                 S_xy = 0.5_rp*(du_y + dv_x)
                 S_xz = 0.5_rp*(du_z + dw_x)
                 S_yz = 0.5_rp*(dv_z + dw_y)

                 sumSS = du_x*du_x + dv_y*dv_y + dw_z*dw_z + &
                         2.0_rp*(S_xy*S_xy + S_xz*S_xz + S_yz*S_yz)

                 DV_xx = du_x*du_x + du_y*dv_x + du_z*dw_x
                 DV_yy = dv_x*du_y + dv_y*dv_y + dv_z*dw_y
                 DV_zz = dw_x*du_z + dw_y*dv_z + dw_z*dw_z

                 SdTrace = -one3*(DV_xx + DV_yy + DV_zz)

                 Sd_xx = DV_xx + SdTrace
                 Sd_yy = DV_yy + SdTrace
                 Sd_zz = DV_zz + SdTrace
                 Sd_xy = (du_x + dv_y)*S_xy + 0.5_rp*(du_z*dw_y + dv_z*dw_x)
                 Sd_xz = (du_x + dw_z)*S_xz + 0.5_rp*(du_y*dv_z + dw_y*dv_x)
                 Sd_yz = (dv_y + dw_z)*S_yz + 0.5_rp*(dv_x*du_z + dw_x*du_y)

                 sumSdSd = Sd_xx*Sd_xx + Sd_yy*Sd_yy + Sd_zz*Sd_zz + &
                           2.0_rp*(Sd_xy*Sd_xy + Sd_xz*Sd_xz + Sd_yz*Sd_yz)

                 SdSd_32 = sumSdSd * sqrt(sumSdSd)

                 if (SdSd_32 > toll) then
                    SdSd_54 = sqrt(SdSd_32 * sumSdSd)
                    S_S__52 = sumSS * sumSS * sqrt(sumSS)

                    delta2  = 1.0_rp/(Jac2D(i,j)*zstep_i(k))**two3
                    mu_turb = Cw2_imu * phi(i,j,k,1) * delta2 &
                                      * SdSd_32 / (S_S__52 + SdSd_54)
                    cp_loc  = GAMMA_T(i,j,k)/(GAMMA_T(i,j,k) - 1.0_rp)

                    VIS(i,j,k) = VIS(i,j,k) + mu_turb
                    LMD(i,j,k) = LMD(i,j,k) + cp_loc * i_Pr_turb * mu_turb
                 endif

              enddo
           enddo
        enddo
        !$acc end parallel

        return
end subroutine compute_wale

end module wale_module
