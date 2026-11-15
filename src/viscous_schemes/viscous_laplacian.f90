module viscous_laplacian_module
use profiling_module
use parameters_module
use mpi_module
use mesh_module
use storage_module

implicit none
private
public :: viscous_flux_3D_laplacian, viscous_heatflux_3D_laplacian

contains

subroutine viscous_heatflux_3D_laplacian
        use parameters_module, only: rp, mu_inf
        use mpi_module       , only: sx,ex,sy,ey,sz,ez
        use mesh_module      , only: dxi_dx, dxi_dy, deta_dx, deta_dy, &
                                     Gamma1_2D, Gamma2_2D,             &
                                     ztastep_i, iz_zta, z_zta2

        implicit none
        real(rp) :: csi_x, csi_y, eta_x, eta_y
        real(rp) :: g11, g12, g22, gam1, gam2
        real(rp) :: i_st_z, i_st2_z, iz_z, d2zeta_dz2
        real(rp) :: dTc,dTe,dTz, dlc,dle,dlz
        real(rp) :: d2Tcc,d2Tee,d2Tce,d2Tzz
        real(rp) :: Tx,Ty,Tz, lx,ly,lz
        real(rp) :: lap_T, rl_, heat_flux, cl1, cl2, clm
        integer  :: l, m, i, j, k

        call StartProfRange("viscous_heatflux_3D_laplacian")

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do       k = sz,ez
           do    j = sy,ey
              do i = sx,ex

                 rl_ = LMD(i,j,k)

                 csi_x = dxi_dx (i,j)
                 csi_y = dxi_dy (i,j)
                 eta_x = deta_dx(i,j)
                 eta_y = deta_dy(i,j)
                 g11   = csi_x*csi_x + csi_y*csi_y
                 g22   = eta_x*eta_x + eta_y*eta_y
                 g12   = csi_x*eta_x + csi_y*eta_y
                 gam1  = Gamma1_2D(i,j)
                 gam2  = Gamma2_2D(i,j)

                 i_st_z     = ztastep_i(k)
                 i_st2_z    = i_st_z*i_st_z
                 iz_z       = iz_zta(k)
                 d2zeta_dz2 = - z_zta2(k)*iz_z**3

                 dTc = 0.0_rp
                 dTe = 0.0_rp
                 dTz = 0.0_rp
                 dlc = 0.0_rp
                 dle = 0.0_rp
                 dlz = 0.0_rp
                 cl2   = central_2(0)
                 d2Tcc = cl2*T(i,j,k)
                 d2Tee = cl2*T(i,j,k)
                 d2Tzz = cl2*T(i,j,k)
                 !$acc loop seq
                 do l = 1,3
                    cl1 = central_1(l)
                    cl2 = central_2(l)
                    dTc = dTc + cl1*(T(i+l,j,k) - T(i-l,j,k))
                    dTe = dTe + cl1*(T(i,j+l,k) - T(i,j-l,k))
                    dTz = dTz + cl1*(T(i,j,k+l) - T(i,j,k-l))
                    dlc = dlc + cl1*(LMD(i+l,j,k) - LMD(i-l,j,k))
                    dle = dle + cl1*(LMD(i,j+l,k) - LMD(i,j-l,k))
                    dlz = dlz + cl1*(LMD(i,j,k+l) - LMD(i,j,k-l))
                    d2Tcc = d2Tcc + cl2*(T(i+l,j,k) + T(i-l,j,k))
                    d2Tee = d2Tee + cl2*(T(i,j+l,k) + T(i,j-l,k))
                    d2Tzz = d2Tzz + cl2*(T(i,j,k+l) + T(i,j,k-l))
                 enddo
                 d2Tce = 0.0_rp
                 !$acc loop seq
                 do m = 1,3
                  !$acc loop seq
                  do l = 1,3
                     clm = central_1(l)*central_1(m)
                     d2Tce = d2Tce + clm*(T(i+l,j+m,k)-T(i-l,j+m,k)-T(i+l,j-m,k)+T(i-l,j-m,k))
                  enddo
                 enddo

                 dTz = i_st_z*dTz
                 dlz = i_st_z*dlz
                 d2Tzz = d2Tzz*i_st2_z*iz_z**2 + dTz*d2zeta_dz2
                 Tz = dTz*iz_z
                 lz = dlz*iz_z
                 Tx = csi_x*dTc + eta_x*dTe
                 Ty = csi_y*dTc + eta_y*dTe
                 lx = csi_x*dlc + eta_x*dle
                 ly = csi_y*dlc + eta_y*dle
                 lap_T = g11*d2Tcc + 2.0_rp*g12*d2Tce + g22*d2Tee + gam1*dTc + gam2*dTe + d2Tzz
                 heat_flux = rl_*lap_T + lx*Tx + ly*Ty + lz*Tz

                 RHS(i,j,k,5) = RHS(i,j,k,5) + mu_inf*heat_flux

              enddo
           enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine viscous_heatflux_3D_laplacian

subroutine viscous_flux_3D_laplacian
        use parameters_module, only: rp, mu_inf
        use mpi_module       , only: sx,ex,sy,ey,sz,ez
        use mesh_module      , only: dxi_dx, dxi_dy, deta_dx, deta_dy, &
                                     Gamma1_2D, Gamma2_2D,             &
                                     ztastep_i, iz_zta, z_zta2

        implicit none
        real(rp), parameter :: two3 = 2.0_rp/3.0_rp

        real(rp) :: csi_x, csi_y, eta_x, eta_y
        real(rp) :: g11, g12, g22, gam1, gam2
        real(rp) :: i_st_z, i_st2_z, iz_z, d2zeta_dz2
        real(rp) :: duc,due,duz, dvc,dve,dvz, dwc,dwe,dwz
        real(rp) :: dmc,dme,dmz
        real(rp) :: d2ucc,d2uee,d2uce,d2uzz
        real(rp) :: d2vcc,d2vee,d2vce,d2vzz
        real(rp) :: d2wcc,d2wee,d2wce,d2wzz
        real(rp) :: ux,uy,uz, vx,vy,vz, wx,wy,wz
        real(rp) :: mx,my,mz
        real(rp) :: lap_u, lap_v, lap_w
        real(rp) :: u_, v_, w_, mu_
        real(rp) :: cl1, cl2, clm
        real(rp) :: div_
        real(rp) :: D11,D12,D13, D22,D23, D33
        real(rp) :: D_dm_1, D_dm_2, D_dm_3
        real(rp) :: div_sigma_1, div_sigma_2, div_sigma_3
        real(rp) :: sum_D_grad, div_sigma_v
        integer  :: l, m, i, j, k

        call StartProfRange("viscous_flux_3D_laplacian")

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do       k = sz,ez
           do    j = sy,ey
              do i = sx,ex

                 u_  = U(i,j,k)
                 v_  = V(i,j,k)
                 w_  = W(i,j,k)
                 mu_ = VIS(i,j,k)

                 csi_x = dxi_dx (i,j)
                 csi_y = dxi_dy (i,j)
                 eta_x = deta_dx(i,j)
                 eta_y = deta_dy(i,j)
                 g11   = csi_x*csi_x + csi_y*csi_y
                 g22   = eta_x*eta_x + eta_y*eta_y
                 g12   = csi_x*eta_x + csi_y*eta_y
                 gam1  = Gamma1_2D(i,j)
                 gam2  = Gamma2_2D(i,j)

                 i_st_z     = ztastep_i(k)
                 i_st2_z    = i_st_z*i_st_z
                 iz_z       = iz_zta(k)
                 d2zeta_dz2 = - z_zta2(k)*iz_z**3

                 duc = 0.0_rp
                 due = 0.0_rp
                 duz = 0.0_rp
                 dvc = 0.0_rp
                 dve = 0.0_rp
                 dvz = 0.0_rp
                 dwc = 0.0_rp
                 dwe = 0.0_rp
                 dwz = 0.0_rp
                 dmc = 0.0_rp
                 dme = 0.0_rp
                 dmz = 0.0_rp
                 cl2   = central_2(0)
                 d2ucc = cl2*U(i,j,k)
                 d2uee = cl2*U(i,j,k)
                 d2uzz = cl2*U(i,j,k)
                 d2vcc = cl2*V(i,j,k)
                 d2vee = cl2*V(i,j,k)
                 d2vzz = cl2*V(i,j,k)
                 d2wcc = cl2*W(i,j,k)
                 d2wee = cl2*W(i,j,k)
                 d2wzz = cl2*W(i,j,k)
                 !$acc loop seq
                 do l = 1,3
                    cl1 = central_1(l)
                    cl2 = central_2(l)

                    duc = duc + cl1*(U(i+l,j,k) - U(i-l,j,k))
                    due = due + cl1*(U(i,j+l,k) - U(i,j-l,k))
                    duz = duz + cl1*(U(i,j,k+l) - U(i,j,k-l))
                    dvc = dvc + cl1*(V(i+l,j,k) - V(i-l,j,k))
                    dve = dve + cl1*(V(i,j+l,k) - V(i,j-l,k))
                    dvz = dvz + cl1*(V(i,j,k+l) - V(i,j,k-l))
                    dwc = dwc + cl1*(W(i+l,j,k) - W(i-l,j,k))
                    dwe = dwe + cl1*(W(i,j+l,k) - W(i,j-l,k))
                    dwz = dwz + cl1*(W(i,j,k+l) - W(i,j,k-l))
                    dmc = dmc + cl1*(VIS(i+l,j,k) - VIS(i-l,j,k))
                    dme = dme + cl1*(VIS(i,j+l,k) - VIS(i,j-l,k))
                    dmz = dmz + cl1*(VIS(i,j,k+l) - VIS(i,j,k-l))

                    d2ucc = d2ucc + cl2*(U(i+l,j,k) + U(i-l,j,k))
                    d2uee = d2uee + cl2*(U(i,j+l,k) + U(i,j-l,k))
                    d2uzz = d2uzz + cl2*(U(i,j,k+l) + U(i,j,k-l))
                    d2vcc = d2vcc + cl2*(V(i+l,j,k) + V(i-l,j,k))
                    d2vee = d2vee + cl2*(V(i,j+l,k) + V(i,j-l,k))
                    d2vzz = d2vzz + cl2*(V(i,j,k+l) + V(i,j,k-l))
                    d2wcc = d2wcc + cl2*(W(i+l,j,k) + W(i-l,j,k))
                    d2wee = d2wee + cl2*(W(i,j+l,k) + W(i,j-l,k))
                    d2wzz = d2wzz + cl2*(W(i,j,k+l) + W(i,j,k-l))
                 enddo

                 d2uce = 0.0_rp
                 d2vce = 0.0_rp
                 d2wce = 0.0_rp
                 !$acc loop seq
                 do m = 1,3
                  !$acc loop seq
                  do l = 1,3
                     clm = central_1(l)*central_1(m)
                     d2uce = d2uce + clm*(U(i+l,j+m,k)-U(i-l,j+m,k)-U(i+l,j-m,k)+U(i-l,j-m,k))
                     d2vce = d2vce + clm*(V(i+l,j+m,k)-V(i-l,j+m,k)-V(i+l,j-m,k)+V(i-l,j-m,k))
                     d2wce = d2wce + clm*(W(i+l,j+m,k)-W(i-l,j+m,k)-W(i+l,j-m,k)+W(i-l,j-m,k))
                  enddo
                 enddo

                 duz = i_st_z*duz
                 dvz = i_st_z*dvz
                 dwz = i_st_z*dwz
                 dmz = i_st_z*dmz
                 d2uzz = d2uzz*i_st2_z*iz_z**2 + duz*d2zeta_dz2
                 d2vzz = d2vzz*i_st2_z*iz_z**2 + dvz*d2zeta_dz2
                 d2wzz = d2wzz*i_st2_z*iz_z**2 + dwz*d2zeta_dz2
                 uz = duz*iz_z
                 vz = dvz*iz_z
                 wz = dwz*iz_z
                 mz = dmz*iz_z

                 ux = csi_x*duc + eta_x*due
                 uy = csi_y*duc + eta_y*due
                 vx = csi_x*dvc + eta_x*dve
                 vy = csi_y*dvc + eta_y*dve
                 wx = csi_x*dwc + eta_x*dwe
                 wy = csi_y*dwc + eta_y*dwe
                 mx = csi_x*dmc + eta_x*dme
                 my = csi_y*dmc + eta_y*dme

                 lap_u = g11*d2ucc + 2.0_rp*g12*d2uce + g22*d2uee + gam1*duc + gam2*due + d2uzz
                 lap_v = g11*d2vcc + 2.0_rp*g12*d2vce + g22*d2vee + gam1*dvc + gam2*dve + d2vzz
                 lap_w = g11*d2wcc + 2.0_rp*g12*d2wce + g22*d2wee + gam1*dwc + gam2*dwe + d2wzz

                 div_ = ux + vy + wz
                 DIV(i,j,k) = div_

                 D11 = ux + ux - two3*div_
                 D12 = uy + vx
                 D13 = uz + wx
                 D22 = vy + vy - two3*div_
                 D23 = vz + wy
                 D33 = wz + wz - two3*div_

                 D_dm_1 = D11*mx + D12*my + D13*mz
                 D_dm_2 = D12*mx + D22*my + D23*mz
                 D_dm_3 = D13*mx + D23*my + D33*mz

                 div_sigma_1 = mu_*lap_u + D_dm_1
                 div_sigma_2 = mu_*lap_v + D_dm_2
                 div_sigma_3 = mu_*lap_w + D_dm_3

                 sum_D_grad = D11*ux + D12*uy + D13*uz + &
                              D12*vx + D22*vy + D23*vz + &
                              D13*wx + D23*wy + D33*wz

                 div_sigma_v = div_sigma_1*u_ + div_sigma_2*v_ + div_sigma_3*w_ + mu_*sum_D_grad

                 RHS(i,j,k,2) = RHS(i,j,k,2) + mu_inf*(div_sigma_1)
                 RHS(i,j,k,3) = RHS(i,j,k,3) + mu_inf*(div_sigma_2)
                 RHS(i,j,k,4) = RHS(i,j,k,4) + mu_inf*(div_sigma_3)
                 RHS(i,j,k,5) = RHS(i,j,k,5) + mu_inf*(div_sigma_v)

              enddo
           enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine viscous_flux_3D_laplacian

end module viscous_laplacian_module
