module viscous_staggered_module
use profiling_module
use parameters_module
use mpi_module
use mesh_module
use storage_module

implicit none
private
public viscous_flux_3D_staggered
public cross_flux_3D_staggered
public heatflux_3D_staggered

contains

subroutine viscous_flux_3D_staggered
        implicit none
        real(rp), parameter :: two3 = 2.0_rp/3.0_rp

        real(rp)               :: DC_xx, DC_xy, DC_xz
        real(rp)               :: DC_yx, DC_yy, DC_yz
        real(rp)               :: DC_zx, DC_zy, DC_zz

        real(rp)               :: DC_dm_x, DC_dm_y, DC_dm_z
        real(rp)               :: div_sigma_x, div_sigma_y, div_sigma_z
        real(rp)               :: Lapl_x, Lapl_y, Lapl_z
        real(rp)               :: Lx_ip, Ly_ip, Lz_ip
        real(rp)               :: Lx_zt, Ly_zt, Lz_zt
        real(rp)               :: du_x, du_y, du_z
        real(rp)               :: dv_x, dv_y, dv_z
        real(rp)               :: dw_x, dw_y, dw_z
        real(rp)               :: dm_x, dm_y, dm_z
        real(rp)               :: c2R_z, c2L_z
        real(rp)               :: vel_x, vel_y, vel_z
        real(rp)               :: i_st_z
        real(rp)               :: i_stR_z, i_stL_z, Jac
        real(rp)               :: csi_x, csi_y, eta_x, eta_y, tmp_c, tmp_e
        real(rp)               :: mu_, muB, muF
        real(rp)               :: mgE, mgW, mgN, mgS
        real(rp)               :: cl1, Lapl_v, muGradVGradV, DC_dm_v, mu_DC_GradV
        real(rp)               :: div_, div_sigma_v
        real(rp)               :: clxR, clxL, clyR, clyL, clzR, clzL, cl
        integer :: i,j,k,l

        call StartProfRange("viscous_flux_3D_staggered")

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do   k = sz,ez
         do  j = sy,ey
          do i = sx,ex

             vel_x  = U(i,j,k)
             vel_y  = V(i,j,k)
             vel_z  = W(i,j,k)

             mu_ = VIS(i,j,k)
             Jac = Jac2D(i,j)
             csi_x = dxi_dx (i,j)
             csi_y = dxi_dy (i,j)
             eta_x = deta_dx(i,j)
             eta_y = deta_dy(i,j)

             mgE = 0.0_rp
             mgW = 0.0_rp
             mgN = 0.0_rp
             mgS = 0.0_rp
             muF = 0.0_rp
             muB = 0.0_rp

             !$acc loop seq
             do l = -2,3

               clxR = mid_point_lele_x(l,i)
               clyR = mid_point_lele_y(l,j)
               clzR = mid_point_lele_z(l,k)

               clxL = mid_point_lele_x(l,i-1)
               clyL = mid_point_lele_y(l,j-1)
               clzL = mid_point_lele_z(l,k-1)

               mgE = mgE + clxR * VIS(i+l,j,k)   * gJ11(i+l,j)
               mgN = mgN + clyR * VIS(i,j+l,k)   * gJ22(i,j+l)
               mgW = mgW + clxL * VIS(i-1+l,j,k) * gJ11(i-1+l,j)
               mgS = mgS + clyL * VIS(i,j-1+l,k) * gJ22(i,j-1+l)

               muF = muF + clzR * VIS(i,j,k+l)
               muB = muB + clzL * VIS(i,j,k-1+l)

             enddo

             i_st_z = zstep_i(k)

             du_x = 0.0_rp
             du_y = 0.0_rp
             du_z = 0.0_rp
             dv_x = 0.0_rp
             dv_y = 0.0_rp
             dv_z = 0.0_rp
             dw_x = 0.0_rp
             dw_y = 0.0_rp
             dw_z = 0.0_rp
             dm_x = 0.0_rp
             dm_y = 0.0_rp
             dm_z = 0.0_rp

             !$acc loop seq
             do l = 1,3
                cl1 = central_1(l)

                du_x = du_x + cl1 * (U(i+l,j,k) - U(i-l,j,k))
                du_y = du_y + cl1 * (U(i,j+l,k) - U(i,j-l,k))
                du_z = du_z + cl1 * (U(i,j,k+l) - U(i,j,k-l))

                dv_x = dv_x + cl1 * (V(i+l,j,k) - V(i-l,j,k))
                dv_y = dv_y + cl1 * (V(i,j+l,k) - V(i,j-l,k))
                dv_z = dv_z + cl1 * (V(i,j,k+l) - V(i,j,k-l))

                dw_x = dw_x + cl1 * (W(i+l,j,k) - W(i-l,j,k))
                dw_y = dw_y + cl1 * (W(i,j+l,k) - W(i,j-l,k))
                dw_z = dw_z + cl1 * (W(i,j,k+l) - W(i,j,k-l))

                dm_x = dm_x + cl1 * (VIS(i+l,j,k) - VIS(i-l,j,k))
                dm_y = dm_y + cl1 * (VIS(i,j+l,k) - VIS(i,j-l,k))
                dm_z = dm_z + cl1 * (VIS(i,j,k+l) - VIS(i,j,k-l))

             enddo

             tmp_c = du_x
             tmp_e = du_y
             du_x = csi_x*tmp_c + eta_x*tmp_e
             du_y = csi_y*tmp_c + eta_y*tmp_e
             du_z = i_st_z*du_z
             tmp_c = dv_x
             tmp_e = dv_y
             dv_x = csi_x*tmp_c + eta_x*tmp_e
             dv_y = csi_y*tmp_c + eta_y*tmp_e
             dv_z = i_st_z*dv_z
             tmp_c = dw_x
             tmp_e = dw_y
             dw_x = csi_x*tmp_c + eta_x*tmp_e
             dw_y = csi_y*tmp_c + eta_y*tmp_e
             dw_z = i_st_z*dw_z
             tmp_c = dm_x
             tmp_e = dm_y
             dm_x = csi_x*tmp_c + eta_x*tmp_e
             dm_y = csi_y*tmp_c + eta_y*tmp_e
             dm_z = i_st_z*dm_z

             div_ = du_x + dv_y + dw_z
             DIV(i,j,k)  = div_

             DC_xx = du_x - two3 * div_
             DC_xy = dv_x
             DC_xz = dw_x
             DC_yx = du_y
             DC_yy = dv_y - two3 * div_
             DC_yz = dw_y
             DC_zx = du_z
             DC_zy = dv_z
             DC_zz = dw_z - two3 * div_

             DC_dm_x = DC_xx*dm_x + DC_xy*dm_y + DC_xz*dm_z
             DC_dm_y = DC_yx*dm_x + DC_yy*dm_y + DC_yz*dm_z
             DC_dm_z = DC_zx*dm_x + DC_zy*dm_y + DC_zz*dm_z

             muGradVGradV = &
             mu_*(du_x*du_x + du_y*du_y + du_z*du_z + &
                  dv_x*dv_x + dv_y*dv_y + dv_z*dv_z + &
                  dw_x*dw_x + dw_y*dw_y + dw_z*dw_z)

             DC_dm_v = DC_dm_x*vel_x + DC_dm_y*vel_y + DC_dm_z*vel_z

             mu_DC_GradV = &
             mu_*(DC_xx*du_x + DC_xy*du_y + DC_xz*du_z + &
                  DC_yx*dv_x + DC_yy*dv_y + DC_yz*dv_z + &
                  DC_zx*dw_x + DC_zy*dw_y + DC_zz*dw_z)
             i_stR_z = izsteph(k)
             i_stL_z = izsteph(k-1)

             Lx_ip = 0.0_rp
             Ly_ip = 0.0_rp
             Lz_ip = 0.0_rp
             Lx_zt = 0.0_rp
             Ly_zt = 0.0_rp
             Lz_zt = 0.0_rp

             !$acc loop seq
             do l = -2, 3

               cl    = central_2_one_half(l)
               c2R_z = i_stR_z * cl
               c2L_z = i_stL_z * cl

               Lx_ip = Lx_ip + cl*(mgE*U(i+l,j,k) - mgW*U(i-1+l,j,k)) &
                             + cl*(mgN*U(i,j+l,k) - mgS*U(i,j-1+l,k))
               Ly_ip = Ly_ip + cl*(mgE*V(i+l,j,k) - mgW*V(i-1+l,j,k)) &
                             + cl*(mgN*V(i,j+l,k) - mgS*V(i,j-1+l,k))
               Lz_ip = Lz_ip + cl*(mgE*W(i+l,j,k) - mgW*W(i-1+l,j,k)) &
                             + cl*(mgN*W(i,j+l,k) - mgS*W(i,j-1+l,k))

               Lx_zt = Lx_zt + i_st_z*(c2R_z*muF*U(i,j,k+l) - c2L_z*muB*U(i,j,k-1+l))
               Ly_zt = Ly_zt + i_st_z*(c2R_z*muF*V(i,j,k+l) - c2L_z*muB*V(i,j,k-1+l))
               Lz_zt = Lz_zt + i_st_z*(c2R_z*muF*W(i,j,k+l) - c2L_z*muB*W(i,j,k-1+l))

             enddo

             lapl_x    = Jac*Lx_ip + Lx_zt
             lapl_y    = Jac*Ly_ip + Ly_zt
             lapl_z    = Jac*Lz_ip + Lz_zt

             div_sigma_x = lapl_x + DC_dm_x
             div_sigma_y = lapl_y + DC_dm_y
             div_sigma_z = lapl_z + DC_dm_z

             Lapl_v = lapl_x*vel_x + lapl_y*vel_y+lapl_z*vel_z

             div_sigma_v = Lapl_v + muGradVGradV + DC_dm_v + mu_DC_GradV

             RHS(i,j,k,2) = RHS(i,j,k,2) + mu_inf * (div_sigma_x)
             RHS(i,j,k,3) = RHS(i,j,k,3) + mu_inf * (div_sigma_y)
             RHS(i,j,k,4) = RHS(i,j,k,4) + mu_inf * (div_sigma_z)
             RHS(i,j,k,5) = RHS(i,j,k,5) + mu_inf * (div_sigma_v)

          enddo
         enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine viscous_flux_3D_staggered

subroutine cross_flux_3D_staggered
        implicit none
        real(rp)               :: Jac, vel_x, vel_y, vel_z
        real(rp)               :: crx_u, crx_v, crx_w, cl1, cn
        real(rp)               :: AxiP, AxiM, AetP, AetM
        real(rp)               :: euP, euM, evP, evM, ewP, ewM
        real(rp)               :: xuP, xuM, xvP, xvM, xwP, xwM
        integer :: i,j,k,l,n

        call StartProfRange("cross_flux_3D_staggered")

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do   k = sz,ez
         do  j = sy,ey
          do i = sx,ex

             Jac   = Jac2D(i,j)
             vel_x = U(i,j,k)
             vel_y = V(i,j,k)
             vel_z = W(i,j,k)

             crx_u = 0.0_rp
             crx_v = 0.0_rp
             crx_w = 0.0_rp

             !$acc loop seq
             do l = 1,3
                cl1  = central_1(l)

                AxiP = VIS(i+l,j,k)*gJ12(i+l,j)
                AxiM = VIS(i-l,j,k)*gJ12(i-l,j)
                AetP = VIS(i,j+l,k)*gJ12(i,j+l)
                AetM = VIS(i,j-l,k)*gJ12(i,j-l)

                euP = 0.0_rp
                euM = 0.0_rp
                evP = 0.0_rp
                evM = 0.0_rp
                ewP = 0.0_rp
                ewM = 0.0_rp
                xuP = 0.0_rp
                xuM = 0.0_rp
                xvP = 0.0_rp
                xvM = 0.0_rp
                xwP = 0.0_rp
                xwM = 0.0_rp

                !$acc loop seq
                do n = 1,3
                   cn = central_1(n)
                   euP = euP + cn*(U(i+l,j+n,k)-U(i+l,j-n,k))
                   euM = euM + cn*(U(i-l,j+n,k)-U(i-l,j-n,k))
                   evP = evP + cn*(V(i+l,j+n,k)-V(i+l,j-n,k))
                   evM = evM + cn*(V(i-l,j+n,k)-V(i-l,j-n,k))
                   ewP = ewP + cn*(W(i+l,j+n,k)-W(i+l,j-n,k))
                   ewM = ewM + cn*(W(i-l,j+n,k)-W(i-l,j-n,k))
                   xuP = xuP + cn*(U(i+n,j+l,k)-U(i-n,j+l,k))
                   xuM = xuM + cn*(U(i+n,j-l,k)-U(i-n,j-l,k))
                   xvP = xvP + cn*(V(i+n,j+l,k)-V(i-n,j+l,k))
                   xvM = xvM + cn*(V(i+n,j-l,k)-V(i-n,j-l,k))
                   xwP = xwP + cn*(W(i+n,j+l,k)-W(i-n,j+l,k))
                   xwM = xwM + cn*(W(i+n,j-l,k)-W(i-n,j-l,k))
                enddo

                crx_u = crx_u + cl1*(AxiP*euP - AxiM*euM) + cl1*(AetP*xuP - AetM*xuM)
                crx_v = crx_v + cl1*(AxiP*evP - AxiM*evM) + cl1*(AetP*xvP - AetM*xvM)
                crx_w = crx_w + cl1*(AxiP*ewP - AxiM*ewM) + cl1*(AetP*xwP - AetM*xwM)
             enddo

             crx_u = Jac*crx_u
             crx_v = Jac*crx_v
             crx_w = Jac*crx_w

             RHS(i,j,k,2) = RHS(i,j,k,2) + mu_inf * (crx_u)
             RHS(i,j,k,3) = RHS(i,j,k,3) + mu_inf * (crx_v)
             RHS(i,j,k,4) = RHS(i,j,k,4) + mu_inf * (crx_w)
             RHS(i,j,k,5) = RHS(i,j,k,5) + mu_inf * (crx_u*vel_x + crx_v*vel_y + crx_w*vel_z)

          enddo
         enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine cross_flux_3D_staggered

subroutine heatflux_3D_staggered
        implicit none
        real(rp)               :: Jac, i_st_z, i_stR_z, i_stL_z
        real(rp)               :: rgE, rgW, rgN, rgS, rlF, rlB
        real(rp)               :: Hf_ip, Hf_zt, heat_flux
        real(rp)               :: cl, cl1, cn, c2R_z, c2L_z
        real(rp)               :: clxR, clxL, clyR, clyL, clzR, clzL
        real(rp)               :: HxiP, HxiM, HetP, HetM
        real(rp)               :: eTP, eTM, xTP, xTM, crx_T
        integer :: i,j,k,l,n

        call StartProfRange("heatflux_3D_staggered")

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do   k = sz,ez
         do  j = sy,ey
          do i = sx,ex

             Jac = Jac2D(i,j)

             rgE = 0.0_rp
             rgW = 0.0_rp
             rgN = 0.0_rp
             rgS = 0.0_rp
             rlF = 0.0_rp
             rlB = 0.0_rp
             !$acc loop seq
             do l = -2,3
               clxR = mid_point_lele_x(l,i)
               clyR = mid_point_lele_y(l,j)
               clzR = mid_point_lele_z(l,k)
               clxL = mid_point_lele_x(l,i-1)
               clyL = mid_point_lele_y(l,j-1)
               clzL = mid_point_lele_z(l,k-1)
               rgE = rgE + clxR * LMD(i+l,j,k)   * gJ11(i+l,j)
               rgN = rgN + clyR * LMD(i,j+l,k)   * gJ22(i,j+l)
               rgW = rgW + clxL * LMD(i-1+l,j,k) * gJ11(i-1+l,j)
               rgS = rgS + clyL * LMD(i,j-1+l,k) * gJ22(i,j-1+l)
               rlF = rlF + clzR * LMD(i,j,k+l)
               rlB = rlB + clzL * LMD(i,j,k-1+l)
             enddo

             i_st_z  = zstep_i(k)
             i_stR_z = izsteph(k)
             i_stL_z = izsteph(k-1)
             Hf_ip = 0.0_rp
             Hf_zt = 0.0_rp
             !$acc loop seq
             do l = -2, 3
               cl    = central_2_one_half(l)
               c2R_z = i_stR_z * cl
               c2L_z = i_stL_z * cl
               Hf_ip = Hf_ip + cl*(rgE*T(i+l,j,k) - rgW*T(i-1+l,j,k)) &
                             + cl*(rgN*T(i,j+l,k) - rgS*T(i,j-1+l,k))
               Hf_zt = Hf_zt + i_st_z*(c2R_z*rlF*T(i,j,k+l) - c2L_z*rlB*T(i,j,k-1+l))
             enddo

             crx_T = 0.0_rp
             if (non_orthogonal) then
             !$acc loop seq
             do l = 1,3
                cl1  = central_1(l)
                HxiP = LMD(i+l,j,k)*gJ12(i+l,j)
                HxiM = LMD(i-l,j,k)*gJ12(i-l,j)
                HetP = LMD(i,j+l,k)*gJ12(i,j+l)
                HetM = LMD(i,j-l,k)*gJ12(i,j-l)
                eTP = 0.0_rp
                eTM = 0.0_rp
                xTP = 0.0_rp
                xTM = 0.0_rp
                !$acc loop seq
                do n = 1,3
                   cn = central_1(n)
                   eTP = eTP + cn*(T(i+l,j+n,k)-T(i+l,j-n,k))
                   eTM = eTM + cn*(T(i-l,j+n,k)-T(i-l,j-n,k))
                   xTP = xTP + cn*(T(i+n,j+l,k)-T(i-n,j+l,k))
                   xTM = xTM + cn*(T(i+n,j-l,k)-T(i-n,j-l,k))
                enddo
                crx_T = crx_T + cl1*(HxiP*eTP - HxiM*eTM) + cl1*(HetP*xTP - HetM*xTM)
             enddo
             endif

             Hf_ip = Hf_ip + crx_T
             heat_flux = Jac*Hf_ip + Hf_zt

             RHS(i,j,k,5) = RHS(i,j,k,5) + mu_inf * (heat_flux)

          enddo
         enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine heatflux_3D_staggered

end module viscous_staggered_module
