module viscous_shared_module
use profiling_module
use parameters_module
use mpi_module
use mpi_comm_module
use mesh_module
use storage_module

implicit none
private
public DivergencyGradient3D, ComputeSGSIsotropicTerm, sgs_isotropic_setup

integer , parameter, private :: fW=1, fE=2, fS=3, fN=4, fB=5, fF=6

contains

subroutine DivergencyGradient3D(U,V,W,VIS,RHS)
        implicit none
        real(rp), dimension(:,:,:)  , allocatable, intent(in)    :: U,V,W,VIS
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: RHS

        if(inplane_axis_aligned) then
           call DivergencyGradient3D_cart(U,V,W,VIS,RHS)
        else
           call DivergencyGradient3D_curv(U,V,W,VIS,RHS)
        endif

        return
end subroutine DivergencyGradient3D

subroutine DivergencyGradient3D_curv(U,V,W,VIS,RHS)
        implicit none
        real(rp), dimension(:,:,:)  , allocatable, intent(in)    :: U,V,W,VIS
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: RHS

        real(rp) :: One3MuInf, mu_, cl1, cl2, clm
        real(rp) :: a, b, c, d
        real(rp) :: cxx, cyy, cxy, exx, eyy, exy
        real(rp) :: i_st_z, i_st2_z, iz_z, d2zeta_dz2, zeta_z
        real(rp) :: uc, ue, vc, ve, wz
        real(rp) :: d2ucc, d2uee, d2uce, d2vcc, d2vee, d2vce, d2wzz
        real(rp) :: wcz, wez, s_z
        real(rp) :: u_xx, u_xy, v_xy, v_yy
        real(rp) :: w_xz, w_yz, w_zz
        real(rp) :: graddiv_x, graddiv_y, graddiv_z
        real(rp) :: RhU_term, RhV_term, RhW_term, RhE_term
        integer  :: i,j,k,l,m

        call StartProfRange("DivergencyGradient3D_curv")

        One3MuInf = 1.0_rp/3.0_rp*mu_inf

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex

                 a   = dxi_dx (i,j);  b   = dxi_dy (i,j)
                 c   = deta_dx(i,j);  d   = deta_dy(i,j)
                 cxx = dxi_dxx(i,j);  cyy = dxi_dyy(i,j);  cxy = dxi_dxy(i,j)
                 exx = deta_dxx(i,j); eyy = deta_dyy(i,j); exy = deta_dxy(i,j)

                 i_st_z     = ztastep_i(k)
                 i_st2_z    = i_st_z*i_st_z
                 iz_z       = iz_zta(k)
                 zeta_z     = i_st_z*iz_z
                 d2zeta_dz2 = - z_zta2(k)*iz_z**3

                 uc = 0.0_rp; ue = 0.0_rp
                 vc = 0.0_rp; ve = 0.0_rp
                 wz = 0.0_rp
                 cl2   = central_2(0)
                 d2ucc = cl2*U(i,j,k); d2uee = cl2*U(i,j,k)
                 d2vcc = cl2*V(i,j,k); d2vee = cl2*V(i,j,k)
                 d2wzz = cl2*W(i,j,k)
                 !$acc loop seq
                 do l = 1,3
                    cl1 = central_1(l)
                    cl2 = central_2(l)

                    uc = uc + cl1*(U(i+l,j,k) - U(i-l,j,k))
                    ue = ue + cl1*(U(i,j+l,k) - U(i,j-l,k))
                    vc = vc + cl1*(V(i+l,j,k) - V(i-l,j,k))
                    ve = ve + cl1*(V(i,j+l,k) - V(i,j-l,k))
                    wz = wz + cl1*(W(i,j,k+l) - W(i,j,k-l))

                    d2ucc = d2ucc + cl2*(U(i+l,j,k) + U(i-l,j,k))
                    d2uee = d2uee + cl2*(U(i,j+l,k) + U(i,j-l,k))
                    d2vcc = d2vcc + cl2*(V(i+l,j,k) + V(i-l,j,k))
                    d2vee = d2vee + cl2*(V(i,j+l,k) + V(i,j-l,k))
                    d2wzz = d2wzz + cl2*(W(i,j,k+l) + W(i,j,k-l))
                 enddo

                 d2uce = 0.0_rp; d2vce = 0.0_rp
                 !$acc loop seq
                 do m = 1,3
                  !$acc loop seq
                  do l = 1,3
                     clm = central_1(l)*central_1(m)
                     d2uce = d2uce + clm*(U(i+l,j+m,k)-U(i-l,j+m,k)-U(i+l,j-m,k)+U(i-l,j-m,k))
                     d2vce = d2vce + clm*(V(i+l,j+m,k)-V(i-l,j+m,k)-V(i+l,j-m,k)+V(i-l,j-m,k))
                  enddo
                 enddo

                 wcz = 0.0_rp; wez = 0.0_rp; s_z = 0.0_rp
                 !$acc loop seq
                 do m = 1,3
                  !$acc loop seq
                  do l = 1,3
                     clm = central_1(l)*central_1(m)
                     wcz = wcz + clm*(W(i+l,j,k+m)-W(i-l,j,k+m)-W(i+l,j,k-m)+W(i-l,j,k-m))
                     wez = wez + clm*(W(i,j+l,k+m)-W(i,j-l,k+m)-W(i,j+l,k-m)+W(i,j-l,k-m))
                     s_z = s_z + clm*( a*(U(i+l,j,k+m)-U(i-l,j,k+m)-U(i+l,j,k-m)+U(i-l,j,k-m))   &
                                     + c*(U(i,j+l,k+m)-U(i,j-l,k+m)-U(i,j+l,k-m)+U(i,j-l,k-m))   &
                                     + b*(V(i+l,j,k+m)-V(i-l,j,k+m)-V(i+l,j,k-m)+V(i-l,j,k-m))   &
                                     + d*(V(i,j+l,k+m)-V(i,j-l,k+m)-V(i,j+l,k-m)+V(i,j-l,k-m)) )
                  enddo
                 enddo

                 u_xx = a*a*d2ucc + 2.0_rp*a*c*d2uce + c*c*d2uee + cxx*uc + exx*ue
                 u_xy = a*b*d2ucc + (a*d+c*b)*d2uce + c*d*d2uee + cxy*uc + exy*ue
                 v_xy = a*b*d2vcc + (a*d+c*b)*d2vce + c*d*d2vee + cxy*vc + exy*ve
                 v_yy = b*b*d2vcc + 2.0_rp*b*d*d2vce + d*d*d2vee + cyy*vc + eyy*ve

                 w_xz = zeta_z*(a*wcz + c*wez)
                 w_yz = zeta_z*(b*wcz + d*wez)
                 w_zz = d2wzz*i_st2_z*iz_z**2 + (i_st_z*wz)*d2zeta_dz2

                 graddiv_x = u_xx + v_xy + w_xz
                 graddiv_y = u_xy + v_yy + w_yz
                 graddiv_z = zeta_z*s_z + w_zz

                 mu_      = One3MuInf*VIS(i,j,k)

                 RhU_term = mu_*graddiv_x
                 RhV_term = mu_*graddiv_y
                 RhW_term = mu_*graddiv_z
                 RhE_term = RhU_term*U(i,j,k) + RhV_term*V(i,j,k) + RhW_term*W(i,j,k)

                 RHS(i,j,k,2) = RHS(i,j,k,2) + RhU_term
                 RHS(i,j,k,3) = RHS(i,j,k,3) + RhV_term
                 RHS(i,j,k,4) = RHS(i,j,k,4) + RhW_term
                 RHS(i,j,k,5) = RHS(i,j,k,5) + RhE_term

              enddo
           enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine DivergencyGradient3D_curv

subroutine DivergencyGradient3D_cart(U,V,W,VIS,RHS)
        implicit none
        real(rp), dimension(:,:,:)  , allocatable, intent(in)    :: U,V,W,VIS
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: RHS

        real(rp) :: One3MuInf, mu_, cl1, cl2, clm
        real(rp) :: a, d
        real(rp) :: cxx, cyy, exx, eyy
        real(rp) :: i_st_z, i_st2_z, iz_z, d2zeta_dz2, zeta_z
        real(rp) :: uc, ue, vc, ve, wz
        real(rp) :: d2ucc, d2uce, d2vce, d2vee, d2wzz
        real(rp) :: wcz, wez, s_z
        real(rp) :: u_xx, u_xy, v_xy, v_yy
        real(rp) :: w_xz, w_yz, w_zz
        real(rp) :: graddiv_x, graddiv_y, graddiv_z
        real(rp) :: RhU_term, RhV_term, RhW_term, RhE_term
        integer  :: i,j,k,l,m

        call StartProfRange("DivergencyGradient3D_cart")

        One3MuInf = 1.0_rp/3.0_rp*mu_inf

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex

                 a   = dxi_dx (i,j);  d   = deta_dy(i,j)
                 cxx = dxi_dxx(i,j);  cyy = dxi_dyy(i,j)
                 exx = deta_dxx(i,j); eyy = deta_dyy(i,j)

                 i_st_z     = ztastep_i(k)
                 i_st2_z    = i_st_z*i_st_z
                 iz_z       = iz_zta(k)
                 zeta_z     = i_st_z*iz_z
                 d2zeta_dz2 = - z_zta2(k)*iz_z**3

                 uc = 0.0_rp; ue = 0.0_rp
                 vc = 0.0_rp; ve = 0.0_rp
                 wz = 0.0_rp
                 cl2   = central_2(0)
                 d2ucc = cl2*U(i,j,k)
                 d2vee = cl2*V(i,j,k)
                 d2wzz = cl2*W(i,j,k)
                 !$acc loop seq
                 do l = 1,3
                    cl1 = central_1(l)
                    cl2 = central_2(l)

                    uc = uc + cl1*(U(i+l,j,k) - U(i-l,j,k))
                    ue = ue + cl1*(U(i,j+l,k) - U(i,j-l,k))
                    vc = vc + cl1*(V(i+l,j,k) - V(i-l,j,k))
                    ve = ve + cl1*(V(i,j+l,k) - V(i,j-l,k))
                    wz = wz + cl1*(W(i,j,k+l) - W(i,j,k-l))

                    d2ucc = d2ucc + cl2*(U(i+l,j,k) + U(i-l,j,k))
                    d2vee = d2vee + cl2*(V(i,j+l,k) + V(i,j-l,k))
                    d2wzz = d2wzz + cl2*(W(i,j,k+l) + W(i,j,k-l))
                 enddo

                 d2uce = 0.0_rp; d2vce = 0.0_rp
                 !$acc loop seq
                 do m = 1,3
                  !$acc loop seq
                  do l = 1,3
                     clm = central_1(l)*central_1(m)
                     d2uce = d2uce + clm*(U(i+l,j+m,k)-U(i-l,j+m,k)-U(i+l,j-m,k)+U(i-l,j-m,k))
                     d2vce = d2vce + clm*(V(i+l,j+m,k)-V(i-l,j+m,k)-V(i+l,j-m,k)+V(i-l,j-m,k))
                  enddo
                 enddo

                 wcz = 0.0_rp; wez = 0.0_rp; s_z = 0.0_rp
                 !$acc loop seq
                 do m = 1,3
                  !$acc loop seq
                  do l = 1,3
                     clm = central_1(l)*central_1(m)
                     wcz = wcz + clm*(W(i+l,j,k+m)-W(i-l,j,k+m)-W(i+l,j,k-m)+W(i-l,j,k-m))
                     wez = wez + clm*(W(i,j+l,k+m)-W(i,j-l,k+m)-W(i,j+l,k-m)+W(i,j-l,k-m))
                     s_z = s_z + clm*( a*(U(i+l,j,k+m)-U(i-l,j,k+m)-U(i+l,j,k-m)+U(i-l,j,k-m))   &
                                     + d*(V(i,j+l,k+m)-V(i,j-l,k+m)-V(i,j+l,k-m)+V(i,j-l,k-m)) )
                  enddo
                 enddo

                 u_xx = a*a*d2ucc + cxx*uc + exx*ue
                 u_xy = (a*d)*d2uce
                 v_xy = (a*d)*d2vce
                 v_yy = d*d*d2vee + cyy*vc + eyy*ve

                 w_xz = zeta_z*(a*wcz)
                 w_yz = zeta_z*(d*wez)
                 w_zz = d2wzz*i_st2_z*iz_z**2 + (i_st_z*wz)*d2zeta_dz2

                 graddiv_x = u_xx + v_xy + w_xz
                 graddiv_y = u_xy + v_yy + w_yz
                 graddiv_z = zeta_z*s_z + w_zz

                 mu_      = One3MuInf*VIS(i,j,k)

                 RhU_term = mu_*graddiv_x
                 RhV_term = mu_*graddiv_y
                 RhW_term = mu_*graddiv_z
                 RhE_term = RhU_term*U(i,j,k) + RhV_term*V(i,j,k) + RhW_term*W(i,j,k)

                 RHS(i,j,k,2) = RHS(i,j,k,2) + RhU_term
                 RHS(i,j,k,3) = RHS(i,j,k,3) + RhV_term
                 RHS(i,j,k,4) = RHS(i,j,k,4) + RhW_term
                 RHS(i,j,k,5) = RHS(i,j,k,5) + RhE_term

              enddo
           enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine DivergencyGradient3D_cart

subroutine ComputeSGSIsotropicTerm()

        use fluid_functions_module, only: laminar_viscosity

        implicit none
        real(rp), parameter :: Cprime     = 45.8_rp
        real(rp), parameter :: two3       = 2.0_rp/3.0_rp
        real(rp), parameter :: two3Cprime = two3*Cprime

        real(rp) :: mu_lam, mu_sgs, idelta2
        real(rp) :: csi_x, csi_y, eta_x, eta_y, zeta_z
        real(rp) :: iso_c, iso_e, iso_z
        real(rp) :: isoSGS_x, isoSGS_y, isoSGS_z, isoEnergy
        integer  :: i,j,k,l,m
        logical  :: phys_W, phys_E, phys_S, phys_N, phys_B, phys_F

        call StartProfRange("ComputeSGSIsotropicTerm")

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do       k = sz,ez
           do    j = sy,ey
              do i = sx,ex

                 idelta2 = (Jac2D(i,j)*zstep_i(k))**two3

                 mu_lam = laminar_viscosity(T(i,j,k),tref,vis_flag)
                 mu_sgs = (VIS(i,j,k) - mu_lam)*mu_inf

                 isoSGS(i,j,k) = - two3Cprime / phi(i,j,k,1) * mu_sgs*mu_sgs*idelta2

              enddo
           enddo
        enddo
        !$acc end parallel

        call mpi_share(mpi_comm_cart,type_send_prim,type_recv_prim,my_neighbour,dims, &
                bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                isoSGS)

        phys_W = (my_neighbour(fW) < 0)
        phys_E = (my_neighbour(fE) < 0)
        phys_S = (my_neighbour(fS) < 0)
        phys_N = (my_neighbour(fN) < 0)
        phys_B = (my_neighbour(fB) < 0) .and. (nz > 1)
        phys_F = (my_neighbour(fF) < 0) .and. (nz > 1)

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do       k = sz,ez
           do    j = sy,ey
              do i = sx,ex

                 iso_c = 0.0_rp
                 if      (phys_W .and. i-sx < 3) then
                    !$acc loop seq
                    do m = 0,6
                       iso_c = iso_c + c1_wbias(i-sx,m) * isoSGS(sx+m,j,k)
                    enddo
                 else if (phys_E .and. ex-i < 3) then
                    !$acc loop seq
                    do m = 0,6
                       iso_c = iso_c + c1_ebias(ex-i,m) * isoSGS(ex-6+m,j,k)
                    enddo
                 else
                    !$acc loop seq
                    do l = -3,3
                       iso_c = iso_c + central_1(l) * isoSGS(i+l,j,k)
                    enddo
                 endif

                 iso_e = 0.0_rp
                 if      (phys_S .and. j-sy < 3) then
                    !$acc loop seq
                    do m = 0,6
                       iso_e = iso_e + c1_wbias(j-sy,m) * isoSGS(i,sy+m,k)
                    enddo
                 else if (phys_N .and. ey-j < 3) then
                    !$acc loop seq
                    do m = 0,6
                       iso_e = iso_e + c1_ebias(ey-j,m) * isoSGS(i,ey-6+m,k)
                    enddo
                 else
                    !$acc loop seq
                    do l = -3,3
                       iso_e = iso_e + central_1(l) * isoSGS(i,j+l,k)
                    enddo
                 endif

                 iso_z = 0.0_rp
                 if      (phys_B .and. k-sz < 3) then
                    !$acc loop seq
                    do m = 0,6
                       iso_z = iso_z + c1_wbias(k-sz,m) * isoSGS(i,j,sz+m)
                    enddo
                 else if (phys_F .and. ez-k < 3) then
                    !$acc loop seq
                    do m = 0,6
                       iso_z = iso_z + c1_ebias(ez-k,m) * isoSGS(i,j,ez-6+m)
                    enddo
                 else
                    !$acc loop seq
                    do l = -3,3
                       iso_z = iso_z + central_1(l) * isoSGS(i,j,k+l)
                    enddo
                 endif

                 csi_x  = dxi_dx (i,j)
                 csi_y  = dxi_dy (i,j)
                 eta_x  = deta_dx(i,j)
                 eta_y  = deta_dy(i,j)
                 zeta_z = ztastep_i(k)*iz_zta(k)

                 isoSGS_x = csi_x*iso_c + eta_x*iso_e
                 isoSGS_y = csi_y*iso_c + eta_y*iso_e
                 isoSGS_z = zeta_z*iso_z

                 isoEnergy = U(i,j,k)*isoSGS_x + &
                             V(i,j,k)*isoSGS_y + &
                             W(i,j,k)*isoSGS_z + &
                             isoSGS(i,j,k)*DIV(i,j,k)

                 RHS(i,j,k,2) = RHS(i,j,k,2) + isoSGS_x
                 RHS(i,j,k,3) = RHS(i,j,k,3) + isoSGS_y
                 RHS(i,j,k,4) = RHS(i,j,k,4) + isoSGS_z
                 RHS(i,j,k,5) = RHS(i,j,k,5) + isoEnergy

              enddo
           enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine ComputeSGSIsotropicTerm

subroutine sgs_isotropic_setup()
        implicit none
        integer :: vlo_x, vhi_x, vlo_y, vhi_y, vlo_z, vhi_z
        logical :: phys_W, phys_E, phys_S, phys_N, phys_B, phys_F, too_narrow

        phys_W = (my_neighbour(fW) < 0)
        phys_E = (my_neighbour(fE) < 0)
        phys_S = (my_neighbour(fS) < 0)
        phys_N = (my_neighbour(fN) < 0)
        phys_B = (my_neighbour(fB) < 0) .and. (nz > 1)
        phys_F = (my_neighbour(fF) < 0) .and. (nz > 1)

        vlo_x = merge(sx, lbx, phys_W) ; vhi_x = merge(ex, ubx, phys_E)
        vlo_y = merge(sy, lby, phys_S) ; vhi_y = merge(ey, uby, phys_N)
        vlo_z = merge(sz, lbz, phys_B) ; vhi_z = merge(ez, ubz, phys_F)

        too_narrow = .false.
        if(phys_W) too_narrow = too_narrow .or. (sx+6 > vhi_x)
        if(phys_E) too_narrow = too_narrow .or. (ex-6 < vlo_x)
        if(phys_S) too_narrow = too_narrow .or. (sy+6 > vhi_y)
        if(phys_N) too_narrow = too_narrow .or. (ey-6 < vlo_y)
        if(phys_B) too_narrow = too_narrow .or. (sz+6 > vhi_z)
        if(phys_F) too_narrow = too_narrow .or. (ez-6 < vlo_z)

        if(too_narrow) then
           print '(1x,a,i0,a,i0,a)', 'ERROR [sgs_isotropic_setup]: rank ', rank, &
                ' subdomain too narrow (< 7 nodes) for the order-', central_fd_order, &
                ' progressively-shifted isoSGS stencil at a physical face -- add points or use fewer ranks there.'
           call secure_stop
        endif

        return
end subroutine sgs_isotropic_setup

end module viscous_shared_module
