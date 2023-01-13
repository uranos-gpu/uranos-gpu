module viscous_module
#ifdef TIME
use performance_module
#endif
use profiling_module
use parameters_module
use mpi_module
use mpi_comm_module
use mesh_module
use storage_module

implicit none
private
public viscous_fluxes, viscous_flux_3D_staggered, viscous_flux_3D

contains
subroutine viscous_fluxes(dims)
        
        implicit none
        integer, intent(in) :: dims

#ifdef TIME
        call mpi_stime(s_vis_time)
#endif
        call StartProfRange("viscous_fluxes")
        
        !!$acc data copyin(mid_point_lele_x,mid_point_lele_y,mid_point_lele_z) &
        !!$acc copyin(xstep_i,ystep_i,zstep_i,xsteph,ysteph,zsteph) &
        !!$acc copyin(U,V,W,T,VIS,LMD) &
        !!$acc copyin(central_2_one_half,central_1,central_2) &
        !!$acc copyin(csistep_i,y_eta,y_eta2,x_csi,x_csi2,etastep_i) &
        !!$acc copy(DIV,RHS) 
        selectcase(dims)

          case(2)
            call viscous_flux_2D(U,V,T,VIS,LMD,DIV,RHS)
            call DivergencyGradient2D(U,V,VIS,DIV,RHS)
          case(3)
            call viscous_flux_3D_staggered(U,V,W,T,VIS,LMD,DIV,RHS)
            call DivergencyGradient3D(U,V,W,VIS,DIV,RHS)

        endselect
        !!$acc end data

        call EndProfRange
#ifdef TIME
        call mpi_etime(s_vis_time,t_vis_calls,t_vis_time)
#endif

        return
end subroutine viscous_fluxes


subroutine viscous_flux_3D

        use parameters_module, only: rp, mu_inf
        use mpi_module       , only: sx,ex,sy,ey,sz,ez
        use mesh_module      , only: csistep_i,etastep_i,ztastep_i, x_csi, y_eta, z_zta, x_csi2, y_eta2, z_zta2
        use storage_module   , only: U, V, W, T, VIS, DIV, RHS, LMD, &
                                     central_fd_order, central_1, central_2

        implicit none

        real(rp), dimension(3)   :: div_D, div_sigma, D_dm
        real(rp), dimension(3)   :: du, dv, dw, dT, dm, dl
        real(rp), dimension(3)   :: d2u, d2v, d2w, d2T
        real(rp), dimension(3)   :: i_st, i_st2
        real(rp), dimension(3)   :: ix_xi, d2xi_dx2
        real(rp), dimension(3)   :: vel
        real(rp), dimension(3,3) :: grad_v, grad_vT, div_vI, D
        real(rp)                 :: cl1, cl2
        real(rp), parameter      :: two3 = 2.0_rp/3.0_rp
        real(rp)                 :: mu_, rl_, heat_flux, div_sigma_v
        real(rp)                 :: div_, sum_D_grad
        integer                  :: fL, fR , l, i,j,k

        fL = + central_fd_order/2
        fR = + central_fd_order/2

        !$omp parallel do collapse(3) default(private), &
        !$omp shared(U,V,W,T,VIS,LMD,DIV,RHS,mu_inf,central_1,central_2), &
        !$omp shared(sx,ex,sy,ey,sz,ez,csistep_i,etastep_i,ztastep_i, x_csi, y_eta, z_zta, x_csi2, y_eta2, z_zta2,fL,fR)
        do k = sz,ez
           do j = sy,ey
              do i = sx,ex

                 i_st(1) = csistep_i(i)
                 i_st(2) = etastep_i(j)
                 i_st(3) = ztastep_i(k)

                 i_st2(:) = i_st(:)*i_st(:)
                 
                 ! --- MOMENTUM VISCOUS THERM ---!
                 vel(1)  = U(i,j,k)
                 vel(2)  = V(i,j,k)
                 vel(3)  = W(i,j,k)
                 mu_     = VIS(i,j,k)
                 rl_     = LMD(i,j,k)
                
                 du = 0.0_rp
                 dv = 0.0_rp
                 dw = 0.0_rp
                 dT = 0.0_rp
                 dm = 0.0_rp
                 dl = 0.0_rp
                 do l = 1,fR
                    cl1 = central_1(l)

                    du(1) = du(1) + cl1 * (U(i+l,j,k) - U(i-l,j,k))
                    du(2) = du(2) + cl1 * (U(i,j+l,k) - U(i,j-l,k))
                    du(3) = du(3) + cl1 * (U(i,j,k+l) - U(i,j,k-l))

                    dv(1) = dv(1) + cl1 * (V(i+l,j,k) - V(i-l,j,k))
                    dv(2) = dv(2) + cl1 * (V(i,j+l,k) - V(i,j-l,k))
                    dv(3) = dv(3) + cl1 * (V(i,j,k+l) - V(i,j,k-l))

                    dw(1) = dw(1) + cl1 * (W(i+l,j,k) - W(i-l,j,k))
                    dw(2) = dw(2) + cl1 * (W(i,j+l,k) - W(i,j-l,k))
                    dw(3) = dw(3) + cl1 * (W(i,j,k+l) - W(i,j,k-l))

                    dT(1) = dT(1) + cl1 * (T(i+l,j,k) - T(i-l,j,k))
                    dT(2) = dT(2) + cl1 * (T(i,j+l,k) - T(i,j-l,k))
                    dT(3) = dT(3) + cl1 * (T(i,j,k+l) - T(i,j,k-l))

                    dm(1) = dm(1) + cl1 * (VIS(i+l,j,k) - VIS(i-l,j,k))
                    dm(2) = dm(2) + cl1 * (VIS(i,j+l,k) - VIS(i,j-l,k))
                    dm(3) = dm(3) + cl1 * (VIS(i,j,k+l) - VIS(i,j,k-l))

                    dl(1) = dl(1) + cl1 * (LMD(i+l,j,k) - LMD(i-l,j,k))
                    dl(2) = dl(2) + cl1 * (LMD(i,j+l,k) - LMD(i,j-l,k))
                    dl(3) = dl(3) + cl1 * (LMD(i,j,k+l) - LMD(i,j,k-l))

                 enddo
                 du = i_st*du
                 dv = i_st*dv
                 dw = i_st*dw
                 dT = i_st*dT
                 dm = i_st*dm
                 dl = i_st*dl
                        
                 d2u = 0.0_rp
                 d2v = 0.0_rp
                 d2w = 0.0_rp
                 d2T = 0.0_rp
                 do l = -fL, fR

                    cl2 = central_2(l)

                    d2u(1) = d2u(1) + cl2 * U(i+l,j,k)
                    d2u(2) = d2u(2) + cl2 * U(i,j+l,k)
                    d2u(3) = d2u(3) + cl2 * U(i,j,k+l)

                    d2v(1) = d2v(1) + cl2 * V(i+l,j,k)
                    d2v(2) = d2v(2) + cl2 * V(i,j+l,k)
                    d2v(3) = d2v(3) + cl2 * V(i,j,k+l)

                    d2w(1) = d2w(1) + cl2 * W(i+l,j,k)
                    d2w(2) = d2w(2) + cl2 * W(i,j+l,k)
                    d2w(3) = d2w(3) + cl2 * W(i,j,k+l)

                    d2T(1) = d2T(1) + cl2 * T(i+l,j,k)
                    d2T(2) = d2T(2) + cl2 * T(i,j+l,k)
                    d2T(3) = d2T(3) + cl2 * T(i,j,k+l)

                 enddo
                 d2u = i_st2*d2u
                 d2v = i_st2*d2v
                 d2w = i_st2*d2w
                 d2T = i_st2*d2T
        
                 ! metrics
                 ix_xi(1) = 1.0_rp/(x_csi(i))
                 ix_xi(2) = 1.0_rp/(y_eta(j))
                 ix_xi(3) = 1.0_rp/(z_zta(k))

                 d2xi_dx2(1) = - x_csi2(i)*ix_xi(1)**3
                 d2xi_dx2(2) = - y_eta2(j)*ix_xi(2)**3
                 d2xi_dx2(3) = - z_zta2(k)*ix_xi(3)**3

                 ! compute second derivative in physical space
                 d2u(:) = d2u(:)*(ix_xi)**2 + du(:) * d2xi_dx2(:)
                 d2v(:) = d2v(:)*(ix_xi)**2 + dv(:) * d2xi_dx2(:)
                 d2w(:) = d2w(:)*(ix_xi)**2 + dw(:) * d2xi_dx2(:)
                 d2T(:) = d2T(:)*(ix_xi)**2 + dT(:) * d2xi_dx2(:)

                 ! compute derivative in physical space
                 du = du*ix_xi
                 dv = dv*ix_xi
                 dw = dw*ix_xi
                 dT = dT*ix_xi
                 dm = dm*ix_xi
                 dl = dl*ix_xi

                 ! calculating D
                 grad_v(1,:) = du(:)
                 grad_v(2,:) = dv(:)
                 grad_v(3,:) = dw(:)

                 grad_vT(:,1) = du(:)
                 grad_vT(:,2) = dv(:)
                 grad_vT(:,3) = dw(:)
                
                 div_ = du(1) + dv(2) + dw(3)
                 div_vI = 0.0_rp
                 div_vI(1,1) = div_
                 div_vI(2,2) = div_
                 div_vI(3,3) = div_
                 DIV(i,j,k)  = div_

                 D = grad_v + grad_vT - two3 * div_vI

                 ! calculating div(D)
                 div_D(1) = d2u(1) + d2u(2) + d2u(3)
                 div_D(2) = d2v(1) + d2v(2) + d2v(3)
                 div_D(3) = d2w(1) + d2w(2) + d2w(3)

                 D_dm(1) = D(1,1)*dm(1) + D(1,2)*dm(2) + D(1,3)*dm(3)
                 D_dm(2) = D(2,1)*dm(1) + D(2,2)*dm(2) + D(2,3)*dm(3)
                 D_dm(3) = D(3,1)*dm(1) + D(3,2)*dm(2) + D(3,3)*dm(3)

                 ! calculating div(sigma)
                 div_sigma = mu_ * div_D + D_dm

                 ! --- ENERGY VISCOUS THERM --- !
                 sum_D_grad = D(1,1)*grad_v(1,1) + D(1,2)*grad_v(1,2) + D(1,3)*grad_v(1,3) + &
                              D(2,1)*grad_v(2,1) + D(2,2)*grad_v(2,2) + D(2,3)*grad_v(2,3) + &
                              D(3,1)*grad_v(3,1) + D(3,2)*grad_v(3,2) + D(3,3)*grad_v(3,3)

                 div_sigma_v = div_sigma(1) * vel(1) + div_sigma(2) * vel(2) + div_sigma(3) * vel(3) + mu_*sum_D_grad
                 
                 heat_flux = rl_ * (d2T(1) + d2T(2) + d2T(3)) + dl(1)*dT(1) + dl(2)*dT(2) + dl(3)*dT(3)

                 RHS(i,j,k,2) = RHS(i,j,k,2) + mu_inf * (div_sigma(1))
                 RHS(i,j,k,3) = RHS(i,j,k,3) + mu_inf * (div_sigma(2))
                 RHS(i,j,k,4) = RHS(i,j,k,4) + mu_inf * (div_sigma(3))
                 RHS(i,j,k,5) = RHS(i,j,k,5) + mu_inf * (div_sigma_v + heat_flux)  
                
              enddo ! i
           enddo ! j
        enddo ! k
        !$omp end parallel do

        return
end subroutine viscous_flux_3D

subroutine viscous_flux_2D(U,V,T,VIS,LMD,DIV,RHS)

        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: RHS
        real(rp), allocatable, dimension(:,:,:)  , intent(inout) :: DIV
        real(rp), allocatable, dimension(:,:,:)  , intent(in)    :: U,V,T
        real(rp), allocatable, dimension(:,:,:)  , intent(in)    :: VIS,LMD

        real(rp)            :: grad_v_xx, grad_v_xy
        real(rp)            :: grad_v_yx, grad_v_yy
        real(rp)            :: grad_vT_xx, grad_vT_xy
        real(rp)            :: grad_vT_yx, grad_vT_yy
        real(rp)            :: div_vI_xx
        real(rp)            :: div_vI_yy
        real(rp)            :: D_xx, D_xy
        real(rp)            :: D_yx, D_yy

        real(rp)            :: div_D_x, div_D_y
        real(rp)            :: div_sigma_x, div_sigma_y
        real(rp)            :: i_st_x, i_st_y
        real(rp)            :: i_st2_x, i_st2_y
        real(rp)            :: vel_x, vel_y
        real(rp)            :: D_dm_x, D_dm_y
        real(rp)            :: cl1_x, cl1_y
        real(rp)            :: du_x, du_y
        real(rp)            :: dv_x, dv_y
        real(rp)            :: dt_x, dt_y
        real(rp)            :: dm_x, dm_y
        real(rp)            :: dl_x, dl_y
        real(rp)            :: d2u_x, d2u_y
        real(rp)            :: d2v_x, d2v_y
        real(rp)            :: d2t_x, d2t_y
        real(rp)            :: ix_xi_x, ix_xi_y
        real(rp)            :: d2xi_dx2_x, d2xi_dx2_y
        real(rp)            :: cl2
        real(rp)            :: mu_, rl_, heat_flux, div_sigma_v, sum_D_grad

        real(rp), parameter :: two3 = 2.0_rp/3.0_rp
        integer             :: i, j, k, fL, fR , l

        call StartProfRange("viscous_flux_2D") 

        fL = central_fd_order/2
        fR = central_fd_order/2

        k = 1

        !$omp parallel do collapse(2) default(private), &
        !$omp shared(U,V,T,VIS,LMD,DIV,RHS,mu_inf,k,central_1,central_2), &
        !$omp shared(sx,ex,sy,ey,csistep_i,etastep_i,x_csi,y_eta,x_csi2,y_eta2,fL,fR)

        !$acc parallel default(present)
        !$acc loop gang,vector collapse(2)
        do j = sy,ey
           do i = sx,ex

              ! grid steps
              i_st_x = csistep_i(i)
              i_st_y = etastep_i(j)

              i_st2_x = i_st_x*i_st_x
              i_st2_y = i_st_y*i_st_y
              
              vel_x =   U(i,j,k)
              vel_y =   V(i,j,k)
              mu_   = VIS(i,j,k)
              rl_   = LMD(i,j,k)
                
              ! === first derivatives
              du_x = 0.0_rp
              du_y = 0.0_rp
              !
              dv_x = 0.0_rp
              dv_y = 0.0_rp
              !
              dt_x = 0.0_rp
              dt_y = 0.0_rp
              !
              dm_x = 0.0_rp
              dm_y = 0.0_rp
              !
              dl_x = 0.0_rp
              dl_y = 0.0_rp

              !$acc loop seq
              do l = 1,fR

                 cl1_x = i_st_x * central_1(l)
                 cl1_y = i_st_y * central_1(l)

                 du_x = du_x + cl1_x * (U(i+l,j,k) - U(i-l,j,k))
                 du_y = du_y + cl1_y * (U(i,j+l,k) - U(i,j-l,k))

                 dv_x = dv_x + cl1_x * (V(i+l,j,k) - V(i-l,j,k))
                 dv_y = dv_y + cl1_y * (V(i,j+l,k) - V(i,j-l,k))

                 dt_x = dt_x + cl1_x * (T(i+l,j,k) - T(i-l,j,k))
                 dt_y = dt_y + cl1_y * (T(i,j+l,k) - T(i,j-l,k))

                 dm_x = dm_x + cl1_x * (VIS(i+l,j,k) - VIS(i-l,j,k))
                 dm_y = dm_y + cl1_y * (VIS(i,j+l,k) - VIS(i,j-l,k))

                 dl_x = dl_x + cl1_x * (LMD(i+l,j,k) - LMD(i-l,j,k))
                 dl_y = dl_y + cl1_y * (LMD(i,j+l,k) - LMD(i,j-l,k))

              enddo

              ! === second derivatives
              d2u_x = 0.0_rp
              d2u_y = 0.0_rp
              !
              d2v_x = 0.0_rp
              d2v_y = 0.0_rp
              !
              d2t_x = 0.0_rp
              d2t_y = 0.0_rp

              !$acc loop seq
              do l = -fL, fR

                 cl2 = central_2(l)

                 d2u_x = d2u_x + cl2 * U(i+l,j,k)
                 d2u_y = d2u_y + cl2 * U(i,j+l,k)

                 d2v_x = d2v_x + cl2 * V(i+l,j,k)
                 d2v_y = d2v_y + cl2 * V(i,j+l,k)

                 d2t_x = d2t_x + cl2 * T(i+l,j,k)
                 d2t_y = d2t_y + cl2 * T(i,j+l,k)

              enddo

              d2u_x = i_st2_x*d2u_x
              d2u_y = i_st2_y*d2u_y
              d2v_x = i_st2_x*d2v_x
              d2v_y = i_st2_y*d2v_y
              d2t_x = i_st2_x*d2t_x
              d2t_y = i_st2_y*d2t_y

              ! metrics
              ix_xi_x = 1.0_rp/(x_csi(i))
              ix_xi_y = 1.0_rp/(y_eta(j))

              d2xi_dx2_x = - x_csi2(i)*ix_xi_x**3
              d2xi_dx2_y = - y_eta2(j)*ix_xi_y**3

              ! compute second derivative in physical space
              d2u_x = d2u_x*(ix_xi_x)**2 + du_x * d2xi_dx2_x
              d2u_y = d2u_y*(ix_xi_y)**2 + du_y * d2xi_dx2_y
              !
              d2v_x = d2v_x*(ix_xi_x)**2 + dv_x * d2xi_dx2_x
              d2v_y = d2v_y*(ix_xi_y)**2 + dv_y * d2xi_dx2_y
              !
              d2t_x = d2t_x*(ix_xi_x)**2 + dt_x * d2xi_dx2_x
              d2t_y = d2t_y*(ix_xi_y)**2 + dt_y * d2xi_dx2_y
                
              ! compute derivative in physical space
              du_x = du_x*ix_xi_x
              du_y = du_y*ix_xi_y
              !
              dv_x = dv_x*ix_xi_x
              dv_y = dv_y*ix_xi_y
              !
              dt_x = dt_x*ix_xi_x
              dt_y = dt_y*ix_xi_y
              !
              dm_x = dm_x*ix_xi_x
              dm_y = dm_y*ix_xi_y
              !
              dl_x = dl_x*ix_xi_x
              dl_y = dl_y*ix_xi_y


              ! --- MOMENTUM VISCOUS THERM ---!

              ! calculating D
              grad_v_xx  = du_x
              grad_v_xy  = du_y
              !
              grad_v_yx  = dv_x
              grad_v_yy  = dv_y

              grad_vT_xx = du_x
              grad_vT_yx = du_y
              !
              grad_vT_xy = dv_x
              grad_vT_yy = dv_y
              !
              div_vI_xx = du_x + dv_y
              div_vI_yy = du_x + dv_y

              DIV(i,j,k) = du_x + dv_y

              D_xx = grad_v_xx + grad_vT_xx - two3 * div_vI_xx
              D_xy = grad_v_xy + grad_vT_xy
              D_yx = grad_v_yx + grad_vT_yx
              D_yy = grad_v_yy + grad_vT_yy - two3 * div_vI_yy

              ! calculating div(D)
              div_D_x = d2u_x + d2u_y
              div_D_y = d2v_x + d2v_y
              
              ! calculating div(sigma)
              D_dm_x = D_xx*dm_x + D_xy*dm_y
              D_dm_y = D_yx*dm_x + D_yy*dm_y

              div_sigma_x = mu_ * div_D_x + D_dm_x
              div_sigma_y = mu_ * div_D_y + D_dm_y

              ! --- ENERGY VISCOUS THERM --- !
              sum_D_grad = D_xx*grad_v_xx + D_xy*grad_v_xy + &
                           D_yx*grad_v_yx + D_yy*grad_v_yy

              div_sigma_v = div_sigma_x*vel_x + div_sigma_y*vel_y + mu_*sum_D_grad
              
              heat_flux = rl_* (d2t_x + d2t_y) + dl_x*dt_x + dl_y*dt_y

              RHS(i,j,k,2) = RHS(i,j,k,2) + mu_inf * (div_sigma_x)
              RHS(i,j,k,3) = RHS(i,j,k,3) + mu_inf * (div_sigma_y)
              RHS(i,j,k,5) = RHS(i,j,k,5) + mu_inf * (div_sigma_v + heat_flux)  

           enddo ! i
        enddo ! j
        !$acc end parallel

        !$omp end parallel do

        call EndProfRange 

        return
end subroutine viscous_flux_2D

subroutine DivergencyGradient2D(U,V,VIS,DIV,RHS)
      
        implicit none
        real(rp), dimension(:,:,:)  , allocatable, intent(in)    :: U,V,VIS
        real(rp), dimension(:,:,:)  , allocatable, intent(inout) :: DIV
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: RHS

        ! local declarations
        real(rp) :: One3MuInf
        real(rp) :: GradDiv_x, GradDiv_y
        real(rp) :: cl1, mu_, rhu_term, rhv_term, rhe_term
        integer  :: i,j,k,l,fR

        call StartProfRange("DivergencyGradient2D")
        
        call mpi_share(mpi_comm_cart,type_send_prim,type_recv_prim,my_neighbour,dims, &
                bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                DIV)

        fR = central_fd_order/2
        One3MuInf = 1.0_rp/3.0_rp*mu_inf

        k = 1
        
        !$acc parallel default(present)
        !$acc loop gang, vector collapse(2)
        do j    = sy,ey
           do i = sx,ex

              GradDiv_x = 0.0_rp
              GradDiv_y = 0.0_rp
              !$acc loop seq
              do l = 1,fR
                 cl1 = central_1(l)

                 GradDiv_x = GradDiv_x + cl1 * (DIV(i+l,j,k) - DIV(i-l,j,k))
                 GradDiv_y = GradDiv_y + cl1 * (DIV(i,j+l,k) - DIV(i,j-l,k))

              enddo

              mu_      = One3MuInf*VIS(i,j,k)
              rhu_term = mu_*xstep_i(i)*GradDiv_x
              rhv_term = mu_*ystep_i(j)*GradDiv_y
              rhe_term = rhu_term*U(i,j,k) + rhv_term*V(i,j,k)

              RHS(i,j,k,2) = RHS(i,j,k,2) + RhU_term
              RHS(i,j,k,3) = RHS(i,j,k,3) + RhV_term
              RHS(i,j,k,5) = RHS(i,j,k,5) + RhE_term
        
           enddo
        enddo
        !$acc end parallel 

        call EndProfRange
        
        return
end subroutine DivergencyGradient2D


subroutine DivergencyGradient3D(U,V,W,VIS,DIV,RHS)
      
        implicit none
        real(rp), dimension(:,:,:)  , allocatable, intent(in)    :: U,V,W,VIS
        real(rp), dimension(:,:,:)  , allocatable, intent(inout) :: DIV
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: RHS

        ! local declarations
        real(rp) :: One3MuInf
        real(rp) :: GradDiv_x, GradDiv_y, GradDiv_z, cl1, mu_
        real(rp) :: RhU_term, RhV_term, RhW_term, rhe_term
        integer  :: i,j,k,l,fR

        call StartProfRange("DivergencyGradient3D_MPI") 

        call mpi_share(mpi_comm_cart,type_send_prim,type_recv_prim,my_neighbour,dims, &
                bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                DIV)

        call EndProfRange 
        call StartProfRange("DivergencyGradient3D") 

        fR = central_fd_order/2
        One3MuInf = 1.0_rp/3.0_rp*mu_inf

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex

                 GradDiv_x = 0.0_rp
                 GradDiv_y = 0.0_rp
                 GradDiv_z = 0.0_rp

                 !$acc loop seq
                 do l = 1,fR
                    cl1 = central_1(l)

                    GradDiv_x = GradDiv_x + cl1 * (DIV(i+l,j,k) - DIV(i-l,j,k))
                    GradDiv_y = GradDiv_y + cl1 * (DIV(i,j+l,k) - DIV(i,j-l,k))
                    GradDiv_z = GradDiv_z + cl1 * (DIV(i,j,k+l) - DIV(i,j,k-l))

                 enddo

                 mu_      = One3MuInf*VIS(i,j,k)

                 RhU_term = mu_*xstep_i(i)*GradDiv_x
                 RhV_term = mu_*ystep_i(j)*GradDiv_y
                 RhW_term = mu_*zstep_i(k)*GradDiv_z
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
end subroutine DivergencyGradient3D




subroutine viscous_flux_3D_staggered(U,V,W,T,VIS,LMD,DIV,RHS)

        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: RHS
        real(rp), allocatable, dimension(:,:,:)  , intent(inout) :: DIV
        real(rp), allocatable, dimension(:,:,:)  , intent(in)    :: U, V, W, T
        real(rp), allocatable, dimension(:,:,:)  , intent(in)    :: VIS, LMD

        real(rp), parameter :: two3 = 2.0_rp/3.0_rp

        real(rp)               :: DC_xx, DC_xy, DC_xz
        real(rp)               :: DC_yx, DC_yy, DC_yz
        real(rp)               :: DC_zx, DC_zy, DC_zz
        real(rp)               :: grad_v_xx, grad_v_xy, grad_v_xz
        real(rp)               :: grad_v_yx, grad_v_yy, grad_v_yz
        real(rp)               :: grad_v_zx, grad_v_zy, grad_v_zz
        real(rp)               :: grad_vT_xx, grad_vT_xy, grad_vT_xz
        real(rp)               :: grad_vT_yx, grad_vT_yy, grad_vT_yz
        real(rp)               :: grad_vT_zx, grad_vT_zy, grad_vT_zz

        real(rp)               :: DC_dm_x, DC_dm_y, DC_dm_z
        real(rp)               :: div_sigma_x, div_sigma_y, div_sigma_z
        real(rp)               :: Lapl_x, Lapl_y, Lapl_z
        real(rp)               :: d_mdu_x, d_mdu_y, d_mdu_z
        real(rp)               :: d_mdv_x, d_mdv_y, d_mdv_z
        real(rp)               :: d_mdw_x, d_mdw_y, d_mdw_z
        real(rp)               :: d_ldT_x, d_ldT_y, d_ldT_z
        real(rp)               :: du_x, du_y, du_z
        real(rp)               :: dv_x, dv_y, dv_z
        real(rp)               :: dw_x, dw_y, dw_z
        real(rp)               :: dm_x, dm_y, dm_z
        real(rp)               :: c2R_x, c2R_y, c2R_z
        real(rp)               :: c2L_x, c2L_y, c2L_z
        real(rp)               :: vel_x, vel_y, vel_z
        real(rp)               :: i_st_x, i_st_y, i_st_z
        real(rp)               :: i_stR_x, i_stR_y, i_stR_z
        real(rp)               :: i_stL_x, i_stL_y, i_stL_z
        real(rp)               :: mu_, muN, muS, muE, muW, muB, muF
        real(rp)               :: rl_, rlN, rlS, rlE, rlW, rlB, rlF
        real(rp)               :: cl1, Lapl_v, muGradVGradV, DC_dm_v, mu_DC_GradV
        real(rp)               :: div_, div_sigma_v, heat_flux
        real(rp)               :: clxR, clxL, clyR, clyL, clzR, clzL, cl
        integer :: i,j,k,l,fR

        call StartProfRange("viscous_flux_3D_staggered") 

        fR = central_fd_order/2

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do   k = sz,ez
         do  j = sy,ey
          do i = sx,ex

             ! velocity
             vel_x  = U(i,j,k)
             vel_y  = V(i,j,k)
             vel_z  = W(i,j,k)
            
             ! viscosity and diffusivity at cell board
             mu_ = VIS(i,j,k)
             muN = 0.0_rp
             muS = 0.0_rp
             muE = 0.0_rp
             muW = 0.0_rp
             muF = 0.0_rp
             muB = 0.0_rp
             rl_ = LMD(i,j,k)
             rlN = 0.0_rp
             rlS = 0.0_rp
             rlE = 0.0_rp
             rlW = 0.0_rp
             rlF = 0.0_rp
             rlB = 0.0_rp

             !$acc loop seq 
             do l = fd_L,fd_R

               clxR = mid_point_lele_x(l,i)
               clyR = mid_point_lele_y(l,j)
               clzR = mid_point_lele_z(l,k)

               clxL = mid_point_lele_x(l,i-1)
               clyL = mid_point_lele_y(l,j-1)
               clzL = mid_point_lele_z(l,k-1)

               muE = muE + clxR * VIS(i+l,j,k)
               muN = muN + clyR * VIS(i,j+l,k)
               muF = muF + clzR * VIS(i,j,k+l)

               muW = muW + clxL * VIS(i-1+l,j,k)
               muS = muS + clyL * VIS(i,j-1+l,k)
               muB = muB + clzL * VIS(i,j,k-1+l)

               rlE = rlE + clxR * LMD(i+l,j,k)
               rlN = rlN + clyR * LMD(i,j+l,k)
               rlF = rlF + clzR * LMD(i,j,k+l)

               rlW = rlW + clxL * LMD(i-1+l,j,k)
               rlS = rlS + clyL * LMD(i,j-1+l,k)
               rlB = rlB + clzL * LMD(i,j,k-1+l)

             enddo

             ! === 1st DERIVATIVE
             i_st_x = xstep_i(i)
             i_st_y = ystep_i(j)
             i_st_z = zstep_i(k)

             du_x = 0.0_rp
             du_y = 0.0_rp
             du_z = 0.0_rp
             !
             dv_x = 0.0_rp
             dv_y = 0.0_rp
             dv_z = 0.0_rp
             !
             dw_x = 0.0_rp
             dw_y = 0.0_rp
             dw_z = 0.0_rp
             !
             dm_x = 0.0_rp
             dm_y = 0.0_rp
             dm_z = 0.0_rp

             !$acc loop seq 
             do l = 1,fR
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

             du_x = i_st_x*du_x
             dv_x = i_st_x*dv_x
             dw_x = i_st_x*dw_x
             dm_x = i_st_x*dm_x
             !
             du_y = i_st_y*du_y
             dv_y = i_st_y*dv_y
             dw_y = i_st_y*dw_y
             dm_y = i_st_y*dm_y
             !
             du_z = i_st_z*du_z
             dv_z = i_st_z*dv_z
             dw_z = i_st_z*dw_z
             dm_z = i_st_z*dm_z

             ! grad V transpose
             grad_v_xx = du_x
             grad_v_xy = du_y
             grad_v_xz = du_z
             !
             grad_v_yx = dv_x
             grad_v_yy = dv_y
             grad_v_yz = dv_z
             !
             grad_v_zx = dw_x
             grad_v_zy = dw_y
             grad_v_zz = dw_z
        
             ! grad V transpose
             grad_vT_xx = du_x
             grad_vT_yx = du_y
             grad_vT_zx = du_z
             !
             grad_vT_xy = dv_x
             grad_vT_yy = dv_y
             grad_vT_zy = dv_z
             !
             grad_vT_xz = dw_x
             grad_vT_yz = dw_y
             grad_vT_zz = dw_z
        
             ! speed divergency
             div_ = du_x + dv_y + dw_z
             DIV(i,j,k)  = div_

             DC_xx = grad_vT_xx - two3 * div_
             DC_xy = grad_vT_xy 
             DC_xz = grad_vT_xz 
             !
             DC_yx = grad_vT_yx 
             DC_yy = grad_vT_yy - two3 * div_
             DC_yz = grad_vT_yz 
             !
             DC_zx = grad_vT_zx
             DC_zy = grad_vT_zy
             DC_zz = grad_vT_zz - two3 * div_
             !
             ! === 2nd DERIVATIVE
             !
             i_stR_x = 1.0_rp/xsteph(i)
             i_stL_x = 1.0_rp/xsteph(i-1)

             i_stR_y = 1.0_rp/ysteph(j)
             i_stL_y = 1.0_rp/ysteph(j-1)

             i_stR_z = 1.0_rp/zsteph(k)
             i_stL_z = 1.0_rp/zsteph(k-1)

             d_mdu_x = 0.0_rp
             d_mdu_y = 0.0_rp
             d_mdu_z = 0.0_rp
             !
             d_mdv_x = 0.0_rp
             d_mdv_y = 0.0_rp
             d_mdv_z = 0.0_rp
             !
             d_mdw_x = 0.0_rp
             d_mdw_y = 0.0_rp
             d_mdw_z = 0.0_rp
             !
             d_ldT_x = 0.0_rp
             d_ldT_y = 0.0_rp
             d_ldT_z = 0.0_rp

             !$acc loop seq 
             do l = fd_L, fd_R

               cl    = central_2_one_half(l) 
               c2R_x = i_stR_x * cl
               c2L_x = i_stL_x * cl
               !
               c2R_y = i_stR_y * cl
               c2L_y = i_stL_y * cl
               !
               c2R_z = i_stR_z * cl
               c2L_z = i_stL_z * cl
               
               d_mdu_x=d_mdu_x + c2R_x*muE*U(i+l,j,k) - c2L_x*muW*U(i-1+l,j,k)
               d_mdu_y=d_mdu_y + c2R_y*muN*U(i,j+l,k) - c2L_y*muS*U(i,j-1+l,k)
               d_mdu_z=d_mdu_z + c2R_z*muF*U(i,j,k+l) - c2L_z*muB*U(i,j,k-1+l)

               d_mdv_x=d_mdv_x + c2R_x*muE*V(i+l,j,k) - c2L_x*muW*V(i-1+l,j,k)
               d_mdv_y=d_mdv_y + c2R_y*muN*V(i,j+l,k) - c2L_y*muS*V(i,j-1+l,k)
               d_mdv_z=d_mdv_z + c2R_z*muF*V(i,j,k+l) - c2L_z*muB*V(i,j,k-1+l)

               d_mdw_x=d_mdw_x + c2R_x*muE*W(i+l,j,k) - c2L_x*muW*W(i-1+l,j,k)
               d_mdw_y=d_mdw_y + c2R_y*muN*W(i,j+l,k) - c2L_y*muS*W(i,j-1+l,k)
               d_mdw_z=d_mdw_z + c2R_z*muF*W(i,j,k+l) - c2L_z*muB*W(i,j,k-1+l)

               d_lDT_x=d_lDT_x + c2R_x*rlE*T(i+l,j,k) - c2L_x*rlW*T(i-1+l,j,k)
               d_lDT_y=d_lDT_y + c2R_y*rlN*T(i,j+l,k) - c2L_y*rlS*T(i,j-1+l,k)
               d_lDT_z=d_lDT_z + c2R_z*rlF*T(i,j,k+l) - c2L_z*rlB*T(i,j,k-1+l)

             enddo

             d_mdu_x = i_st_x*d_mdu_x
             d_mdv_x = i_st_x*d_mdv_x
             d_mdw_x = i_st_x*d_mdw_x
             d_lDT_x = i_st_x*d_lDT_x
             !
             d_mdu_y = i_st_y*d_mdu_y
             d_mdv_y = i_st_y*d_mdv_y
             d_mdw_y = i_st_y*d_mdw_y
             d_lDT_y = i_st_y*d_lDT_y
             !
             d_mdu_z = i_st_z*d_mdu_z
             d_mdv_z = i_st_z*d_mdv_z
             d_mdw_z = i_st_z*d_mdw_z
             d_lDT_z = i_st_z*d_lDT_z

             !
             ! === MOMENTUM COMPONENTS
             !
             ! 1) u Laplacian
             lapl_x = d_mdu_x + d_mdu_y + d_mdu_z
             lapl_y = d_mdv_x + d_mdv_y + d_mdv_z
             lapl_z = d_mdw_x + d_mdw_y + d_mdw_z
        
             ! 2) DCompressible x Gradient(mu)
             DC_dm_x = DC_xx*dm_x + DC_xy*dm_y + DC_xz*dm_z
             DC_dm_y = DC_yx*dm_x + DC_yy*dm_y + DC_yz*dm_z
             DC_dm_z = DC_zx*dm_x + DC_zy*dm_y + DC_zz*dm_z

             div_sigma_x = lapl_x + DC_dm_x
             div_sigma_y = lapl_y + DC_dm_y
             div_sigma_z = lapl_z + DC_dm_z

             !
             ! === ENERGY COMPONENTS
             !
             ! 1) Laplacian(u) x Velocity
             Lapl_v = lapl_x*vel_x + lapl_y*vel_y+lapl_z*vel_z

             ! 2) mu Gradient(V) x Gradient(V)
             muGradVGradV = &
             mu_*(Grad_V_xx**2 + Grad_V_xy**2 + Grad_V_xz**2 + &
                  Grad_V_yx**2 + Grad_V_yy**2 + Grad_V_yz**2 + &
                  Grad_V_zx**2 + Grad_V_zy**2 + Grad_V_zz**2)

             ! 3) Gradient(mu) x DC x Velocity
             DC_dm_v = DC_dm_x*vel_x + DC_dm_y*vel_y + DC_dm_z*vel_z

             ! 4) mu DC x Gradient(V)
             mu_DC_GradV = &
             mu_*(DC_xx*grad_v_xx + DC_xy*grad_v_xy + DC_xz*grad_v_xz + &
                  DC_yx*grad_v_yx + DC_yy*grad_v_yy + DC_yz*grad_v_yz + &
                  DC_zx*grad_v_zx + DC_zy*grad_v_zy + DC_zz*grad_v_zz)

             div_sigma_v = Lapl_v + muGradVGradV + DC_dm_v + mu_DC_GradV       
             heat_flux   = d_lDT_x + d_lDT_y + d_lDT_z

             RHS(i,j,k,2) = RHS(i,j,k,2) + mu_inf * (div_sigma_x)
             RHS(i,j,k,3) = RHS(i,j,k,3) + mu_inf * (div_sigma_y)
             RHS(i,j,k,4) = RHS(i,j,k,4) + mu_inf * (div_sigma_z)
             RHS(i,j,k,5) = RHS(i,j,k,5) + mu_inf * (div_sigma_v + heat_flux)  
        
          enddo
         enddo
        enddo
        !$acc end parallel

        call EndProfRange 

        return
end subroutine viscous_flux_3D_staggered





end module viscous_module
