module sgs_module
use parameters_module, only: rp
use mpi_module
use mpi_comm_module
use profiling_module

implicit none
private
real(rp), parameter :: Pr_turb   = 0.9_rp
real(rp), parameter :: i_Pr_turb = 1.0_rp/Pr_turb

public compute_subgrid_model

contains
subroutine compute_subgrid_model
! -------------------------------------------------------------------------
!
!       Selection of the subgrid model for LES computation
!
!       1. 'Smagorinsky' classical Smagorinsky model
!       2. 'WALE'        Wall adaptive Local Eddy-viscosity model
!       
! -------------------------------------------------------------------------
        use parameters_module, only: sgs_model, dims
        use storage_module   , only: phi, VIS, LMD

        implicit none

        call StartProfRange("compute_subgrid_model")
        call StartProfRange("sgs_model")

        if(dims /= 3) stop ' LES requires 3D!'
        !
        ! === selection of th e subgrid scale model
        !
        selectcase(trim(sgs_model))

          case('Smagorinsky')
            call compute_smagorinsky(VIS,LMD)

          case('DynamicalSmagorinsky')
            call compute_dynamical_smagorinsky(VIS,LMD)

          case('WALE')
            call compute_wale(VIS,LMD)

          case('sigma_model')
            call compute_sigma_model(VIS,LMD)

          case('MixedTimeScale')
            call compute_MixedTimeScale(VIS,LMD)

          case default
            write(*,'(a)') ' FATAL ERROR!'
            write(*,'(a)') trim(sgs_model), ' is not implemented. The program will be stopped'
            call secure_stop

        endselect
        call EndProfRange

        !
        ! === commucate halos via MPI
        !
        call StartProfRange("sgs_mpi_calls")
        if(mpi_flag) then
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
        !
        ! === apply physical boundary conditions to turbulent viscosity
        !
        call TrbVisBC(phi,VIS,LMD)

        call EndProfRange

        return
end subroutine compute_subgrid_model


subroutine compute_smagorinsky(VIS,LMD)
        
        use parameters_module, only: cp, central_1, central_fd_order, mu_inf
        use mpi_module       , only: sx,ex, sy,ey, sz,ez
        use storage_module   , only: phi, U, V, W
        use mesh_module      , only: xstep, ystep, zstep, xstep_i, ystep_i, zstep_i
        use bc_module, only: face_type

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: VIS
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: LMD

        ! local declarations
        real(rp)                 :: sqrtS_S_
        real(rp)                 :: i_st_x, i_st_y, i_st_z
        real(rp)                 :: cl1_x, cl1_y, cl1_z
        real(rp)                 :: du_x, du_y, du_z
        real(rp)                 :: dv_x, dv_y, dv_z
        real(rp)                 :: dw_x, dw_y, dw_z
        real(rp)                 :: S_xx, S_xy, S_xz
        real(rp)                 :: S_yx, S_yy, S_yz
        real(rp)                 :: S_zx, S_zy, S_zz
        real(rp), parameter      :: two3 = 2.0_rp/3.0_rp, Cs2 = 0.12_rp**2
        real(rp)                 :: delta2
        real(rp)                 :: mu_turb, lb_turb, i_mu_inf
        integer                  :: fR
        integer                  :: i,j,k,l

        fR       = central_fd_order/2
        i_mu_inf = 1.0_rp/mu_inf

        !$omp parallel do collapse(3), default(private), &
        !$omp shared(sx,sy,sz,ex,ey,ez,xstep,ystep,zstep,xstep_i,ystep_i,zstep_i),&
        !$omp shared(phi,U,V,W,central_1,VIS,LMD,i_mu_inf,fR)
        
        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex
                        
                 ! === compute grid effective step
                 delta2 = (xstep(i)*ystep(j)*zstep(k))**two3

                 ! === compute derivatives
                 i_st_x = xstep_i(i)
                 i_st_y = ystep_i(j)
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
                 
                 !$acc loop seq
                 do l = 1,fR

                    cl1_x = i_st_x * central_1(l)
                    cl1_y = i_st_y * central_1(l)
                    cl1_z = i_st_z * central_1(l)

                    du_x = du_x + cl1_x * (U(i+l,j,k) - U(i-l,j,k))
                    du_y = du_y + cl1_y * (U(i,j+l,k) - U(i,j-l,k))
                    du_z = du_z + cl1_z * (U(i,j,k+l) - U(i,j,k-l))

                    dv_x = dv_x + cl1_x * (V(i+l,j,k) - V(i-l,j,k))
                    dv_y = dv_y + cl1_y * (V(i,j+l,k) - V(i,j-l,k))
                    dv_z = dv_z + cl1_z * (V(i,j,k+l) - V(i,j,k-l))

                    dw_x = dw_x + cl1_x * (W(i+l,j,k) - W(i-l,j,k))
                    dw_y = dw_y + cl1_y * (W(i,j+l,k) - W(i,j-l,k))
                    dw_z = dw_z + cl1_z * (W(i,j,k+l) - W(i,j,k-l))

                 enddo

                 ! === compute strain tensor
                 S_xx = 0.5_rp*(du_x + du_x)
                 S_xy = 0.5_rp*(du_y + dv_x)
                 S_xz = 0.5_rp*(du_z + dw_x)
                 S_yx = 0.5_rp*(dv_x + du_y)
                 S_yy = 0.5_rp*(dv_y + dv_y)
                 S_yz = 0.5_rp*(dv_z + dw_y)
                 S_zx = 0.5_rp*(dw_x + du_z)
                 S_zy = 0.5_rp*(dw_y + dv_z)
                 S_zz = 0.5_rp*(dw_z + dw_z)


                 sqrtS_S_ = sqrt(2.0_rp*(S_xx*S_xx+S_xy*S_xy+S_xz*S_xz + &
                                         S_yx*S_yx+S_yy*S_yy+S_yz*S_yz + &
                                         S_zx*S_zx+S_zy*S_zy+S_zz*S_zz))
                 
                 
                 ! === turbulent viscosity
                 mu_turb = phi(i,j,k,1) * Cs2 * delta2 * sqrtS_S_

                 ! === turbulent diffusivity
                 lb_turb = cp*mu_turb*i_Pr_turb

                 ! === update variables
                 VIS(i,j,k) = VIS(i,j,k) + i_mu_inf*mu_turb
                 LMD(i,j,k) = LMD(i,j,k) + i_mu_inf*lb_turb

              enddo
           enddo
        enddo
        !$acc end parallel

        !$omp end parallel do

        return
end subroutine compute_smagorinsky


subroutine compute_dynamical_smagorinsky(VIS,LMD)
        
        use parameters_module, only: central_1, central_fd_order
        use mpi_module       , only: sx,ex, sy,ey, sz,ez
        use storage_module   , only: phi, U, V, W
        use mesh_module      , only: xstep, ystep, zstep

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: VIS
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: LMD

        ! local declarations
        real(rp), parameter        :: one2 = 1.0_rp/2.0_rp
        real(rp), parameter        :: two3 = 2.0_rp/3.0_rp
        real(rp), parameter        :: alph = 8.0_rp**(two3)

        integer                    :: i,j,k,ii,jj,kk
        integer                    :: fR, l, n, m, q
        integer , dimension(3,7)   :: id
        real(rp), dimension(7)     :: rho, absS
        real(rp), dimension(3,7)   :: rui
        real(rp), dimension(3)     :: mean_rui
        real(rp)                   :: mean_rho, mean_abS
        real(rp), dimension(3,7)   :: du, dv, dw
        real(rp), dimension(3)     :: mean_du, mean_dv, mean_dw
        real(rp), dimension(3)     :: st, ist, cl1
        real(rp)                   :: i_mu_inf, CsDelta2, mu_turb, lb_turb
        real(rp)                   :: rho_, i_rho, ihRho

        real(rp), dimension(3,3)   :: meanRuiMeanRuj, meanRuiRuj
        real(rp), dimension(3,3)   :: mean_Str, mean_SS
        real(rp), dimension(3,3)   :: LG, MG
        real(rp), dimension(3,3,7) :: gradV, gradVT, Str
        real(rp), dimension(7)     :: cf

        fR = central_fd_order/2
        i_mu_inf = 1.0_rp/mu_inf

        cf = (/2.0_rp, 1._rp,1._rp,1._rp,1._rp,1._rp,1._rp/)/8.0_rp

        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex

                 ! averaging position
                 id(:,1) = (/i,j,k/)
                 id(:,2) = (/i-1,j,k/)
                 id(:,3) = (/i+1,j,k/)
                 id(:,4) = (/i,j-1,k/)
                 id(:,5) = (/i,j+1,k/)
                 id(:,6) = (/i,j,k-1/)
                 id(:,7) = (/i,j,k+1/)

                 mean_du = 0.0_rp
                 mean_dv = 0.0_rp
                 mean_dw = 0.0_rp

                 mean_rho = 0.0_rp
                 mean_rui = 0.0_rp
                 meanRuiRuj = 0.0_rp
                 mean_Str   = 0.0_rp
                 mean_AbS   = 0.0_rp
                 mean_SS    = 0.0_rp
                 do n = 1,7
                        
                    ii = id(1,n)
                    jj = id(2,n)
                    kk = id(3,n)

                    st(1) = xstep(ii)
                    st(2) = ystep(jj)
                    st(3) = zstep(kk)
                    ist   = 1.0_rp/st

                    du(:,n) = 0.0_rp
                    dv(:,n) = 0.0_rp
                    dw(:,n) = 0.0_rp
                    do l = 1, fR
                       cl1 = ist * central_1(l)

                       du(1,n) = du(1,n) + cl1(1) * (U(ii+l,jj,kk) - U(ii-l,jj,kk))
                       du(2,n) = du(2,n) + cl1(2) * (U(ii,jj+l,kk) - U(ii,jj-l,kk))
                       du(3,n) = du(3,n) + cl1(3) * (U(ii,jj,kk+l) - U(ii,jj,kk-l))

                       dv(1,n) = dv(1,n) + cl1(1) * (V(ii+l,jj,kk) - V(ii-l,jj,kk))
                       dv(2,n) = dv(2,n) + cl1(2) * (V(ii,jj+l,kk) - V(ii,jj-l,kk))
                       dv(3,n) = dv(3,n) + cl1(3) * (V(ii,jj,kk+l) - V(ii,jj,kk-l))

                       dw(1,n) = dw(1,n) + cl1(1) * (W(ii+l,jj,kk) - W(ii-l,jj,kk))
                       dw(2,n) = dw(2,n) + cl1(2) * (W(ii,jj+l,kk) - W(ii,jj-l,kk))
                       dw(3,n) = dw(3,n) + cl1(3) * (W(ii,jj,kk+l) - W(ii,jj,kk-l))

                    enddo
                
                    ! n-th location quantities
                    rho(n)   =   phi(ii,jj,kk,1)
                    rui(:,n) = (/phi(ii,jj,kk,2),phi(ii,jj,kk,3),phi(ii,jj,kk,4)/)

                    gradV (1,:,n) = du(:,n)
                    gradV (2,:,n) = dv(:,n)
                    gradV (3,:,n) = dw(:,n)

                    gradVT(:,1,n) = gradV (1,:,n)
                    gradVT(:,2,n) = gradV (2,:,n)
                    gradVT(:,3,n) = gradV (3,:,n)

                    Str(:,:,n) = 0.5_rp*(gradV(:,:,n) + gradVT(:,:,n))
                    absS(n)  = sqrt(2.0_rp*sum(Str(:,:,n)*Str(:,:,n)))

                    ! mean quantities
                    mean_rho = mean_rho + cf(n)*rho  (n)
                    mean_rui = mean_rui + cf(n)*rui(:,n)

                    mean_Str = mean_Str + cf(n)*Str(:,:,n)
                    mean_abS = mean_abS + cf(n)*absS(n)   
                    mean_SS  = mean_SS  + cf(n)*absS(n)*Str(:,:,n)

                    ! tilde(rhui*rhuj)
                    do m    = 1,3
                       do q = 1,3
                          meanRuiRuj(q,m) = meanRuiruj(q,m) + cf(n)*rui(q,n)*rui(m,n)
                       enddo
                    enddo

                 enddo
                
                 ! tilde(rhoui) * tilde(rhouj)
                 do m = 1,3
                    do q = 1,3
                       MeanRuiMeanRuj(q,m) = mean_rui(q) * mean_rui(m)
                    enddo
                 enddo
               
                 rho_  = phi(i,j,k,1)
                 i_Rho = 1.0_rp/rho_
                 ihRho = 1.0_rp/mean_rho
                
                 LG = i_Rho * meanRuiRuj - ihRho * meanRuiMeanRuj

                 MG = alph * mean_rho * mean_abS * mean_Str - &
                             rho_     * mean_SS

                 ! Dynamical Smagorinsky constant according to Lilly
                 CsDelta2 = -one2 * sum(LG*MG)/(sum(MG*MG))

                 mu_turb = i_mu_inf*phi(i,j,k,1) * CsDelta2 * absS(1)
                 if(mu_turb < -0.3_rp*VIS(i,j,k)) mu_turb = -0.3_rp*VIS(i,j,k)

                 lb_turb = cp*mu_turb*i_Pr_turb

                 VIS(i,j,k) = VIS(i,j,k) + mu_turb
                 LMD(i,j,k) = LMD(i,j,k) + lb_turb


              enddo
           enddo
        enddo
        

        return
end subroutine compute_dynamical_smagorinsky




subroutine compute_wale(VIS,LMD)

        use parameters_module, only: cp, central_1, mu_inf, bc
        use mpi_module       , only: sx,ex, sy,ey, sz,ez
        use storage_module   , only: phi, U, V, W
        use mesh_module      , only: xstep, ystep, zstep, xstep_i, ystep_i, zstep_i

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: VIS
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: LMD

        ! local declarations
        real(rp), parameter :: one3 = 1.0_rp/3.0_rp, two3 = 2.0_rp/3.0_rp
        real(rp)            :: Cw2
        real(rp)            :: delta2
        real(rp)            :: sqrtS_S_, sqrtSdSd, sumSS, sumSdSd, SdTrace
        real(rp)            :: S_S__52
        real(rp)            :: SdSd_32
        real(rp)            :: SdSd_54
        real(rp)            :: mu_turb, lb_turb, i_mu_inf

        real(rp)            :: i_st_x, i_st_y, i_st_z, cl
        real(rp)            :: cl1_x, cl1_y, cl1_z
        real(rp)            :: du_x, du_y, du_z
        real(rp)            :: dv_x, dv_y, dv_z
        real(rp)            :: dw_x, dw_y, dw_z
        real(rp)            :: S_xx, S_xy, S_xz
        real(rp)            :: S_yx, S_yy, S_yz
        real(rp)            :: S_zx, S_zy, S_zz
        real(rp)            :: Sd_xx, Sd_xy, Sd_xz
        real(rp)            :: Sd_yx, Sd_yy, Sd_yz
        real(rp)            :: Sd_zx, Sd_zy, Sd_zz
        real(rp)            :: DV_ij_2_xx, DV_ij_2_xy, DV_ij_2_xz
        real(rp)            :: DV_ij_2_yx, DV_ij_2_yy, DV_ij_2_yz
        real(rp)            :: DV_ij_2_zx, DV_ij_2_zy, DV_ij_2_zz

        integer             :: i,j,k,l
        
        i_mu_inf = 1.0_rp/mu_inf

        Cw2 = 0.325_rp**2
        !Cw2 = 0.5_rp**2 ! wmles good value

        !$omp parallel do collapse(3), default(private), &
        !$omp shared(sx,sy,sz,ex,ey,ez,xstep,ystep,zstep,xstep_i,ystep_i,zstep_i),&
        !$omp shared(phi,U,V,W,central_1,VIS,LMD,i_mu_inf,fR,Cw2)

        !$acc parallel default(present)
        !$acc loop gang,vector collapse(3)
        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex
                        
                 ! === compute grid effective step
                 delta2 = (xstep(i)*ystep(j)*zstep(k))**two3

                 ! === compute derivatives
                 i_st_x = xstep_i(i)
                 i_st_y = ystep_i(j)
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
                 !$acc loop seq
                 do l = 1,3
        
                    cl = central_1(l)
                    cl1_x = i_st_x * cl
                    cl1_y = i_st_y * cl
                    cl1_z = i_st_z * cl

                    du_x = du_x + cl1_x * (U(i+l,j,k) - U(i-l,j,k))
                    du_y = du_y + cl1_y * (U(i,j+l,k) - U(i,j-l,k))
                    du_z = du_z + cl1_z * (U(i,j,k+l) - U(i,j,k-l))

                    dv_x = dv_x + cl1_x * (V(i+l,j,k) - V(i-l,j,k))
                    dv_y = dv_y + cl1_y * (V(i,j+l,k) - V(i,j-l,k))
                    dv_z = dv_z + cl1_z * (V(i,j,k+l) - V(i,j,k-l))

                    dw_x = dw_x + cl1_x * (W(i+l,j,k) - W(i-l,j,k))
                    dw_y = dw_y + cl1_y * (W(i,j+l,k) - W(i,j-l,k))
                    dw_z = dw_z + cl1_z * (W(i,j,k+l) - W(i,j,k-l))

                 enddo

                 ! Strain tensor
                 S_xx = 0.5_rp*(du_x + du_x)
                 S_xy = 0.5_rp*(du_y + dv_x)
                 S_xz = 0.5_rp*(du_z + dw_x)
                 S_yx = 0.5_rp*(dv_x + du_y)
                 S_yy = 0.5_rp*(dv_y + dv_y)
                 S_yz = 0.5_rp*(dv_z + dw_y)
                 S_zx = 0.5_rp*(dw_x + du_z)
                 S_zy = 0.5_rp*(dw_y + dv_z)
                 S_zz = 0.5_rp*(dw_z + dw_z)

                 !
                 ! === traceless symmatric part of the square of the velocity gradient
                 !
                 ! square of the velocity gradients
                 !DV_ij_2(:,:) = matmul(grad_v , grad_v )
                 DV_ij_2_xx = du_x * du_x + du_y * dv_x + du_z * dw_x
                 DV_ij_2_xy = du_x * du_y + du_y * dv_y + du_z * dw_y
                 DV_ij_2_xz = du_x * du_z + du_y * dv_z + du_z * dw_z

                 DV_ij_2_yx = dv_x * du_x + dv_y * dv_x + dv_z * dw_x
                 DV_ij_2_yy = dv_x * du_y + dv_y * dv_y + dv_z * dw_y
                 DV_ij_2_yz = dv_x * du_z + dv_y * dv_z + dv_z * dw_z
                                                                                      
                 DV_ij_2_zx = dw_x * du_x + dw_y * dv_x + dw_z * dw_x
                 DV_ij_2_zy = dw_x * du_y + dw_y * dv_y + dw_z * dw_y
                 DV_ij_2_zz = dw_x * du_z + dw_y * dv_z + dw_z * dw_z

                 !DV_ji_2(:,:) = matmul(grad_vT, grad_vT)
                 !remember ... A^T*A^T = (A*A)^T

                 SdTrace = - one3 * (  du_x*du_x +   dv_y*dv_y +   dw_z*dw_z + &
                                     2*du_y*dv_x + 2*du_z*dw_x + 2*dv_z*dw_y)

                 ! Sd tensor
                 Sd_xx = DV_ij_2_xx + SdTrace
                 Sd_xy = 0.5_rp*(DV_ij_2_xy + DV_ij_2_yx) 
                 Sd_xz = 0.5_rp*(DV_ij_2_xz + DV_ij_2_zx) 
                 Sd_yx = 0.5_rp*(DV_ij_2_yx + DV_ij_2_xy) 
                 Sd_yy = DV_ij_2_yy + SdTrace
                 Sd_yz = 0.5_rp*(DV_ij_2_yz + DV_ij_2_zy) 
                 Sd_zx = 0.5_rp*(DV_ij_2_zx + DV_ij_2_xz) 
                 Sd_zy = 0.5_rp*(DV_ij_2_zy + DV_ij_2_yz) 
                 Sd_zz = DV_ij_2_zz + SdTrace


                 ! compute tensor norms
                 sumSS   = S_xx*S_xx+S_xy*S_xy+S_xz*S_xz + &
                           S_yx*S_yx+S_yy*S_yy+S_yz*S_yz + &
                           S_zx*S_zx+S_zy*S_zy+S_zz*S_zz

                 sumSdSd = Sd_xx*Sd_xx+Sd_xy*Sd_xy+Sd_xz*Sd_xz + &
                           Sd_yx*Sd_yx+Sd_yy*Sd_yy+Sd_yz*Sd_yz + &
                           Sd_zx*Sd_zx+Sd_zy*Sd_zy+Sd_zz*Sd_zz

                 sqrtS_S_ = sqrt(sumSS)
                 sqrtSdSd = sqrt(sumSdSd)

                 SdSd_32 = sqrtSdSd**3
                 SdSd_54 = sqrt(sqrtSdSd**5)
                 S_S__52 = sqrtS_S_**5
                
                 ! === turbulent viscosity
                 mu_turb = phi(i,j,k,1) * Cw2 * delta2 * SdSd_32/(S_S__52 + SdSd_54)
                 if(SdSd_32 < 1.0E-14_rp) mu_turb = 0.0_rp

                 ! === turbulent diffusivity
                 lb_turb = cp*mu_turb*i_Pr_turb

                 ! === update variables
                 VIS(i,j,k) = VIS(i,j,k) + i_mu_inf*mu_turb
                 LMD(i,j,k) = LMD(i,j,k) + i_mu_inf*lb_turb

            enddo
          enddo
        enddo
        !$acc end parallel

        !$omp end parallel do



        return
end subroutine compute_wale



subroutine compute_sigma_model(VIS,LMD)
        
        use parameters_module, only: cp, central_1, central_fd_order, mu_inf
        use mpi_module       , only: sx,ex, sy,ey, sz,ez
        use storage_module   , only: phi, U, V, W
        use mesh_module      , only: xstep, ystep, zstep, xstep_i, ystep_i, zstep_i

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: VIS
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: LMD

        ! local declarations
        real(rp)            :: i_st_x, i_st_y, i_st_z
        real(rp)            :: cl1_x, cl1_y, cl1_z
        real(rp)            :: du_x, du_y, du_z
        real(rp)            :: dv_x, dv_y, dv_z
        real(rp)            :: dw_x, dw_y, dw_z
        real(rp)            :: G_xx, G_xy, G_xz
        real(rp)            :: G_yx, G_yy, G_yz
        real(rp)            :: G_zx, G_zy, G_zz
        real(rp)            :: lambda_x, lambda_y, lambda_z
        real(rp)            :: s_x, s_y, s_z
        real(rp), parameter      :: two3 = 2.0_rp/3.0_rp, Cs2 = 1.0_rp**2 !Cs2 = 1.5_rp**2 (secondo teoria!)
        real(rp)                 :: delta2
        real(rp)                 :: mu_turb, lb_turb, i_mu_inf
        integer                  :: fR
        integer                  :: i,j,k,l
        
        ! declarations for eigenvalues computations
        real(rp), parameter :: one3 = 1.0_rp/3.0_rp, one6 = 1.0_rp/6.0_rp
        real(rp), parameter :: twoPithird = 2*pi*one3
        real(rp), parameter :: toll = 1.0E-14_rp
        real(rp)            :: diag_xx, diag_yy, diag_zz
        real(rp)            :: B_xx, B_xy, B_xz
        real(rp)            :: B_yx, B_yy, B_yz
        real(rp)            :: B_zx, B_zy, B_zz
        real(rp)            :: trac, q, p, p1, p2, phi_, r, ip

        fR       = central_fd_order/2
        i_mu_inf = 1.0_rp/mu_inf

        !$omp parallel do collapse(3), default(private), &
        !$omp shared(sx,sy,sz,ex,ey,ez,xstep,ystep,zstep,xstep_i,ystep_i,zstep_i),&
        !$omp shared(phi,U,V,W,central_1,VIS,LMD,i_mu_inf,fR)
        !$acc parallel default(present)
        !$acc loop gang,vector collapse(3)
        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex
                        
                 ! === compute grid effective step
                 delta2 = (xstep(i)*ystep(j)*zstep(k))**two3

                 ! === compute derivatives
                 i_st_x = xstep_i(i)
                 i_st_y = ystep_i(j)
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
                 !$acc loop seq
                 do l = 1,fR

                    cl1_x = i_st_x * central_1(l)
                    cl1_y = i_st_y * central_1(l)
                    cl1_z = i_st_z * central_1(l)

                    du_x = du_x + cl1_x * (U(i+l,j,k) - U(i-l,j,k))
                    du_y = du_y + cl1_y * (U(i,j+l,k) - U(i,j-l,k))
                    du_z = du_z + cl1_z * (U(i,j,k+l) - U(i,j,k-l))

                    dv_x = dv_x + cl1_x * (V(i+l,j,k) - V(i-l,j,k))
                    dv_y = dv_y + cl1_y * (V(i,j+l,k) - V(i,j-l,k))
                    dv_z = dv_z + cl1_z * (V(i,j,k+l) - V(i,j,k-l))

                    dw_x = dw_x + cl1_x * (W(i+l,j,k) - W(i-l,j,k))
                    dw_y = dw_y + cl1_y * (W(i,j+l,k) - W(i,j-l,k))
                    dw_z = dw_z + cl1_z * (W(i,j,k+l) - W(i,j,k-l))

                 enddo

                 ! === compute G tensor
                 G_xx = du_x*du_x+dv_x*dv_x+dw_x*dw_x
                 G_yx = du_y*du_x+dv_y*dv_x+dw_y*dw_x
                 G_zx = du_z*du_x+dv_z*dv_x+dw_z*dw_x
                 G_xy = du_x*du_y+dv_x*dv_y+dw_x*dw_y
                 G_yy = du_y*du_y+dv_y*dv_y+dw_y*dw_y
                 G_zy = du_z*du_y+dv_z*dv_y+dw_z*dw_y
                 G_xz = du_x*du_z+dv_x*dv_z+dw_x*dw_z
                 G_yz = du_y*du_z+dv_y*dv_z+dw_y*dw_z
                 G_zz = du_z*du_z+dv_z*dv_z+dw_z*dw_z

                 ! === compute G eigenvalues and singular values of DU
                 diag_xx = G_xx
                 diag_yy = G_yy
                 diag_zz = G_zz
                 trac    = diag_xx + diag_yy + diag_zz
                   p1 = G_xy*G_xy + G_xz*G_xz + G_yz*G_yz

                 if(p1 < toll) then ! G is diagonal

                   lambda_x = max(diag_xx, diag_yy, diag_zz)
                   lambda_z = min(diag_xx, diag_yy, diag_zz)
                   lambda_y = trac - (lambda_x + lambda_z) 

                 else

                   q  = trac*one3
                   p2 = diag_xx*diag_xx + diag_yy*diag_yy + diag_zz*diag_zz - 3*q*q + 2*p1
                   p  = sqrt(p2*one6)
                   
                   ! B matrix
                   ip   = 1.0_rp/p
                   B_xx = ip*(G_xx - q)
                   B_xy = ip*(G_xy)
                   B_xz = ip*(G_xz)
                   !
                   B_yx = ip*(G_yx)
                   B_yy = ip*(G_yy - q)
                   B_yz = ip*(G_yz)
                   !
                   B_zx = ip*(G_zx)
                   B_zy = ip*(G_zy)
                   B_zz = ip*(G_zz - q)

                   r = 0.5_rp*(  B_xx*(B_yy*B_zz - B_yz*B_zy)  &
                               - B_xy*(B_yx*B_zz - B_yz*B_zx)  &
                               + B_xz*(B_yx*B_zy - B_yy*B_zx) )

                   if(r < -1.0_rp) then
                     phi_ = pi*one3
                   elseif(r > 1.0_rp) then
                     phi_ = 0.0_rp
                   else
                     phi_ = acos(r)*one3
                   endif

                   lambda_x = q + 2*p*cos(phi_)
                   lambda_z = q + 2*p*cos(phi_ + twoPithird)
                   lambda_y = trac - (lambda_x + lambda_z)

                 endif
                 ! eigenvalues computed
        
                 !> abs is just to avoid lambda to negative (even if very very small)
                 s_x = sqrt(abs(lambda_x)) 
                 s_y = sqrt(abs(lambda_y)) 
                 s_z = sqrt(abs(lambda_z)) 

                 ! === turbulent viscosity
                 mu_turb = phi(i,j,k,1) * Cs2 * delta2 * s_z* (s_x-s_y) * (s_y-s_z)/(s_x**2)

                 ! === turbulent diffusivity
                 lb_turb = cp*mu_turb*i_Pr_turb

                 ! === update variables
                 VIS(i,j,k) = VIS(i,j,k) + i_mu_inf*mu_turb
                 LMD(i,j,k) = LMD(i,j,k) + i_mu_inf*lb_turb

              enddo
           enddo
        enddo
        !$acc end parallel
        !$omp end parallel do

        return
end subroutine compute_sigma_model



subroutine compute_MixedTimeScale(VIS,LMD)
        
        use parameters_module, only: cp, central_1, central_fd_order, mu_inf
        use mpi_module       , only: sx,ex, sy,ey, sz,ez
        use storage_module   , only: phi, U, V, W
        use mesh_module      , only: xstep, ystep, zstep, xstep_i, ystep_i, zstep_i

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: VIS
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: LMD

        ! local declarations
        real(rp), parameter   :: two3 = 2.0_rp/3.0_rp
        real(rp), parameter   :: one7 = 1.0_rp/7.0_rp
        real(rp), parameter   :: Cmts = 0.03_rp
        real(rp), parameter   :: C_t  = 10.0_rp
        real(rp)              :: i_st_x, i_st_y, i_st_z
        real(rp)              :: cl1_x, cl1_y, cl1_z
        real(rp)              :: du_x, du_y, du_z
        real(rp)              :: dv_x, dv_y, dv_z
        real(rp)              :: dw_x, dw_y, dw_z
        real(rp)              :: S_xx, S_xy, S_xz
        real(rp)              :: S_yx, S_yy, S_yz
        real(rp)              :: S_zx, S_zy, S_zz
        real(rp)              :: mean_phi_1,mean_phi_2,mean_phi_3,mean_phi_4
        real(rp)              :: imeanr, mean_u, mean_v, mean_w
        real(rp)              :: asS, kes, TsA, TsB, T_S
        real(rp)              :: delta2, sumSS
        real(rp)              :: mu_turb, lb_turb, i_mu_inf
        integer               :: fR
        integer               :: i,j,k,l
        
        fR       = central_fd_order/2
        i_mu_inf = 1.0_rp/mu_inf

        !$acc parallel default(present)
        !$acc loop gang,vector collapse(3)
        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex
                        
                 ! === compute grid effective step
                 delta2 = (xstep(i)*ystep(j)*zstep(k))**two3

                 ! === compute derivatives
                 i_st_x = xstep_i(i)
                 i_st_y = ystep_i(j)
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
                 !$acc loop seq
                 do l = 1,fR

                    cl1_x = i_st_x * central_1(l)
                    cl1_y = i_st_y * central_1(l)
                    cl1_z = i_st_z * central_1(l)

                    du_x = du_x + cl1_x * (U(i+l,j,k) - U(i-l,j,k))
                    du_y = du_y + cl1_y * (U(i,j+l,k) - U(i,j-l,k))
                    du_z = du_z + cl1_z * (U(i,j,k+l) - U(i,j,k-l))

                    dv_x = dv_x + cl1_x * (V(i+l,j,k) - V(i-l,j,k))
                    dv_y = dv_y + cl1_y * (V(i,j+l,k) - V(i,j-l,k))
                    dv_z = dv_z + cl1_z * (V(i,j,k+l) - V(i,j,k-l))

                    dw_x = dw_x + cl1_x * (W(i+l,j,k) - W(i-l,j,k))
                    dw_y = dw_y + cl1_y * (W(i,j+l,k) - W(i,j-l,k))
                    dw_z = dw_z + cl1_z * (W(i,j,k+l) - W(i,j,k-l))

                 enddo

                 ! === compute strain tensor
                 S_xx = 0.5_rp*(du_x + du_x)
                 S_xy = 0.5_rp*(du_y + dv_x)
                 S_xz = 0.5_rp*(du_z + dw_x)
                 S_yx = 0.5_rp*(dv_x + du_y)
                 S_yy = 0.5_rp*(dv_y + dv_y)
                 S_yz = 0.5_rp*(dv_z + dw_y)
                 S_zx = 0.5_rp*(dw_x + du_z)
                 S_zy = 0.5_rp*(dw_y + dv_z)
                 S_zz = 0.5_rp*(dw_z + dw_z)
                
                 ! compute double filtered speed
                 !mean_phi = (phi(i,j,k,:) + phi(i+1,j,k,:) + phi(i-1,j,k,:) + &
                 !                           phi(i,j+1,k,:) + phi(i,j-1,k,:) + &
                 !                           phi(i,j,k+1,:) + phi(i,j,k-1,:))*one7

                 mean_phi_1 = (phi(i,j,k,1) + phi(i+1,j,k,1) + phi(i-1,j,k,1) + &
                                              phi(i,j+1,k,1) + phi(i,j-1,k,1) + &
                                              phi(i,j,k+1,1) + phi(i,j,k-1,1))*one7
                 mean_phi_2 = (phi(i,j,k,2) + phi(i+1,j,k,2) + phi(i-1,j,k,2) + &
                                              phi(i,j+1,k,2) + phi(i,j-1,k,2) + &
                                              phi(i,j,k+1,2) + phi(i,j,k-1,2))*one7
                 mean_phi_3 = (phi(i,j,k,3) + phi(i+1,j,k,3) + phi(i-1,j,k,3) + &
                                              phi(i,j+1,k,3) + phi(i,j-1,k,3) + &
                                              phi(i,j,k+1,3) + phi(i,j,k-1,3))*one7
                 mean_phi_4 = (phi(i,j,k,4) + phi(i+1,j,k,4) + phi(i-1,j,k,4) + &
                                              phi(i,j+1,k,4) + phi(i,j-1,k,4) + &
                                              phi(i,j,k+1,4) + phi(i,j,k-1,4))*one7

                 imeanr = 1.0_rp/mean_phi_1
                 mean_u = mean_phi_2*imeanr
                 mean_v = mean_phi_3*imeanr
                 mean_w = mean_phi_4*imeanr
        
                 ! compute turbulent viscosity
                 kes = (U(i,j,k) - mean_u)**2+&
                       (V(i,j,k) - mean_v)**2+&
                       (W(i,j,k) - mean_w)**2

                 sumSS = S_xx*S_xx+S_xy*S_xy+S_xz*S_xz + &
                         S_yx*S_yx+S_yy*S_yy+S_yz*S_yz + &
                         S_zx*S_zx+S_zy*S_zy+S_zz*S_zz
                
                 asS = sqrt(2.0_rp*sumSS)
                 TsA = sqrt(kes/delta2)
                 TsB = asS/C_t
                 T_s = 1.0_rp/(TsA + TsB)

                 mu_turb = phi(i,j,k,1) * Cmts*kes*T_s
                 lb_turb = cp*mu_turb*i_Pr_turb

                 ! === update variables
                 VIS(i,j,k) = VIS(i,j,k) + i_mu_inf*mu_turb
                 LMD(i,j,k) = LMD(i,j,k) + i_mu_inf*lb_turb

              enddo
           enddo
        enddo
        !$acc end parallel

        return
end subroutine compute_MixedTimeScale

subroutine TrbVisBC(phi,VIS,LMD)

        use bc_module, only: face_type

        implicit none
        real(rp), dimension(:,:,:)  , allocatable, intent(inout) :: VIS, LMD
        real(rp), dimension(:,:,:,:), allocatable, intent(in)    :: phi

        type(face_type), dimension(6) :: all_bound
        type(face_type)               :: bound
        real(rp)                      :: tWall
        integer                       :: f

        call StartProfRange("TrbVisBC")

        if(.not.mpi_flag) my_neighbour = MPI_PROC_NULL

        all_bound%node = (/ sx , ex , sy , ey , sz ,  ez/)
        all_bound%norm = (/ -1 ,  1 , -1 ,  1 , -1 ,  1 /)
        all_bound%face = (/ 'W', 'E', 'S', 'N', 'B', 'F'/)

        do f = 1, 6
           if(my_neighbour(f) == MPI_PROC_NULL) then

             bound%node = all_bound(f)%node
             bound%norm = all_bound(f)%norm
             bound%face = all_bound(f)%face

             selectcase(trim(bc(f)))
               case('adiabatic_wall','nscbc_adiabatic_wall')
                 tWall = Trat*(1.0_rp + Prandtl**(1.0_rp/3.0_rp)*0.5_rp*(gamma0-1.0_rp)*Mach**2)
                 call WallResolvedTrbVis(bound,tWall,VIS,LMD)

               case('isothermal_wall','nscbc_isothermal_wall')
                 tWall = 1.0_rp
                 call WallResolvedTrbVis(bound,tWall,VIS,LMD)

               case('idws_isothermal')
                 tWall = 1.0_rp
                 call iDifferentialWallModelledTrbVis(bound,tWall,phi,VIS,LMD)

               case('idws_adiabatic')
                 tWall = Trat*(1.0_rp + Prandtl**(1.0_rp/3.0_rp)*0.5_rp*(gamma0-1.0_rp)*Mach**2)
                 call iDifferentialWallModelledTrbVis(bound,tWall,phi,VIS,LMD)

               case('static_dws_isothermal')
                 tWall = 1.0_rp
                 call StaticDifferentialWallModelledTrbVis(bound,tWall,phi,VIS,LMD)

               case('static_dws_adiabatic')
                 tWall = Trat*(1.0_rp + Prandtl**(1.0_rp/3.0_rp)*0.5_rp*(gamma0-1.0_rp)*Mach**2)
                 call StaticDifferentialWallModelledTrbVis(bound,tWall,phi,VIS,LMD)

               case('istatic_dws_adiabatic')
                 tWall = Trat*(1.0_rp + Prandtl**(1.0_rp/3.0_rp)*0.5_rp*(gamma0-1.0_rp)*Mach**2)
                 call iStaticDifferentialWallModelledTrbVis(bound,tWall,phi,VIS,LMD)

               case default
                 call ZeroGradientTrbVis(bound,VIS,LMD)

             endselect


           endif
        enddo

        call EndProfRange


        return
end subroutine TrbVisBC

subroutine ZeroGradientTrbVis(b,VIS,LMD)

        use bc_module, only: face_type

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: VIS,LMD
        type(face_type)                        , intent(in)    :: b

        ! local declarations
        integer :: i,j,k
        integer :: ii,ig,ji,jg,ki,kg,bnode,bnorm

        bnode = b%node
        bnorm = b%norm

        selectcase(b%face)
          case('E','W')
            !$acc parallel default(present)
            !$acc loop gang,vector collapse(2)
            do       k = sz,ez
               do    j = sy,ey
                  !$acc loop seq
                  do i = 1,GN

                     ig = bnode + bnorm*i
                     ii = bnode - bnorm*(i-1)

                     VIS(ig,j,k) = VIS(ii,j,k)
                     LMD(ig,j,k) = LMD(ii,j,k)

                  enddo
               enddo
            enddo
            !$acc end parallel

          case('S','N')
            !$acc parallel default(present)
            !$acc loop gang,vector collapse(2)
            do       k = sz,ez
               do    i = sx,ex
                  !$acc loop seq
                  do j = 1,GN

                     jg = bnode + bnorm*j
                     ji = bnode - bnorm*(j-1)

                     VIS(i,jg,k) = VIS(i,ji,k)
                     LMD(i,jg,k) = LMD(i,ji,k)

                  enddo
               enddo
            enddo
            !$acc end parallel

          case('B','F')
            !$acc parallel default(present)
            !$acc loop gang,vector collapse(2)
            do       i = sx,ex
               do    j = sy,ey
                  !$acc loop seq
                  do k = 1,GN

                     kg = bnode + bnorm*k
                     ki = bnode - bnorm*(k-1)

                     VIS(i,j,kg) = VIS(i,j,ki)
                     LMD(i,j,kg) = LMD(i,j,ki)

                  enddo
               enddo
            enddo
            !$acc end parallel

        endselect


        return
end subroutine ZeroGradientTrbVis


subroutine WallResolvedTrbVis(b,tWall,VIS,LMD)

        use fluid_functions_module, only: laminar_viscosity
        use bc_module             , only: face_type
        use parameters_module     , only: k_inf 

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: VIS,LMD
        real(rp)                               , intent(in)    :: tWall
        type(face_type)                        , intent(in)    :: b

        ! local declaration
        real(rp) :: mWall, kWall
        integer  :: i,j,k,ji,jg,bnode,bnorm

        mWall = laminar_viscosity(tWall,Tref,vis_flag)
        kWall = k_inf * mWall

        bnode = b%node
        bnorm = b%norm

        selectcase(b%face)

          case('S','N')
            !$acc parallel default(present)
            !$acc loop gang,vector collapse(2)
            do       k = sz,ez
               do    i = sx,ex
                  !$acc loop seq
                  do j = 1,GN

                     jg = bnode + bnorm*j
                     ji = bnode - bnorm*(j-1)

                     VIS(i,jg,k) = 2*mWall - VIS(i,ji,k)
                     LMD(i,jg,k) = 2*kWall - LMD(i,ji,k)

                  enddo
               enddo
            enddo
            !$acc end parallel

          case default
            print*, ' WallResolvedTrbVis is not implemented for face ', &
                      trim(b%face)
            stop

        endselect


        return
end subroutine WallResolvedTrbVis



subroutine iDifferentialWallModelledTrbVis(b,tWall,phi,VIS,LMD)

        use parameters_module      , only: gamma0, Prandtl
        use bc_module              , only: face_type
        use mesh_module            , only: y, xstep, zstep
        use fluid_functions_module , only: laminar_viscosity, getmuT
        use storage_module         , only: WMLES_DATA_LW, WMLES_DATA_UW, WMLES_DATA, &
                                           nvwmlesdata
        use wmles_module

        implicit none
        real(rp), dimension(:,:,:)  , allocatable, intent(inout) :: VIS,LMD
        real(rp), dimension(:,:,:,:), allocatable, intent(in)    :: phi
        type(face_type)                          , intent(in)    :: b
        real(rp)                                 , intent(in)    :: tWall

        ! local declaration
        real(rp) :: yw , hw , u_w, T_w
        real(rp) :: r_1, ir1, u_1, v_1, w_1, ek1, p_1, u_p
        real(rp) :: r_h, irh, u_h, v_h, w_h, ekh, p_h, T_1, T_h
        integer  :: j0, j1, i,j,k,ji, jg
        
        real(rp), parameter :: ClipPar = 100.0_rp
        real(rp)            :: mW, lW, wm_sensor
        real(rp)            :: tauW_WM, qauW_WM, tauW_WR, qauW_WR
        real(rp)            :: mu_R, lm_R

        real(rp) :: irho0
        real(rp) :: utau0, inuw0
        real(rp) :: xPl0, yPl0, zPl0, yPlWalln
        integer  :: jInt, jl, b_norm, b_node, b_face, l
        real(rp) :: dyh, idyh, iy1

        call StartProfRange("iDifferentialWallModelledTrbVis")

        ! set wall quantities
        T_w = tWall
        u_w = 0.0_rp

        mW = laminar_viscosity(tWall,Tref,vis_flag)
        lW = k_inf*mW

        b_node = b%node
        b_norm = b%norm
        if(b%face == 'S') b_face = 0
        if(b%face == 'N') b_face = 1

        selectcase(b%face)
        
          case('S','N')
        
          j1 = b_node
          j0 = b_node + b_norm
          yW = 0.5_rp*(y(j1) + y(j0))
          jl = b_node - b_norm*(wmles_indx-1)
          dyh  = 0.5_rp*abs(y(j1) - y(j0))
          idyh = 1.0_rp/dyh
          iy1  = 1.0_rp/(y(j1)-yW)

          !$acc parallel default(present)
          !$acc loop gang, vector collapse(2)
          do k    = sz,ez
             do i = sx,ex
                !
                ! get the numerical shear stress and heat flux
                !
                r_1 = phi(i,j1,k,1)
                ir1 = 1.0_rp/r_1
                u_1 = phi(i,j1,k,2)*ir1
                v_1 = phi(i,j1,k,3)*ir1
                w_1 = phi(i,j1,k,4)*ir1
                ek1 = 0.5_rp*(u_1*u_1 + v_1*v_1 + w_1*w_1)
                p_1 = (gamma0-1.0_rp) * (phi(i,j1,k,5) - r_1*ek1)
                T_1 = p_1*ir1

                tauW_WR = mu_inf*mW*(u_1+w_1)*idyh
                qauW_WR = mu_inf*lW*(T_1-T_w)*idyh
                !
                ! get the viscous length
                !
                irho0 = T_w/p_1
                utau0 = sqrt(abs(tauW_WR)*irho0)
                inuw0 = utau0/(irho0*mW*mu_inf)

                xPl0  = xstep(i)*inuw0
                yPl0  = dyh     *inuw0
                zPl0  = zstep(k)*inuw0

                if(xPl0 < xpt_wm .and. yPl0 < ypt_wm .and. zPl0 < zpt_wm) then
                  ! WALL RESOLVED
        
                  mu_R = 1.0_rp
                  lm_R = 1.0_rp

                  wm_sensor = 0.0_rp

                else
                  ! WALL MODELLED
                  !
                  ! === looking at the grid point where yPlus>wmles_interface
                  !
                  ji = 1
                  do j = wmles_strI,ey
                    ji = b_node - b_norm*(j-1)
                    yPlWalln = yPl0*(y(ji)-yW)*iy1
                    if(yPlWalln > wmles_intr) exit
                  enddo
                  jInt = ji
                  hw = abs(y(jInt) - yW)

                  !
                  ! === get the LES field at the interface location
                  !
                  r_h = phi(i,jInt,k,1)
                  irh = 1.0_rp/r_h
                  u_h = phi(i,jInt,k,2)*irh
                  v_h = phi(i,jInt,k,3)*irh
                  w_h = phi(i,jInt,k,4)*irh
                  ekh = 0.5_rp*(u_h*u_h + v_h*v_h + w_h*w_h)
                  p_h = (gamma0-1.0_rp) * (phi(i,jInt,k,5) - r_h*ekh)
                  T_h = p_h*irh

                  u_p = sqrt(u_h*u_h + w_h*w_h)

                  !
                  ! === get the correct tauWall and qWall from ODE model
                  !
                  call OdeWMLES2(gamma0,Prandtl,mu_inf,Tref,vis_flag,hw,inuw0,u_w,u_p,T_w,T_h,p_h,tauW_WM,qauW_WM)

                  mu_R = abs(tauW_WM/tauW_WR)
                  lm_R = abs(qauW_WM/qauW_WR)

                  wm_sensor = 1.0_rp

                endif

                !
                ! enforce the effective viscosity
                !
                if(mu_R >  clipPar) mu_R =  ClipPar
                if(lm_R >  clipPar) lm_R =  ClipPar

                do j = 1,GN
        
                   jg = b_node + b_norm*j
                   ji = b_node - b_norm*(j-1)

                   VIS(i,jg,k) = 2*mu_R*mW - VIS(i,ji,k)
                   LMD(i,jg,k) = 2*lm_R*lW - LMD(i,ji,k)

                enddo
                !
                ! save results for post treatments
                !
                if(b_face == 0) then
                    WMLES_DATA_LW(i,k,1) = mu_R
                    WMLES_DATA_LW(i,k,2) = lm_R
                    WMLES_DATA_LW(i,k,3) = mu_R*tauW_WR
                    WMLES_DATA_LW(i,k,4) = lm_R*qauW_WR
                    WMLES_DATA_LW(i,k,5) = wm_sensor
                elseif(b_face == 1) then
                    WMLES_DATA_UW(i,k,1) = mu_R
                    WMLES_DATA_UW(i,k,2) = lm_R
                    WMLES_DATA_UW(i,k,3) = mu_R*tauW_WR
                    WMLES_DATA_UW(i,k,4) = lm_R*qauW_WR
                    WMLES_DATA_UW(i,k,5) = 1.0_rp
                endif

                if(ic == 'turbulent_channel') then
                  do l = 1,nvWMLESData
                     WMLES_DATA(i,k,l) = 0.5_rp*(WMLES_DATA_LW(i,k,l) + WMLES_DATA_UW(i,k,l))
                  enddo
                endif


             enddo
          enddo
          !$acc end parallel

          case default
          print*, 'DifferentialWallStressModel not implemented for face ', trim(b%face)
          stop
        
        endselect


        call EndProfRange
        
        return
end subroutine iDifferentialWallModelledTrbVis









subroutine iStaticDifferentialWallModelledTrbVis(b,tWall,phi,VIS,LMD)

        use parameters_module      , only: gamma0, Prandtl
        use bc_module              , only: face_type
        use mesh_module            , only: y, xstep, zstep
        use fluid_functions_module , only: laminar_viscosity, getmuT
        use storage_module         , only: WMLES_DATA_LW
        use wmles_module

        implicit none
        real(rp), dimension(:,:,:)  , allocatable, intent(inout) :: VIS,LMD
        real(rp), dimension(:,:,:,:), allocatable, intent(in)    :: phi
        type(face_type)                          , intent(in)    :: b
        real(rp)                                 , intent(in)    :: tWall

        ! local declaration
        integer , parameter       :: nw = wmles_npts
        real(rp), dimension(1:nw) :: u_wm,T_wm
        
        real(rp) :: yw , hw , u_w, T_w
        real(rp) :: r_1, ir1, u_1, v_1, w_1, ek1, p_1, u_p
        real(rp) :: r_h, irh, u_h, v_h, w_h, ekh, p_h, T_1, T_h
        integer  :: j0, j1, i,j,k,ji, jg
        
        real(rp), parameter :: ClipPar = 100.0_rp
        real(rp)            :: mW, lW, wm_sensor
        real(rp)            :: tauW_WM, qauW_WM, tauW_WR, qauW_WR
        real(rp)            :: mu_R, lm_R

        real(rp) :: rhoWall0, tauWall0
        real(rp) :: utaWall0, dltWall0, inuWall0
        real(rp) :: xPlWall0, yPlWall0, zPlWall0
        integer  :: jl

        real(rp) :: dyh, xpR, ypR, zpR

        ! set wall quantities
        T_w = tWall
        u_w = 0.0_rp

        mW = laminar_viscosity(tWall,Tref,vis_flag)
        lW = k_inf*mW

        selectcase(b%face)
        
          case('S')
        
          j1 = b%node
          j0 = b%node - 1
          yW = 0.5_rp*(y(j1) + y(j0))
          jl = b%node - b%norm*(wmles_indx-1)

          dyh  = 0.5_rp*(y(j1) - y(j0))

          do k    = sz,ez
             do i = sx,ex
                !
                ! get the numerical shear stress and heat flux
                !
                r_1 = phi(i,j1,k,1)
                ir1 = 1.0_rp/r_1
                u_1 = phi(i,j1,k,2)*ir1
                v_1 = phi(i,j1,k,3)*ir1
                w_1 = phi(i,j1,k,4)*ir1
                ek1 = 0.5_rp*(u_1*u_1 + v_1*v_1 + w_1*w_1)
                p_1 = (gamma0-1.0_rp) * (phi(i,j1,k,5) - r_1*ek1)
                T_1 = p_1*ir1

                tauW_WR = mu_inf*mW*(u_1+w_1)/dyh
                qauW_WR = mu_inf*lW*(T_1-T_w)/dyh
                !
                ! === get the LES field at the jl location
                !
                r_h = phi(i,jl,k,1)
                irh = 1.0_rp/r_h
                u_h = phi(i,jl,k,2)*irh
                v_h = phi(i,jl,k,3)*irh
                w_h = phi(i,jl,k,4)*irh
                ekh = 0.5_rp*(u_h*u_h + v_h*v_h + w_h*w_h)
                p_h = (gamma0-1.0_rp) * (phi(i,jl,k,5) - r_h*ekh)
                T_h = p_h*irh

                u_p = sqrt(u_h*u_h + w_h*w_h)
                !
                ! estimate viscous length based on resolved shear stress
                !
                tauWall0 = tauW_WR
                rhoWall0 = p_1/T_w
                utaWall0 = sqrt(abs(tauWall0)/rhoWall0)
                inuWall0 = rhoWall0*utaWall0/(mW*mu_inf)
                !
                ! get internal resolutions estimation
                !
                xPlWall0 = xstep(i)*inuWall0
                yPlWall0 = dyh     *inuWall0
                zPlWall0 = zstep(k)*inuWall0
                dltWall0 = yPlWall0/dyh

                xpR = xPlWall0/xpt_wm
                ypR = yPlWall0/ypt_wm
                zpR = zPlWall0/zpt_wm
                
                if(xpR < 1.0_rp .and. ypR < 1.0_rp .and. zpR < 1.0_rp) then
                  ! WALL RESOLVED
        
                  mu_R = 1.0_rp
                  lm_R = 1.0_rp

                  wm_sensor = 0.0_rp

                else
                  ! WALL MODELLED
                  !
                  ! === get the grid
                  !
                  hw = abs(y(jl) - yW)
                  call LogLawWMLES(mu_inf,Tref,vis_flag,hw,T_w,u_p,p_h,tauWall0)
                  rhoWall0 = p_h/T_w
                  utaWall0 = sqrt(abs(tauWall0)/rhoWall0)
                  dltWall0 = rhoWall0*utaWall0/(mW*mu_inf)

                  !
                  ! === get the correct tauWall and qWall from ODE model
                  !
                  call OdeWMLES(gamma0,Prandtl,mu_inf,Tref,vis_flag,hw,dltWall0,u_w,u_p,T_w,T_h,p_h,u_wm,T_wm,tauW_WM,qauW_WM)

                  mu_R = abs(tauW_WM/tauW_WR)
                  lm_R = abs(qauW_WM/qauW_WR)

                  wm_sensor = 1.0_rp

                endif

                !
                ! enforce the effective viscosity
                !
                if(mu_R >  clipPar) mu_R =  ClipPar
                if(lm_R >  clipPar) lm_R =  ClipPar

                do j = 1,GN
        
                   jg = b%node + b%norm*j
                   ji = b%node - b%norm*(j-1)

                   VIS(i,jg,k) = 2*mu_R*mW - VIS(i,ji,k)
                   LMD(i,jg,k) = 2*lm_R*lW - LMD(i,ji,k)

                enddo
                !
                ! save results for post treatments
                !
                WMLES_DATA_LW(i,k,1) = mu_R
                WMLES_DATA_LW(i,k,2) = lm_R
                WMLES_DATA_LW(i,k,3) = mu_R*tauW_WR
                WMLES_DATA_LW(i,k,4) = lm_R*qauW_WR
                WMLES_DATA_LW(i,k,5) = wm_sensor

             enddo
          enddo

          case default
          print*, 'iStaticWallStressModel not implemented for face ', trim(b%face)
          stop
        
        endselect
        


        
        return
end subroutine iStaticDifferentialWallModelledTrbVis























































subroutine StaticDifferentialWallModelledTrbVis(b,tWall,phi,VIS,LMD)

        use parameters_module      , only: gamma0, k_inf, mu_inf
        use storage_module         , only: WMLES_DATA_LW,WMLES_DATA_UW,WMLES_DATA,nvwmlesdata
        use bc_module              , only: face_type
        use mesh_module            , only: y
        use fluid_functions_module , only: laminar_viscosity, getmuT
        use real_to_integer_module , only: locate
        use fileModule             , only: str
        use wmles_module

        implicit none
        real(rp), dimension(:,:,:)  , allocatable, intent(inout) :: VIS,LMD
        real(rp), dimension(:,:,:,:), allocatable, intent(in)    :: phi
        type(face_type)                          , intent(in)    :: b
        real(rp)                                 , intent(in)    :: tWall

        ! local declaration
        real(rp), parameter       :: ClipPar = 100.0_rp
        integer , parameter       :: nw = wmles_npts
        real(rp), dimension(1:nw) :: u_wm,T_wm
        
        real(rp) :: yw , u_w, T_w, hw
        real(rp) :: r_h, irh, u_h, v_h, w_h, ekh, p_h, T_h, u_p
        real(rp) :: r_1, ir1, u_1, v_1, w_1, ek1, p_1, T_1
        real(rp) :: rhoWall0, tauWall0, utaWall0, dltWall0

        real(rp) :: dyh, mW, lW
        real(rp) :: tauW_WM, tauW_WR, mu_R
        real(rp) :: qauW_WM, qauW_WR, lm_R
        integer  :: j0, j1, i,j,k, ji, jg, jInt, jl, l
        integer  :: b_node, b_norm
        integer  :: b_face
        
        ! set wall quantities
        T_w = tWall
        u_w = 0.0_rp

        mW = laminar_viscosity(tWall,Tref,vis_flag)
        lW = k_inf*mW

        b_norm = b%norm
        b_node = b%node
        if(b%face == 'S') b_face = 0
        if(b%face == 'N') b_face = 1

        selectcase(b%face)
        
          case('S','N')
        
          j0 = b%node + b%norm
          j1 = b%node
          yW = 0.5_rp*(y(j1) + y(j0))
          jInt = wmles_strI
          jl = b%node - b%norm*(jInt-1)

          dyh  = 0.5_rp*abs(y(j1) - y(j0))

          !$acc parallel default(present)
          !$acc loop collapse(2) private(u_wm, T_wm)
          do k    = sz,ez
             do i = sx,ex
                !
                ! get the numerical shear stress and heat flux
                !
                r_1 = phi(i,j1,k,1)
                ir1 = 1.0_rp/r_1
                u_1 = phi(i,j1,k,2)*ir1
                v_1 = phi(i,j1,k,3)*ir1
                w_1 = phi(i,j1,k,4)*ir1
                ek1 = 0.5_rp*(u_1*u_1 + v_1*v_1 + w_1*w_1)
                p_1 = (gamma0-1.0_rp) * (phi(i,j1,k,5) - r_1*ek1)
                T_1 = p_1*ir1

                tauW_WR = mu_inf*mW*(u_1+w_1)/dyh
                qauW_WR = mu_inf*lW*(T_1-T_w)/dyh
                ! 
                ! get les quantities at the jInt node
                !
                r_h = phi(i,jl,k,1)
                irh = 1.0_rp/r_h
                u_h = phi(i,jl,k,2)*irh
                v_h = phi(i,jl,k,3)*irh
                w_h = phi(i,jl,k,4)*irh
                ekh = 0.5_rp*(u_h*u_h + v_h*v_h + w_h*w_h)
                p_h = (gamma0-1.0_rp) * (phi(i,jl,k,5) - r_h*ekh)
                T_h = p_h*irh

                u_p = sqrt(u_h*u_h + w_h*w_h)
                !
                ! estimate delta_nu with Reichardt's law and get the grid
                !
                hw = abs(y(jl)-yW)
                call LogLawWMLES(mu_inf,Tref,vis_flag,hw,T_w,u_p,p_h,tauWall0)
                rhoWall0 = p_h/T_w
                utaWall0 = sqrt(abs(tauWall0)/rhoWall0)
                dltWall0 = rhoWall0*utaWall0/(mW*mu_inf)

                !
                ! get the correct tauWall and qWall from ODE model
                !
                call OdeWMLES(gamma0,Prandtl,mu_inf,Tref,vis_flag,hw,dltWall0,u_w,u_p,T_w,T_h,p_h,u_wm,T_wm,&
                              tauW_WM,qauW_WM)

                !
                ! enforce the effective viscosity and diffusivity
                !
                mu_R = abs(tauW_WM/tauW_WR)
                lm_R = abs(qauW_WM/qauW_WR)

                if(mu_R >  clipPar) mu_R =  ClipPar
                if(lm_R >  clipPar) lm_R =  ClipPar

                do j = 1,GN
        
                   jg = b_node + b_norm*j
                   ji = b_node - b_norm*(j-1)

                   VIS(i,jg,k) = 2*mu_R*mW - VIS(i,ji,k)
                   LMD(i,jg,k) = 2*lm_R*lW - LMD(i,ji,k)

                enddo

                if(b_face==0) then
                    WMLES_DATA_LW(i,k,1) = mu_R
                    WMLES_DATA_LW(i,k,2) = lm_R
                    WMLES_DATA_LW(i,k,3) = mu_R*tauW_WR
                    WMLES_DATA_LW(i,k,4) = lm_R*qauW_WR
                    WMLES_DATA_LW(i,k,5) = 1.0_rp
                elseif(b_face==1) then
                    WMLES_DATA_UW(i,k,1) = mu_R
                    WMLES_DATA_UW(i,k,2) = lm_R
                    WMLES_DATA_UW(i,k,3) = mu_R*tauW_WR
                    WMLES_DATA_UW(i,k,4) = lm_R*qauW_WR
                    WMLES_DATA_UW(i,k,5) = 1.0_rp
                endif

                if(ic == 'turbulent_channel') then
                do l = 1,nvWMLESData
                   WMLES_DATA(i,k,l) = 0.5_rp*(WMLES_DATA_LW(i,k,l) + WMLES_DATA_UW(i,k,l))
                enddo
                endif

             enddo
          enddo
          !$acc end parallel

          case default
          print*, 'StaticWallStressModel not implemented for face ', trim(b%face)
          stop
        
        endselect
        


        
        return
end subroutine StaticDifferentialWallModelledTrbVis







end module sgs_module
