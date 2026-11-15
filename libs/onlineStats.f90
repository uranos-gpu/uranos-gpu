module onlineStats

use mpi
use statistics_module
use profiling_module, only: StartProfRange, EndProfRange
use thermodynamics_module, only: entropy_t, heat_capacity_t, gamma_Tref

implicit none

contains

subroutine compute_1DStats(phi,T,GAMMA_T,VIS,nx,nz,nVar,itStat,&
                          sx,ex,sy,ey,sz,ez,comm,vmean)
        use parameters_module     , only: rp, MPI_RP, Tref, vis_flag, cuda_aware
        use fluid_functions_module, only: laminar_viscosity
        use storage_module        , only: central_1, &
                                          stat_pln1D, stat_lcl1D, stat_gbl1D
        use mesh_module           , only: ystep_i
        use mpi_module            , only: lby,uby

        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(in)    :: phi
        real(rp), dimension(:,:,:)  , allocatable, intent(in)    :: T
        real(rp), dimension(:,:,:)  , allocatable, intent(in)    :: GAMMA_T
        real(rp), dimension(:,:,:)  , allocatable, intent(in)    :: VIS
        real(rp), dimension(:,:)    , allocatable, intent(inout) :: vmean

        integer , intent(in) :: comm
        integer , intent(in) :: nx,nz,nVar,sx,ex,sy,ey,sz,ez, itStat

        integer  :: err, npt
        real(rp) :: r_, u_, v_, w_, p_, T_, re, ir, mL, mT, g_, M2
        real(rp) :: cl, r_l, rul, u_l, r_y, u_y, ruy
        real(rp) :: irNxNz, r_itStat, r_itStatm1
        integer  :: i,j,k, l

        call StartProfRange('compute_1DStats')

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(2)
        do    j = sy,ey
           do i = sx,ex

              do l = 1,nVar
                 stat_pln1D(i,j,l) = 0.0_rp
              enddo

              do k = sz,ez

                 r_ = phi(i,j,k,1)
                 ir = 1.0_rp/r_
                 u_ = phi(i,j,k,2)*ir
                 v_ = phi(i,j,k,3)*ir
                 w_ = phi(i,j,k,4)*ir
                 re = phi(i,j,k,5)
                 T_ = T(i,j,k)
                 p_ = r_*T_

                 mL = laminar_viscosity(T_,Tref,vis_flag)
                 mT = VIS(i,j,k) - mL

                 stat_pln1D(i,j,1) = stat_pln1D(i,j,1) + r_
                 stat_pln1D(i,j,2) = stat_pln1D(i,j,2) + r_*u_
                 stat_pln1D(i,j,3) = stat_pln1D(i,j,3) + r_*v_
                 stat_pln1D(i,j,4) = stat_pln1D(i,j,4) + r_*w_
                 stat_pln1D(i,j,5) = stat_pln1D(i,j,5) + re

                 stat_pln1D(i,j,6) = stat_pln1D(i,j,6) + (r_*u_)*(r_*u_)*ir
                 stat_pln1D(i,j,7) = stat_pln1D(i,j,7) + (r_*v_)*(r_*v_)*ir
                 stat_pln1D(i,j,8) = stat_pln1D(i,j,8) + (r_*w_)*(r_*w_)*ir
                 stat_pln1D(i,j,9) = stat_pln1D(i,j,9) + (r_*u_)*(r_*v_)*ir

                 stat_pln1D(i,j,10) = stat_pln1D(i,j,10) + r_*r_

                 stat_pln1D(i,j,11) = stat_pln1D(i,j,11) + p_
                 stat_pln1D(i,j,12) = stat_pln1D(i,j,12) + p_*p_

                 stat_pln1D(i,j,13) = stat_pln1D(i,j,13) + T_
                 stat_pln1D(i,j,14) = stat_pln1D(i,j,14) + r_*T_*T_

                 stat_pln1D(i,j,15) = stat_pln1D(i,j,15) + mL
                 stat_pln1D(i,j,16) = stat_pln1D(i,j,16) + mT

                 r_y = 0.0_rp
                 ruy = 0.0_rp
                 u_y = 0.0_rp
                 do l = -3,3

                    cl   = central_1(l)

                    r_l = phi(i,j+l,k,1)
                    rul = phi(i,j+l,k,2)
                    u_l = rul/r_l

                    r_y = r_y + cl * r_l
                    u_y = u_y + cl * u_l
                    ruy = ruy + cl * rul

                 enddo
                 r_y = r_y*ystep_i(j)
                 u_y = u_y*ystep_i(j)
                 ruy = ruy*ystep_i(j)

                 stat_pln1D(i,j,17) = stat_pln1D(i,j,17) + r_y
                 stat_pln1D(i,j,18) = stat_pln1D(i,j,18) + ruy

                 stat_pln1D(i,j,19) = stat_pln1D(i,j,19) + mL*u_y
                 stat_pln1D(i,j,20) = stat_pln1D(i,j,20) + mT*u_y

                 stat_pln1D(i,j,21) = stat_pln1D(i,j,21) + (r_*u_*u_*u_)
                 stat_pln1D(i,j,22) = stat_pln1D(i,j,22) + (r_*v_*v_*v_)
                 stat_pln1D(i,j,23) = stat_pln1D(i,j,23) + (r_*w_*w_*w_)
                 stat_pln1D(i,j,24) = stat_pln1D(i,j,24) + (r_*u_*u_*v_)

                 stat_pln1D(i,j,25) = stat_pln1D(i,j,25) + (r_*u_*u_*u_*u_)
                 stat_pln1D(i,j,26) = stat_pln1D(i,j,26) + (r_*v_*v_*v_*v_)
                 stat_pln1D(i,j,27) = stat_pln1D(i,j,27) + (r_*w_*w_*w_*w_)

                 g_ = GAMMA_T(i,j,k)
                 stat_pln1D(i,j,28) = stat_pln1D(i,j,28) + g_/(g_ - 1.0_rp)
                 stat_pln1D(i,j,29) = stat_pln1D(i,j,29) + g_
                 stat_pln1D(i,j,30) = stat_pln1D(i,j,30) + mL*ir

                 stat_pln1D(i,j,31) = stat_pln1D(i,j,31) + r_*u_*T_
                 stat_pln1D(i,j,32) = stat_pln1D(i,j,32) + r_*v_*T_
                 stat_pln1D(i,j,33) = stat_pln1D(i,j,33) + T_*T_

                 M2 = (u_*u_ + v_*v_ + w_*w_)/(g_*T_)
                 stat_pln1D(i,j,34) = stat_pln1D(i,j,34) + sqrt(M2)
                 stat_pln1D(i,j,35) = stat_pln1D(i,j,35) + M2

              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(2)
        do    l = 1,nVar
           do j = sy,ey
              stat_lcl1D(j,l) = 0.0_rp
              do i = sx,ex
                 stat_lcl1D(j,l) = stat_lcl1D(j,l) + stat_pln1D(i,j,l)
              enddo
           enddo
        enddo
        !$acc end parallel

        npt = nVar*(uby-lby+1)
        if(cuda_aware) then
          !$acc host_data use_device(stat_lcl1D, stat_gbl1D)
          call mpi_allreduce(stat_lcl1D,stat_gbl1D,npt,MPI_RP,mpi_sum,comm,err)
          !$acc end host_data
        else
          !$acc update host(stat_lcl1D)
          call mpi_allreduce(stat_lcl1D,stat_gbl1D,npt,MPI_RP,mpi_sum,comm,err)
          !$acc update device(stat_gbl1D)
        endif

        irNxNz     = 1.0_rp/real(nx*nz,rp)
        r_itStat   = real(itStat,rp)
        r_itStatm1 = real(itStat-1,rp)
        !$acc parallel default(present)
        !$acc loop gang, vector collapse(2)
        do    l = 1,nVar
           do j = sy,ey
              vmean(j,l) = (vmean(j,l)*r_itStatm1 + stat_gbl1D(j,l)*irNxNz)/r_itStat
           enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine compute_1DStats

subroutine compute_2DStats(phi,U,V,W,T,GAMMA_T,VIS,LMD,isoSGS,nz,nVar,itStat, &
                           sx,ex,sy,ey,sz,ez,comm,vmean)
        use parameters_module     , only: rp, MPI_RP, gamma0, Prandtl, &
                                          hybrid_weno, Tref, vis_flag, thermally_perfect_gas, cuda_aware
        use fluid_functions_module, only: laminar_viscosity
        use mpi_module            , only: lbx,ubx, lby,uby
        use storage_module        , only: central_1, central_2, SSENSOR, weno_flag, &
                                          stat_lcl2D, stat_gbl2D
        use mesh_module           , only: zstep_i, dxi_dx, dxi_dy, deta_dx, deta_dy, &
                                          Gamma1_2D, Gamma2_2D, ztastep_i, iz_zta, z_zta2

        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(in)    :: phi
        real(rp), dimension(:,:,:)  , allocatable, intent(in)    :: U,V,W,T,GAMMA_T
        real(rp), dimension(:,:,:)  , allocatable, intent(in)    :: VIS, LMD, isoSGS
        real(rp), dimension(:,:,:)  , allocatable, intent(inout) :: vmean

        integer , intent(in) :: comm, itStat
        integer , intent(in) :: nz,nVar,sx,ex,sy,ey,sz,ez

        real(rp), parameter :: two3 = 2.0_rp/3.0_rp
        integer             :: err, npt

        real(rp) :: r_, u_, v_, w_, p_, T_, re, ir
        real(rp) :: gamma_, c_, ve
        real(rp) :: mL, mT

        real(rp) :: irNz, r_itStat, r_itStatm1

        real(rp) :: cl, cl2
        real(rp) :: idz
        real(rp) :: csi_x, csi_y, eta_x, eta_y
        real(rp) :: g11, g12, g22, gam1, gam2
        real(rp) :: dc, de
        real(rp) :: Tc, Te, Tz_b
        real(rp) :: d2Tce, clm, iz_z
        real(rp) :: r_x, r_y, r_z
        real(rp) :: rux, ruy, ruz
        real(rp) :: rvx, rvy, rvz
        real(rp) :: rwx, rwy, rwz
        real(rp) :: u_x, u_y, u_z
        real(rp) :: v_x, v_y, v_z
        real(rp) :: w_x, w_y, w_z
        real(rp) :: T_x, T_y, T_z
        real(rp) :: T2x, T2y, T2z
        real(rp) :: l_x, l_y, l_z

        real(rp) :: tau_11, tau_12, tau_13
        real(rp) :: tau_21, tau_22, tau_23
        real(rp) :: tau_31, tau_32, tau_33
        real(rp) :: div_

        real(rp) :: sensor, weno_flag_

        real(rp) :: entropy0
        real(rp) :: phi_visL, phi_visT
        real(rp) :: d2T

        integer  :: i, j, k, l, m

        call StartProfRange('compute_2DStats')

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(2)
        do    j = sy,ey
           do i = sx,ex

              do l = 1,nVar
                 stat_lcl2D(i,j,l) = 0.0_rp
              enddo

              do k = sz,ez

                 r_     = phi(i,j,k,1)
                 ir     = 1.0_rp/r_
                 u_     = phi(i,j,k,2)*ir
                 v_     = phi(i,j,k,3)*ir
                 w_     = phi(i,j,k,4)*ir
                 re     = phi(i,j,k,5)
                 T_     = T(i,j,k)
                 p_     = r_*T_
                 gamma_ = GAMMA_T(i,j,k)
                 c_     = sqrt(gamma_*T_)
                 ve     = sqrt(u_*u_ + v_*v_ + w_*w_)

                 mL = laminar_viscosity(T_,Tref,vis_flag)
                 mT = VIS(i,j,k) - mL

                 stat_lcl2D(i,j,1) = stat_lcl2D(i,j,1) + r_
                 stat_lcl2D(i,j,2) = stat_lcl2D(i,j,2) + r_*u_
                 stat_lcl2D(i,j,3) = stat_lcl2D(i,j,3) + r_*v_
                 stat_lcl2D(i,j,4) = stat_lcl2D(i,j,4) + r_*w_
                 stat_lcl2D(i,j,5) = stat_lcl2D(i,j,5) + re

                 stat_lcl2D(i,j,6) = stat_lcl2D(i,j,6) + (r_*u_)*(r_*u_)*ir
                 stat_lcl2D(i,j,7) = stat_lcl2D(i,j,7) + (r_*v_)*(r_*v_)*ir
                 stat_lcl2D(i,j,8) = stat_lcl2D(i,j,8) + (r_*w_)*(r_*w_)*ir
                 stat_lcl2D(i,j,9) = stat_lcl2D(i,j,9) + (r_*u_)*(r_*v_)*ir

                 stat_lcl2D(i,j,10) = stat_lcl2D(i,j,10) + r_*r_

                 stat_lcl2D(i,j,11) = stat_lcl2D(i,j,11) + p_
                 stat_lcl2D(i,j,12) = stat_lcl2D(i,j,12) + p_*p_

                 stat_lcl2D(i,j,13) = stat_lcl2D(i,j,13) + T_
                 stat_lcl2D(i,j,14) = stat_lcl2D(i,j,14) + r_*T_*T_

                 stat_lcl2D(i,j,15) = stat_lcl2D(i,j,15) + mL
                 stat_lcl2D(i,j,16) = stat_lcl2D(i,j,16) + mT

                 stat_lcl2D(i,j,17) = stat_lcl2D(i,j,17) + r_*T_

                 stat_lcl2D(i,j,18) = stat_lcl2D(i,j,18) + ve/c_
                 stat_lcl2D(i,j,19) = stat_lcl2D(i,j,19) + r_*ve/c_

                 csi_x = dxi_dx (i,j) ; csi_y = dxi_dy (i,j)
                 eta_x = deta_dx(i,j) ; eta_y = deta_dy(i,j)
                 idz   = zstep_i(k)

                 r_x = 0.0_rp ; rux = 0.0_rp ; rvx = 0.0_rp ; rwx = 0.0_rp
                 r_y = 0.0_rp ; ruy = 0.0_rp ; rvy = 0.0_rp ; rwy = 0.0_rp
                 r_z = 0.0_rp ; ruz = 0.0_rp ; rvz = 0.0_rp ; rwz = 0.0_rp
                 u_x = 0.0_rp ; u_y = 0.0_rp ; u_z = 0.0_rp
                 v_x = 0.0_rp ; v_y = 0.0_rp ; v_z = 0.0_rp
                 w_x = 0.0_rp ; w_y = 0.0_rp ; w_z = 0.0_rp
                 T_x = 0.0_rp ; T_y = 0.0_rp ; T_z = 0.0_rp
                 T2x = 0.0_rp ; T2y = 0.0_rp ; T2z = 0.0_rp
                 l_x = 0.0_rp ; l_y = 0.0_rp ; l_z = 0.0_rp

                 do l = -3,3

                    cl  = central_1(l)
                    cl2 = central_2(l)

                    r_x = r_x + cl*phi(i+l,j,k,1)
                    rux = rux + cl*phi(i+l,j,k,2)
                    rvx = rvx + cl*phi(i+l,j,k,3)
                    rwx = rwx + cl*phi(i+l,j,k,4)

                    r_y = r_y + cl*phi(i,j+l,k,1)
                    ruy = ruy + cl*phi(i,j+l,k,2)
                    rvy = rvy + cl*phi(i,j+l,k,3)
                    rwy = rwy + cl*phi(i,j+l,k,4)

                    r_z = r_z + cl*phi(i,j,k+l,1)
                    ruz = ruz + cl*phi(i,j,k+l,2)
                    rvz = rvz + cl*phi(i,j,k+l,3)
                    rwz = rwz + cl*phi(i,j,k+l,4)

                    u_x = u_x + cl*U(i+l,j,k)
                    u_y = u_y + cl*U(i,j+l,k)
                    u_z = u_z + cl*U(i,j,k+l)

                    v_x = v_x + cl*V(i+l,j,k)
                    v_y = v_y + cl*V(i,j+l,k)
                    v_z = v_z + cl*V(i,j,k+l)

                    w_x = w_x + cl*W(i+l,j,k)
                    w_y = w_y + cl*W(i,j+l,k)
                    w_z = w_z + cl*W(i,j,k+l)

                    T_x = T_x + cl*T(i+l,j,k)
                    T_y = T_y + cl*T(i,j+l,k)
                    T_z = T_z + cl*T(i,j,k+l)

                    T2x = T2x + cl2*T(i+l,j,k)
                    T2y = T2y + cl2*T(i,j+l,k)
                    T2z = T2z + cl2*T(i,j,k+l)

                    l_x = l_x + cl*LMD(i+l,j,k)
                    l_y = l_y + cl*LMD(i,j+l,k)
                    l_z = l_z + cl*LMD(i,j,k+l)

                 enddo

                 d2Tce = 0.0_rp
                 do m = 1,3
                    do l = 1,3
                       clm = central_1(l)*central_1(m)
                       d2Tce = d2Tce + clm*(T(i+l,j+m,k) - T(i-l,j+m,k) &
                                          -  T(i+l,j-m,k) + T(i-l,j-m,k))
                    enddo
                 enddo

                 dc = r_x  ; de = r_y
                 r_x = csi_x*dc + eta_x*de ; r_y = csi_y*dc + eta_y*de ; r_z = r_z*idz
                 dc = rux  ; de = ruy
                 rux = csi_x*dc + eta_x*de ; ruy = csi_y*dc + eta_y*de ; ruz = ruz*idz
                 dc = rvx  ; de = rvy
                 rvx = csi_x*dc + eta_x*de ; rvy = csi_y*dc + eta_y*de ; rvz = rvz*idz
                 dc = rwx  ; de = rwy
                 rwx = csi_x*dc + eta_x*de ; rwy = csi_y*dc + eta_y*de ; rwz = rwz*idz
                 dc = u_x  ; de = u_y
                 u_x = csi_x*dc + eta_x*de ; u_y = csi_y*dc + eta_y*de ; u_z = u_z*idz
                 dc = v_x  ; de = v_y
                 v_x = csi_x*dc + eta_x*de ; v_y = csi_y*dc + eta_y*de ; v_z = v_z*idz
                 dc = w_x  ; de = w_y
                 w_x = csi_x*dc + eta_x*de ; w_y = csi_y*dc + eta_y*de ; w_z = w_z*idz
                 dc = l_x  ; de = l_y
                 l_x = csi_x*dc + eta_x*de ; l_y = csi_y*dc + eta_y*de ; l_z = l_z*idz

                 Tc = T_x ; Te = T_y ; Tz_b = T_z
                 T_x = csi_x*Tc + eta_x*Te
                 T_y = csi_y*Tc + eta_y*Te
                 T_z = Tz_b*idz

                 g11  = csi_x*csi_x + csi_y*csi_y
                 g22  = eta_x*eta_x + eta_y*eta_y
                 g12  = csi_x*eta_x + csi_y*eta_y
                 gam1 = Gamma1_2D(i,j)
                 gam2 = Gamma2_2D(i,j)
                 iz_z = iz_zta(k)
                 T2z  = T2z*idz*idz - (ztastep_i(k)*Tz_b)*z_zta2(k)*iz_z**3

                 stat_lcl2D(i,j,20) = stat_lcl2D(i,j,20) + r_x
                 stat_lcl2D(i,j,21) = stat_lcl2D(i,j,21) + r_y
                 stat_lcl2D(i,j,22) = stat_lcl2D(i,j,22) + r_z

                 stat_lcl2D(i,j,23) = stat_lcl2D(i,j,23) + rux
                 stat_lcl2D(i,j,24) = stat_lcl2D(i,j,24) + ruy
                 stat_lcl2D(i,j,25) = stat_lcl2D(i,j,25) + ruz

                 stat_lcl2D(i,j,26) = stat_lcl2D(i,j,26) + rvx
                 stat_lcl2D(i,j,27) = stat_lcl2D(i,j,27) + rvy
                 stat_lcl2D(i,j,28) = stat_lcl2D(i,j,28) + rvz

                 stat_lcl2D(i,j,29) = stat_lcl2D(i,j,29) + rwx
                 stat_lcl2D(i,j,30) = stat_lcl2D(i,j,30) + rwy
                 stat_lcl2D(i,j,31) = stat_lcl2D(i,j,31) + rwz

                 div_ = u_x + v_y + w_z

                 tau_11 = mL*(u_x + u_x - two3 * div_)
                 tau_12 = mL*(u_y + v_x)
                 tau_13 = mL*(u_z + w_x)

                 tau_21 = mL*(v_x + u_y)
                 tau_22 = mL*(v_y + v_y - two3 * div_)
                 tau_23 = mL*(v_z + w_y)

                 tau_31 = mL*(w_x + u_z)
                 tau_32 = mL*(w_y + v_z)
                 tau_33 = mL*(w_z + w_z - two3 * div_)

                 phi_visL = tau_11 * u_x + tau_12 * u_y + tau_13 * u_z + &
                            tau_21 * v_x + tau_22 * v_y + tau_23 * v_z + &
                            tau_31 * w_x + tau_32 * w_y + tau_33 * w_z

                 phi_visT = phi_visL * mT/mL

                 stat_lcl2D(i,j,37) = stat_lcl2D(i,j,37) + phi_visL
                 stat_lcl2D(i,j,38) = stat_lcl2D(i,j,38) + phi_visT

                 stat_lcl2D(i,j,39) = stat_lcl2D(i,j,39) + phi_visL/T_
                 stat_lcl2D(i,j,40) = stat_lcl2D(i,j,40) + phi_visT/T_

                 d2T = g11*T2x + 2.0_rp*g12*d2Tce + g22*T2y + gam1*Tc + gam2*Te + T2z
                 stat_lcl2D(i,j,44) = stat_lcl2D(i,j,44) + LMD(i,j,k)
                 stat_lcl2D(i,j,45) = stat_lcl2D(i,j,45) + d2T

                 stat_lcl2D(i,j,46) = stat_lcl2D(i,j,46) + l_x
                 stat_lcl2D(i,j,47) = stat_lcl2D(i,j,47) + l_y
                 stat_lcl2D(i,j,48) = stat_lcl2D(i,j,48) + l_z

                 stat_lcl2D(i,j,49) = stat_lcl2D(i,j,49) + T_x
                 stat_lcl2D(i,j,50) = stat_lcl2D(i,j,50) + T_y
                 stat_lcl2D(i,j,51) = stat_lcl2D(i,j,51) + T_z

                 stat_lcl2D(i,j,52) = stat_lcl2D(i,j,52) + LMD(i,j,k) * d2T / T_
                 stat_lcl2D(i,j,53) = stat_lcl2D(i,j,53) + (l_x * T_x + l_y * T_y + l_z * T_z) / T_

                 stat_lcl2D(i,j, 95) = stat_lcl2D(i,j, 95) + LMD(i,j,k)*(T_x*T_x + T_y*T_y + T_z*T_z)/(T_*T_)
                 stat_lcl2D(i,j, 96) = stat_lcl2D(i,j, 96) + LMD(i,j,k)*T_x/T_
                 stat_lcl2D(i,j, 97) = stat_lcl2D(i,j, 97) + LMD(i,j,k)*T_y/T_
                 stat_lcl2D(i,j, 98) = stat_lcl2D(i,j, 98) + LMD(i,j,k)*T_z/T_

                 stat_lcl2D(i,j,54) = stat_lcl2D(i,j,54) + isoSGS(i,j,k)

                 stat_lcl2D(i,j,55) = stat_lcl2D(i,j,55) + tau_11
                 stat_lcl2D(i,j,56) = stat_lcl2D(i,j,56) + tau_12
                 stat_lcl2D(i,j,57) = stat_lcl2D(i,j,57) + tau_13
                 stat_lcl2D(i,j,58) = stat_lcl2D(i,j,58) + tau_21
                 stat_lcl2D(i,j,59) = stat_lcl2D(i,j,59) + tau_22
                 stat_lcl2D(i,j,60) = stat_lcl2D(i,j,60) + tau_23
                 stat_lcl2D(i,j,61) = stat_lcl2D(i,j,61) + tau_31
                 stat_lcl2D(i,j,62) = stat_lcl2D(i,j,62) + tau_32
                 stat_lcl2D(i,j,63) = stat_lcl2D(i,j,63) + tau_33

                 stat_lcl2D(i,j,64) = stat_lcl2D(i,j,64) + u_
                 stat_lcl2D(i,j,65) = stat_lcl2D(i,j,65) + v_
                 stat_lcl2D(i,j,66) = stat_lcl2D(i,j,66) + w_

                 stat_lcl2D(i,j,67) = stat_lcl2D(i,j,67) + tau_11*u_ + tau_21*v_ + tau_31*w_
                 stat_lcl2D(i,j,68) = stat_lcl2D(i,j,68) + tau_12*u_ + tau_22*v_ + tau_32*w_

                 stat_lcl2D(i,j,69) = stat_lcl2D(i,j,69) + tau_11*u_x + tau_12*u_y + tau_13*u_z &
                                                   + tau_21*v_x + tau_22*v_y + tau_23*v_z &
                                                   + tau_31*w_x + tau_32*w_y + tau_33*w_z

                 stat_lcl2D(i,j,70) = stat_lcl2D(i,j,70) + u_x
                 stat_lcl2D(i,j,71) = stat_lcl2D(i,j,71) + u_y
                 stat_lcl2D(i,j,72) = stat_lcl2D(i,j,72) + u_z
                 stat_lcl2D(i,j,73) = stat_lcl2D(i,j,73) + v_x
                 stat_lcl2D(i,j,74) = stat_lcl2D(i,j,74) + v_y
                 stat_lcl2D(i,j,75) = stat_lcl2D(i,j,75) + v_z
                 stat_lcl2D(i,j,76) = stat_lcl2D(i,j,76) + w_x
                 stat_lcl2D(i,j,77) = stat_lcl2D(i,j,77) + w_y
                 stat_lcl2D(i,j,78) = stat_lcl2D(i,j,78) + w_z

                 stat_lcl2D(i,j,79) = stat_lcl2D(i,j,79) + p_*div_
                 stat_lcl2D(i,j,80) = stat_lcl2D(i,j,80) + p_*u_
                 stat_lcl2D(i,j,81) = stat_lcl2D(i,j,81) + p_*v_

                 stat_lcl2D(i,j,82) = stat_lcl2D(i,j,82) + r_*u_*u_*u_
                 stat_lcl2D(i,j,83) = stat_lcl2D(i,j,83) + r_*v_*v_*v_
                 stat_lcl2D(i,j,84) = stat_lcl2D(i,j,84) + r_*w_*w_*u_
                 stat_lcl2D(i,j,85) = stat_lcl2D(i,j,85) + r_*w_*w_*v_
                 stat_lcl2D(i,j,86) = stat_lcl2D(i,j,86) + r_*u_*u_*v_
                 stat_lcl2D(i,j,87) = stat_lcl2D(i,j,87) + r_*v_*v_*u_

                 stat_lcl2D(i,j,88) = stat_lcl2D(i,j,88) + (r_*u_)*(r_*w_)*ir
                 stat_lcl2D(i,j,89) = stat_lcl2D(i,j,89) + (r_*v_)*(r_*w_)*ir

                 stat_lcl2D(i,j,90) = stat_lcl2D(i,j,90) + r_*T_*T_*T_
                 stat_lcl2D(i,j,91) = stat_lcl2D(i,j,91) + r_*T_*T_*T_*T_
                 stat_lcl2D(i,j,92) = stat_lcl2D(i,j,92) + r_*u_*u_*u_*u_

                 stat_lcl2D(i,j,93) = stat_lcl2D(i,j,93) + gamma_
                 stat_lcl2D(i,j,94) = stat_lcl2D(i,j,94) + r_ * gamma_

              enddo
           enddo
        enddo
        !$acc end parallel

        if(hybrid_weno) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(2)
          do    j = sy,ey
             do i = sx,ex
                do k = sz,ez
                   sensor = SSENSOR(i,j,k)
                   stat_lcl2D(i,j,32) = stat_lcl2D(i,j,32) + sensor
                   stat_lcl2D(i,j,33) = stat_lcl2D(i,j,33) + sensor*sensor

                   weno_flag_ = real(weno_flag(i,j,k),rp)
                   stat_lcl2D(i,j,34) = stat_lcl2D(i,j,34) + weno_flag_
                   stat_lcl2D(i,j,35) = stat_lcl2D(i,j,35) + weno_flag_*weno_flag_
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(thermally_perfect_gas > 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(2)
          do    j = sy,ey
             do i = sx,ex
                do k = sz,ez
                   r_ = phi(i,j,k,1)
                   ir = 1.0_rp/r_
                   T_ = T(i,j,k)
                   p_ = r_*T_
                   entropy0 = entropy_t(T_) - log(p_)
                   stat_lcl2D(i,j, 36) = stat_lcl2D(i,j, 36) + r_*entropy0
                   stat_lcl2D(i,j, 41) = stat_lcl2D(i,j, 41) + r_ * entropy0 * (phi(i,j,k,2)*ir)
                   stat_lcl2D(i,j, 42) = stat_lcl2D(i,j, 42) + r_ * entropy0 * (phi(i,j,k,3)*ir)
                   stat_lcl2D(i,j, 43) = stat_lcl2D(i,j, 43) + r_ * entropy0 * (phi(i,j,k,4)*ir)
                   stat_lcl2D(i,j, 99) = stat_lcl2D(i,j, 99) + r_*entropy0*entropy0
                   stat_lcl2D(i,j,100) = stat_lcl2D(i,j,100) + entropy0
                enddo
             enddo
          enddo
          !$acc end parallel
        else
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(2)
          do    j = sy,ey
             do i = sx,ex
                do k = sz,ez
                   r_ = phi(i,j,k,1)
                   ir = 1.0_rp/r_
                   T_ = T(i,j,k)
                   p_ = r_*T_
                   entropy0 = 1.0_rp/(gamma0 - 1.0_rp) * log(p_*ir**(gamma0))
                   stat_lcl2D(i,j, 36) = stat_lcl2D(i,j, 36) + r_*entropy0
                   stat_lcl2D(i,j, 41) = stat_lcl2D(i,j, 41) + r_ * entropy0 * (phi(i,j,k,2)*ir)
                   stat_lcl2D(i,j, 42) = stat_lcl2D(i,j, 42) + r_ * entropy0 * (phi(i,j,k,3)*ir)
                   stat_lcl2D(i,j, 43) = stat_lcl2D(i,j, 43) + r_ * entropy0 * (phi(i,j,k,4)*ir)
                   stat_lcl2D(i,j, 99) = stat_lcl2D(i,j, 99) + r_*entropy0*entropy0
                   stat_lcl2D(i,j,100) = stat_lcl2D(i,j,100) + entropy0
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        npt = nVar*(ubx-lbx+1)*(uby-lby+1)
        if(cuda_aware) then
          !$acc host_data use_device(stat_lcl2D, stat_gbl2D)
          call mpi_allreduce(stat_lcl2D,stat_gbl2D,npt,MPI_RP,mpi_sum,comm,err)
          !$acc end host_data
        else
          !$acc update host(stat_lcl2D)
          call mpi_allreduce(stat_lcl2D,stat_gbl2D,npt,MPI_RP,mpi_sum,comm,err)
          !$acc update device(stat_gbl2D)
        endif

        irNz       = 1.0_rp/real(nz,rp)
        r_itStat   = real(itStat,rp)
        r_itStatm1 = real(itStat-1,rp)
        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do       l = 1,nVar
           do    j = sy,ey
              do i = sx,ex
                 vmean(i,j,l) = (vmean(i,j,l)*r_itStatm1 + stat_gbl2D(i,j,l)*irNz)/r_itStat
              enddo
           enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine compute_2DStats

subroutine compute_3DStats(phi,U,V,W,T,sx,ex,sy,ey,sz,ez,comm,mpi_flag,vol,Ek,Ew,Mmax,rhoS,rhoE)
        use parameters_module    , only: rp, gamma0, cv, thermally_perfect_gas, MPI_RP
        use storage_module       , only: central_1, central_fd_order
        use mesh_module          , only: Jac2D, dxi_dx, dxi_dy, deta_dx, deta_dy, zstep, zstep_i

        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(in)  :: phi
        real(rp), dimension(:,:,:)  , allocatable, intent(in)  :: U,V,W,T

        integer , intent(in)  :: sx,ex,sy,ey,sz,ez,comm
        logical , intent(in)  :: mpi_flag
        real(rp), intent(out) :: vol, Ek, Ew, Mmax, rhoS, rhoE

        real(rp), dimension(5) :: lcl, gbl
        real(rp)               :: s_vol, s_ke, s_ens, s_rs, s_re, lcl_M

        real(rp) :: r_, u_, v_, w_, T_, p_, s_, ve2, c2, dv, cp_

        real(rp) :: du_xi, du_et, du_zt
        real(rp) :: dv_xi, dv_et, dv_zt
        real(rp) :: dw_xi, dw_et

        real(rp) :: uy, uz, vx, vz, wx, wy
        real(rp) :: dxidx, dxidy, detadx, detady, izt, cl
        real(rp) :: ox, oy, oz

        integer  :: i, j, k, l, err = 0

        s_vol = 0.0_rp
        s_ke  = 0.0_rp
        s_ens = 0.0_rp
        s_rs  = 0.0_rp
        s_re  = 0.0_rp
        lcl_M = 0.0_rp
        !$acc parallel default(present) copy(s_vol,s_ke,s_ens,s_rs,s_re,lcl_M)
        !$acc loop gang, vector collapse(3) &
        !$acc   reduction(+:s_vol,s_ke,s_ens,s_rs,s_re) reduction(max:lcl_M)
        do       k = sz,ez
           do    j = sy,ey
              do i = sx,ex

                 r_  = phi(i,j,k,1)
                 u_  = U(i,j,k)
                 v_  = V(i,j,k)
                 w_  = W(i,j,k)
                 T_  = T(i,j,k)
                 p_  = r_*T_
                 ve2 = u_*u_ + v_*v_ + w_*w_
                 if(thermally_perfect_gas > 0) then
                    cp_ = heat_capacity_t(T_)
                    c2  = cp_/(cp_ - 1.0_rp)*T_
                 else
                    c2  = gamma_Tref*T_
                 endif

                 if(thermally_perfect_gas > 0) then
                    s_ = entropy_t(T_) - log(p_)
                 else
                    s_ = cv*log(p_/(r_**gamma0))
                 endif

                 du_xi = 0.0_rp; du_et = 0.0_rp; du_zt = 0.0_rp
                 dv_xi = 0.0_rp; dv_et = 0.0_rp; dv_zt = 0.0_rp
                 dw_xi = 0.0_rp; dw_et = 0.0_rp
                 do l = 1, central_fd_order/2
                    cl    = central_1(l)
                    du_xi = du_xi + cl*(U(i+l,j,k) - U(i-l,j,k))
                    du_et = du_et + cl*(U(i,j+l,k) - U(i,j-l,k))
                    du_zt = du_zt + cl*(U(i,j,k+l) - U(i,j,k-l))
                    dv_xi = dv_xi + cl*(V(i+l,j,k) - V(i-l,j,k))
                    dv_et = dv_et + cl*(V(i,j+l,k) - V(i,j-l,k))
                    dv_zt = dv_zt + cl*(V(i,j,k+l) - V(i,j,k-l))
                    dw_xi = dw_xi + cl*(W(i+l,j,k) - W(i-l,j,k))
                    dw_et = dw_et + cl*(W(i,j+l,k) - W(i,j-l,k))
                 enddo

                 dxidx  = dxi_dx (i,j); dxidy  = dxi_dy (i,j)
                 detadx = deta_dx(i,j); detady = deta_dy(i,j)
                 izt    = zstep_i(k)

                 uy = dxidy*du_xi + detady*du_et
                 uz = izt  *du_zt
                 vx = dxidx*dv_xi + detadx*dv_et
                 vz = izt  *dv_zt
                 wx = dxidx*dw_xi + detadx*dw_et
                 wy = dxidy*dw_xi + detady*dw_et

                 ox = wy - vz
                 oy = uz - wx
                 oz = vx - uy

                 dv = zstep(k)/Jac2D(i,j)

                 s_vol = s_vol + dv
                 s_ke  = s_ke  + 0.5_rp*r_*ve2*dv
                 s_ens = s_ens + 0.5_rp*r_*(ox*ox + oy*oy + oz*oz)*dv
                 s_rs  = s_rs  + r_*s_*dv
                 s_re  = s_re  + phi(i,j,k,5)*dv
                 lcl_M = max(lcl_M, sqrt(ve2/c2))

              enddo
           enddo
        enddo
        !$acc end parallel
        lcl = (/s_vol, s_ke, s_ens, s_rs, s_re/)

        if(mpi_flag) then
          call mpi_allreduce(lcl,   gbl,  5, MPI_RP, mpi_sum, comm, err)
          call mpi_allreduce(lcl_M, Mmax, 1, MPI_RP, mpi_max, comm, err)
        else
          gbl  = lcl
          Mmax = lcl_M
        endif

        vol = gbl(1)
        Ek  = gbl(2)/vol
        Ew  = gbl(3)/vol
        rhoS = gbl(4)/vol
        rhoE = gbl(5)/vol

        return
end subroutine compute_3DStats

end module onlineStats
