module hybrid_wenoEP_PEPTP_module
use parameters_module
use mpi_module
use storage_module
use profiling_module
use mesh_module               , only: dx_deta, dy_deta, dx_dxi, dy_dxi, Jac2D
use euler_matrices_module     , only: euler_system_n
use weno_reconstruction_module, only: weno_reconstruction
use thermodynamics_module     , only: cv_Tref, T_switch, T_cold, cp_coef, aL_adim__0, aH_adim__0
implicit none
private
public :: hybrid_wenoEP_PEPTP_X, hybrid_wenoEP_PEPTP_Y, hybrid_wenoEP_PEPTP_Z

contains

pure function logmean_pep(a,b) result(lm)
!$acc routine seq
        implicit none
        real(rp), intent(in) :: a,b
        real(rp) :: lm
        real(rp) :: d, f, s, u

        s = a + b
        d = b - a
        f = d/s
        u = f*f

        if(u < 1.0e-2_rp) then
          lm = 0.5_rp*s/(1.0_rp + u*(1.0_rp/3.0_rp + u*(1.0_rp/5.0_rp + u*(1.0_rp/7.0_rp &
             + u*(1.0_rp/9.0_rp + u*(1.0_rp/11.0_rp + u*(1.0_rp/13.0_rp + u/15.0_rp)))))))
        else
          lm = d/(log(b) - log(a))
        endif
end function logmean_pep

pure subroutine pep_tp_pair_core(rho_l,rho_r,u_l,u_r,v_l,v_r,w_l,w_r, &
                                 p_l,p_r,t_l,t_r,eint_l,eint_r,vel_l,vel_r, &
                                 t_l2,t_r2,t_l3,t_r3,t_l4,t_r4, &
                                 frho,frhoe,prhoE,k_m)
!$acc routine seq
        implicit none
        real(rp), intent(in)  :: rho_l, rho_r
        real(rp), intent(in)  :: u_l, u_r, v_l, v_r, w_l, w_r
        real(rp), intent(in)  :: p_l, p_r, t_l, t_r
        real(rp), intent(in)  :: eint_l, eint_r, vel_l, vel_r
        real(rp), intent(in)  :: t_l2, t_r2, t_l3, t_r3, t_l4, t_r4
        real(rp), intent(out) :: frho, frhoe, prhoE, k_m

        real(rp) :: b_m2, b_m1, b0, b1, b2, b3, b4
        real(rp) :: t_lr, t_log
        real(rp) :: cv_l, cv_r, n_fact, d_fact, n_inside, d_inside, inv_t
        real(rp) :: rho_harm, vel_m, alpha_m, lambda_l, lambda_r, psi_m
        integer  :: kb

        if(aL_adim__0 == 0.0_rp .and. aH_adim__0 == 0.0_rp) then
          b_m2 = 0.0_rp
          b_m1 = 0.0_rp
          b0   = cv_Tref
          b1   = 0.0_rp
          b2   = 0.0_rp
          b3   = 0.0_rp
          b4   = 0.0_rp
        else
          if(sqrt(t_l*t_r) <= T_cold) then
            kb = 0
          else if(sqrt(t_l*t_r) <= T_switch) then
            kb = 1
          else
            kb = 2
          endif
          b_m2 = cp_coef(-2,kb)
          b_m1 = cp_coef(-1,kb)
          b0   = cp_coef( 0,kb) - 1.0_rp
          b1   = cp_coef( 1,kb)
          b2   = cp_coef( 2,kb)
          b3   = cp_coef( 3,kb)
          b4   = cp_coef( 4,kb)
        endif

        t_lr = t_l*t_r
        t_log = logmean_pep(t_l,t_r)

        cv_l = b0 + b1*t_l + b2*t_l2 + b3*t_l3 + b4*t_l4 + b_m1/t_l + b_m2/t_l2
        cv_r = b0 + b1*t_r + b2*t_r2 + b3*t_r3 + b4*t_r4 + b_m1/t_r + b_m2/t_r2

        vel_m = 0.5_rp*(vel_l + vel_r)
        prhoE = 0.5_rp*(vel_l*p_r + vel_r*p_l)

        n_inside = b1 + b2*(t_r + t_l) + b3*(t_r2 + t_r*t_l + t_l2) &
                 + b4*(t_r3 + t_r2*t_l + t_r*t_l2 + t_l3) &
                 - b_m1/t_lr - b_m2*(t_r + t_l)/(t_r2*t_l2)

        d_inside = -2.0_rp*b_m2/t_lr - b_m1/t_log &
                 + 0.5_rp*b1*(t_r + t_l) &
                 + (2.0_rp/3.0_rp)*b2*(t_r2 + t_r*t_l + t_l2) &
                 + 0.75_rp*b3*(t_r3 + t_r2*t_l + t_r*t_l2 + t_l3) &
                 + 0.8_rp*b4*(t_r4 + t_r3*t_l + t_r2*t_l2 + t_r*t_l3 + t_l4)

        n_fact = -n_inside
        d_fact = -d_inside

        rho_harm = 2.0_rp*rho_l*rho_r/(rho_l + rho_r)
        if(d_fact /= 0.0_rp) then
          inv_t = n_fact/d_fact
          inv_t = min(max(inv_t, 1.0_rp/max(t_l,t_r)), 1.0_rp/min(t_l,t_r))
          frho  = (prhoE + 0.25_rp*(vel_r - vel_l)*(p_r - p_l))*inv_t
        else
          frho = rho_harm*vel_m
        endif

        alpha_m = 0.5_rp*((eint_l - t_l*cv_l) + (eint_r - t_r*cv_r))
        lambda_l = -cv_l*t_l/rho_l
        lambda_r = -cv_r*t_r/rho_r
        psi_m = vel_m*0.5_rp*(rho_r*rho_r*lambda_r + rho_l*rho_l*lambda_l)
        frhoe = alpha_m*frho - psi_m

        k_m = 0.5_rp*(u_l*u_r + v_l*v_r + w_l*w_r)
end subroutine pep_tp_pair_core

subroutine hybrid_wenoEP_PEPTP_X(weno_num,lbx,ubx,istr,iend,shock_recon, flag)

        implicit none

        integer, parameter :: vec_size  = 64
        integer, parameter :: chunkSize = vec_size - 2*GN

        integer, intent(in) :: lbx,ubx
        integer, intent(in) :: weno_num
        integer, intent(in) :: istr,iend
        integer, intent(in) :: shock_recon
        integer, intent(in) :: flag

        real(rp),   dimension(3,1-GN:chunkSize+GN,5) :: tilde_op
        real(rp),   dimension(1-GN:chunkSize+GN,13)  :: pri_1D
        real(rp),   dimension(1-GN:chunkSize+GN,5)   :: phi_arr, flx_arr, flx
        integer(1), dimension(1-GN:chunkSize+GN)     :: weno_flag_chunk

        real(rp) :: gm
        real(rp) :: dye, dxe, dye_f, dxe_f, km, ikm, nrm_x, nrm_y, Ucon_roe

        real(rp), dimension(5,-3:4) :: fp, fm
        real(rp), dimension(5)      :: fhatp, fhatm
        real(rp), dimension(5,5)    :: right, left
        real(rp), dimension(-3:3)   :: central_1_p

        real(rp) :: phi1, phi2, phi3, phi4, phi5, inner_sum
        real(rp) :: sqrt_rho0, sqrt_rho1, hFn
        real(rp) :: uroe, vroe, wroe, hroe, croe, den_roe, groe, troe
        real(rp) :: f1, f2, f3, f4, f5, w1, w2, w3, w4, w5, wc, gc

        real(rp), parameter :: one8 = 1.0_rp/8.0_rp
        real(rp) :: ir, u_, v_, w_, ek, p_, idx, t_, htot_
        real(rp) :: weight, lambda_max
        integer  :: i,j,k,l,m,il, n, s
        integer  :: str_x, end_x, igbl, istr_x, iend_x
        integer(1), parameter :: one = 1
        real(rp) :: frho, frhoe, prhoE, k_m, u_m, v_m, w_m
        real(rp) :: tl2, tr2, tl3, tr3, tl4, tr4

        call StartProfRange("hybrid_wenoEP_PEPTP_X")

        !$acc parallel vector_length(vec_size) default(present) &
        !$acc private(tilde_op,pri_1D,phi_arr) &
        !$acc private(flx_arr,flx,weno_flag_chunk) &
        !$acc create(central_1_p)
        !$acc cache(central_1_p,tilde_op,pri_1D,phi_arr,flx_arr,flx,weno_flag_chunk)
        !$acc loop gang collapse(2)
        do k   = sz,ez
          do j = sy,ey
             central_1_p = 2*central_1

             do str_x = lbx,ubx,chunkSize
                end_x = min(str_x+ChunkSize, ubx-2*GN+1) - str_x

                istr_x = max(str_x+GN       , istr     ) - str_x - GN + 1
                iend_x = min(str_x+ChunkSize, iend-GN+1) - str_x
                if(iend_x - istr_x > -1) then
                !$acc loop vector
                do i = 1-GN,end_x+GN
                   igbl = i + str_x-1+GN

                   phi1 = phi(igbl,j,k,1)
                   phi2 = phi(igbl,j,k,2)
                   phi3 = phi(igbl,j,k,3)
                   phi4 = phi(igbl,j,k,4)
                   phi5 = phi(igbl,j,k,5)

                   phi_arr(i,1) = phi1
                   phi_arr(i,2) = phi2
                   phi_arr(i,3) = phi3
                   phi_arr(i,4) = phi4
                   phi_arr(i,5) = phi5

                   gm  = GAMMA_T(igbl,j,k)

                   ir = 1.0_rp/phi1
                   u_ = phi2*ir
                   v_ = phi3*ir
                   w_ = phi4*ir
                   ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)

                   p_    = P(igbl,j,k)
                   T_    = T(igbl,j,k)
                   htot_ = ENTHALPY(igbl,j,k) + ek

                   dye = dy_deta(igbl,j)
                   dxe = dx_deta(igbl,j)

                   flx_arr(i,1) = dye* phi2           - dxe* phi3
                   flx_arr(i,2) = dye*(phi2*u_ + p_)  - dxe*(phi3*u_)
                   flx_arr(i,3) = dye*(phi2*v_)       - dxe*(phi3*v_ + p_)
                   flx_arr(i,4) = dye*(phi2*w_)       - dxe*(phi3*w_)
                   flx_arr(i,5) = dye*((phi5+p_)*u_)  - dxe*((phi5+p_)*v_)

                   pri_1D(i,1)  = p_
                   pri_1D(i,2)  = u_
                   pri_1D(i,3)  = v_
                   pri_1D(i,4)  = w_
                   pri_1D(i,5)  = htot_
                   pri_1D(i,6)  = phi1
                   pri_1D(i,7)  = gm
                   pri_1D(i,8)  = T_
                   pri_1D(i,9)  = dye*u_ - dxe*v_
                   pri_1D(i,10) = p_*dye
                   pri_1D(i,11) = p_*dxe
                   pri_1D(i,12) = p_/phi1
                   pri_1D(i,13) = phi5*ir - ek

                   weno_flag_chunk(i) = iand(weno_flag_xyz(igbl,j,k),one)

                enddo

             !$acc loop vector collapse(2)
             do i = 1-GN,end_x+GN-3
                do l = 1,3

                   il = i+l

                   tl2 = pri_1D(i,12)*pri_1D(i,12)
                   tr2 = pri_1D(il,12)*pri_1D(il,12)
                   tl3 = tl2*pri_1D(i,12)
                   tr3 = tr2*pri_1D(il,12)
                   tl4 = tl3*pri_1D(i,12)
                   tr4 = tr3*pri_1D(il,12)
                   call pep_tp_pair_core(pri_1D(i,6),pri_1D(il,6), &
                                         pri_1D(i,2),pri_1D(il,2), &
                                         pri_1D(i,3),pri_1D(il,3), &
                                         pri_1D(i,4),pri_1D(il,4), &
                                         pri_1D(i,1),pri_1D(il,1), &
                                         pri_1D(i,12),pri_1D(il,12), &
                                         pri_1D(i,13),pri_1D(il,13), &
                                         pri_1D(i,9),pri_1D(il,9), &
                                         tl2,tr2,tl3,tr3,tl4,tr4, &
                                         frho,frhoe,prhoE,k_m)

                   u_m = 0.5_rp*(pri_1D(i,2) + pri_1D(il,2))
                   v_m = 0.5_rp*(pri_1D(i,3) + pri_1D(il,3))
                   w_m = 0.5_rp*(pri_1D(i,4) + pri_1D(il,4))

                   tilde_op(l,i,1) = frho
                   tilde_op(l,i,2) = frho*u_m + 0.5_rp*(pri_1D(i,10) + pri_1D(il,10))
                   tilde_op(l,i,3) = frho*v_m - 0.5_rp*(pri_1D(i,11) + pri_1D(il,11))
                   tilde_op(l,i,4) = frho*w_m
                   tilde_op(l,i,5) = frhoe + frho*k_m + prhoE

                enddo
             enddo

             !$acc loop vector
             do i = 0,end_x
             if (weno_flag_chunk(i) == weno_smooth) then

                do n = 1,5
                   hFn = 0.0_rp
                   do l = 1, 3
                      inner_sum = 0.0_rp
                      do m = 0, l-1
                         inner_sum = inner_sum + tilde_op(l,i-m,n)
                      enddo

                      hFn = hFn + central_1_p(l)*inner_sum
                   enddo
                   flx(i,n) = hFn
                enddo
             endif
             enddo

             !$acc loop vector private(right,left,fp,fm,fhatp,fhatm)
             do i = 0,end_x
             if (weno_flag_chunk(i) == weno_shock) then

                igbl = i + str_x-1+GN

                sqrt_rho0 = sqrt(pri_1D(i,6))
                sqrt_rho1 = sqrt(pri_1D(i+1,6))
                den_roe   = 1.0_rp/(sqrt_rho0 + sqrt_rho1)

                uroe = den_roe*(pri_1D(i,2)*sqrt_rho0+pri_1D(i+1,2)*sqrt_rho1)
                vroe = den_roe*(pri_1D(i,3)*sqrt_rho0+pri_1D(i+1,3)*sqrt_rho1)
                wroe = den_roe*(pri_1D(i,4)*sqrt_rho0+pri_1D(i+1,4)*sqrt_rho1)
                hroe = den_roe*(pri_1D(i,5)*sqrt_rho0+pri_1D(i+1,5)*sqrt_rho1)
                groe = den_roe*(pri_1D(i,7)*sqrt_rho0+pri_1D(i+1,7)*sqrt_rho1)
                troe = den_roe*(pri_1D(i,8)*sqrt_rho0+pri_1D(i+1,8)*sqrt_rho1)

                croe = sqrt(groe * troe)

                dye_f = 0.5_rp*(dy_deta(igbl,j) + dy_deta(igbl+1,j))
                dxe_f = 0.5_rp*(dx_deta(igbl,j) + dx_deta(igbl+1,j))
                km    = sqrt(dye_f*dye_f + dxe_f*dxe_f)
                ikm   = 1.0_rp/km
                nrm_x =  dye_f*ikm
                nrm_y = -dxe_f*ikm
                Ucon_roe = dye_f*uroe - dxe_f*vroe

                call euler_system_n(nrm_x,nrm_y,uroe,vroe,wroe,troe,groe,hroe,right,left)

               lambda_max = max(abs(Ucon_roe-croe*km),abs(Ucon_roe+croe*km))

               do s = -weno_num+1,weno_num
                f1 = flx_arr(i+s,1)
                f2 = flx_arr(i+s,2)
                f3 = flx_arr(i+s,3)
                f4 = flx_arr(i+s,4)
                f5 = flx_arr(i+s,5)
                w1 = phi_arr(i+s,1)
                w2 = phi_arr(i+s,2)
                w3 = phi_arr(i+s,3)
                w4 = phi_arr(i+s,4)
                w5 = phi_arr(i+s,5)
                do m = 1,5
                 wc = 0._rp
                 gc = 0._rp

                 wc = wc + left(m,1) * w1
                 wc = wc + left(m,2) * w2
                 wc = wc + left(m,3) * w3
                 wc = wc + left(m,4) * w4
                 wc = wc + left(m,5) * w5
                 gc = gc + left(m,1) * f1
                 gc = gc + left(m,2) * f2
                 gc = gc + left(m,3) * f3
                 gc = gc + left(m,4) * f4
                 gc = gc + left(m,5) * f5

                 fp(m,s) = 0.5_rp * (gc+wc*lambda_max)
                 fm(m,s) = gc - 0.5_rp * (gc+wc*lambda_max)
                 enddo
                enddo

                call weno_reconstruction(weno_num,shock_recon,fp,fm,fhatp,fhatm)

                do l = 1,5
                   flx(i,l) = 0.0_rp
                   do m = 1,5
                      flx(i,l) = flx(i,l) + right(l,m) * (fhatp(m) + fhatm(m))
                   enddo
                enddo
             endif
             enddo

             !$acc loop vector
             do i = istr_x,iend_x

                igbl = i + str_x-1+GN

                idx = Jac2D(igbl,j)
                RHS(igbl,j,k,1) = - idx * (flx(i,1) - flx(i-1,1))
                RHS(igbl,j,k,2) = - idx * (flx(i,2) - flx(i-1,2))
                RHS(igbl,j,k,3) = - idx * (flx(i,3) - flx(i-1,3))
                RHS(igbl,j,k,4) = - idx * (flx(i,4) - flx(i-1,4))
                RHS(igbl,j,k,5) = - idx * (flx(i,5) - flx(i-1,5))

             enddo
             endif
             enddo
          enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine hybrid_wenoEP_PEPTP_X

subroutine hybrid_wenoEP_PEPTP_Y(weno_num,lby,uby,jstr,jend,shock_recon)

        implicit none

        integer, parameter :: vec_size  = 64
        integer, parameter :: chunkSize = vec_size - 2*GN

        integer, intent(in) :: lby,uby
        integer, intent(in) :: weno_num
        integer, intent(in) :: jstr,jend
        integer, intent(in) :: shock_recon

        real(rp),   dimension(3,1-GN:chunkSize+GN,5) :: tilde_op
        real(rp),   dimension(1-GN:chunkSize+GN,13)  :: pri_1D
        real(rp),   dimension(1-GN:chunkSize+GN,5)   :: flx_arr, flx
        real(rp),   dimension(0:3,1-GN:chunkSize+GN,5)   :: phi_arr
        integer(1), dimension(1-GN:chunkSize+GN)     :: weno_flag_chunk

        real(rp) :: gm
        real(rp) :: dyx, dxx, dyx_f, dxx_f, km, ikm, nrm_x, nrm_y, Vcon_roe

        real(rp), dimension(5,-3:4) :: fp, fm
        real(rp), dimension(5)      :: fhatp, fhatm
        real(rp), dimension(5,5)    :: right, left
        real(rp), dimension(-3:3)   :: central_1_p

        real(rp) :: phi1, phi2, phi3, phi4, phi5, inner_sum
        real(rp) :: sqrt_rho0, sqrt_rho1, hGn
        real(rp) :: uroe, vroe, wroe, hroe, croe, den_roe, groe, troe
        real(rp) :: f1, f2, f3, f4, f5, w1, w2, w3, w4, w5, wc, gc

        real(rp), parameter :: one8 = 1.0_rp/8.0_rp
        real(rp) :: ir, u_, v_, w_, ek, p_, idy, t_, htot_
        real(rp) :: weight, lambda_max
        integer  :: i,ii,j,k,l,m,jl, n, s,iimax
        integer  :: str_y, end_y, jgbl, jstr_y, jend_y
        integer(1), parameter :: two = 2
        real(rp) :: frho, frhoe, prhoE, k_m, u_m, v_m, w_m
        real(rp) :: tl2, tr2, tl3, tr3, tl4, tr4
        integer(1), parameter :: two_weno_smooth = 2*weno_smooth
        integer(1), parameter :: two_weno_shock  = 2*weno_shock

        call StartProfRange("hybrid_wenoEP_PEPTP_Y")

        !$acc parallel vector_length(vec_size) default(present) &
        !$acc private(tilde_op,pri_1D,phi_arr) &
        !$acc private(flx_arr,flx,weno_flag_chunk) &
        !$acc create(central_1_p)
        !$acc cache(central_1_p,tilde_op,pri_1D,phi_arr,flx_arr,flx,weno_flag_chunk)
        !$acc loop gang collapse(2)
        do k   = sz,ez
          do i = sx,ex,4
             central_1_p = 2*central_1

             iimax = min(3,ex-i)

             do str_y = lby,uby,chunkSize
                end_y = min(str_y+ChunkSize, uby-2*GN+1) - str_y

                jstr_y = max(str_y+GN       , jstr     ) - str_y - GN + 1
                jend_y = min(str_y+ChunkSize, jend-GN+1) - str_y
                if(jend_y - jstr_y > -1) then
                !$acc loop vector collapse(2)
                do j = 1-GN,end_y+GN
                   do ii = 0,iimax
                      jgbl = j + str_y-1+GN
                      phi_arr(ii,j,1) = phi(i+ii,jgbl,k,1)
                      phi_arr(ii,j,2) = phi(i+ii,jgbl,k,2)
                      phi_arr(ii,j,3) = phi(i+ii,jgbl,k,3)
                      phi_arr(ii,j,4) = phi(i+ii,jgbl,k,4)
                      phi_arr(ii,j,5) = phi(i+ii,jgbl,k,5)
                   enddo
                enddo
                do ii=0,iimax
                   !$acc loop vector
                   do j = 1-GN,end_y+GN
                      jgbl = j + str_y-1+GN

                      phi1 = phi_arr(ii,j,1)
                      phi2 = phi_arr(ii,j,2)
                      phi3 = phi_arr(ii,j,3)
                      phi4 = phi_arr(ii,j,4)
                      phi5 = phi_arr(ii,j,5)

                      gm   = GAMMA_T(i+ii,jgbl,k)

                      ir = 1.0_rp/phi1
                      u_ = phi2*ir
                      v_ = phi3*ir
                      w_ = phi4*ir
                      ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)

                      p_      = P(i+ii,jgbl,k)
                      T_      = T(i+ii,jgbl,k)
                      htot_   = ENTHALPY(i+ii,jgbl,k) + ek

                      dyx = dy_dxi(i+ii,jgbl)
                      dxx = dx_dxi(i+ii,jgbl)

                      flx_arr(j,1) = -dyx* phi2           + dxx* phi3
                      flx_arr(j,2) = -dyx*(phi2*u_ + p_)  + dxx*(phi3*u_)
                      flx_arr(j,3) = -dyx*(phi2*v_)       + dxx*(phi3*v_ + p_)
                      flx_arr(j,4) = -dyx*(phi2*w_)       + dxx*(phi3*w_)
                      flx_arr(j,5) = -dyx*((phi5+p_)*u_)  + dxx*((phi5+p_)*v_)

                      pri_1D(j,1)  = p_
                      pri_1D(j,2)  = u_
                      pri_1D(j,3)  = v_
                      pri_1D(j,4)  = w_
                      pri_1D(j,5)  = htot_
                      pri_1D(j,6)  = phi1
                      pri_1D(j,7)  = gm
                      pri_1D(j,8)  = T_
                      pri_1D(j,9)  = -dyx*u_ + dxx*v_
                      pri_1D(j,10) = p_*dxx
                      pri_1D(j,11) = p_*dyx
                      pri_1D(j,12) = p_/phi1
                      pri_1D(j,13) = phi5*ir - ek

                      weno_flag_chunk(j) = iand(weno_flag_xyz(i+ii,jgbl,k),two)

                   enddo

                   !$acc loop vector collapse(2)
                   do j = 1-GN,end_y+GN-3
                      do l = 1,3

                         jl = j+l

                         tl2 = pri_1D(j,12)*pri_1D(j,12)
                         tr2 = pri_1D(jl,12)*pri_1D(jl,12)
                         tl3 = tl2*pri_1D(j,12)
                         tr3 = tr2*pri_1D(jl,12)
                         tl4 = tl3*pri_1D(j,12)
                         tr4 = tr3*pri_1D(jl,12)
                         call pep_tp_pair_core(pri_1D(j,6),pri_1D(jl,6), &
                                               pri_1D(j,2),pri_1D(jl,2), &
                                               pri_1D(j,3),pri_1D(jl,3), &
                                               pri_1D(j,4),pri_1D(jl,4), &
                                               pri_1D(j,1),pri_1D(jl,1), &
                                               pri_1D(j,12),pri_1D(jl,12), &
                                               pri_1D(j,13),pri_1D(jl,13), &
                                               pri_1D(j,9),pri_1D(jl,9), &
                                               tl2,tr2,tl3,tr3,tl4,tr4, &
                                               frho,frhoe,prhoE,k_m)

                         u_m = 0.5_rp*(pri_1D(j,2) + pri_1D(jl,2))
                         v_m = 0.5_rp*(pri_1D(j,3) + pri_1D(jl,3))
                         w_m = 0.5_rp*(pri_1D(j,4) + pri_1D(jl,4))

                         tilde_op(l,j,1) = frho
                         tilde_op(l,j,2) = frho*u_m - 0.5_rp*(pri_1D(j,11) + pri_1D(jl,11))
                         tilde_op(l,j,3) = frho*v_m + 0.5_rp*(pri_1D(j,10) + pri_1D(jl,10))
                         tilde_op(l,j,4) = frho*w_m
                         tilde_op(l,j,5) = frhoe + frho*k_m + prhoE

                      enddo
                   enddo

                   !$acc loop vector
                   do j = 0,end_y
                      if(weno_flag_chunk(j) == two_weno_smooth) then

                        do n = 1,5
                           hGn = 0.0_rp
                           do l = 1, 3
                              inner_sum = 0.0_rp
                              do m = 0, l-1
                                 inner_sum = inner_sum + tilde_op(l,j-m,n)
                              enddo

                              hGn = hGn + central_1_p(l)*inner_sum
                           enddo
                           flx(j,n) = hGn
                        enddo
                      endif
                   enddo

                   !$acc loop vector private(right,left,fp,fm,fhatp,fhatm)
                   do j = 0,end_y
                      if(weno_flag_chunk(j) == two_weno_shock) then

                      jgbl = j + str_y-1+GN

                      sqrt_rho0 = sqrt(pri_1D(j,6))
                      sqrt_rho1 = sqrt(pri_1D(j+1,6))
                      den_roe   = 1.0_rp/(sqrt_rho0 + sqrt_rho1)

                      uroe = den_roe*(pri_1D(j,2)*sqrt_rho0+pri_1D(j+1,2)*sqrt_rho1)
                      vroe = den_roe*(pri_1D(j,3)*sqrt_rho0+pri_1D(j+1,3)*sqrt_rho1)
                      wroe = den_roe*(pri_1D(j,4)*sqrt_rho0+pri_1D(j+1,4)*sqrt_rho1)
                      hroe = den_roe*(pri_1D(j,5)*sqrt_rho0+pri_1D(j+1,5)*sqrt_rho1)
                      groe = den_roe*(pri_1D(j,7)*sqrt_rho0+pri_1D(j+1,7)*sqrt_rho1)
                      troe = den_roe*(pri_1D(j,8)*sqrt_rho0+pri_1D(j+1,8)*sqrt_rho1)

                      croe = sqrt(groe * troe)

                      dyx_f = 0.5_rp*(dy_dxi(i+ii,jgbl) + dy_dxi(i+ii,jgbl+1))
                      dxx_f = 0.5_rp*(dx_dxi(i+ii,jgbl) + dx_dxi(i+ii,jgbl+1))
                      km    = sqrt(dyx_f*dyx_f + dxx_f*dxx_f)
                      ikm   = 1.0_rp/km
                      nrm_x = -dyx_f*ikm
                      nrm_y =  dxx_f*ikm
                      Vcon_roe = -dyx_f*uroe + dxx_f*vroe

                      call euler_system_n(nrm_x,nrm_y,uroe,vroe,wroe,troe,groe,hroe,right,left)

                      lambda_max = max(abs(Vcon_roe-croe*km),abs(Vcon_roe+croe*km))

                      do s = -weno_num+1,weno_num
                       f1 = flx_arr(j+s,1)
                       f2 = flx_arr(j+s,2)
                       f3 = flx_arr(j+s,3)
                       f4 = flx_arr(j+s,4)
                       f5 = flx_arr(j+s,5)
                       w1 = phi_arr(ii,j+s,1)
                       w2 = phi_arr(ii,j+s,2)
                       w3 = phi_arr(ii,j+s,3)
                       w4 = phi_arr(ii,j+s,4)
                       w5 = phi_arr(ii,j+s,5)
                       do m = 1,5
                        wc = 0._rp
                        gc = 0._rp

                        wc = wc + left(m,1) * w1
                        wc = wc + left(m,2) * w2
                        wc = wc + left(m,3) * w3
                        wc = wc + left(m,4) * w4
                        wc = wc + left(m,5) * w5
                        gc = gc + left(m,1) * f1
                        gc = gc + left(m,2) * f2
                        gc = gc + left(m,3) * f3
                        gc = gc + left(m,4) * f4
                        gc = gc + left(m,5) * f5

                        fp(m,s) = 0.5_rp * (gc+wc*lambda_max)
                        fm(m,s) = gc - 0.5_rp * (gc+wc*lambda_max)
                        enddo
                       enddo

                      call weno_reconstruction(weno_num,shock_recon,fp,fm,fhatp,fhatm)

                      do l = 1,5
                         flx(j,l) = 0.0_rp
                         do m = 1,5
                            flx(j,l) = flx(j,l) + right(l,m) * (fhatp(m) + fhatm(m))
                         enddo
                      enddo
                      endif
                   enddo

                   !$acc loop vector
                   do j = jstr_y,jend_y
                      jgbl = j + str_y-1+GN
                      idy = Jac2D(i+ii,jgbl)
                      phi_arr(ii,j,1) = -idy * (flx(j,1) - flx(j-1,1))
                      phi_arr(ii,j,2) = -idy * (flx(j,2) - flx(j-1,2))
                      phi_arr(ii,j,3) = -idy * (flx(j,3) - flx(j-1,3))
                      phi_arr(ii,j,4) = -idy * (flx(j,4) - flx(j-1,4))
                      phi_arr(ii,j,5) = -idy * (flx(j,5) - flx(j-1,5))
                   enddo
                enddo

                !$acc loop vector collapse(2)
                do j = jstr_y,jend_y
                   do ii = 0,iimax
                      jgbl = j + str_y-1+GN
                      RHS(i+ii,jgbl,k,1) = RHS(i+ii,jgbl,k,1) + phi_arr(ii,j,1)
                      RHS(i+ii,jgbl,k,2) = RHS(i+ii,jgbl,k,2) + phi_arr(ii,j,2)
                      RHS(i+ii,jgbl,k,3) = RHS(i+ii,jgbl,k,3) + phi_arr(ii,j,3)
                      RHS(i+ii,jgbl,k,4) = RHS(i+ii,jgbl,k,4) + phi_arr(ii,j,4)
                      RHS(i+ii,jgbl,k,5) = RHS(i+ii,jgbl,k,5) + phi_arr(ii,j,5)
                   enddo
                enddo
                endif
             enddo
          enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine hybrid_wenoEP_PEPTP_Y

subroutine hybrid_wenoEP_PEPTP_Z(weno_num,lbz,ubz,kstr,kend,shock_recon)

        implicit none
        integer, parameter :: vec_size = 64
        integer, parameter :: chunkSize = vec_size - 2*GN

        integer, intent(in) :: lbz,ubz
        integer, intent(in) :: weno_num
        integer, intent(in) :: kstr,kend
        integer, intent(in) :: shock_recon

        real(rp),   dimension(3,1-GN:chunkSize+GN,5) :: tilde_op
        real(rp),   dimension(1-GN:chunkSize+GN,10)  :: pri_1D
        real(rp),   dimension(1-GN:chunkSize+GN,5)   :: flx_arr, flx
        real(rp),   dimension(0:3,1-GN:chunkSize+GN,5)   :: phi_arr

        integer(1), dimension(1-GN:chunkSize+GN)     :: weno_flag_chunk

        real(rp):: gm, gm1, hgm

        real(rp), dimension(5,-3:4) :: fp, fm
        real(rp), dimension(5)      :: fhatp, fhatm
        real(rp), dimension(5,5)    :: right, left
        real(rp), dimension(-3:3)   :: central_1_p

        real(rp) :: phi1, phi2, phi3, phi4, phi5, inner_sum
        real(rp) :: sqrt_rho0, sqrt_rho1, hHn
        real(rp) :: uroe, vroe, wroe, hroe, ekroe, croe, den_roe, groe, proe, rroe, troe
        real(rp) :: f1, f2, f3, f4, f5, w1, w2, w3, w4, w5, wc, gc

        real(rp), parameter :: one8 = 1.0_rp/8.0_rp
        real(rp) :: ir, u_, v_, w_, ek, p_, idz, h_, t_, htot_, e_, a_
        real(rp) :: weight, lambda_max
        real(rp) :: vnc, hgU, hgV, hgW, hgE, cc, un_c2
        integer  :: i,ii,j,k,l,m,kl, n, s, iimax
        integer  :: str_z, end_z, kgbl, kstr_z, kend_z
        integer(1), parameter :: four = 4
        integer(1), parameter :: four_weno_smooth = four*weno_smooth
        integer(1), parameter :: four_weno_shock  = four*weno_shock
        real(rp) :: frho, frhoe, prhoE, k_m, u_m, v_m, w_m
        real(rp) :: tl2, tr2, tl3, tr3, tl4, tr4
        real(rp) :: b1, b2, b3, dp_de, ci

        call StartProfRange("hybrid_wenoEP_PEPTP_Z")

        !$acc parallel vector_length(vec_size) default(present) &
        !$acc private(tilde_op,pri_1D,phi_arr) &
        !$acc private(flx_arr,flx,weno_flag_chunk) &
        !$acc create(central_1_p)
        !$acc cache(central_1_p,tilde_op,pri_1D,phi_arr,flx_arr,flx,weno_flag_chunk)
        !$acc loop gang collapse(2)
        do j    = sy,ey
           do i = sx,ex,4
              central_1_p = 2*central_1

              iimax = min(3,ex-i)

              do str_z = lbz,ubz,chunkSize
                 end_z = min(str_z+ChunkSize, ubz-2*GN+1) - str_z

                 kstr_z = max(str_z+GN       , kstr     ) - str_z - GN + 1
                 kend_z = min(str_z+ChunkSize, kend-GN+1) - str_z
                 if(kend_z - kstr_z > -1) then
                 !$acc loop vector collapse(2)
                 do k = 1-GN,end_z+GN
                    do ii=0,iimax
                       kgbl = k + str_z-1+GN
                       phi_arr(ii,k,1) = phi(i+ii,j,kgbl,1)
                       phi_arr(ii,k,2) = phi(i+ii,j,kgbl,2)
                       phi_arr(ii,k,3) = phi(i+ii,j,kgbl,3)
                       phi_arr(ii,k,4) = phi(i+ii,j,kgbl,4)
                       phi_arr(ii,k,5) = phi(i+ii,j,kgbl,5)
                    enddo
                 enddo
                 do ii=0,iimax
                    !$acc loop vector
                    do k = 1-GN,end_z+GN
                       kgbl = k + str_z-1+GN

                       phi1 = phi_arr(ii,k,1)
                       phi2 = phi_arr(ii,k,2)
                       phi3 = phi_arr(ii,k,3)
                       phi4 = phi_arr(ii,k,4)
                       phi5 = phi_arr(ii,k,5)

                       gm   = GAMMA_T(i+ii,j,kgbl)
                       gm1  = gm - 1.0_rp
                       hgm  = gm / gm1

                       ir = 1.0_rp/phi1
                       u_ = phi2*ir
                       v_ = phi3*ir
                       w_ = phi4*ir
                       ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
                       e_ = phi5 * ir - ek

                       p_     = P(i+ii,j,kgbl)
                       T_     = T(i+ii,j,kgbl)
                       htot_  = ENTHALPY(i+ii,j,kgbl) + ek

                       flx_arr(k,1) =  phi4
                       flx_arr(k,2) =  phi4    * u_
                       flx_arr(k,3) =  phi4    * v_
                       flx_arr(k,4) =  phi4    * w_ + p_
                       flx_arr(k,5) = (phi5+p_)* w_

                       pri_1D(k,1) = p_
                       pri_1D(k,2) = u_
                       pri_1D(k,3) = v_
                       pri_1D(k,4) = w_
                       pri_1D(k,5) = htot_
                       pri_1D(k,6) = phi1
                       pri_1D(k,7) = gm
                       pri_1D(k,8) = T_
                       pri_1D(k,9)  = p_/phi1
                       pri_1D(k,10) = e_

                       weno_flag_chunk(k) = iand(weno_flag_xyz(i+ii,j,kgbl),four)

                    enddo

                    !$acc loop vector collapse(2)
                    do k = 1-GN,end_z+GN-3
                       do l = 1,3

                          kl = k+l

                          tl2 = pri_1D(k,9)*pri_1D(k,9)
                          tr2 = pri_1D(kl,9)*pri_1D(kl,9)
                          tl3 = tl2*pri_1D(k,9)
                          tr3 = tr2*pri_1D(kl,9)
                          tl4 = tl3*pri_1D(k,9)
                          tr4 = tr3*pri_1D(kl,9)
                          call pep_tp_pair_core(pri_1D(k,6),pri_1D(kl,6), &
                                                pri_1D(k,2),pri_1D(kl,2), &
                                                pri_1D(k,3),pri_1D(kl,3), &
                                                pri_1D(k,4),pri_1D(kl,4), &
                                                pri_1D(k,1),pri_1D(kl,1), &
                                                pri_1D(k,9),pri_1D(kl,9), &
                                                pri_1D(k,10),pri_1D(kl,10), &
                                                pri_1D(k,4),pri_1D(kl,4), &
                                                tl2,tr2,tl3,tr3,tl4,tr4, &
                                                frho,frhoe,prhoE,k_m)

                          u_m = 0.5_rp*(pri_1D(k,2) + pri_1D(kl,2))
                          v_m = 0.5_rp*(pri_1D(k,3) + pri_1D(kl,3))
                          w_m = 0.5_rp*(pri_1D(k,4) + pri_1D(kl,4))

                          tilde_op(l,k,1) = frho
                          tilde_op(l,k,2) = frho*u_m
                          tilde_op(l,k,3) = frho*v_m
                          tilde_op(l,k,4) = frho*w_m + 0.5_rp*(pri_1D(k,1) + pri_1D(kl,1))
                          tilde_op(l,k,5) = frhoe + frho*k_m + prhoE

                       enddo
                    enddo

                    !$acc loop vector
                    do k = 0,end_z
                       if(weno_flag_chunk(k) == four_weno_smooth) then

                          do n = 1,5
                             hHn = 0.0_rp
                             do l = 1, 3
                                inner_sum = 0.0_rp
                                do m = 0, l-1
                                   inner_sum = inner_sum + tilde_op(l,k-m,n)
                                enddo

                                hHn = hHn + central_1_p(l)*inner_sum
                             enddo
                             flx(k,n) = hHn
                          enddo
                       endif
                    enddo

                    !$acc loop vector private(right,left,fp,fm,fhatp,fhatm)
                    do k = 0,end_z
                       if (weno_flag_chunk(k) == four_weno_shock) then
                          sqrt_rho0 = sqrt(pri_1D(k,6))
                          sqrt_rho1 = sqrt(pri_1D(k+1,6))
                          den_roe   = 1.0_rp/(sqrt_rho0 + sqrt_rho1)

                          rroe = sqrt_rho0 * sqrt_rho1
                          proe = den_roe*(pri_1D(k,1)*sqrt_rho0+pri_1D(k+1,1)*sqrt_rho1)
                          uroe = den_roe*(pri_1D(k,2)*sqrt_rho0+pri_1D(k+1,2)*sqrt_rho1)
                          vroe = den_roe*(pri_1D(k,3)*sqrt_rho0+pri_1D(k+1,3)*sqrt_rho1)
                          wroe = den_roe*(pri_1D(k,4)*sqrt_rho0+pri_1D(k+1,4)*sqrt_rho1)
                          hroe = den_roe*(pri_1D(k,5)*sqrt_rho0+pri_1D(k+1,5)*sqrt_rho1)
                          groe = den_roe*(pri_1D(k,7)*sqrt_rho0+pri_1D(k+1,7)*sqrt_rho1)
                          troe = den_roe*(pri_1D(k,8)*sqrt_rho0+pri_1D(k+1,8)*sqrt_rho1)

                          ekroe= 0.5_rp*(uroe*uroe + vroe*vroe + wroe*wroe)
                          croe = sqrt(groe * troe)

                          dp_de   = rroe * (groe - 1.0_rp)
                          b3      = hroe - troe - rroe * troe / dp_de
                          b2      = dp_de / (rroe * croe * croe)
                          b1      = troe / (croe *  croe) - b2 * (hroe - troe - 2.0_rp * ekroe)

                          right(1,1) = 1.0_rp
                          right(2,1) = uroe
                          right(3,1) = vroe
                          right(4,1) = wroe - croe
                          right(5,1) = hroe - wroe*croe

                          right(1,2) = 1.0_rp
                          right(2,2) = uroe
                          right(3,2) = vroe
                          right(4,2) = wroe
                          right(5,2) = b3

                          right(1,3) = 0.0_rp
                          right(2,3) = 1.0_rp
                          right(3,3) = 0.0_rp
                          right(4,3) = 0.0_rp
                          right(5,3) = uroe

                          right(1,4) = 0.0_rp
                          right(2,4) = 0.0_rp
                          right(3,4) = 1.0_rp
                          right(4,4) = 0.0_rp
                          right(5,4) = vroe

                          right(1,5) = 1.0_rp
                          right(2,5) = uroe
                          right(3,5) = vroe
                          right(4,5) = wroe + croe
                          right(5,5) = hroe + wroe*croe

                          ci = 1.0_rp / croe

                          left(1,1) = 0.5_rp*(b1 + wroe*ci)
                          left(2,1) = 1.0_rp - b1
                          left(3,1) = -uroe
                          left(4,1) = -vroe
                          left(5,1) = 0.5_rp*(b1 - wroe*ci)

                          left(1,2) = -0.5_rp*(b2*uroe)
                          left(2,2) =  b2*uroe
                          left(3,2) =  1.0_rp
                          left(4,2) =  0.0_rp
                          left(5,2) = -0.5_rp*(b2*uroe)

                          left(1,3) = -0.5_rp*(b2*vroe)
                          left(2,3) =  b2*vroe
                          left(3,3) =  0.0_rp
                          left(4,3) =  1.0_rp
                          left(5,3) = -0.5_rp*(b2*vroe)

                          left(1,4) = -0.5_rp*(b2*wroe + ci)
                          left(2,4) =  b2*wroe
                          left(3,4) =  0.0_rp
                          left(4,4) =  0.0_rp
                          left(5,4) = -0.5_rp*(b2*wroe - ci)

                          left(1,5) =  0.5_rp*b2
                          left(2,5) = -b2
                          left(3,5) =  0.0_rp
                          left(4,5) =  0.0_rp
                          left(5,5) =  0.5_rp*b2

                          lambda_max = max(abs(wroe-croe),abs(wroe+croe))

                          do s = -weno_num+1,weno_num
                           f1 = flx_arr(k+s,1)
                           f2 = flx_arr(k+s,2)
                           f3 = flx_arr(k+s,3)
                           f4 = flx_arr(k+s,4)
                           f5 = flx_arr(k+s,5)
                           w1 = phi_arr(ii,k+s,1)
                           w2 = phi_arr(ii,k+s,2)
                           w3 = phi_arr(ii,k+s,3)
                           w4 = phi_arr(ii,k+s,4)
                           w5 = phi_arr(ii,k+s,5)
                           do m = 1,5
                            wc = 0._rp
                            gc = 0._rp

                            wc = wc + left(m,1) * w1
                            wc = wc + left(m,2) * w2
                            wc = wc + left(m,3) * w3
                            wc = wc + left(m,4) * w4
                            wc = wc + left(m,5) * w5
                            gc = gc + left(m,1) * f1
                            gc = gc + left(m,2) * f2
                            gc = gc + left(m,3) * f3
                            gc = gc + left(m,4) * f4
                            gc = gc + left(m,5) * f5

                            fp(m,s) = 0.5_rp * (gc+wc*lambda_max)
                            fm(m,s) = gc - 0.5_rp * (gc+wc*lambda_max)
                            enddo
                           enddo

                          call weno_reconstruction(weno_num,shock_recon,fp,fm,fhatp,fhatm)

                          do l = 1,5
                             flx(k,l) = 0.0_rp
                             do m = 1,5
                                flx(k,l) = flx(k,l) + right(l,m) * (fhatp(m) + fhatm(m))
                             enddo
                          enddo
                       endif
                    enddo

                    !$acc loop vector
                    do k = kstr_z,kend_z
                       kgbl = k + str_z-1+GN
                       idz = zstep_i(kgbl)
                       phi_arr(ii,k,1) =  - idz * (flx(k,1) - flx(k-1,1))
                       phi_arr(ii,k,2) =  - idz * (flx(k,2) - flx(k-1,2))
                       phi_arr(ii,k,3) =  - idz * (flx(k,3) - flx(k-1,3))
                       phi_arr(ii,k,4) =  - idz * (flx(k,4) - flx(k-1,4))
                       phi_arr(ii,k,5) =  - idz * (flx(k,5) - flx(k-1,5))
                    enddo
                 enddo

                 !$acc loop vector collapse(2)
                 do k = kstr_z,kend_z
                    do ii=0,iimax
                       kgbl = k + str_z-1+GN
                       RHS(i+ii,j,kgbl,1) = RHS(i+ii,j,kgbl,1) + phi_arr(ii,k,1)
                       RHS(i+ii,j,kgbl,2) = RHS(i+ii,j,kgbl,2) + phi_arr(ii,k,2)
                       RHS(i+ii,j,kgbl,3) = RHS(i+ii,j,kgbl,3) + phi_arr(ii,k,3)
                       RHS(i+ii,j,kgbl,4) = RHS(i+ii,j,kgbl,4) + phi_arr(ii,k,4)
                       RHS(i+ii,j,kgbl,5) = RHS(i+ii,j,kgbl,5) + phi_arr(ii,k,5)
                    enddo
                 enddo
              endif
              enddo
           enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine hybrid_wenoEP_PEPTP_Z

end module hybrid_wenoEP_PEPTP_module
