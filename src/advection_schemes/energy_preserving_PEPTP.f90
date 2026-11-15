module energy_preserving_PEPTP_module
use parameters_module
use storage_module
use mesh_module, only: dx_dxi, dx_deta, dy_dxi, dy_deta, Jac2D
use thermodynamics_module, only: cv_Tref, T_switch, T_cold, cp_coef, aL_adim__0, aH_adim__0
use profiling_module
implicit none
private
public :: energy_preserving_PEPTP_X, energy_preserving_PEPTP_Y, energy_preserving_PEPTP_Z

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

subroutine energy_preserving_PEPTP_X(lbx,ubx,istr,iend)

        implicit none
        integer, intent(in) :: lbx,ubx
        integer, intent(in) :: istr,iend

        integer, parameter :: vec_size  = 32
        integer, parameter :: chunkSize = 2*vec_size - 2*GN

        real(rp), dimension(3,1-GN:chunkSize+GN,5) :: tilde_op
        real(rp), dimension(1-GN:chunkSize+GN,10) :: pri_1D
        real(rp), dimension(1-GN:chunkSize+GN,5)  :: flx
        real(rp), dimension(-3:3)                 :: central_1_p

        integer :: i, j, k, l, il, m, n
        integer :: istr_x, iend_x, igbl, str_x, end_x

        real(rp) :: ir, u_, v_, w_, ek, p_, inner_sum, hFn, idx
        real(rp) :: tl2, tr2, tl3, tr3, tl4, tr4
        real(rp) :: phi1, phi2, phi3, phi4, phi5
        real(rp) :: frho, frhoe, prhoE, k_m, u_m, v_m, w_m

        call StartProfRange("energy_preserving_PEPTP_X")

        !$acc parallel vector_length(vec_size) default(present) &
        !$acc private(tilde_op,pri_1D,flx) &
        !$acc create(central_1_p)
        !$acc cache(central_1_p,tilde_op,pri_1D,flx)
        !$acc loop gang collapse(2)
        do k = sz,ez
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

                     ir = 1.0_rp/phi1
                     u_ = phi2*ir
                     v_ = phi3*ir
                     w_ = phi4*ir
                     ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)

                     p_ = P(igbl,j,k)
                     pri_1D(i,1)  = p_
                     pri_1D(i,2)  = u_
                     pri_1D(i,3)  = v_
                     pri_1D(i,4)  = w_
                     pri_1D(i,5)  = phi5*ir - ek
                     pri_1D(i,6)  = phi1
                     pri_1D(i,7)  = p_/phi1
                     pri_1D(i,8)  = dy_deta(igbl,j)*u_ - dx_deta(igbl,j)*v_
                     pri_1D(i,9)  = p_*dy_deta(igbl,j)
                     pri_1D(i,10) = p_*dx_deta(igbl,j)
                  enddo

                  !$acc loop vector collapse(2)
                  do i = 1-GN,end_x+GN-3
                     do l = 1,3

                        il = i+l
                        tl2 = pri_1D(i,7)*pri_1D(i,7)
                        tr2 = pri_1D(il,7)*pri_1D(il,7)
                        tl3 = tl2*pri_1D(i,7)
                        tr3 = tr2*pri_1D(il,7)
                        tl4 = tl3*pri_1D(i,7)
                        tr4 = tr3*pri_1D(il,7)
                        call pep_tp_pair_core(pri_1D(i,6),pri_1D(il,6), &
                                              pri_1D(i,2),pri_1D(il,2), &
                                              pri_1D(i,3),pri_1D(il,3), &
                                              pri_1D(i,4),pri_1D(il,4), &
                                              pri_1D(i,1),pri_1D(il,1), &
                                              pri_1D(i,7),pri_1D(il,7), &
                                              pri_1D(i,5),pri_1D(il,5), &
                                              pri_1D(i,8),pri_1D(il,8), &
                                              tl2,tr2,tl3,tr3,tl4,tr4, &
                                              frho,frhoe,prhoE,k_m)

                        u_m = 0.5_rp*(pri_1D(i,2) + pri_1D(il,2))
                        v_m = 0.5_rp*(pri_1D(i,3) + pri_1D(il,3))
                        w_m = 0.5_rp*(pri_1D(i,4) + pri_1D(il,4))

                        tilde_op(l,i,1) = frho
                        tilde_op(l,i,2) = frho*u_m + 0.5_rp*(pri_1D(i,9) + pri_1D(il,9))
                        tilde_op(l,i,3) = frho*v_m - 0.5_rp*(pri_1D(i,10) + pri_1D(il,10))
                        tilde_op(l,i,4) = frho*w_m
                        tilde_op(l,i,5) = frhoe + frho*k_m + prhoE

                     enddo
                  enddo

                  !$acc loop vector
                  do i = 0,end_x
                     do n = 1,5
                        hFn = 0.0_rp
                        do l = 1,3
                           inner_sum = 0.0_rp
                           do m = 0,l-1
                              inner_sum = inner_sum + tilde_op(l,i-m,n)
                           enddo
                           hFn = hFn + central_1_p(l)*inner_sum
                        enddo
                        flx(i,n) = hFn
                     enddo
                  enddo

                  !$acc loop vector
                  do i = istr_x,iend_x
                     igbl = i + str_x-1+GN
                     idx = Jac2D(igbl,j)
                     RHS(igbl,j,k,1) = -idx*(flx(i,1) - flx(i-1,1))
                     RHS(igbl,j,k,2) = -idx*(flx(i,2) - flx(i-1,2))
                     RHS(igbl,j,k,3) = -idx*(flx(i,3) - flx(i-1,3))
                     RHS(igbl,j,k,4) = -idx*(flx(i,4) - flx(i-1,4))
                     RHS(igbl,j,k,5) = -idx*(flx(i,5) - flx(i-1,5))
                  enddo
                endif
             enddo
           enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine energy_preserving_PEPTP_X

subroutine energy_preserving_PEPTP_Y(lby,uby,jstr,jend)

        implicit none
        integer, parameter :: vec_size  = 64
        integer, parameter :: chunkSize = vec_size - 2*GN

        real(rp), dimension(3,1-GN:chunkSize+GN,5)   :: tilde_op
        real(rp), dimension(1-GN:chunkSize+GN,13)    :: pri_1D
        real(rp), dimension(1-GN:chunkSize+GN,5)     :: flx
        real(rp), dimension(0:3,1-GN:chunkSize+GN,5) :: phi_arr
        real(rp), dimension(-3:3)                    :: central_1_p

        integer, intent(in) :: lby, uby
        integer, intent(in) :: jstr,jend

        integer :: i, j, k, l, m, jl, n, str_y, end_y, ii, iimax, jstr_y, jend_y, jgbl

        real(rp) :: ir, u_, v_, w_, ek, p_, inner_sum, hGn, idy
        real(rp) :: phi1, phi2, phi3, phi4, phi5
        real(rp) :: frho, frhoe, prhoE, k_m, u_m, v_m, w_m

        call StartProfRange("energy_preserving_PEPTP_Y")

        !$acc parallel vector_length(vec_size) default(present) &
        !$acc private(tilde_op,pri_1D,phi_arr,flx) &
        !$acc create(central_1_p)
        !$acc cache(central_1_p,tilde_op,pri_1D,phi_arr,flx)
        !$acc loop gang collapse(2)
        do k = sz,ez
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

                  do ii = 0,iimax
                     !$acc loop vector
                     do j = 1-GN,end_y+GN
                        jgbl = j + str_y-1+GN

                        phi1 = phi_arr(ii,j,1)
                        phi2 = phi_arr(ii,j,2)
                        phi3 = phi_arr(ii,j,3)
                        phi4 = phi_arr(ii,j,4)
                        phi5 = phi_arr(ii,j,5)

                        ir = 1.0_rp/phi1
                        u_ = phi2*ir
                        v_ = phi3*ir
                        w_ = phi4*ir
                        ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)

                        p_ = P(i+ii,jgbl,k)
                        pri_1D(j,1)  = p_
                        pri_1D(j,2)  = u_
                        pri_1D(j,3)  = v_
                        pri_1D(j,4)  = w_
                        pri_1D(j,5)  = phi5*ir - ek
                        pri_1D(j,6)  = phi1
                        pri_1D(j,7)  = p_/phi1
                        pri_1D(j,8)  = -dy_dxi(i+ii,jgbl)*u_ + dx_dxi(i+ii,jgbl)*v_
                        pri_1D(j,9)  = p_*dx_dxi(i+ii,jgbl)
                        pri_1D(j,10) = p_*dy_dxi(i+ii,jgbl)
                        pri_1D(j,11) = pri_1D(j,7)*pri_1D(j,7)
                        pri_1D(j,12) = pri_1D(j,11)*pri_1D(j,7)
                        pri_1D(j,13) = pri_1D(j,12)*pri_1D(j,7)
                     enddo

                     !$acc loop vector collapse(2)
                     do j = 1-GN,end_y+GN-3
                        do l = 1,3

                           jl = j+l
                           call pep_tp_pair_core(pri_1D(j,6),pri_1D(jl,6), &
                                                 pri_1D(j,2),pri_1D(jl,2), &
                                                 pri_1D(j,3),pri_1D(jl,3), &
                                                 pri_1D(j,4),pri_1D(jl,4), &
                                                 pri_1D(j,1),pri_1D(jl,1), &
                                                 pri_1D(j,7),pri_1D(jl,7), &
                                                 pri_1D(j,5),pri_1D(jl,5), &
                                                 pri_1D(j,8),pri_1D(jl,8), &
                                                 pri_1D(j,11),pri_1D(jl,11), &
                                                 pri_1D(j,12),pri_1D(jl,12), &
                                                 pri_1D(j,13),pri_1D(jl,13), &
                                                 frho,frhoe,prhoE,k_m)

                           u_m = 0.5_rp*(pri_1D(j,2) + pri_1D(jl,2))
                           v_m = 0.5_rp*(pri_1D(j,3) + pri_1D(jl,3))
                           w_m = 0.5_rp*(pri_1D(j,4) + pri_1D(jl,4))

                           tilde_op(l,j,1) = frho
                           tilde_op(l,j,2) = frho*u_m - 0.5_rp*(pri_1D(j,10) + pri_1D(jl,10))
                           tilde_op(l,j,3) = frho*v_m + 0.5_rp*(pri_1D(j,9) + pri_1D(jl,9))
                           tilde_op(l,j,4) = frho*w_m
                           tilde_op(l,j,5) = frhoe + frho*k_m + prhoE

                        enddo
                     enddo

                     !$acc loop vector
                     do j = 0,end_y
                        do n = 1,5
                           hGn = 0.0_rp
                           do l = 1,3
                              inner_sum = 0.0_rp
                              do m = 0,l-1
                                 inner_sum = inner_sum + tilde_op(l,j-m,n)
                              enddo
                              hGn = hGn + central_1_p(l)*inner_sum
                           enddo
                           flx(j,n) = hGn
                        enddo
                     enddo

                     !$acc loop vector
                     do j = jstr_y,jend_y
                        jgbl = j + str_y-1+GN
                        idy = Jac2D(i+ii,jgbl)
                        phi_arr(ii,j,1) = -idy*(flx(j,1) - flx(j-1,1))
                        phi_arr(ii,j,2) = -idy*(flx(j,2) - flx(j-1,2))
                        phi_arr(ii,j,3) = -idy*(flx(j,3) - flx(j-1,3))
                        phi_arr(ii,j,4) = -idy*(flx(j,4) - flx(j-1,4))
                        phi_arr(ii,j,5) = -idy*(flx(j,5) - flx(j-1,5))
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
end subroutine energy_preserving_PEPTP_Y

subroutine energy_preserving_PEPTP_Z(lbz,ubz,kstr,kend)

        implicit none
        integer, parameter :: vec_size = 64
        integer, parameter :: chunkSize = vec_size - 2*GN

        integer, intent(in) :: lbz,ubz
        integer, intent(in) :: kstr,kend

        real(rp), dimension(3,1-GN:chunkSize+GN,5)   :: tilde_op
        real(rp), dimension(1-GN:chunkSize+GN,7)     :: pri_1D
        real(rp), dimension(1-GN:chunkSize+GN,5)     :: flx
        real(rp), dimension(0:3,1-GN:chunkSize+GN,5) :: phi_arr
        real(rp), dimension(-3:3)                    :: central_1_p

        integer :: i, j, k, l, m, kl, n, ii, iimax, str_z, end_z, kstr_z, kend_z, kgbl

        real(rp) :: ir, u_, v_, w_, ek, p_, inner_sum, hHn, idz
        real(rp) :: tl2, tr2, tl3, tr3, tl4, tr4
        real(rp) :: phi1, phi2, phi3, phi4, phi5
        real(rp) :: frho, frhoe, prhoE, k_m, u_m, v_m, w_m

        call StartProfRange("energy_preserving_PEPTP_Z")

        !$acc parallel vector_length(vec_size) default(present) &
        !$acc private(tilde_op,pri_1D,phi_arr,flx) &
        !$acc create(central_1_p)
        !$acc cache(central_1_p,tilde_op,pri_1D,phi_arr,flx)
        !$acc loop gang collapse(2)
        do j = sy,ey
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
                    do ii = 0,iimax
                       kgbl = k + str_z-1+GN
                       phi_arr(ii,k,1) = phi(i+ii,j,kgbl,1)
                       phi_arr(ii,k,2) = phi(i+ii,j,kgbl,2)
                       phi_arr(ii,k,3) = phi(i+ii,j,kgbl,3)
                       phi_arr(ii,k,4) = phi(i+ii,j,kgbl,4)
                       phi_arr(ii,k,5) = phi(i+ii,j,kgbl,5)
                    enddo
                 enddo

                 do ii = 0,iimax
                    !$acc loop vector
                    do k = 1-GN,end_z+GN
                       kgbl = k + str_z-1+GN

                       phi1 = phi_arr(ii,k,1)
                       phi2 = phi_arr(ii,k,2)
                       phi3 = phi_arr(ii,k,3)
                       phi4 = phi_arr(ii,k,4)
                       phi5 = phi_arr(ii,k,5)

                       ir = 1.0_rp/phi1
                       u_ = phi2*ir
                       v_ = phi3*ir
                       w_ = phi4*ir
                       ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)

                       p_ = P(i+ii,j,kgbl)
                       pri_1D(k,1) = p_
                       pri_1D(k,2) = u_
                       pri_1D(k,3) = v_
                       pri_1D(k,4) = w_
                       pri_1D(k,5) = phi5*ir - ek
                       pri_1D(k,6) = phi1
                       pri_1D(k,7) = p_/phi1
                    enddo

                    !$acc loop vector collapse(2)
                    do k = 1-GN,end_z+GN-3
                       do l = 1,3

                          kl = k+l
                          tl2 = pri_1D(k,7)*pri_1D(k,7)
                          tr2 = pri_1D(kl,7)*pri_1D(kl,7)
                          tl3 = tl2*pri_1D(k,7)
                          tr3 = tr2*pri_1D(kl,7)
                          tl4 = tl3*pri_1D(k,7)
                          tr4 = tr3*pri_1D(kl,7)
                          call pep_tp_pair_core(pri_1D(k,6),pri_1D(kl,6), &
                                                pri_1D(k,2),pri_1D(kl,2), &
                                                pri_1D(k,3),pri_1D(kl,3), &
                                                pri_1D(k,4),pri_1D(kl,4), &
                                                pri_1D(k,1),pri_1D(kl,1), &
                                                pri_1D(k,7),pri_1D(kl,7), &
                                                pri_1D(k,5),pri_1D(kl,5), &
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
                       do n = 1,5
                          hHn = 0.0_rp
                          do l = 1,3
                             inner_sum = 0.0_rp
                             do m = 0,l-1
                                inner_sum = inner_sum + tilde_op(l,k-m,n)
                             enddo
                             hHn = hHn + central_1_p(l)*inner_sum
                          enddo
                          flx(k,n) = hHn
                       enddo
                    enddo

                    !$acc loop vector
                    do k = kstr_z,kend_z
                       kgbl = k + str_z-1+GN
                       idz = zstep_i(kgbl)
                       phi_arr(ii,k,1) = -idz*(flx(k,1) - flx(k-1,1))
                       phi_arr(ii,k,2) = -idz*(flx(k,2) - flx(k-1,2))
                       phi_arr(ii,k,3) = -idz*(flx(k,3) - flx(k-1,3))
                       phi_arr(ii,k,4) = -idz*(flx(k,4) - flx(k-1,4))
                       phi_arr(ii,k,5) = -idz*(flx(k,5) - flx(k-1,5))
                    enddo
                 enddo

                 !$acc loop vector collapse(2)
                 do k = kstr_z,kend_z
                    do ii = 0,iimax
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
end subroutine energy_preserving_PEPTP_Z

end module energy_preserving_PEPTP_module
