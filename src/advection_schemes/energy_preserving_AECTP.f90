module energy_preserving_AECTP_module
use parameters_module
use storage_module
use mesh_module, only: dx_dxi, dx_deta, dy_dxi, dy_deta, Jac2D
use thermodynamics_module, only: cv_Tref, T_switch, T_cold, cp_coef, h_coef, aL_adim__0, aH_adim__0
use profiling_module
implicit none
private
public :: energy_preserving_AECTP_X, energy_preserving_AECTP_Y, energy_preserving_AECTP_Z

contains

pure function logmean_tp(a,b) result(lm)
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
end function logmean_tp

pure function aec_tp_internal_energy_mean(tl,tr,inv_tl,inv_tr,logtl,logtr) result(e_mean)
!$acc routine seq
        implicit none
        real(rp), intent(in) :: tl,tr
        real(rp), intent(in) :: inv_tl, inv_tr
        real(rp), intent(in) :: logtl, logtr
        real(rp) :: e_mean
        real(rp) :: c_m2, c_m1, c_0, c_p1, c_p2, c_p3, c_p4, h_ref
        real(rp) :: invt_log, meanlogt
        real(rp) :: tbar, sumlow, sumup, sumref, ee_star, h_ref_value
        integer  :: kb

        invt_log = logmean_tp(inv_tl,inv_tr)

        if(aL_adim__0 == 0.0_rp .and. aH_adim__0 == 0.0_rp) then
          e_mean = cv_Tref/invt_log
          return
        endif

        if(sqrt(tl*tr) <= T_cold) then
          kb = 0
        else if(sqrt(tl*tr) <= T_switch) then
          kb = 1
        else
          kb = 2
        endif
        c_m2 = cp_coef(-2,kb); c_m1 = cp_coef(-1,kb); c_0 = cp_coef(0,kb)
        c_p1 = cp_coef( 1,kb); c_p2 = cp_coef( 2,kb)
        c_p3 = cp_coef( 3,kb); c_p4 = cp_coef( 4,kb)
        h_ref = h_coef(kb)

        meanlogt = 0.5_rp*(logtl + logtr)

        tbar   = 0.5_rp*(inv_tl + inv_tr)
        sumlow = c_m2/(-1.0_rp)*(tbar - 1.0_rp)

        tbar  = tr*tl
        sumup = c_p1/2.0_rp*(tbar - 1.0_rp)

        tbar  = 0.5_rp*(tr*tr*tl + tr*tl*tl)
        sumup = sumup + c_p2/3.0_rp*(tbar - 1.0_rp)

        tbar  = (tr*tr*tr*tl + tr*tr*tl*tl + tr*tl*tl*tl)/3.0_rp
        sumup = sumup + c_p3/4.0_rp*(tbar - 1.0_rp)

        tbar  = (tr**4*tl + tr**3*tl*tl + tr*tr*tl**3 + tr*tl**4)/4.0_rp
        sumup = sumup + c_p4/5.0_rp*(tbar - 1.0_rp)

        sumref = c_0
        h_ref_value = h_ref - c_m2 + c_0 + c_p1/2.0_rp + c_p2/3.0_rp + &
                      c_p3/4.0_rp + c_p4/5.0_rp

        ee_star = c_m1*(-0.5_rp*(inv_tl + inv_tr)/invt_log + 1.0_rp + meanlogt) + &
                  (c_0 - 1.0_rp)/invt_log + sumlow + sumup

        e_mean = h_ref_value - sumref + ee_star
end function aec_tp_internal_energy_mean

subroutine energy_preserving_AECTP_X(lbx,ubx,istr,iend)

        implicit none
        integer, intent(in) :: lbx,ubx
        integer, intent(in) :: istr,iend

        integer, parameter :: vec_size  = 64
        integer, parameter :: chunkSize = vec_size - 2*GN

        real(rp), dimension(3,1-GN:chunkSize+GN,5) :: tilde_op
        real(rp), dimension(1-GN:chunkSize+GN,11) :: pri_1D
        real(rp), dimension(1-GN:chunkSize+GN,5)   :: flx
        real(rp), dimension(-3:3)                  :: central_1_p

        integer :: i, j, k, l, il, m, n
        integer :: istr_x, iend_x, igbl, str_x, end_x

        real(rp) :: ir, u_, v_, w_, p_, weight, inner_sum, hFn, Ucon
        real(rp) :: phi1, phi2, phi3, phi4, idx
        real(rp) :: ui_uil, rhom_log, eintm_log
        real(rp), parameter :: one4 = 1.0_rp/4.0_rp

        call StartProfRange("energy_preserving_AECTP_X")

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

                     ir = 1.0_rp/phi1
                     u_ = phi2*ir
                     v_ = phi3*ir
                     w_ = phi4*ir

                     p_    = P(igbl,j,k)
                      pri_1D(i,1) = p_
                      pri_1D(i,2) = u_
                      pri_1D(i,3) = v_
                      pri_1D(i,4) = w_
                      pri_1D(i,5) = phi1

                      pri_1D(i,6) = dy_deta(igbl,j)*u_ - dx_deta(igbl,j)*v_
                      pri_1D(i,7) = p_*dy_deta(igbl,j)
                      pri_1D(i,8) = p_*dx_deta(igbl,j)
                      pri_1D(i,10) = pri_1D(i,1)/pri_1D(i,5)
                      pri_1D(i,11) = 1.0_rp/pri_1D(i,10)
                      pri_1D(i,9) = log(pri_1D(i,10))

                  enddo

                  !$acc loop vector collapse(2)
                  do i = 1-GN,end_x+GN-3
                     do l = 1,3

                        il = i+l
                        ui_uil   = pri_1D(i,2) + pri_1D(il,2)
                        Ucon     = pri_1D(i,6) + pri_1D(il,6)
                        rhom_log = logmean_tp(pri_1D(i,5),pri_1D(il,5))
                        eintm_log = aec_tp_internal_energy_mean(pri_1D(i,10), pri_1D(il,10), &
                                                                pri_1D(i,11), pri_1D(il,11), &
                                                                pri_1D(i,9), pri_1D(il,9))
                        weight = one4*rhom_log*Ucon

                        tilde_op(l,i,1) = 2.0_rp*weight
                        tilde_op(l,i,2) = weight*ui_uil + 0.5_rp*(pri_1D(i,7) + pri_1D(il,7))
                        tilde_op(l,i,3) = weight*(pri_1D(i,3) + pri_1D(il,3)) &
                                        - 0.5_rp*(pri_1D(i,8) + pri_1D(il,8))
                        tilde_op(l,i,4) = weight*(pri_1D(i,4) + pri_1D(il,4))
                        tilde_op(l,i,5) = weight*(pri_1D(i,2)*pri_1D(il,2) &
                                        + pri_1D(i,3)*pri_1D(il,3) &
                                        + pri_1D(i,4)*pri_1D(il,4) &
                                        + 2.0_rp*eintm_log) &
                                        + 0.5_rp*(pri_1D(i,1)*pri_1D(il,6) &
                                        + pri_1D(il,1)*pri_1D(i,6))

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
end subroutine energy_preserving_AECTP_X

subroutine energy_preserving_AECTP_Y(lby,uby,jstr,jend)

        implicit none
        integer, parameter :: vec_size  = 64
        integer, parameter :: chunkSize = vec_size - 2*GN

        real(rp), dimension(3,1-GN:chunkSize+GN,5)   :: tilde_op
        real(rp), dimension(1-GN:chunkSize+GN,11)   :: pri_1D
        real(rp), dimension(1-GN:chunkSize+GN,5)     :: flx
        real(rp), dimension(0:3,1-GN:chunkSize+GN,5) :: phi_arr
        real(rp), dimension(-3:3)                    :: central_1_p

        integer, intent(in) :: lby, uby
        integer, intent(in) :: jstr,jend

        integer :: i, j, k, l, m, jl, n, str_y, end_y, ii, iimax, jstr_y, jend_y, jgbl

        real(rp) :: ir, u_, v_, w_, p_, weight, inner_sum, hGn, Vcon
        real(rp) :: phi1, phi2, phi3, phi4, idy
        real(rp) :: vj_vjl, rhom_log, eintm_log
        real(rp), parameter :: one4 = 1.0_rp/4.0_rp

        call StartProfRange("energy_preserving_AECTP_Y")

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

                  do ii=0,iimax
                     !$acc loop vector
                     do j = 1-GN,end_y+GN
                        jgbl = j + str_y-1+GN

                        phi1 = phi_arr(ii,j,1)
                        phi2 = phi_arr(ii,j,2)
                        phi3 = phi_arr(ii,j,3)
                        phi4 = phi_arr(ii,j,4)

                        ir = 1.0_rp/phi1
                        u_ = phi2*ir
                        v_ = phi3*ir
                        w_ = phi4*ir

                        p_    = P(i+ii,jgbl,k)
                        pri_1D(j,1) = p_
                        pri_1D(j,2) = u_
                        pri_1D(j,3) = v_
                        pri_1D(j,4) = w_
                        pri_1D(j,5) = phi1

                        pri_1D(j,6) = -dy_dxi(i+ii,jgbl)*u_ + dx_dxi(i+ii,jgbl)*v_
                        pri_1D(j,7) = p_*dx_dxi(i+ii,jgbl)
                        pri_1D(j,8) = p_*dy_dxi(i+ii,jgbl)
                        pri_1D(j,10) = pri_1D(j,1)/pri_1D(j,5)
                        pri_1D(j,11) = 1.0_rp/pri_1D(j,10)
                        pri_1D(j,9) = log(pri_1D(j,10))

                     enddo

                     !$acc loop vector collapse(2)
                     do j = 1-GN,end_y+GN-3
                        do l = 1,3

                           jl = j+l
                           vj_vjl   = pri_1D(j,3) + pri_1D(jl,3)
                           Vcon     = pri_1D(j,6) + pri_1D(jl,6)
                           rhom_log = logmean_tp(pri_1D(j,5),pri_1D(jl,5))
                           eintm_log = aec_tp_internal_energy_mean(pri_1D(j,10), pri_1D(jl,10), &
                                                                   pri_1D(j,11), pri_1D(jl,11), &
                                                                   pri_1D(j,9), pri_1D(jl,9))
                           weight = one4*rhom_log*Vcon

                           tilde_op(l,j,1) = 2.0_rp*weight
                           tilde_op(l,j,2) = weight*(pri_1D(j,2) + pri_1D(jl,2)) &
                                           - 0.5_rp*(pri_1D(j,8) + pri_1D(jl,8))
                           tilde_op(l,j,3) = weight*vj_vjl + 0.5_rp*(pri_1D(j,7) + pri_1D(jl,7))
                           tilde_op(l,j,4) = weight*(pri_1D(j,4) + pri_1D(jl,4))
                           tilde_op(l,j,5) = weight*(pri_1D(j,2)*pri_1D(jl,2) &
                                           + pri_1D(j,3)*pri_1D(jl,3) &
                                           + pri_1D(j,4)*pri_1D(jl,4) &
                                           + 2.0_rp*eintm_log) &
                                           + 0.5_rp*(pri_1D(j,1)*pri_1D(jl,6) &
                                           + pri_1D(jl,1)*pri_1D(j,6))

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
end subroutine energy_preserving_AECTP_Y

subroutine energy_preserving_AECTP_Z(lbz,ubz,kstr,kend)

        implicit none

        integer, parameter :: vec_size = 64
        integer, parameter :: chunkSize = vec_size - 2*GN

        integer, intent(in) :: lbz,ubz
        integer, intent(in) :: kstr,kend

        real(rp), dimension(3,1-GN:chunkSize+GN,5)   :: tilde_op
        real(rp), dimension(1-GN:chunkSize+GN,8)     :: pri_1D
        real(rp), dimension(0:3,1-GN:chunkSize+GN,5) :: phi_arr
        real(rp), dimension(-3:3)   :: central_1_p

        integer :: i, j, k, l, m, kl, n, ii, iimax, str_z, end_z, kstr_z, kend_z, kgbl

        real(rp) :: ir, u_, v_, w_, p_, weight, inner_sum, hHn
        real(rp) :: hHnm1, inner_sum_m
        real(rp) :: phi1, phi2, phi3, phi4, idz
        real(rp) :: wk_wkl, rhom_log, eintm_log
        real(rp), parameter :: one4 = 1.0_rp/4.0_rp

        call StartProfRange("energy_preserving_AECTP_Z")

        !$acc parallel vector_length(vec_size) default(present) &
        !$acc private(tilde_op,pri_1D,phi_arr) &
        !$acc create(central_1_p)
        !$acc cache(central_1_p,tilde_op,pri_1D,phi_arr)
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

                       ir = 1.0_rp/phi1
                       u_ = phi2*ir
                       v_ = phi3*ir
                       w_ = phi4*ir

                       p_    = P(i+ii,j,kgbl)
                       pri_1D(k,1) = p_
                       pri_1D(k,2) = u_
                       pri_1D(k,3) = v_
                       pri_1D(k,4) = w_
                       pri_1D(k,5) = phi1
                       pri_1D(k,7) = pri_1D(k,1)/pri_1D(k,5)
                       pri_1D(k,8) = 1.0_rp/pri_1D(k,7)
                       pri_1D(k,6) = log(pri_1D(k,7))

                    enddo

                  !$acc loop vector collapse(2)
                  do k = 1-GN,end_z+GN-3
                     do l = 1,3

                        kl = k+l
                        wk_wkl   = pri_1D(k,4) + pri_1D(kl,4)
                        rhom_log = logmean_tp(pri_1D(k,5),pri_1D(kl,5))
                        eintm_log = aec_tp_internal_energy_mean(pri_1D(k,7), pri_1D(kl,7), &
                                                                pri_1D(k,8), pri_1D(kl,8), &
                                                                pri_1D(k,6), pri_1D(kl,6))
                        weight = one4*rhom_log*wk_wkl

                        tilde_op(l,k,1) = 2.0_rp*weight
                        tilde_op(l,k,2) = weight*(pri_1D(k,2) + pri_1D(kl,2))
                        tilde_op(l,k,3) = weight*(pri_1D(k,3) + pri_1D(kl,3))
                        tilde_op(l,k,4) = weight*wk_wkl + 0.5_rp*(pri_1D(k,1) + pri_1D(kl,1))
                        tilde_op(l,k,5) = weight*(pri_1D(k,2)*pri_1D(kl,2) &
                                          + pri_1D(k,3)*pri_1D(kl,3) &
                                          + pri_1D(k,4)*pri_1D(kl,4) &
                                          + 2.0_rp*eintm_log) &
                                          + 0.5_rp*(pri_1D(k,1)*pri_1D(kl,4) &
                                          + pri_1D(kl,1)*pri_1D(k,4))

                     enddo
                  enddo

                    !$acc loop vector
                    do k = kstr_z,kend_z
                       kgbl = k + str_z-1+GN
                       idz = zstep_i(kgbl)
                       do n = 1,5
                          hHn   = 0.0_rp
                          hHnm1 = 0.0_rp
                          do l = 1,3
                             inner_sum   = 0.0_rp
                             inner_sum_m = 0.0_rp
                             do m = 0,l-1
                                inner_sum   = inner_sum   + tilde_op(l,k-m  ,n)
                                inner_sum_m = inner_sum_m + tilde_op(l,k-1-m,n)
                             enddo
                             hHn   = hHn   + central_1_p(l)*inner_sum
                             hHnm1 = hHnm1 + central_1_p(l)*inner_sum_m
                          enddo
                          phi_arr(ii,k,n) =  - idz * (hHn - hHnm1)
                       enddo
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
end subroutine energy_preserving_AECTP_Z

end module energy_preserving_AECTP_module
