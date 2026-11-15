module energy_preserving_KEEP_module
use parameters_module
use storage_module
use mesh_module, only: dx_dxi, dx_deta, dy_dxi, dy_deta, Jac2D
use profiling_module
implicit none
private
public :: energy_preserving_KEEP_X, energy_preserving_KEEP_Y, energy_preserving_KEEP_Z

contains

subroutine energy_preserving_KEEP_X(lbx,ubx,istr,iend)

        implicit none
        integer, intent(in) :: lbx,ubx
        integer, intent(in) :: istr,iend

        integer, parameter :: vec_size  = 64
        integer, parameter :: chunkSize = vec_size - 2*GN

        real(rp), dimension(3,1-GN:chunkSize+GN,5) :: tilde_op
        real(rp), dimension(1-GN:chunkSize+GN,9)  :: pri_1D
        real(rp), dimension(-3:3)                  :: central_1_p

        integer :: i, j, k, l, il, m, n
        integer :: istr_x, iend_x, igbl, str_x, end_x

        real(rp) :: ir, u_, v_, w_, ek, p_, weight, inner_sum, hFn, Ucon
        real(rp) :: phi1, phi2, phi3, phi4, phi5, idx
        real(rp) :: rhoi_rhoil, ui_uil, einti_eintil
        real(rp) :: drho, deint, eint_fac
        real(rp) :: rho2, drho2, eint2, deint2
        real(rp) :: inm, hFm
        real(rp), parameter :: one8= 1.0_rp/8.0_rp
        real(rp), parameter :: one3= 1.0_rp/3.0_rp
        real(rp), parameter :: one5= 1.0_rp/5.0_rp

        call StartProfRange("energy_preserving_KEEP_X")

        !$acc parallel vector_length(vec_size) default(present) &
        !$acc private(tilde_op,pri_1D) &
        !$acc create(central_1_p)
        !$acc cache(central_1_p,tilde_op,pri_1D)
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

                     p_    = P(igbl,j,k)

                     pri_1D(i,1) = p_
                     pri_1D(i,2) = u_
                     pri_1D(i,3) = v_
                     pri_1D(i,4) = w_
                     pri_1D(i,5) = phi1
                     pri_1D(i,6) = phi5*ir - ek

                     pri_1D(i,7) = dy_deta(igbl,j)*u_ - dx_deta(igbl,j)*v_
                     pri_1D(i,8) = p_*dy_deta(igbl,j)
                     pri_1D(i,9) = p_*dx_deta(igbl,j)

                  enddo

                  !$acc loop vector collapse(2)
                  do i = 1-GN,end_x+GN-3
                     do l = 1,3

                        il = i+l
                        ui_uil = pri_1D(i,2) + pri_1D(il,2)
                        rhoi_rhoil = (pri_1D(i,5) + pri_1D(il,5))
                        drho = pri_1D(il,5) - pri_1D(i,5)
                        einti_eintil =  pri_1D(il,6) + pri_1D(i,6)
                        deint = pri_1D(il,6) - pri_1D(i,6)
                        Ucon = pri_1D(i,7) + pri_1D(il,7)
                        rho2   = rhoi_rhoil*rhoi_rhoil
                        drho2  = drho*drho
                        eint2  = einti_eintil*einti_eintil
                        deint2 = deint*deint
                        weight = one8 * rho2*rho2*rhoi_rhoil * Ucon &
                               / (rho2*(rho2 + one3*drho2) + one5*drho2*drho2)
                        eint_fac = einti_eintil * (eint2*(eint2 + one3*deint2) + one5*deint2*deint2) &
                               / (eint2*(eint2 + deint2) + deint2*deint2)

                        tilde_op(l,i,1) = 2*weight
                        tilde_op(l,i,2) = weight*ui_uil + 0.5_rp*(pri_1D(i,8) + pri_1D(il,8))
                        tilde_op(l,i,3) = weight*(pri_1D(i,3) + pri_1D(il,3)) &
                                          - 0.5_rp*(pri_1D(i,9) + pri_1D(il,9))
                        tilde_op(l,i,4) = weight*(pri_1D(i,4) + pri_1D(il,4))
                        tilde_op(l,i,5) = weight*(pri_1D(i,2) * pri_1D(il,2) &
                                          + pri_1D(i,3) * pri_1D(il,3) &
                                          + pri_1D(i,4) * pri_1D(il,4) &
                                          + eint_fac) &
                                          + 0.5_rp * (pri_1D(i,1) * pri_1D(il,7) + pri_1D(il,1) * pri_1D(i,7))

                     enddo
                  enddo

                  !$acc loop vector
                  do i = istr_x,iend_x
                     igbl = i + str_x-1+GN
                     idx = Jac2D(igbl,j)
                     do n = 1,5
                        hFn = 0.0_rp
                        hFm = 0.0_rp
                        do l = 1,3
                           inner_sum = 0.0_rp
                           inm       = 0.0_rp
                           do m = 0,l-1
                              inner_sum = inner_sum + tilde_op(l,i-m,n)
                              inm       = inm       + tilde_op(l,i-1-m,n)
                           enddo
                           hFn = hFn + central_1_p(l)*inner_sum
                           hFm = hFm + central_1_p(l)*inm
                        enddo
                        RHS(igbl,j,k,n) = -idx*(hFn - hFm)
                     enddo
                  enddo
               endif
             enddo
           enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine energy_preserving_KEEP_X

subroutine energy_preserving_KEEP_Y(lby,uby,jstr,jend)

        implicit none
        integer, parameter :: vec_size  = 64
        integer, parameter :: chunkSize = vec_size - 2*GN

        real(rp), dimension(3,1-GN:chunkSize+GN,5)   :: tilde_op
        real(rp), dimension(1-GN:chunkSize+GN,9)    :: pri_1D
        real(rp), dimension(0:3,1-GN:chunkSize+GN,5) :: phi_arr
        real(rp), dimension(-3:3)                    :: central_1_p

        integer, intent(in) :: lby, uby
        integer, intent(in) :: jstr,jend

        integer :: i, j, k, l, m, jl, n, str_y, end_y, ii, iimax, jstr_y, jend_y, jgbl

        real(rp) :: ir, u_, v_, w_, ek, p_, weight, inner_sum, hGn, Vcon
        real(rp) :: phi1, phi2, phi3, phi4, phi5, idy
        real(rp) :: rhoj_rhojl, vj_vjl, eintj_eintjl
        real(rp) :: drho, deint, eint_fac
        real(rp) :: rho2, drho2, eint2, deint2
        real(rp) :: inm, hGm
        real(rp), parameter :: one8= 1.0_rp/8.0_rp
        real(rp), parameter :: one3= 1.0_rp/3.0_rp
        real(rp), parameter :: one5= 1.0_rp/5.0_rp

        call StartProfRange("energy_preserving_KEEP_Y")

        !$acc parallel vector_length(vec_size) default(present) &
        !$acc private(tilde_op,pri_1D,phi_arr) &
        !$acc create(central_1_p)
        !$acc cache(central_1_p,tilde_op,pri_1D,phi_arr)
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
                        phi5 = phi_arr(ii,j,5)

                        ir = 1.0_rp/phi1
                        u_ = phi2*ir
                        v_ = phi3*ir
                        w_ = phi4*ir
                        ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)

                        p_    = P(i+ii,jgbl,k)

                        pri_1D(j,1) = p_
                        pri_1D(j,2) = u_
                        pri_1D(j,3) = v_
                        pri_1D(j,4) = w_
                        pri_1D(j,5) = phi1
                        pri_1D(j,6) = phi5*ir - ek

                        pri_1D(j,7) = -dy_dxi(i+ii,jgbl)*u_ + dx_dxi(i+ii,jgbl)*v_
                        pri_1D(j,8) = p_*dx_dxi(i+ii,jgbl)
                        pri_1D(j,9) = p_*dy_dxi(i+ii,jgbl)

                     enddo

                  !$acc loop vector collapse(2)
                  do j = 1-GN,end_y+GN-3
                     do l = 1,3

                           jl = j+l
                           vj_vjl = pri_1D(j,3) + pri_1D(jl,3)
                           rhoj_rhojl = pri_1D(j,5) + pri_1D(jl,5)
                           drho = pri_1D(jl,5) - pri_1D(j,5)
                           eintj_eintjl =  pri_1D(jl,6) + pri_1D(j,6)
                           deint = pri_1D(jl,6) - pri_1D(j,6)
                           Vcon = pri_1D(j,7) + pri_1D(jl,7)
                           rho2   = rhoj_rhojl*rhoj_rhojl
                           drho2  = drho*drho
                           eint2  = eintj_eintjl*eintj_eintjl
                           deint2 = deint*deint
                           weight = one8 * rho2*rho2*rhoj_rhojl * Vcon &
                                  / (rho2*(rho2 + one3*drho2) + one5*drho2*drho2)
                           eint_fac = eintj_eintjl * (eint2*(eint2 + one3*deint2) + one5*deint2*deint2) &
                                  / (eint2*(eint2 + deint2) + deint2*deint2)

                           tilde_op(l,j,1) = 2*weight
                           tilde_op(l,j,2) = weight*(pri_1D(j,2) + pri_1D(jl,2)) &
                                             - 0.5_rp*(pri_1D(j,9) + pri_1D(jl,9))
                           tilde_op(l,j,3) = weight*vj_vjl + 0.5_rp*(pri_1D(j,8) + pri_1D(jl,8))
                           tilde_op(l,j,4) = weight*(pri_1D(j,4) + pri_1D(jl,4))
                           tilde_op(l,j,5) = weight*(pri_1D(j,2) * pri_1D(jl,2) &
                                             + pri_1D(j,3) * pri_1D(jl,3) &
                                             + pri_1D(j,4) * pri_1D(jl,4) &
                                             + eint_fac) &
                                             + 0.5_rp * (pri_1D(j,1) * pri_1D(jl,7) + pri_1D(jl,1) * pri_1D(j,7))

                     enddo
                  enddo

                     !$acc loop vector
                     do j = jstr_y,jend_y
                        jgbl = j + str_y-1+GN
                        idy = Jac2D(i+ii,jgbl)
                        do n = 1,5
                           hGn = 0.0_rp
                           hGm = 0.0_rp
                           do l = 1,3
                              inner_sum = 0.0_rp
                              inm       = 0.0_rp
                              do m = 0,l-1
                                 inner_sum = inner_sum + tilde_op(l,j-m,n)
                                 inm       = inm       + tilde_op(l,j-1-m,n)
                              enddo
                              hGn = hGn + central_1_p(l)*inner_sum
                              hGm = hGm + central_1_p(l)*inm
                           enddo
                           phi_arr(ii,j,n) = -idy*(hGn - hGm)
                        enddo
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
end subroutine energy_preserving_KEEP_Y

subroutine energy_preserving_KEEP_Z(lbz,ubz,kstr,kend)

        implicit none

        integer, parameter :: vec_size = 64
        integer, parameter :: chunkSize = vec_size

        integer, intent(in) :: lbz,ubz
        integer, intent(in) :: kstr,kend

        real(rp), dimension(3,1-GN:chunkSize+GN,5)   :: tilde_op
        real(rp), dimension(1-GN:chunkSize+GN,6)     :: pri_1D
        real(rp), dimension(0:3,1-GN:chunkSize+GN,6) :: phi_arr
        real(rp), dimension(-3:3)   :: central_1_p

        integer :: i, j, k, l, m, kl, n, ii, iimax, str_z, end_z, kstr_z, kend_z, kgbl

        real(rp) :: ir, u_, v_, w_, ek, p_, weight, inner_sum, hHn
        real(rp) :: phi1, phi2, phi3, phi4, phi5, idz
        real(rp) :: rhok_rhokl, wk_wkl, eintk_eintkl
        real(rp) :: drho, deint, eint_fac
        real(rp) :: rho2, drho2, eint2, deint2
        real(rp) :: inm, hHm
        real(rp), parameter :: one8= 1.0_rp/8.0_rp
        real(rp), parameter :: one3= 1.0_rp/3.0_rp
        real(rp), parameter :: one5= 1.0_rp/5.0_rp

        call StartProfRange("energy_preserving_KEEP_Z")

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
                       phi_arr(ii,k,6) = P(i+ii,j,kgbl)
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

                       ir = 1.0_rp/phi1
                       u_ = phi2*ir
                       v_ = phi3*ir
                       w_ = phi4*ir
                       ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)

                       p_    = phi_arr(ii,k,6)

                       pri_1D(k,1) = p_
                       pri_1D(k,2) = u_
                       pri_1D(k,3) = v_
                       pri_1D(k,4) = w_
                       pri_1D(k,5) = phi1
                       pri_1D(k,6) = phi5*ir - ek

                    enddo

                  !$acc loop vector collapse(2)
                  do k = 1-GN,end_z+GN-3
                     do l = 1,3

                        kl = k+l
                        wk_wkl = pri_1D(k,4) + pri_1D(kl,4)
                        rhok_rhokl = pri_1D(k,5) + pri_1D(kl,5)
                        drho = pri_1D(kl,5) - pri_1D(k,5)
                        eintk_eintkl = pri_1D(kl,6) + pri_1D(k,6)
                        deint = pri_1D(kl,6) - pri_1D(k,6)
                        rho2   = rhok_rhokl*rhok_rhokl
                        drho2  = drho*drho
                        eint2  = eintk_eintkl*eintk_eintkl
                        deint2 = deint*deint
                        weight = one8 * rho2*rho2*rhok_rhokl * wk_wkl &
                               / (rho2*(rho2 + one3*drho2) + one5*drho2*drho2)
                        eint_fac = eintk_eintkl * (eint2*(eint2 + one3*deint2) + one5*deint2*deint2) &
                               / (eint2*(eint2 + deint2) + deint2*deint2)

                        tilde_op(l,k,1) = 2.0_rp * weight
                        tilde_op(l,k,2) = weight  * (pri_1D(k,2) + pri_1D(kl,2))
                        tilde_op(l,k,3) = weight  * (pri_1D(k,3) + pri_1D(kl,3))
                        tilde_op(l,k,4) = weight  * wk_wkl &
                                + 0.5_rp * (pri_1D(k,1) + pri_1D(kl,1))
                        tilde_op(l,k,5) = weight * ((pri_1D(k,2) * pri_1D(kl,2)) &
                                + (pri_1D(k,3) * pri_1D(kl,3)) &
                                + (pri_1D(k,4) * pri_1D(kl,4)) &
                                +  eint_fac) &
                                + 0.5_rp * (pri_1D(kl,4) * pri_1D(k,1) &
                                + pri_1D(k,4) * pri_1D(kl,1))

                     enddo
                  enddo

                    !$acc loop vector
                    do k = kstr_z,kend_z
                       kgbl = k + str_z-1+GN
                       idz = zstep_i(kgbl)
                       do n = 1,5
                          hHn = 0.0_rp
                          hHm = 0.0_rp
                          do l = 1,3
                             inner_sum = 0.0_rp
                             inm       = 0.0_rp
                             do m = 0,l-1
                                inner_sum = inner_sum + tilde_op(l,k-m,n)
                                inm       = inm       + tilde_op(l,k-1-m,n)
                             enddo
                             hHn = hHn + central_1_p(l)*inner_sum
                             hHm = hHm + central_1_p(l)*inm
                          enddo
                          phi_arr(ii,k,n) = -idz*(hHn - hHm)
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
end subroutine energy_preserving_KEEP_Z

end module energy_preserving_KEEP_module
