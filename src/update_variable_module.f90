module update_variable_module
use parameters_module
use storage_module
use profiling_module

implicit none
private
public update_all, limit_positivity, report_positivity_limiter, finalize_positivity_limiter

integer(8), save :: lim_count(1) = 0_8
!$acc declare create(lim_count)
logical , save :: lim_count_init = .false.
real(rp), save :: lim_calls = 0.0_rp

real(rp), save :: lim_snd = 0.0_rp, lim_rcv = 0.0_rp
real(rp), save :: lim_calls_posted = 0.0_rp
integer , save :: lim_req = 0, lim_it_posted = 0
logical , save :: lim_pending = .false.

contains
subroutine update_all

        use thermodynamics_module

        implicit none
        real(rp) :: r_, ir, u_, v_, w_, p_, T_, mu_, e_, cp_
        real(rp) :: tnr, dT_, h_t, invT, r2
        real(rp) :: sutA, sutB
        integer  :: i, j, k, iter

        integer , parameter :: NR_MAXIT = 50
        real(rp), parameter :: NR_TOLL2 = 1.0e-12_rp
        real(rp), parameter :: NR_TMIN  = 1.0e-2_rp

        call StartProfRange("update_all")

        sutA = suth_S / Tref
        sutB = 1.0_rp + sutA

        if(thermally_perfect_gas > 0) then

          !$acc parallel default(present)
          !$acc loop gang,vector collapse(3)
          do k       = lbz,ubz
             do j    = lby,uby
                do i = lbx,ubx

                   r_ = phi(i,j,k,1)
                   ir = 1.0_rp/r_
                   u_ = phi(i,j,k,2)*ir
                   v_ = phi(i,j,k,3)*ir
                   w_ = phi(i,j,k,4)*ir

                   e_ = phi(i,j,k,5)*ir - 0.5_rp*(u_*u_ + v_*v_ + w_*w_)

                   U(i,j,k) = u_
                   V(i,j,k) = v_
                   W(i,j,k) = w_

                   tnr = e_ / cv_Tref
                   do iter = 1, NR_MAXIT
                      invT = 1.0_rp / tnr
                      if (tnr <= T_cold) then
                         r2  = aC_adim_p4*tnr**4
                         cp_ = aC_adim__0 + r2
                         h_t = hC_Tref + tnr*(aC_adim__0 + 0.2_rp*r2)
                      else if (tnr <= T_switch) then
                         r2  = aL_adim_m2 * invT
                         cp_ = invT*(aL_adim_m1 + r2) + aL_adim__0 + &
                               tnr*(aL_adim_p1 + tnr*(aL_adim_p2 + tnr*(aL_adim_p3 + tnr*aL_adim_p4)))
                         h_t = hL_Tref - r2 + aL_adim_m1*log(tnr) + &
                               tnr*(aL_adim__0 + tnr*(hL_adim_p1 + tnr*(hL_adim_p2 + tnr*(hL_adim_p3 + tnr*hL_adim_p4))))
                      else
                         r2  = aH_adim_m2 * invT
                         cp_ = invT*(aH_adim_m1 + r2) + aH_adim__0 + &
                               tnr*(aH_adim_p1 + tnr*(aH_adim_p2 + tnr*(aH_adim_p3 + tnr*aH_adim_p4)))
                         h_t = hH_Tref - r2 + aH_adim_m1*log(tnr) + &
                               tnr*(aH_adim__0 + tnr*(hH_adim_p1 + tnr*(hH_adim_p2 + tnr*(hH_adim_p3 + tnr*hH_adim_p4))))
                      end if
                      dT_ = -(e_ - h_t + tnr) / (1.0_rp - cp_)
                      tnr = max(tnr + dT_, NR_TMIN)
                      if (dT_*dT_ <= NR_TOLL2) exit
                   end do
                   T_ = tnr
                   p_ = r_ * T_

                   P(i,j,k) = p_
                   T(i,j,k) = T_

                   invT = 1.0_rp / T_
                   if (T_ <= T_cold) then
                      r2  = aC_adim_p4*T_**4
                      cp_ = aC_adim__0 + r2
                      ENTHALPY(i,j,k) = hC_Tref + T_*(aC_adim__0 + 0.2_rp*r2)
                   else if (T_ <= T_switch) then
                      cp_ = invT*(aL_adim_m1 + invT*aL_adim_m2) + aL_adim__0 + &
                            T_*(aL_adim_p1 + T_*(aL_adim_p2 + T_*(aL_adim_p3 + T_*aL_adim_p4)))
                      ENTHALPY(i,j,k) = hL_Tref + invT*hL_adim_m2 + aL_adim_m1*log(T_) + &
                            T_*(aL_adim__0 + T_*(hL_adim_p1 + T_*(hL_adim_p2 + T_*(hL_adim_p3 + T_*hL_adim_p4))))
                   else
                      cp_ = invT*(aH_adim_m1 + invT*aH_adim_m2) + aH_adim__0 + &
                            T_*(aH_adim_p1 + T_*(aH_adim_p2 + T_*(aH_adim_p3 + T_*aH_adim_p4)))
                      ENTHALPY(i,j,k) = hH_Tref + invT*hH_adim_m2 + aH_adim_m1*log(T_) + &
                            T_*(aH_adim__0 + T_*(hH_adim_p1 + T_*(hH_adim_p2 + T_*(hH_adim_p3 + T_*hH_adim_p4))))
                   end if
                   GAMMA_T(i,j,k) = cp_ / (cp_ - 1.0_rp)

                   if (vis_flag == 0) then
                      mu_ = sqrt(T_)
                      mu_ = mu_*sqrt(mu_)
                   else if (vis_flag == 1) then
                      mu_ = T_*sqrt(T_) * sutB/(T_ + sutA)
                   else if (vis_flag == 2) then
                      mu_ = 1.0_rp
                   else
                      mu_ = 0.0_rp
                   end if

                   VIS(i,j,k) = mu_
                   LMD(i,j,k) = cp_ * iPrandtl * mu_

                enddo
             enddo
          enddo
          !$acc end parallel
        else

          !$acc parallel default(present)
          !$acc loop gang,vector collapse(3)
          do k       = lbz,ubz
             do j    = lby,uby
                do i = lbx,ubx

                   r_ = phi(i,j,k,1)
                   ir = 1.0_rp/r_
                   u_ = phi(i,j,k,2)*ir
                   v_ = phi(i,j,k,3)*ir
                   w_ = phi(i,j,k,4)*ir

                   p_ = cv_i*(phi(i,j,k,5) - 0.5_rp*r_*(u_*u_ + v_*v_ + w_*w_))

                   U(i,j,k) = u_
                   V(i,j,k) = v_
                   W(i,j,k) = w_

                   T_ = p_*ir
                   P(i,j,k) = p_
                   T(i,j,k) = T_

                   ENTHALPY(i,j,k) = cp0 * T_

                   if (vis_flag == 0) then
                      mu_ = sqrt(T_)
                      mu_ = mu_*sqrt(mu_)
                   else if (vis_flag == 1) then
                      mu_ = T_*sqrt(T_) * sutB/(T_ + sutA)
                   else if (vis_flag == 2) then
                      mu_ = 1.0_rp
                   else
                      mu_ = 0.0_rp
                   end if

                   VIS(i,j,k) = mu_
                   LMD(i,j,k) = k_inf * mu_

                enddo
             enddo
          enddo
          !$acc end parallel

        endif

        call EndProfRange

        return
end subroutine update_all

subroutine limit_positivity
        implicit none
        real(rp) :: r_, ke, ket, f, rfl, efl
        integer  :: i, j, k, l

        call StartProfRange("limit_positivity")

        rfl = limiter_rho_floor
        efl = limiter_e_floor

        if(.not. lim_count_init) then
          lim_count(1) = 0_8
          !$acc update device(lim_count)
          lim_count_init = .true.
        endif

        !$acc parallel default(present)
        !$acc loop gang,vector collapse(3)
        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex

                 l  = 0
                 r_ = phi(i,j,k,1)
                 if(.not.(r_ >= rfl)) then
                    r_ = rfl
                    phi(i,j,k,1) = r_
                    l  = 1
                 endif

                 ke  = 0.5_rp*( phi(i,j,k,2)*phi(i,j,k,2)   &
                              + phi(i,j,k,3)*phi(i,j,k,3)   &
                              + phi(i,j,k,4)*phi(i,j,k,4) )/r_
                 ket = phi(i,j,k,5) - r_*efl

                 if(.not.(ket >= ke)) then
                    l = 1
                    if(ket >= 0.0_rp .and. ke > 0.0_rp) then
                       f = sqrt(ket/ke)
                       phi(i,j,k,2) = f*phi(i,j,k,2)
                       phi(i,j,k,3) = f*phi(i,j,k,3)
                       phi(i,j,k,4) = f*phi(i,j,k,4)
                    else
                       if(.not.(ke >= 0.0_rp .and. ke < huge(1.0_rp))) then
                          phi(i,j,k,2) = 0.0_rp
                          phi(i,j,k,3) = 0.0_rp
                          phi(i,j,k,4) = 0.0_rp
                          ke = 0.0_rp
                       endif
                       phi(i,j,k,5) = r_*efl + ke
                    endif
                 endif

                 if(l == 1) then
                    !$acc atomic update
                    lim_count(1) = lim_count(1) + 1_8
                 endif

              enddo
           enddo
        enddo
        !$acc end parallel

        lim_calls = lim_calls + 1.0_rp

        call EndProfRange

        return
end subroutine limit_positivity

subroutine report_positivity_limiter
        implicit none
        integer :: ierr

        call harvest_limiter_window(allow_stop = .true.)

        !$acc update host(lim_count)
        lim_snd      = real(lim_count(1), rp)
        lim_count(1) = 0_8
        !$acc update device(lim_count)
        lim_calls_posted = lim_calls
        lim_calls        = 0.0_rp
        lim_it_posted    = it
        call MPI_Iallreduce(lim_snd, lim_rcv, 1, MPI_RP, MPI_SUM, mpi_comm_cart, lim_req, ierr)
        lim_pending = .true.

        return
end subroutine report_positivity_limiter

subroutine harvest_limiter_window(allow_stop)
        implicit none
        logical, intent(in) :: allow_stop
        real(rp) :: frac
        integer  :: ierr

        if(.not. lim_pending) return

        call MPI_Wait(lim_req, MPI_STATUS_IGNORE, ierr)
        lim_pending = .false.

        if(lim_rcv > 0.5_rp) then

          frac = lim_rcv/(real(nx,rp)*real(ny,rp)*real(nz,rp)*max(lim_calls_posted,1.0_rp))

          if(rank == root) write(*,'(a,i8,a,i12,a,es10.3)') &
             ' [limiter] it=', lim_it_posted, '  limited cells =', int(lim_rcv,8), '  frac/stage =', frac

          if(allow_stop .and. limiter_abort_frac > 0.0_rp .and. frac > limiter_abort_frac) then
            if(rank == root) then
              write(*,'(a)')                ' URANOS: POSITIVITY LIMITER RUNAWAY -- cooperative stop'
              write(*,'(a,es10.3,a,es10.3)') '   limited per-stage cell fraction ', frac, &
                                             '  exceeds limiter_abort_frac ', limiter_abort_frac
              write(*,'(a)')                '   the limiter is masking a real problem (under-resolution or blow-up).'
            endif
            istop = 1
          endif

        endif

        return
end subroutine harvest_limiter_window

subroutine finalize_positivity_limiter
        implicit none

        call harvest_limiter_window(allow_stop = .false.)

        return
end subroutine finalize_positivity_limiter

end module update_variable_module
