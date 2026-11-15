module time_module
use parameters_module
use storage_module
use mpi_module
use shock_detection_module
use profiling_module
use output_probes     , only: flush_probes
use mesh_module, only: dxi_dx, dxi_dy, deta_dx, deta_dy, m_xi, m_eta, iz_zta, ztastep_i

implicit none
private
public init_runge_kutta, runge_kutta, compute_dt, last_iteration

real(rp) :: machineTmax        = -1.0_rp
logical  :: wallclock_resolved = .false.

contains

subroutine init_runge_kutta
        implicit none
        integer :: i, j, k, l

        call StartProfRange("init_runge_kutta")

        time_n = time

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = lbx,ubx
                    phi_n(i,j,k,l) = phi(i,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        if(hybrid_weno) call compute_hybrid_weno_flag

        call EndProfRange

        return
end subroutine init_runge_kutta

subroutine runge_kutta
        implicit none
        real(rp) :: aik, bik
        integer  :: i,j,k,l

        call StartProfRange("runge_kutta")

        time = time_n + Dt * c_rk(ik)

        aik = a_rk(ik)
        bik = b_rk(ik)

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = sz,ez
              do    j = sy,ey
                 do i = sx,ex
                    phi(i,j,k,l) = aik * phi_n(i,j,k,l) &
                                 + bik * (phi(i,j,k,l) + dt * RHS(i,j,k,l))
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        if(ik == n_step) call check_finite_time_step

        call EndProfRange

        return
end subroutine runge_kutta

subroutine compute_dt
        implicit none
        real(rp), parameter :: four_thirds = 4.0_rp/3.0_rp
        real(rp) :: my_dt, dt_conv, dt_diff
        real(rp) :: max_conv, max_diff
        real(rp) :: gbl_conv, gbl_diff
        real(rp) :: c_, uu, vv, ww, ucon, vcon, wcon, sloc, m2, irho, dloc, m_zta
        integer  :: i,j,k

        call StartProfRange("compute_dt")

          if(viscous) then

            max_conv = 0.0_rp
            max_diff = 0.0_rp
            !$acc parallel default(present) copy(max_conv,max_diff)
            !$acc loop gang, vector collapse(3) &
            !$acc   reduction(max:max_conv,max_diff)
            do       k = sz,ez
               do    j = sy,ey
                  do i = sx,ex
                     c_   = sqrt(GAMMA_T(i,j,k)*T(i,j,k))
                     uu   = U(i,j,k)
                     vv   = V(i,j,k)
                     ww   = W(i,j,k)
                     m_zta = ztastep_i(k)*iz_zta(k)
                     ucon = uu*dxi_dx (i,j) + vv*dxi_dy (i,j)
                     vcon = uu*deta_dx(i,j) + vv*deta_dy(i,j)
                     wcon = ww*m_zta
                     sloc =      abs(ucon) + m_xi (i,j)*c_   &
                            +       abs(vcon) + m_eta(i,j)*c_   &
                            +       abs(wcon) + m_zta      *c_

                     irho = 1.0_rp/phi(i,j,k,1)
                     m2   = m_xi(i,j)**2 + m_eta(i,j)**2 + m_zta**2
                     dloc = max(four_thirds*VIS(i,j,k),                     &
                                LMD(i,j,k)*(GAMMA_T(i,j,k) - 1.0_rp)) * irho

                     max_conv = max(max_conv, sloc)
                     max_diff = max(max_diff, dloc*m2)
                  enddo
               enddo
            enddo
            !$acc end parallel

          else

            max_conv = 0.0_rp
            !$acc parallel default(present) copy(max_conv)
            !$acc loop gang, vector collapse(3) reduction(max:max_conv)
            do       k = sz,ez
               do    j = sy,ey
                  do i = sx,ex
                     c_   = sqrt(GAMMA_T(i,j,k)*T(i,j,k))
                     uu   = U(i,j,k)
                     vv   = V(i,j,k)
                     ww   = W(i,j,k)
                     m_zta = ztastep_i(k)*iz_zta(k)
                     ucon = uu*dxi_dx (i,j) + vv*dxi_dy (i,j)
                     vcon = uu*deta_dx(i,j) + vv*deta_dy(i,j)
                     wcon = ww*m_zta
                     sloc =      abs(ucon) + m_xi (i,j)*c_   &
                            +       abs(vcon) + m_eta(i,j)*c_   &
                            +       abs(wcon) + m_zta      *c_
                     max_conv = max(max_conv, sloc)
                  enddo
               enddo
            enddo
            !$acc end parallel

            max_diff = 0.0_rp

          endif

        call MPI_allreduce(max_conv, gbl_conv, 1, MPI_RP, MPI_MAX, mpi_comm_cart, mpi_err)
        call MPI_allreduce(max_diff, gbl_diff, 1, MPI_RP, MPI_MAX, mpi_comm_cart, mpi_err)

        dt_conv = CFL/gbl_conv
        if(gbl_diff > 0.0_rp) then
          dt_diff = Fourier/(mu_inf*gbl_diff)
          dt      = 1.0_rp/(1.0_rp/dt_conv + 1.0_rp/dt_diff)
        else
          dt      = dt_conv
        endif
        my_dt = dt

        if(dt >= 0.0_rp .and. dt <= tmax .and. time + dt >= tmax) then
           dt    = tmax - time
           istop = 1
        endif

        call EndProfRange

        return
end subroutine compute_dt

subroutine check_finite_time_step
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use, intrinsic :: iso_fortran_env, only: error_unit
        implicit none
        logical :: sane

        sane = ieee_is_finite(time) .and. ieee_is_finite(dt) .and. &
               dt > 0.0_rp .and. dt <= tmax
        if(logical_CFL) then
           sane = sane .and. time <= tmax + 1.0e3_rp*spacing(tmax)
        else
           sane = sane .and. time <= tmax + dt
        endif
        if(sane) return

        if(rank == root) then
           write(error_unit,'(A)')        ''
           write(error_unit,'(A)')        ' URANOS: SIMULATION DIVERGED -- non-finite or unphysical time/dt (blow-up)'
           write(error_unit,'(A,I0)')     '   iteration : ', it
           write(error_unit,'(A,ES12.5)') '   time      : ', time
           write(error_unit,'(A,ES12.5)') '   dt        : ', dt
           flush(error_unit)
        endif

        call flush_probes
        call MPI_Barrier(mpi_comm_cart, mpi_err)
        call secure_stop

end subroutine check_finite_time_step

subroutine last_iteration(sTime, istop)
        use mpi
        implicit none
        real(rp), intent(in)  :: sTime
        integer , intent(out) :: istop

        real(rp)      :: actualTime, budget
        integer       :: ierr
        logical       :: from_batch

        call StartProfRange("last_iteration")

        istop = 0

        if(.not. wallclock_resolved) then

           wallclock_resolved = .true.

           if(rank == root) then

              call batch_wallclock_budget(budget, from_batch)

              if(from_batch) then

                 machineTmax = (MPI_WTIME() - sTime) + budget

                 write(*,'(A,F10.2)') '[last_iteration] batch deadline, stopping after [s] : ', machineTmax

              else

                 machineTmax = -1.0_rp

              endif

           endif

           call MPI_BCAST(machineTmax, 1, MPI_RP, root, mpi_comm_cart, ierr)

        endif

        if(machineTmax > 0.0_rp) then
           actualTime = MPI_WTIME()
           if(actualTime - sTime >= machineTmax) istop = 1
           call MPI_allreduce(MPI_IN_PLACE, istop, 1, MPI_INTEGER, MPI_MAX, mpi_comm_cart, ierr)
           if(istop == 1) then
             if(rank == root) then
               write(*,'(A)')        '[last_iteration] Wall-clock limit reached: stopping at next safe point.'
               write(*,'(A, F10.2)') '  Elapsed time [s] : ', actualTime - sTime
               write(*,'(A, F10.2)') '  Limit       [s]  : ', machineTmax
             endif
           endif
        endif

        call EndProfRange

        return
end subroutine last_iteration

subroutine batch_wallclock_budget(budget, ok)
        implicit none
        real(rp), intent(out) :: budget
        logical , intent(out) :: ok

        real(8), parameter :: margin_frac = 0.1_8
        real(8), parameter :: margin_max  = 1800.0_8

        character(64) :: buffer
        real(8)       :: t_end, t_start, t_now, span, margin
        integer       :: buflen, status, ios

        budget = -1.0_rp
        ok     = .false.

        call get_environment_variable('SLURM_JOB_END_TIME', buffer, buflen, status)
        if(status /= 0 .or. buflen == 0) return
        read(buffer, *, iostat = ios) t_end
        if(ios /= 0) return

        t_now = epoch_seconds()
        if(t_end - t_now <= 0.0_rp) return

        span = t_end - t_now
        call get_environment_variable('SLURM_JOB_START_TIME', buffer, buflen, status)
        if(status == 0 .and. buflen > 0) then
          read(buffer, *, iostat = ios) t_start
          if(ios == 0 .and. t_end - t_start > 0.0_rp) span = t_end - t_start
        endif

        margin = min(margin_frac * span, margin_max)
        budget = (t_end - t_now) - margin

        if(budget <= 0.0_rp) return

        ok = .true.

        return
end subroutine batch_wallclock_budget

function epoch_seconds() result(secs)
        implicit none
        real(rp) :: secs

        integer :: v(8), y, m, era, yoe, doy, doe, days

        call date_and_time(values = v)

        y = v(1)
        m = v(2)
        if(m <= 2) y = y - 1

        era = (y - modulo(y, 400)) / 400
        yoe = y - era*400
        if(v(2) > 2) then
          doy = (153*(v(2) - 3) + 2)/5 + v(3) - 1
        else
          doy = (153*(v(2) + 9) + 2)/5 + v(3) - 1
        endif
        doe  = yoe*365 + yoe/4 - yoe/100 + doy
        days = era*146097 + doe - 719468

        secs =   real(days, rp) * 86400.0_rp &
             & + real(v(5), rp) * 3600.0_rp  &
             & + real(v(6), rp) * 60.0_rp    &
             & + real(v(7), rp)              &
             & + real(v(8), rp) * 1.0e-3_rp  &
             & - real(v(4), rp) * 60.0_rp

        return
end function epoch_seconds

end module time_module
