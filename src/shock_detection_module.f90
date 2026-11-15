module shock_detection_module

        use parameters_module
        use mpi_comm_module
        use storage_module
        use omp_lib
        use profiling_module
        use mesh_module, only: dxi_dx, dxi_dy, deta_dx, deta_dy, zstep_i

        implicit none
        private
        public :: compute_hybrid_weno_flag, spread_hybrid_weno_flag

contains

subroutine compute_hybrid_weno_flag
        implicit none

        call StartProfRange("compute_hybrid_weno_flag")

        selectcase(trim(sensor))

          case('density')
            call compute_density_shock_sensor

          case('density_jump')
            call compute_density_jump_shock_sensor

          case('ducros')
            call compute_ducros_sensor

          case('weno')
            call compute_pure_weno_flag

          case default
            write(*,'(a)')      ' FATAL ERROR!'
            write(*,'(a,a,a)')  ' Shock detection ', trim(sensor), &
                                ' is not implemented. The program will be stopped'
            call secure_stop

        endselect

        call force_weno_extrapolation_boundary

        call EndProfRange
        call StartProfRange("compute_hybrid_weno_flag_MPI")

        if (mpi_flag) then
           call mpi_share_int1(mpi_comm_cart,type_send_flag,type_recv_flag, &
                   my_neighbour, &
                   i13D_bfr_send_E, i13D_bfr_send_W, i13D_bfr_recv_E, i13D_bfr_recv_W, &
                   i13D_bfr_send_N, i13D_bfr_send_S, i13D_bfr_recv_N, i13D_bfr_recv_S, &
                   i13D_bfr_send_B, i13D_bfr_send_F, i13D_bfr_recv_B, i13D_bfr_recv_F, &
                   weno_flag)
        endif

        call spread_hybrid_weno_flag

        if (mpi_flag) then
           call mpi_share_int1(mpi_comm_cart,type_send_flag,type_recv_flag, &
                   my_neighbour, &
                   i13D_bfr_send_E, i13D_bfr_send_W, i13D_bfr_recv_E, i13D_bfr_recv_W, &
                   i13D_bfr_send_N, i13D_bfr_send_S, i13D_bfr_recv_N, i13D_bfr_recv_S, &
                   i13D_bfr_send_B, i13D_bfr_send_F, i13D_bfr_recv_B, i13D_bfr_recv_F, &
                   weno_flag_xyz)
        endif

        call EndProfRange

        return
end subroutine compute_hybrid_weno_flag

subroutine compute_density_shock_sensor
        implicit none
        real(rp), parameter :: hw_toll = 0.3_rp

        real(rp) :: r_csi, r_eta, r_zta, cl1
        real(rp) :: drho_x, drho_y, drho_z, sshock
        integer  :: i, j, k, l

        selectcase(dims)

          case(2)
            k = 1
            !$acc parallel default(present)
            !$acc loop gang, vector collapse(2)
            do j    = sy-1,ey+1
               do i = sx-1,ex+1

                  r_csi = 0.0_rp
                  r_eta = 0.0_rp
                  !$acc loop seq
                  do l = -3, 3
                     cl1   = central_1(l)
                     r_csi = r_csi + cl1 * phi(i+l,j,k,1)
                     r_eta = r_eta + cl1 * phi(i,j+l,k,1)
                  enddo

                  drho_x = r_csi*dxi_dx(i,j) + r_eta*deta_dx(i,j)
                  drho_y = r_csi*dxi_dy(i,j) + r_eta*deta_dy(i,j)

                  sshock = max(abs(drho_x), abs(drho_y))

                  if (sshock >= hw_toll) then
                     weno_flag(i,j,k) = weno_shock
                  else
                     weno_flag(i,j,k) = weno_smooth
                  endif

               enddo
            enddo
            !$acc end parallel

          case(3)
            !$acc parallel default(present)
            !$acc loop gang, vector collapse(3)
            do k       = sz-1,ez+1
               do j    = sy-1,ey+1
                  do i = sx-1,ex+1

                     r_csi = 0.0_rp
                     r_eta = 0.0_rp
                     r_zta = 0.0_rp
                     !$acc loop seq
                     do l = -3, 3
                        cl1   = central_1(l)
                        r_csi = r_csi + cl1 * phi(i+l,j,k,1)
                        r_eta = r_eta + cl1 * phi(i,j+l,k,1)
                        r_zta = r_zta + cl1 * phi(i,j,k+l,1)
                     enddo

                     drho_x = r_csi*dxi_dx(i,j) + r_eta*deta_dx(i,j)
                     drho_y = r_csi*dxi_dy(i,j) + r_eta*deta_dy(i,j)
                     drho_z = r_zta*zstep_i(k)

                     sshock = max(abs(drho_x), abs(drho_y), abs(drho_z))

                     if (sshock >= hw_toll) then
                        weno_flag(i,j,k) = weno_shock
                     else
                        weno_flag(i,j,k) = weno_smooth
                     endif

                  enddo
               enddo
            enddo
            !$acc end parallel

        endselect

        return
end subroutine compute_density_shock_sensor

subroutine compute_density_jump_shock_sensor
        implicit none
        real(rp), parameter :: eps1 = 0.01_rp, eps2 = 0.05_rp

        real(rp) :: M_eps1, M_deps
        real(rp) :: r_, u_, v_, w_, M, hw_toll
        real(rp) :: drho_x, drho_y, drho_z, sshock
        integer  :: i, j, k

        M_eps1 = Mach*eps1
        M_deps = Mach*(eps2 - eps1)

        selectcase(dims)

          case(2)
            k = 1
            !$acc parallel default(present)
            !$acc loop gang, vector collapse(2)
            do j    = sy-1,ey+1
               do i = sx-1,ex+1

                  r_ = phi(i,j,k,1)
                  u_ = U(i,j,k)
                  v_ = V(i,j,k)
                  w_ = W(i,j,k)

                  M       = sqrt((u_*u_ + v_*v_ + w_*w_)/(GAMMA_T(i,j,k)*T(i,j,k)))
                  hw_toll = M_eps1 + M_deps/(1.0_rp + (M + 0.2_rp)**20)

                  drho_x = phi(i+1,j,k,1) - r_
                  drho_y = phi(i,j+1,k,1) - r_

                  sshock = max(abs(drho_x), abs(drho_y))

                  if (sshock > hw_toll) then
                     weno_flag(i,j,k) = weno_shock
                  else
                     weno_flag(i,j,k) = weno_smooth
                  endif

               enddo
            enddo
            !$acc end parallel

          case(3)
            !$acc parallel default(present)
            !$acc loop gang, vector collapse(3)
            do k       = sz-1,ez+1
               do j    = sy-1,ey+1
                  do i = sx-1,ex+1

                     r_ = phi(i,j,k,1)
                     u_ = U(i,j,k)
                     v_ = V(i,j,k)
                     w_ = W(i,j,k)

                     M       = sqrt((u_*u_ + v_*v_ + w_*w_)/(GAMMA_T(i,j,k)*T(i,j,k)))
                     hw_toll = M_eps1 + M_deps/(1.0_rp + (M + 0.2_rp)**20)

                     drho_x = phi(i+1,j,k,1) - r_
                     drho_y = phi(i,j+1,k,1) - r_
                     drho_z = phi(i,j,k+1,1) - r_

                     sshock = max(abs(drho_x), abs(drho_y), abs(drho_z))

                     if (sshock > hw_toll) then
                        weno_flag(i,j,k) = weno_shock
                     else
                        weno_flag(i,j,k) = weno_smooth
                     endif

                  enddo
               enddo
            enddo
            !$acc end parallel

        endselect

        return
end subroutine compute_density_jump_shock_sensor

subroutine compute_ducros_sensor
        implicit none
        real(rp) :: u_csi, u_eta, u_zta, v_csi, v_eta, v_zta, w_csi, w_eta, w_zta, cl1
        real(rp) :: csi_x, csi_y, eta_x, eta_y, zta_z
        real(rp) :: ux, uy, uz, vx, vy, vz, wx, wy, wz
        real(rp) :: div_, vor_x, vor_y, vor_z
        real(rp) :: div2, vor2, inv_norm, ducros, eps
        integer  :: i, j, k, l

        if (dims /= 3) then
           write(*,'(a)') ' FATAL ERROR!'
           write(*,'(a)') ' Ducros shock sensor is implemented for 3D only.'
           call secure_stop
        endif

        eps = u_inf*u_inf

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex

                 u_csi = 0.0_rp
                 u_eta = 0.0_rp
                 u_zta = 0.0_rp

                 v_csi = 0.0_rp
                 v_eta = 0.0_rp
                 v_zta = 0.0_rp

                 w_csi = 0.0_rp
                 w_eta = 0.0_rp
                 w_zta = 0.0_rp

                 !$acc loop seq
                 do l = -3, 3
                    cl1 = central_1(l)

                    u_csi = u_csi + cl1 * U(i+l,j,k)
                    u_eta = u_eta + cl1 * U(i,j+l,k)
                    u_zta = u_zta + cl1 * U(i,j,k+l)

                    v_csi = v_csi + cl1 * V(i+l,j,k)
                    v_eta = v_eta + cl1 * V(i,j+l,k)
                    v_zta = v_zta + cl1 * V(i,j,k+l)

                    w_csi = w_csi + cl1 * W(i+l,j,k)
                    w_eta = w_eta + cl1 * W(i,j+l,k)
                    w_zta = w_zta + cl1 * W(i,j,k+l)
                 enddo

                 csi_x = dxi_dx (i,j)
                 csi_y = dxi_dy (i,j)
                 eta_x = deta_dx(i,j)
                 eta_y = deta_dy(i,j)
                 zta_z = zstep_i(k)

                 ux = u_csi*csi_x + u_eta*eta_x
                 uy = u_csi*csi_y + u_eta*eta_y
                 uz = u_zta*zta_z

                 vx = v_csi*csi_x + v_eta*eta_x
                 vy = v_csi*csi_y + v_eta*eta_y
                 vz = v_zta*zta_z

                 wx = w_csi*csi_x + w_eta*eta_x
                 wy = w_csi*csi_y + w_eta*eta_y
                 wz = w_zta*zta_z

                 div_  = ux + vy + wz
                 vor_x = wy - vz
                 vor_y = uz - wx
                 vor_z = vx - uy

                 div2 = div_*div_
                 vor2 = vor_x*vor_x + vor_y*vor_y + vor_z*vor_z

                 inv_norm = 1.0_rp/sqrt(div2 + vor2 + eps)
                 ducros   = max(-div_*inv_norm, 0.0_rp)

                 SSENSOR(i,j,k) = ducros
                 if (ducros > ducrosToll) then
                    weno_flag(i,j,k) = weno_shock
                 else
                    weno_flag(i,j,k) = weno_smooth
                 endif

              enddo
           enddo
        enddo
        !$acc end parallel

        return
end subroutine compute_ducros_sensor

subroutine compute_pure_weno_flag
        implicit none
        integer :: i, j, k

        if (dims == 2) then
           k = 1
           !$acc parallel default(present)
           !$acc loop gang, vector collapse(2)
           do j    = sy-1, ey+1
              do i = sx-1, ex+1
                 weno_flag(i,j,k) = weno_shock
              enddo
           enddo
           !$acc end parallel
        else
           !$acc parallel default(present)
           !$acc loop gang, vector collapse(3)
           do k       = sz-1, ez+1
              do j    = sy-1, ey+1
                 do i = sx-1, ex+1
                    weno_flag(i,j,k) = weno_shock
                 enddo
              enddo
           enddo
           !$acc end parallel
        endif

        return
end subroutine compute_pure_weno_flag

subroutine force_weno_extrapolation_boundary
        implicit none
        integer :: i, j, k, n, i0, j0, k0
        integer :: klo, khi

        klo = sz-1 ; khi = ez+1
        if (dims == 2) then
           klo = 1 ; khi = 1
        endif

        do n = 1, 2
           if (.not. extrapolation_bc(n)) cycle
           i0 = merge(sx, ex, n == 1)
           !$acc parallel default(present)
           !$acc loop gang, vector collapse(2)
           do k    = klo, khi
              do j = sy-1, ey+1
                 weno_flag(i0,j,k) = weno_shock
              enddo
           enddo
           !$acc end parallel
        enddo

        do n = 1, 2
           if (.not. extrapolation_bc(n+2)) cycle
           j0 = merge(sy, ey, n == 1)
           !$acc parallel default(present)
           !$acc loop gang, vector collapse(2)
           do k    = klo, khi
              do i = sx-1, ex+1
                 weno_flag(i,j0,k) = weno_shock
              enddo
           enddo
           !$acc end parallel
        enddo

        do n = 1, 2
           if (.not. extrapolation_bc(n+4)) cycle
           k0 = merge(sz, ez, n == 1)
           !$acc parallel default(present)
           !$acc loop gang, vector collapse(2)
           do j    = sy-1, ey+1
              do i = sx-1, ex+1
                 weno_flag(i,j,k0) = weno_shock
              enddo
           enddo
           !$acc end parallel
        enddo

        return
end subroutine force_weno_extrapolation_boundary

subroutine spread_hybrid_weno_flag
        implicit none
        integer :: i, j, k
        integer :: bL, bR

        bL = weno_num - 1
        bR = weno_num

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do k       = sz-1,ez
           do j    = sy-1,ey
              do i = sx-1,ex

                 weno_flag_xyz(i,j,k) = minval(weno_flag(i-bL:i+bR,j,k))
                 weno_flag_xyz(i,j,k) = weno_flag_xyz(i,j,k) + 2*minval(weno_flag(i,j-bL:j+bR,k))
                 weno_flag_xyz(i,j,k) = weno_flag_xyz(i,j,k) + 4*minval(weno_flag(i,j,k-bL:k+bR))

              enddo
           enddo
        enddo
        !$acc end parallel

        return
end subroutine spread_hybrid_weno_flag

end module shock_detection_module
