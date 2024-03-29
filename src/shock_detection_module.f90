module shock_detection_module
use parameters_module
use mpi_comm_module
use storage_module
use omp_lib
use profiling_module

implicit none
private
public compute_hybrid_weno_flag, spread_hybrid_weno_flag,  compute_ducros_sensor_lower_wall

contains
subroutine compute_hybrid_weno_flag
        implicit none

        call StartProfRange("compute_hybrid_weno_flag") 

        if    (trim(sensor) == 'density') then
          call compute_density_shock_sensor

        elseif(trim(sensor) == 'density_jump') then
          call compute_density_jump_shock_sensor

        elseif(trim(sensor) == 'ducros' ) then
          call compute_ducros_sensor

        else
          print*, ' Shock detection', trim(sensor), ' is not implemented!'
          stop

        endif

        call EndProfRange
        call StartProfRange("compute_hybrid_weno_flag_MPI") 
        
        if(mpi_flag) then
        call mpi_share_int1(mpi_comm_cart,type_send_flag,type_recv_flag,&
                my_neighbour, &
                i13D_bfr_send_E, i13D_bfr_send_W, i13D_bfr_recv_E, i13D_bfr_recv_W, &
                i13D_bfr_send_N, i13D_bfr_send_S, i13D_bfr_recv_N, i13D_bfr_recv_S, &
                i13D_bfr_send_B, i13D_bfr_send_F, i13D_bfr_recv_B, i13D_bfr_recv_F, &
                weno_flag)
        endif

        call spread_hybrid_weno_flag

        if(mpi_flag) then
        call mpi_share_int1(mpi_comm_cart,type_send_flag,type_recv_flag,&
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
! -----------------------------------------------------------------------------------
!
!       This subroutine compute the hybrid weno flag with a density based shock sensor
!
! -----------------------------------------------------------------------------------

        implicit none
        real(rp), parameter :: hw_toll = 0.3_rp
        real(rp)            :: delta_rho1, delta_rho2, delta_rho3
        real(rp)            :: rho_, shock_sensor
        integer             :: i,j,k

        selectcase(dims)
        
          !
          ! === 2D CASES
          !
          case(2)
          k = 1
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(2)
          do j    = sy-1,ey+1
             do i = sx-1,ex+1
                
                 rho_ = phi(i,j,k,1)
                 delta_rho1 = (phi(i+1,j,k,1) - rho_)*xstep_i(i)
                 delta_rho2 = (phi(i,j+1,k,1) - rho_)*ystep_i(j)

                 shock_sensor = max(abs(delta_rho1),abs(delta_rho2))

                 weno_flag(i,j,k) = weno_smooth
                 if(shock_sensor >= hw_toll) weno_flag(i,j,k) = weno_shock

              enddo
          enddo
          !$acc end parallel

          !
          ! === 3D CASES
          !
          case(3)
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do k       = sz-1,ez+1
             do j    = sy-1,ey+1
                do i = sx-1,ex+1
                  
                   rho_ = phi(i,j,k,1)
                   delta_rho1 = (phi(i+1,j,k,1) - rho_)*xstep_i(i)
                   delta_rho2 = (phi(i,j+1,k,1) - rho_)*ystep_i(j)
                   delta_rho3 = (phi(i,j,k+1,1) - rho_)*zstep_i(k)

                   shock_sensor = max(abs(delta_rho1),abs(delta_rho2),abs(delta_rho3))
                
                   weno_flag(i,j,k) = weno_smooth
                   if(shock_sensor >= hw_toll) weno_flag(i,j,k) = weno_shock

                enddo
             enddo
          enddo
          !$acc end parallel

        endselect

        return
end subroutine compute_density_shock_sensor


subroutine compute_density_jump_shock_sensor
! -----------------------------------------------------------------------------------
!
!       This subroutine compute the hybrid weno flag with a density based shock sensor
!
! -----------------------------------------------------------------------------------
        implicit none
        !real(rp), parameter       :: hw_toll = 0.5_rp ! Double mach reflection
        !real(rp), parameter       :: hw_toll = 0.05_rp ! Four quadrant A, B
        real(rp), parameter :: eps1 = 0.01_rp, eps2 = 0.05_rp
        real(rp)            :: hw_toll, shock_sensor
        real(rp)            :: delta_rho1, delta_rho2, delta_rho3
        real(rp)            :: r_, ru, rv, rw, ir, rek, p_, M
        integer             :: i,j,k

        selectcase(dims)

          !
          ! === 2D CASES
          !
          case(2)
          k = 1
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(2)
          do j    = sy-1,ey+1
             do i = sx-1,ex+1

                 r_ = phi(i,j,k,1)
                 ru = phi(i,j,k,2)
                 rv = phi(i,j,k,3)
                 ir = 1.0_rp/r_

                 rek= 0.5_rp*ir*(ru*ru + rv*rv)
                 p_ = (gamma0-1._rp)*(phi(i,j,k,5) - rek)
                 M  = sqrt(2*rek/(gamma0*p_))
                      
                 hw_toll = Mach*eps1 + Mach*(eps2-eps1)/(1.0_rp+(M+0.2_rp)**20)
                
                 delta_rho1 = phi(i+1,j,k,1) - r_
                 delta_rho2 = phi(i,j+1,k,1) - r_

                 shock_sensor = max(abs(delta_rho1),abs(delta_rho2))
                
                 weno_flag(i,j,k) = weno_smooth
                 if(shock_sensor > hw_toll) weno_flag(i,j,k) = weno_shock

              enddo
          enddo
          !$acc end parallel

          !
          ! === 3D CASES
          !
          case(3)
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do k       = sz-1,ez+1
             do j    = sy-1,ey+1
                do i = sx-1,ex+1
                  
                   r_ = phi(i,j,k,1)
                   ru = phi(i,j,k,2)
                   rv = phi(i,j,k,3)
                   rw = phi(i,j,k,4)
                   ir = 1.0_rp/r_

                   rek= 0.5_rp*ir*(ru*ru + rv*rv + rw*rw)
                   p_ = (gamma0-1._rp)*(phi(i,j,k,5) - rek)
                   M  = sqrt(2*rek/(gamma0*p_))
                        
                   hw_toll = Mach*eps1 + Mach*(eps2-eps1)/(1.0_rp+(M+0.2_rp)**20)

                   delta_rho1 = phi(i+1,j,k,1) - r_
                   delta_rho2 = phi(i,j+1,k,1) - r_
                   delta_rho3 = phi(i,j,k+1,1) - r_

                   shock_sensor = max(abs(delta_rho1),abs(delta_rho2),abs(delta_rho3))

                   weno_flag(i,j,k) = weno_smooth
                   if(shock_sensor > hw_toll) weno_flag(i,j,k) = weno_shock

                enddo
            enddo
          enddo
          !$acc end parallel

        endselect

        return
end subroutine compute_density_jump_shock_sensor



subroutine compute_ducros_sensor
        
        implicit none
        real(rp) :: ducros, div_, div2, vor2, eps
        real(rp) :: vor_x, vor_y, vor_z
        real(rp) :: ux, uy, uz, vx, vy, vz, wx, wy, wz
        real(rp) :: idx, idy, idz
        real(rp) :: cl1
        integer  :: i,j,k, l

        eps = u_inf**2

        selectcase(dims)

          !
          ! === 2D CASES
          !
          case(2)
          k = 1
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(2)
          do j    = sy-1,ey+1
             do i = sx-1,ex+1

                idx = xstep_i(i)
                idy = ystep_i(j)

                ux = 0.0_rp
                uy = 0.0_rp
                !
                vx = 0.0_rp
                vy = 0.0_rp
                !$acc loop seq
                do l = -3, 3

                   cl1 = central_1(l)
                
                   !
                   ux = ux + cl1 * U(i+l,j,k)
                   uy = uy + cl1 * U(i,j+l,k)
                   !
                   vx = vx + cl1 * V(i+l,j,k)
                   vy = vy + cl1 * V(i,j+l,k)

                enddo
                ux = ux*idx
                vx = vx*idx
                !
                uy = uy*idy
                vy = vy*idy
                
                ! divergency
                div_ = ux + vy
                   
                ! vorticity
                vor_z = vx - uy

                div2 = div_*div_
                vor2 = vor_z*vor_z

                ducros = max(-div_/sqrt((div2 + vor2 + eps)), 0.0_rp)
        
                weno_flag(i,j,k) = weno_smooth
                if (ducros > ducrosToll) weno_flag(i,j,k) = weno_shock

             enddo
          enddo
          !$acc end parallel
       
          !
          ! === 3D CASES
          !
          case(3)
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do k       = sz,ez
             do j    = sy,ey
                do i = sx,ex

                   idx = xstep_i(i)
                   idy = ystep_i(j)
                   idz = zstep_i(k)

                   ux = 0.0_rp
                   uy = 0.0_rp
                   uz = 0.0_rp
                   !
                   vx = 0.0_rp
                   vy = 0.0_rp
                   vz = 0.0_rp
                   !
                   wx = 0.0_rp
                   wy = 0.0_rp
                   wz = 0.0_rp
                    
                   !$acc loop seq
                   do l = -3, 3

                      cl1 = central_1(l)
                      !
                      ux = ux + cl1 * U(i+l,j,k)
                      uy = uy + cl1 * U(i,j+l,k)
                      uz = uz + cl1 * U(i,j,k+l)
                      !
                      vx = vx + cl1 * V(i+l,j,k)
                      vy = vy + cl1 * V(i,j+l,k)
                      vz = vz + cl1 * V(i,j,k+l)
                      !
                      wx = wx + cl1 * W(i+l,j,k)
                      wy = wy + cl1 * W(i,j+l,k)
                      wz = wz + cl1 * W(i,j,k+l)

                   enddo
                   ux = ux*idx
                   vx = vx*idx
                   wx = wx*idx
                   !
                   uy = uy*idy
                   vy = vy*idy
                   wy = wy*idy
                   !
                   uz = uz*idz
                   vz = vz*idz
                   wz = wz*idz
                
                   ! divergency
                   div_ = ux + vy + wz
                        
                   ! vorticity
                   vor_x = wy - vz
                   vor_y = uz - wx
                   vor_z = vx - uy

                   div2 = div_*div_
                   vor2 = vor_x*vor_x + vor_y*vor_y + vor_z*vor_z

                   ducros = max(-div_/sqrt((div2 + vor2 + eps)), 0.0_rp)

                   SSENSOR(i,j,k) = ducros

                   weno_flag(i,j,k) = weno_smooth
                   if (ducros > ducrosToll) weno_flag(i,j,k) = weno_shock

                enddo
             enddo
          enddo
          !$acc end parallel

        endselect

        return
end subroutine compute_ducros_sensor





subroutine compute_ducros_sensor_lower_wall

        implicit none
        integer, parameter :: jpos = 3

        real(rp) :: ducros, div_, div2, vor2, eps
        real(rp) :: vor_x, vor_y, vor_z
        real(rp) :: ux, uy, uz, vx, vy, vz, wx, wy, wz
        real(rp) :: irx, iry, irz
        real(rp) :: idx, idy, idz, cl1
        integer  :: i,j,k,l,fdR

        fdR = central_fd_order/2
        eps = u_inf**2

        selectcase(dims)

          !
          ! === 2D CASES
          !
          case(2)
          write(*,*) "Ducros lower wall is not implemented in 2D"
          stop
       
          !
          ! === 3D CASES
          !
          case(3)

          !$omp parallel do collapse(3) default(private), &
          !$omp shared(phi,xstep_i,ystep_i,zstep_i,SSENSOR), &
          !$omp shared(weno,sx,ex,sy,ey,sz,ez,ducrosToll), &
          !$omp shared(central_1,fdR,eps)
          do k       = sz,ez
             do j    = sy,ey
                do i = sx,ex

                   if(j <= jpos) then

                     SSENSOR(i,j,k) = 0.0_rp
                     weno_flag(i,j,k) = weno_smooth

                   else

                     idx = xstep_i(i)
                     idy = ystep_i(j)
                     idz = zstep_i(k)

                     ux = 0.0_rp
                     uy = 0.0_rp
                     uz = 0.0_rp
                     !
                     vx = 0.0_rp
                     vy = 0.0_rp
                     vz = 0.0_rp
                     !
                     wx = 0.0_rp
                     wy = 0.0_rp
                     wz = 0.0_rp
                     do l = -fdR, fdR

                        cl1 = central_1(l)
                
                        irx = 1.0_rp/phi(i+l,j,k,1)
                        iry = 1.0_rp/phi(i,j+l,k,1)
                        irz = 1.0_rp/phi(i,j,k+l,1)
                        !
                        ux = ux + cl1 * phi(i+l,j,k,2)*irx
                        uy = uy + cl1 * phi(i,j+l,k,2)*iry
                        uz = uz + cl1 * phi(i,j,k+l,2)*irz
                        !
                        vx = vx + cl1 * phi(i+l,j,k,3)*irx
                        vy = vy + cl1 * phi(i,j+l,k,3)*iry
                        vz = vz + cl1 * phi(i,j,k+l,3)*irz
                        !
                        wx = wx + cl1 * phi(i+l,j,k,4)*irx
                        wy = wy + cl1 * phi(i,j+l,k,4)*iry
                        wz = wz + cl1 * phi(i,j,k+l,4)*irz

                     enddo
                     ux = ux*idx
                     vx = vx*idx
                     wx = wx*idx
                     !
                     uy = uy*idy
                     vy = vy*idy
                     wy = wy*idy
                     !
                     uz = uz*idz
                     vz = vz*idz
                     wz = wz*idz
                
                     ! divergency
                     div_ = ux + vy + wz
                          
                     ! vorticity
                     vor_x = wy - vz
                     vor_y = uz - wx
                     vor_z = vx - uy

                     div2 = div_*div_
                     vor2 = vor_x*vor_x + vor_y*vor_y + vor_z*vor_z

                     ducros = max(-div_/sqrt((div2 + vor2 + eps)), 0.0_rp)

                     SSENSOR(i,j,k) = ducros

                     weno_flag(i,j,k) = weno_smooth
                     if (ducros > ducrosToll) weno_flag(i,j,k) = weno_shock

                   endif

                enddo
             enddo
          enddo
          !$omp end parallel do

        endselect

        return
end subroutine compute_ducros_sensor_lower_wall









subroutine spread_hybrid_weno_flag
        implicit none
        integer :: i,j,k
        integer :: bL, bR

        bL = weno_num - 1
        bR = weno_num 
       
        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do k       = sz-1,ez
           do j    = sy-1,ey
              do i = sx-1,ex
                 
                 weno_flag_xyz(i,j,k) = minval(weno_flag(i-bL:i+BR,j,k))
                 weno_flag_xyz(i,j,k) = weno_flag_xyz(i,j,k) + 2*minval(weno_flag(i,j-bL:j+bR,k))
                 weno_flag_xyz(i,j,k) = weno_flag_xyz(i,j,k) + 4*minval(weno_flag(i,j,k-bL:k+bR))

              enddo
           enddo
        enddo
        !$acc end parallel


        return
end subroutine spread_hybrid_weno_flag



end module shock_detection_module
