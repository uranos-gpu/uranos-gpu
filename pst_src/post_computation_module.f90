module post_computation_module
use storage_module
use post_storage_module
use fluid_functions_module
use shock_detection_module
use post_solutions_module
use integration_module
implicit none

private
public compute_fields, compute_urms, compute_vrms, compute_wrms, compute_ekt

contains
subroutine compute_fields
        implicit none
        
        if(velocity)            call compute_velocity
        if(pressure)            call compute_pressure
        if(temperature)         call compute_temperature
        if(mach_)               call compute_mach_number
        if(vorticity)           call compute_vorticity
        if(speed_div)           call compute_divergency
        if(hybrid_weno)         call compute_hybrid_weno_flag
        if(vorticity_magnitude) call compute_vorticity_magnitude
        if(qcriterion)          call compute_q_criterion
        if(sdensity)            call compute_schlieren_density
        if(charactheristic)     call compute_charactheristic_x
        if(viscous)             call compute_molecular_viscosity
        if(MixedVelocity)       call compute_mixed_velocities

        return
end subroutine compute_fields

subroutine compute_velocity
        implicit none
        real(rp) :: i_rho
        
        !$omp parallel do collapse(3) default(private), shared(phi,U,V,W)
        do k = lbz,ubz
           do j = lby,uby
              do i = lbx, ubx
                          
                 i_rho = 1._rp/phi(i,j,k,1)

                 U(i,j,k) = phi(i,j,k,2) * i_rho
                 V(i,j,k) = phi(i,j,k,3) * i_rho
                 W(i,j,k) = phi(i,j,k,4) * i_rho

              enddo
           enddo
        enddo
        !$omp end parallel do
        return
end subroutine compute_velocity

subroutine compute_pressure
        implicit none
        real(rp) :: rho_, i_rho, u_, v_, w_

        !$omp parallel do collapse(3) default(private), shared(phi,P)
        do k = lbz,ubz
           do j = lby,uby
              do i = lbx, ubx

                 rho_  = phi(i,j,k,1)
                 i_rho = 1._rp/rho_

                 u_ = phi(i,j,k,2)*i_rho
                 v_ = phi(i,j,k,3)*i_rho
                 w_ = phi(i,j,k,4)*i_rho

                 P(i,j,k)=(gamma0-1._rp)*(phi(i,j,k,5)-0.5_rp*rho_*(u_**2 + v_**2 + w_**2))

              enddo
           enddo
        enddo
        !$omp end parallel do
        
        return
end subroutine compute_pressure

subroutine compute_temperature
        implicit none
        real(rp) :: rho_, i_rho, u_, v_, w_, p_

        !$omp parallel do collapse(3) default(private), shared(phi,T)
        do k = lbz,ubz
           do j = lby,uby
              do i = lbx, ubx

                 rho_  = phi(i,j,k,1)
                 i_rho = 1._rp/rho_

                 u_ = phi(i,j,k,2)*i_rho
                 v_ = phi(i,j,k,3)*i_rho
                 w_ = phi(i,j,k,4)*i_rho

                 p_ = (gamma0-1._rp)*(phi(i,j,k,5)-0.5_rp*rho_*(u_**2 + v_**2 + w_**2))

                 T(i,j,k) = p_ * i_rho

              enddo
           enddo
        enddo
        !$omp end parallel do

        return
end subroutine compute_temperature

subroutine compute_mach_number
        implicit none
        real(rp) :: u_, v_, w_, rho_, i_rho, p_, T_, i_c

        !$omp parallel do collapse(3) default(private), shared(phi,Mach_x,Mach_y,Mach_z)
        do k = lbz,ubz
           do j = lby,uby
              do i = lbx, ubx
                 
                 rho_  = phi(i,j,k,1)
                 i_rho = 1._rp/rho_

                 u_ = phi(i,j,k,2) * i_rho
                 v_ = phi(i,j,k,3) * i_rho
                 w_ = phi(i,j,k,4) * i_rho

                 p_ = (gamma0-1._rp)*(phi(i,j,k,5)-0.5_rp*rho_*(u_**2 + v_**2 + w_**2))

                 T_ = p_ * i_rho

                 i_c = 1._rp/(sqrt(gamma0*T_))

                 Mach_x(i,j,k) = u_ * i_c
                 Mach_y(i,j,k) = v_ * i_c
                 Mach_z(i,j,k) = w_ * i_c

              enddo
           enddo
        enddo
        !$omp end parallel do
        return
end subroutine compute_mach_number


subroutine compute_vorticity
        implicit none
        real(rp), dimension(3) :: vor

        !$omp parallel do collapse(3) default(private), shared(vor_x,vor_y,vor_z), &
        !$omp shared(sx,ex,sy,ey,sz,ez)
        do k = sz-1,ez+1
           do j = sy-1,ey+1
              do i = sx-1,ex+1

                 vor = vort(i,j,k)

                 vor_x(i,j,k) = vor(1)
                 vor_y(i,j,k) = vor(2)
                 vor_z(i,j,k) = vor(3)

              enddo
           enddo
        enddo
        !$omp end parallel do

        return
end subroutine compute_vorticity


subroutine compute_divergency
        implicit none

        !$omp parallel do collapse(3) default(private), shared(div), &
        !$omp shared(sx,ex,sy,ey,sz,ez)
        do k = sz-1,ez+1
           do j = sy-1,ey+1
              do i = sx-1,ex+1

                 div(i,j,k) = divergency(i,j,k)

              enddo
           enddo
        enddo
        !$omp end parallel do
        return
end subroutine compute_divergency


subroutine compute_schlieren_density
! ---------------------------------------------------------------------
!
!       This subroutine compute the numerical schlieren of a variable
!       Ref: Boukharfane "A combined ghost-point-forcing / direct-forcing 
!                         immersed boundary method for compressible
!                         flow simulations", Computer and Fluids 2018
!
! ---------------------------------------------------------------------
        implicit none
        real(rp), dimension(3)    :: grad_rho, i_step
        integer                   :: s
        real(rp)                  :: norm_d, max_norm, cs, alpha
        
        alpha = 5.0_rp
        ! --- Cycling all cell to get the max((sqrt(var_x**2 + var_y**2)) --- !
        max_norm = 0._rp

        do k = sz-1,ez+1
           do j = sy-1,ey+1
              do i = sx-1,ex+1

                 i_step(1) = xstep_i(i)
                 i_step(2) = ystep_i(j)
                 i_step(3) = zstep_i(k)

                 grad_rho(:) = 0.0_rp
                 do s = -central_fd_order/2, central_fd_order/2
                    
                    cs = central_1(s)

                    grad_rho(1) = grad_rho(1) + i_step(1) * cs * phi(i+s,j,k,1)
                    grad_rho(2) = grad_rho(2) + i_step(2) * cs * phi(i,j+s,k,1)
                    grad_rho(3) = grad_rho(3) + i_step(3) * cs * phi(i,j,k+s,1)

                 enddo

                 norm_d    = sqrt(grad_rho(1)**2 + grad_rho(2)**2 + grad_rho(3)**2)
                 max_norm  = max(norm_d, max_norm)

              enddo
           enddo
        enddo
        
        
        ! --- computing the numerical Schlieren --- !
        do k = sz-1,ez+1
           do j = sy-1,ey+1
              do i = sx-1,ex+1

                 i_step(1) = xstep_i(i)
                 i_step(2) = ystep_i(j)
                 i_step(3) = zstep_i(k)

                 grad_rho(:) = 0._rp
                 do s = -central_fd_order/2, central_fd_order/2
                    
                    cs = central_1(s)

                    grad_rho(1) = grad_rho(1) + i_step(1) * cs * phi(i+s,j,k,1)
                    grad_rho(2) = grad_rho(2) + i_step(2) * cs * phi(i,j+s,k,1)
                    grad_rho(3) = grad_rho(3) + i_step(3) * cs * phi(i,j,k+s,1)

                 enddo

                 norm_d    = sqrt(grad_rho(1)**2 + grad_rho(2)**2 + grad_rho(3)**2)

                 s_density(i,j,k) = exp(-alpha* norm_d / max_norm)
              enddo
           enddo
        enddo
        return
endsubroutine compute_schlieren_density


subroutine compute_charactheristic_x
        implicit none
        real(rp), dimension(5) :: lambda
        real(rp)               :: i_dx, r_, u_, c
        real(rp)               :: r_x, p_x, u_x, v_x, w_x, cs
        integer                :: s, is
        
        if(velocity    .eqv. .false. .or. &
           pressure    .eqv. .false. .or. &
           temperature .eqv. .false. ) then
           print*, 'velocity, pressure and temperature must be computed for characteristics.'
        endif

        do k = sz-1,ez+1
           do j = sy-1,ey+1
              do i = sx-1,ex+1

                 i_dx = xstep_i(i)
                
                 r_ = phi(i,j,k,1)
                 u_ =   U(i,j,k)
                 c = sqrt(gamma0*T(i,j,k))

                 lambda(1) = u_ - c
                 lambda(2) = u_
                 lambda(3) = u_
                 lambda(4) = u_
                 lambda(5) = u_ + c

                 r_x = 0.0_rp
                 p_x = 0.0_rp
                 u_x = 0.0_rp
                 v_x = 0.0_rp
                 w_x = 0.0_rp
                 do s = -central_fd_order/2, central_fd_order/2
                    
                    cs = central_1(s)
                    is = i+s

                    r_x = r_x + i_dx * cs * phi(is,j,k,1)
                    p_x = p_x + i_dx * cs *   P(is,j,k)
                    u_x = u_x + i_dx * cs *   U(is,j,k)
                    v_x = v_x + i_dx * cs *   V(is,j,k)
                    w_x = w_x + i_dx * cs *   W(is,j,k)

                 enddo

                 L1(i,j,k) = lambda(1) * (p_x - r_*c*u_x)
                 L2(i,j,k) = lambda(2) * (c**2*r_x - p_x)
                 L3(i,j,k) = lambda(3) * v_x
                 L4(i,j,k) = lambda(4) * w_x
                 L5(i,j,k) = lambda(5) * (p_x + r_*c*u_x)

              enddo
           enddo
        enddo


        return
end subroutine compute_charactheristic_x


subroutine compute_vorticity_magnitude
! -------------------------------------------------------
!       Computation of the vorticity magnitude
! -------------------------------------------------------

        use fluid_functions_module, only: vort

        implicit none
        real(rp), dimension(3) :: vor
       
        !$omp parallel do collapse(3) default(private), &
        !$omp shared(sx,ex,sy,ey,sz,ez,abs_vor)
        do k = sz-1,ez+1
           do j = sy-1,ey+1
              do i = sx-1,ex+1

                 vor = vort(i,j,k)

                 abs_vor(i,j,k) = sqrt(vor(1)*vor(1) + vor(2)*vor(2) + vor(3)*vor(3))

              enddo
           enddo
        enddo
        !$omp end parallel do

        return
end subroutine compute_vorticity_magnitude




subroutine compute_q_criterion

        use storage_module     , only: U,V,W,central_1, central_fd_order
        use post_storage_module, only: Q_CRITERION
        use mpi_module         , only: sx,ex, sy,ey, sz,ez

        implicit none
        real(rp), dimension(3)   :: i_st, cl1
        real(rp), dimension(3)   :: du, dv, dw
        real(rp), dimension(3,3) :: grad_v, grad_vT
        real(rp), dimension(3,3) :: omega_tensor, strai_tensor
        real(rp)                 :: norm_omega2, norm_strai2, q_value
        integer                  :: i,j,k,l,fR

        fR = central_fd_order/2
        
        !$omp workshare
        q_criterion(:,:,:) = 0.0_rp
        !$omp end workshare
        
        !$omp parallel do collapse(3) default(private), shared(sx,ex,sy,ey,sz,ez,u_inf), &
        !$omp shared(xstep_i, ystep_i, zstep_i,U,V,W,Q_criterion, central_1,fR)
        do       k = sz,ez
           do    j = sy,ey
              do i = sx,ex

                 ! grid steps
                 i_st(1) = xstep_i(i)
                 i_st(2) = ystep_i(j)
                 i_st(3) = zstep_i(k)


                 ! === first derivatives
                 du = 0.0_rp
                 dv = 0.0_rp
                 dw = 0.0_rp
                 do l = 1,fR

                    cl1(:) = i_st(:) * central_1(l)

                    du(1) = du(1) + cl1(1) * (U(i+l,j,k) - U(i-l,j,k))
                    du(2) = du(2) + cl1(2) * (U(i,j+l,k) - U(i,j-l,k))
                    du(3) = du(3) + cl1(3) * (U(i,j,k+l) - U(i,j,k-l))

                    dv(1) = dv(1) + cl1(1) * (V(i+l,j,k) - V(i-l,j,k))
                    dv(2) = dv(2) + cl1(2) * (V(i,j+l,k) - V(i,j-l,k))
                    dv(3) = dv(3) + cl1(3) * (V(i,j,k+l) - V(i,j,k-l))

                    dw(1) = dw(1) + cl1(1) * (W(i+l,j,k) - W(i-l,j,k))
                    dw(2) = dw(2) + cl1(2) * (W(i,j+l,k) - W(i,j-l,k))
                    dw(3) = dw(3) + cl1(3) * (W(i,j,k+l) - W(i,j,k-l))

                 enddo

                 ! === calculating u grandient
                 grad_v(1,:)  = du
                 grad_v(2,:)  = dv
                 grad_v(3,:)  = dw

                 grad_vT(:,1) = du
                 grad_vT(:,2) = dv
                 grad_vT(:,3) = dw

                 omega_tensor = 0.5_rp*(grad_v - grad_vT)
                 strai_tensor = 0.5_rp*(grad_v + grad_vT)
        
                 ! compute q invariant
                 norm_omega2 = sum(abs(omega_tensor(:,:))**2)
                 norm_strai2 = sum(abs(strai_tensor(:,:))**2)

                 q_value = 0.5_rp*(norm_omega2 - norm_strai2)/(u_inf**2)

                 Q_CRITERION(i,j,k) = q_value

              enddo
           enddo
        enddo
        !$omp end parallel do



        return
end subroutine compute_q_criterion

subroutine compute_molecular_viscosity
        implicit none
        integer  :: i,j,k
        real(rp) :: T_
        
        if(.not. allocated(T)) then
          stop ' Molecular viscosity computation requires temperature'
        endif


        !$omp parallel do collapse(3) default(private), shared(lbx,ubx, lby,uby, lbz, ubz, VIS,T)
        do k       = lbz,ubz
           do j    = lby,uby
              do i = lbx,ubx

                 T_  = T(i,j,k)
                 VIS(i,j,k) = laminar_viscosity(T_,tref,vis_flag)

              enddo
           enddo
        enddo
        !$omp end parallel do


        return
end subroutine compute_molecular_viscosity


subroutine compute_mixed_velocities
        implicit  none

        real(rp) :: u_
        integer  :: i,j,k

        if(.not.velocity) stop 'Mixed velocity computation requires velocity'

        do k       = lbz,ubz
           do j    = lby,uby
              do i = lbx,ubx
                 
                 u_ = U(i,j,k)

                 UV(i,j,k) = u_*V(i,j,k)
                 UW(i,j,k) = u_*W(i,j,k)

              enddo
           enddo
        enddo

        return
end subroutine compute_mixed_velocities


subroutine compute_urms(wrk)
        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: wrk

        ! local
        integer  :: i,j,k
        real(rp) :: mean_rhuu, meanR_u_u, meantau11, urms

        if(.not. allocated(wrk)) stop ' wrk array is not allocate for compute_urms'
        
        k = 1
        do    j = lby,uby
           do i = lbx,ubx

              mean_rhuu = Re_av%ruu(i,j,k)
              meanR_u_u = Re_av%rhu(i,j,k)*Re_av%rhu(i,j,k)/Re_av%rho(i,j,k)
              meantau11 = mean_rhuu - meanR_u_u
             
              urms = sqrt(meantau11)/u_inf

              wrk(i,j,k) = urms

           enddo
        enddo

        return
end subroutine compute_urms


subroutine compute_vrms(wrk)
        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: wrk

        ! local
        integer  :: i,j,k
        real(rp) :: mean_rhvv, meanR_v_v, meantau22, vrms

        if(.not. allocated(wrk)) stop ' wrk array is not allocate for compute_urms'
        
        k = 1
        do    j = lby,uby
           do i = lbx,ubx

              mean_rhvv = Re_av%rvv(i,j,k)
              meanR_v_v = Re_av%rhv(i,j,k)*Re_av%rhv(i,j,k)/Re_av%rho(i,j,k)
              meantau22 = mean_rhvv - meanR_v_v
             
              vrms = sqrt(meantau22)/u_inf

              wrk(i,j,k) = vrms

           enddo
        enddo

        return
end subroutine compute_vrms


subroutine compute_wrms(wrk)
        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: wrk

        ! local
        integer  :: i,j,k
        real(rp) :: mean_rhww, meanR_w_w, meantau33, wrms

        if(.not. allocated(wrk)) stop ' wrk array is not allocate for compute_urms'
        
        k = 1
        do    j = lby,uby
           do i = lbx,ubx

              mean_rhww = Re_av%rww(i,j,k)
              meanR_w_w = Re_av%rhw(i,j,k)*Re_av%rhw(i,j,k)/Re_av%rho(i,j,k)
              meantau33 = mean_rhww - meanR_w_w
             
              wrms = sqrt(meantau33)/u_inf

              wrk(i,j,k) = wrms

           enddo
        enddo

        return
end subroutine compute_wrms


subroutine compute_ekt(wrk)
        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: wrk

        ! local
        integer  :: i,j,k
        real(rp) :: mean_rhuu, mean_rhvv, mean_rhww
        real(rp) :: meanR_u_u, meanR_v_v, meanR_w_w
        real(rp) :: ekt, r_, ir, ruu, rvv, rww

        if(.not. allocated(wrk)) stop ' wrk array is not allocate for compute_urms'
        
        k = 1
        do    j = lby,uby
           do i = lbx,ubx

              r_ = Re_av%rho(i,j,k)
              ir = 1.0_rp/r_

              mean_rhuu = Re_av%ruu(i,j,k)
              mean_rhvv = Re_av%rvv(i,j,k)
              mean_rhww = Re_av%rww(i,j,k)

              meanR_u_u = Re_av%rhu(i,j,k)*Re_av%rhu(i,j,k)*ir
              meanR_v_v = Re_av%rhv(i,j,k)*Re_av%rhv(i,j,k)*ir
              meanR_w_w = Re_av%rhw(i,j,k)*Re_av%rhw(i,j,k)*ir

              ruu = mean_rhuu - meanR_u_u
              rvv = mean_rhvv - meanR_v_v
              rww = mean_rhww - meanR_w_w

              ekt = (ruu + rvv + rww)/u_inf**2

              wrk(i,j,k) = ekt

           enddo
        enddo

        return
end subroutine compute_ekt
















end module post_computation_module
