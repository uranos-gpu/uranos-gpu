module post_solutions_module
! ---------------------------------------------------------
!
!       This module contains some analytical solution of 
!       Navier-Stokes equations
!       
! ---------------------------------------------------------
use parameters_module
use storage_module
use post_storage_module
implicit none

private
public compute_exact_fields

contains
subroutine compute_exact_fields
        implicit none

        ! periodic euler exact solution
        if(    ic == 'periodic_euler_x') then
                call exact_periodic_euler_x(exact%rho,error%rho)
        elseif(ic == 'periodic_euler_y') then
                call exact_periodic_euler_y(exact%rho,error%rho)
        elseif(ic == 'periodic_euler_z') then
                call exact_periodic_euler_z(exact%rho,error%rho)

        elseif(ic == 'shock_tube_x') then
                call exact_shock_tube_x(exact%rho)

        elseif(ic == 'linear_ode') then
                call exact_linear_ode(exact%rho,error%rho)
        elseif(ic == 'linear_advection') then
                call exact_linear_advection(exact%rho,error%rho)

        ! isentropic vortex exact solution
        elseif(ic == 'isentropic_vortex_x') then
                call exact_isentropic_vortex_x(exact%rho,error%rho)

        ! poiseuille exact solutions
        elseif(ic == 'poiseuille_x' .or. ic == 'inflow_poiseuille') then
                call exact_steady_poiseuille('x',exact%vel,exact%tmp,error%vel,error%tmp)
        elseif(ic == 'poiseuille_y') then
                call exact_steady_poiseuille('y',exact%vel,exact%tmp,error%vel,error%tmp)
        elseif(ic == 'poiseuille_z') then
                call exact_steady_poiseuille('z',exact%vel,exact%tmp,error%vel,error%tmp)

        ! couette exact solution
        elseif(ic == 'couette_x') then
                call exact_startup_couette('x',exact%vel,error%vel)
        elseif(ic == 'couette_y') then
                call exact_startup_couette('y',exact%vel,error%vel)
        elseif(ic == 'couette_z') then
                call exact_startup_couette('z',exact%vel,error%vel)
        
        ! steady couette
        elseif(ic == 'steady_couette_x') then
                call exact_steady_couette('x',exact%vel,error%vel)
        elseif(ic == 'steady_couette_y') then
                call exact_steady_couette('y',exact%vel,error%vel)
        elseif(ic == 'steady_couette_z') then
                call exact_steady_couette('z',exact%vel,error%vel)

        ! I stokes problem exact solution
        elseif(ic == 'I_stokes_problem_x') then
                call exact_I_stokes_problem('x',exact%vel,error%vel)
        elseif(ic == 'I_stokes_problem_y') then
                call exact_I_stokes_problem('y',exact%vel,error%vel)
        elseif(ic == 'I_stokes_problem_z') then
                call exact_I_stokes_problem('z',exact%vel,error%vel)

        ! II stokes problem exact solution
        elseif(ic == 'II_stokes_problem_x') then
                call exact_II_stokes_problem('x',exact%vel,error%vel)
        elseif(ic == 'II_stokes_problem_y') then
                call exact_II_stokes_problem('y',exact%vel,error%vel)
        elseif(ic == 'II_stokes_problem_z') then
                call exact_II_stokes_problem('z',exact%vel,error%vel)

        elseif(ic == 'KolmogorovFlow') then
                call exact_KolmogorovFlow(exact%vel,error%vel)

        ! Blasius
        elseif(ic == 'blasius' .and. f == 1) then
                call exact_blasius
        endif

        return
end subroutine compute_exact_fields


subroutine exact_periodic_euler_x(exact,error)
! -------------------------------------------------
! Exact solution of periodic euler equations
! -------------------------------------------------
        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: error, exact
        real(rp)                                               :: x_, period

        period = 2.0_rp

        do k = sz,ez
           do j = sy,ey
              do i = sx,ex

                 x_ = mod(x(i) - Mach*time - 0.5_rp*period, period)
                 x_ = mod(x_ + period, period) - 0.5_rp*period

                 exact(i,j,k) = 1._rp + 0.2_rp*sin(pi*x_)

                 error(i,j,k) = phi(i,j,k,1) - exact(i,j,k)

              enddo
          enddo
        enddo
        
        return
end subroutine exact_periodic_euler_x

subroutine exact_periodic_euler_y(exact,error)
! -------------------------------------------------
! Exact solution of periodic euler equations
! -------------------------------------------------
        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: error, exact
        real(rp)                                               :: y_, period

        period = 2.0_rp

        do k = sz,ez
           do j = sy,ey
              do i = sx,ex

                 y_ = mod(y(j) - Mach*time - 0.5_rp*period, period)
                 y_ = mod(y_ + period, period) - 0.5_rp*period

                 exact(i,j,k) = 1._rp + 0.2_rp*sin(pi*y_)

                 error(i,j,k) = phi(i,j,k,1) - exact(i,j,k)

              enddo
          enddo
        enddo

        return
end subroutine exact_periodic_euler_y


subroutine exact_periodic_euler_z(exact,error)
! -------------------------------------------------
! Exact solution of periodic euler equations
! -------------------------------------------------
        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: error, exact
        real(rp)                                               :: z_, period

        period = 2.0_rp

        do k = sz,ez
           do j = sy,ey
              do i = sx,ex

                 z_ = mod(z(k) - Mach*time - 0.5_rp*period, period)
                 z_ = mod(z_ + period, period) - 0.5_rp*period

                 exact(i,j,k) = 1._rp + 0.2_rp*sin(pi*z_)

                 error(i,j,k) = phi(i,j,k,1) - exact(i,j,k)

              enddo
          enddo
        enddo

        return
end subroutine exact_periodic_euler_z


subroutine exact_linear_ode(exact,error)
! -----------------------------------------------
! Exact solution of y'+y=0
! -----------------------------------------------
        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: error, exact
        real(rp)                                               :: exp_t
        
        exp_t = exp(-time)

        do k = sz,ez
           do j = sy,ey
              do i = sx,ex

                exact(i,j,k) = exp_t
                error(i,j,k) = phi(i,j,k,1) - exact(i,j,k)

              enddo
          enddo
        enddo

        return
end subroutine exact_linear_ode


subroutine exact_linear_advection(exact,error)
! -----------------------------------------------
! Exact solution of u_t + u_x = 0
! -----------------------------------------------
        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: error, exact
        real(rp)                                               :: period, x_

        period = 2._rp

        do k = sz,ez
           do j = sy,ey
              do i = sx,ex

                 x_ = mod(x(i) - time - 0.5_rp*period, period)
                 x_ = mod(x_ + period, period) - 0.5_rp*period

                 exact(i,j,k) = 1._rp + 0.2_rp*sin(pi*x_)

                 error(i,j,k) = phi(i,j,k,1) - exact(i,j,k)

              enddo
          enddo
        enddo

        
        return
end subroutine exact_linear_advection


subroutine exact_isentropic_vortex_x(exact,error)
        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: error, exact
        real(rp)                                               :: x0, y0, beta, alpha, period
        real(rp)                                               :: x_, y_, r2, exp_

        x0   = 0._rp
        y0   = 0._rp
        beta = 5._rp
        alpha = (beta**2)*(gamma0-1._rp)/(8._rp*gamma0*pi**2)
        period = 10._rp
        
        do k = sz,ez
           do j = sy,ey

              y_ = y(j) - y0

              do i = sx,ex

                 x_ = mod(x(i) - x0 - Mach*time - 0.5_rp*period, period)
                 x_ = mod(x_ + period, period) - 0.5_rp*period

                 r2 = x_**2 + y_**2
                 exp_ = exp(1._rp-r2)
                 
                 exact(i,j,k) = (1._rp - alpha * exp_)**(1._rp/(gamma0-1._rp))

                 error(i,j,k) = phi(i,j,k,1) - exact(i,j,k)

              enddo
           enddo
        enddo

        return
end subroutine exact_isentropic_vortex_x


subroutine exact_I_stokes_problem(dir,exact,error)
        implicit none
        character(1)                           , intent(in)    :: dir
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: error, exact
        real(rp)                                               :: erf_
        
        selectcase(dir)

                case('x')
                do j = sy-1,ey+1

                         erf_ = erf((y(j))/(2*sqrt(mu_inf*time)))

                    do k = sz-1,ez+1
                      do i = sx-1,ex-1

                         exact(i,j,k) = 0.5_rp*u_inf * (1.0_rp - erf_)

                         error(i,j,k) = phi(i,j,k,2)/phi(i,j,k,1) - exact(i,j,k)

                      enddo
                   enddo
                enddo

                case('y')
                do i = sx-1,ex+1

                         erf_ = erf((x(i))/(2*sqrt(mu_inf*time)))

                   do k = sz-1,ez+1
                      do j = sy-1,ey+1

                         exact(i,j,k) = 0.5_rp*u_inf * (1.0_rp - erf_)

                         error(i,j,k) = phi(i,j,k,3)/phi(i,j,k,1) - exact(i,j,k)

                      enddo
                   enddo
                enddo

                case('z')
                do k = sz-1,ez+1

                         erf_ = erf((z(k))/(2*sqrt(mu_inf*time)))

                   do j = sy-1,ey+1
                      do i = sx-1,ex+1

                         exact(i,j,k) = 0.5_rp*u_inf * (1.0_rp - erf_)

                         error(i,j,k) = phi(i,j,k,2)/phi(i,j,k,1) - exact(i,j,k)

                      enddo
                   enddo
                enddo


        endselect

        return
end subroutine exact_I_stokes_problem


subroutine exact_II_stokes_problem(dir,exact,error)
        implicit none
        character(1)                           , intent(in)    :: dir
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: error, exact

        ! local declarations
        real(rp) :: omega, omega_c

        omega = 0.25_rp*pi
        omega_c = sqrt((omega/(2*mu_inf)))
        selectcase(dir)

                case('x')
                do k = sz-1,ez+1
                   do j = sy-1,ey+1
                      do i = sx-1,ex+1

                        exact(i,j,k) = u_inf*exp(-omega_c*y(j)) * cos(omega*time - omega_c*y(j))

                        error(i,j,k) = phi(i,j,k,2)/phi(i,j,k,1) - exact(i,j,k)

                      enddo
                   enddo
                enddo

                case('y')
                do k = sz-1,ez+1
                   do j = sy-1,ey+1
                      do i = sx-1,ex+1

                        exact(i,j,k) = u_inf*exp(-omega_c*x(i)) * cos(omega*time - omega_c*x(i))

                        error(i,j,k) = phi(i,j,k,3)/phi(i,j,k,1) - exact(i,j,k)

                      enddo
                   enddo
                enddo

                case('z')
                do k = sz-1,ez+1
                   do j = sy-1,ey+1
                      do i = sx-1,ex+1

                        exact(i,j,k) = u_inf*exp(-omega_c*z(k)) * cos(omega*time - omega_c*z(k))

                        error(i,j,k) = phi(i,j,k,2)/phi(i,j,k,1) - exact(i,j,k)

                      enddo
                   enddo
                enddo

        endselect
        return
end subroutine exact_II_stokes_problem



subroutine exact_steady_poiseuille(dir,exact_u,exact_T,error_u,error_t)
        implicit none
        character(1)                           , intent(in)    :: dir
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: error_u, exact_u
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: error_T, exact_T

        ! local declarations
        real(rp) :: dp_dx, dp_dy, c_1, c_2
        
        select case(dir)

                case('x')
                dp_dx = -2._rp*u_inf
                c_1   = 0.0_rp!-0.5_rp * dp_dx * (ymin+ymax)
                c_2   = - 0.5_rp*dp_dx

                do k = sz-1,ez+1
                   do j = sy-1,ey+1
                      do i = sx-1,ex+1

                         ! exact Poseuille speed
                         exact_u(i,j,k) = 0.5_rp * dp_dx * y(j)**2 + c_1 * y(j) + c_2

                         error_u(i,j,k) = phi(i,j,k,2)/phi(i,j,k,1) - exact_u(i,j,k)

                         ! Exact Poiseuille temperature
                         exact_T(i,j,k) = 1._rp - 0.5_rp* (gamma0-1._rp)/gamma0 * Prandtl * exact_u(i,j,k)**2

                         error_T(i,j,k) = T(i,j,k) - exact_T(i,j,k)

                      enddo
                   enddo
                enddo

                case('y')
                dp_dy = -2._rp*u_inf
                c_1 = -0.5_rp * dp_dy * (xmin+xmax)
                c_2 =  0.5_rp * dp_dy * xmin * (xmin+xmax) - 0.5_rp*dp_dy*xmin**2

                do k = sz-1,ez+1
                   do j = sy-1,ey+1
                      do i = sx-1,ex+1

                         ! exact Poseuille speed
                         exact_u(i,j,k) = 0.5_rp * dp_dy * x(i)**2 + c_1 * x(i) + c_2

                         error_u(i,j,k) = phi(i,j,k,3)/phi(i,j,k,1) - exact_u(i,j,k)

                         ! Exact Poiseuille temperature
                         exact_T(i,j,k) = 1._rp - 0.5_rp* (gamma0-1._rp)/gamma0 * Prandtl * exact_u(i,j,k)**2

                         error_T(i,j,k) = T(i,j,k) - exact_T(i,j,k)

                      enddo
                   enddo
                enddo

                case('z')
                dp_dx = -2._rp*u_inf
                c_1 = -0.5_rp * dp_dx * (zmin+zmax)
                c_2 =  0.5_rp * dp_dx * zmin * (zmin+zmax) - 0.5_rp*dp_dx*zmin**2

                do k = sz-1,ez+1
                   do j = sy-1,ey+1
                      do i = sx-1,ex+1

                         ! exact Poseuille speed
                         exact_u(i,j,k) = 0.5_rp * dp_dx * z(k)**2 + c_1 * z(k) + c_2

                         error_u(i,j,k) = phi(i,j,k,2)/phi(i,j,k,1) - exact_u(i,j,k)

                         ! Exact Poiseuille temperature
                         exact_T(i,j,k) = 1._rp - 0.5_rp* (gamma0-1._rp)/gamma0 * Prandtl * exact_u(i,j,k)**2

                         error_T(i,j,k) = T(i,j,k) - exact_T(i,j,k)

                      enddo
                   enddo
                enddo


        end select

        return
end subroutine exact_steady_poiseuille


subroutine exact_KolmogorovFlow(exact_u,error_u)
        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: error_u, exact_u

        do k = sz-1,ez+1
           do j = sy-1,ey+1
              do i = sx-1,ex+1

                 exact_u(i,j,k) = 0.1_rp*u_inf/mu_inf*(Ly/(8*pi))**2*sin(8*pi*y(j)/Ly)

                 error_u(i,j,k) = phi(i,j,k,2)/phi(i,j,k,1) - exact_u(i,j,k)

              enddo
           enddo
        enddo

        return
end subroutine exact_KolmogorovFlow



subroutine exact_startup_couette(dir,exact,error)
        implicit none
        character(1)                           , intent(in)    :: dir
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: error, exact

        ! local declaration
        real(rp) :: x_, y_, z_                
        real(rp) :: h               !< channel span
        real(rp) :: w               !< viscous frequency
        real(rp) :: ut              !< unisteady part of the solution
        integer  :: n,ntot = 1000   !< integer to compute the unsteady part of the solution

        select case(dir)

                case('x')
                h = ymax - ymin
                w = pi**2 * mu_inf * time/ h**2

                do k = sz-1,ez+1
                   do j = sy-1,ey+1

                      ! unsteady part of the solution
                      y_ = y(j)
                      ut = 0.0_rp
                      do n = 1, ntot

                         ut = ut + 1.0_rp/real(n,rp) * exp(-n**2 * w) * sin(n*pi*(1.0_rp - y_/h))

                      enddo
                      ut = ut* (2*u_inf/pi)

                      do i = sx-1,ex+1
                      
                         exact(i,j,k) = u_inf * y_/h - ut

                         error(i,j,k) = phi(i,j,k,2)/phi(i,j,k,1) - exact(i,j,k)

                      enddo
                   enddo
                enddo

                case('y')
                
                h = xmax - xmin
                w = pi**2 * mu_inf * time/ h**2

                do k = sz-1,ez+1
                   do i = sx-1,ex+1

                      ! unsteady part of the solution
                      x_ = x(i)
                      ut = 0.0_rp
                      do n = 1, ntot

                         ut = ut + 1.0_rp/real(n,rp) * exp(-n**2 * w) * sin(n*pi*(1.0_rp - x_/h))

                      enddo
                      ut = ut* (2*u_inf/pi)

                      do j = sy-1,ey+1
                      
                         exact(i,j,k) = u_inf * x_/h - ut

                         error(i,j,k) = phi(i,j,k,3)/phi(i,j,k,1) - exact(i,j,k)

                      enddo
                   enddo
                enddo

                case('z')

                h = zmax - zmin
                w = pi**2 * mu_inf * time/ h**2

                do k = sz-1,ez+1

                   ! unsteady part of the solution
                   z_ = z(k)
                   ut = 0.0_rp
                   do n = 1, ntot

                      ut = ut + 1.0_rp/real(n,rp) * exp(-n**2 * w) * sin(n*pi*(1.0_rp - z_/h))

                   enddo
                   ut = ut* (2*u_inf/pi)

                   do i = sx-1,ex+1
                      do j = sy-1,ey+1

                         exact(i,j,k) = u_inf * z_/h - ut

                         error(i,j,k) = phi(i,j,k,2)/phi(i,j,k,1) - exact(i,j,k)

                      enddo
                   enddo
                enddo

        endselect

        return
end subroutine exact_startup_couette


subroutine exact_steady_couette(dir,exact,error)
        implicit none
        character(1)                           , intent(in)    :: dir
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: error, exact

        ! local declaration
        real(rp) :: h               !< channel span

        select case(dir)

                case('x')
                h = ymax - ymin

                do k = sz-1,ez+1
                   do j = sy-1,ey+1
                      do i = sx-1, ex+1

                         exact(i,j,k) = u_inf * y(j)/h

                         error(i,j,k) = phi(i,j,k,2)/phi(i,j,k,1) - exact(i,j,k)

                      enddo
                   enddo
                enddo

                case('y')
                h = ymax - ymin

                do k = sz-1,ez+1
                   do j = sy-1,ey+1
                      do i = sx-1, ex+1

                         exact(i,j,k) = u_inf * x(i)/h

                         error(i,j,k) = phi(i,j,k,3)/phi(i,j,k,1) - exact(i,j,k)

                      enddo
                   enddo
                enddo

                case('z')
                h = zmax - zmin

                do k = sz-1,ez+1
                   do j = sy-1,ey+1
                      do i = sx-1, ex+1

                         exact(i,j,k) = u_inf * z(k)/h

                         error(i,j,k) = phi(i,j,k,2)/phi(i,j,k,1) - exact(i,j,k)

                      enddo
                   enddo
                enddo

        end select
        return
end subroutine exact_steady_couette


subroutine exact_shock_tube_x(R)

        use math_tools_module, only: newton_raphson

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: R

        ! local declaration
        integer , parameter :: NR_itmax = 1000

        real(rp), parameter :: G1 = 2*gamma0/(gamma0+1.0_rp)
        real(rp), parameter :: G2 = (gamma0-1.0_rp)/(gamma0+1.0_rp)
        real(rp), parameter :: G4 = 2.0_rp/(gamma0+1.0_rp)
        real(rp), parameter :: G5 = 2*gamma0/(gamma0-1.0_rp)

        real(rp), parameter :: r_L = 1.0_rp, r_R = 0.125_rp
        real(rp), parameter :: p_L = 1.0_rp, p_R = 0.1_rp
        real(rp)            :: T_L      , T_R
        real(rp)            :: a_L      , a_R

        real(rp) :: Ms
        real(rp) :: phi_p, phi_r, phi_T
        real(rp) :: p1, r1, u1
        real(rp) :: p2, r2, u2, a2

        real(rp) :: ui, pi, ai

        real(rp), parameter :: x0 = 0.0_rp   !< shock initial position
        real(rp)            :: xi         !< grid position
        real(rp)            :: x1         !< left fan position
        real(rp)            :: x2         !< right fan positon
        real(rp)            :: x3         !< contact discontinuity position
        real(rp)            :: x4         !< shock position
        
        T_r = p_R/r_R
        T_L = p_L/r_L
        a_R = sqrt(gamma0*T_R)
        a_L = sqrt(gamma0*T_L)

        ! compute pressure ratio witheen the shock
        call newton_raphson(Mach_number,NR_itmax, toll_equality, 1.7_rp, Ms)

        ! ====== RH quotients
        phi_p = G1*Ms**2 - G2
        phi_r = (2.0_rp/(gamma0+1)*(1.0_rp/Ms**2) + G2)**(-1)
        phi_T = phi_p / phi_r

        ! ====== STATE 1
        p1 = p_R * phi_p
        r1 = r_R * phi_r
        u1 = (2.0_rp/(gamma0+1._rp)) *(Ms- 1._rp/Ms)

        ! ====== STATE 2
        p2 = p1
        u2 = u1
        r2 = r_L*(p2/p_L)**(1.0_rp/gamma0)
        a2 = sqrt(gamma0*p2/r2)

        ! ====== shock positions
        x1 = x0 - a_l*time
        x2 = x0 + (u2-a2)*time
        x3 = x0 + u2*time
        x4 = x0 + Ms*time

        do k = sz,ez
           do j = sy,ey
              do i = sx,ex
                
                 xi = x(i)
                 if    (xi >= xmin .and. xi < x1) then
                   R(i,j,k) = r_L

                 elseif(xi >= x1   .and. xi < x2) then
                   ui       = G4 * (a_l + (xi-x0)/time)
                   ai       = a_l - 0.5_rp*(gamma0-1.0_rp) * ui
                   pi       = p_l*(ai/a_l)**G5
                   R(i,j,k) = (gamma0*pi)/(ai**2)

                 elseif(xi >= x2   .and. xi < x3) then
                   R(i,j,k) = r2

                 elseif(xi >= x3   .and. xi < x4) then
                   R(i,j,k) = r1

                 elseif(xi >= x4   .and. xi <= xmax) then
                   R(i,j,k) = r_R

                 endif

              enddo
           enddo
        enddo


        return
end subroutine exact_shock_tube_x


function Mach_number(Ms) result(fM)
        implicit none
        real(rp), intent(in) :: Ms
        real(rp)             :: fM

        ! local declarations
        real(rp), parameter :: r_L = 1.0_rp
        real(rp), parameter :: p_L = 1.0_rp, p_R = 0.1_rp
        real(rp)            :: a_L
        real(rp)            :: aa,bb,cc,dd, inner

        ! right/left state
        a_L = sqrt(gamma0*p_L/r_L)

        aa = (gamma0+1)/(gamma0-1)
        bb = (gamma0-1)/(gamma0+1)
        cc = (2*gamma0)/(gamma0+1)
        dd = (gamma0-1)/(2*gamma0)
       
        inner = 1.0_rp - (p_R/p_L * (cc*Ms**2 - bb)) **(dd)

        fM = Ms -1.0_rp/Ms - a_l * aa * inner

        return
end function Mach_number



subroutine exact_blasius

        use input_module

        implicit none
        real(rp), dimension(:), allocatable :: eta
        real(rp), dimension(:), allocatable :: f0, f1, f2
        real(rp)                            :: y_, d_eta, eta_, lb, ub, errb, toll
        real(rp), parameter                 :: x0 = 1.0_rp
        integer , parameter                 :: n_huge = 10000, itermax = 100
        integer                             :: j, err = 0, ntot = 100, funit, iter
        character(dl)                       :: filename

        !
        ! ==== find were eta > 10
        !
        loop1:do j = 0,n_huge

           y_   = (j-0.5_rp)/real(n_huge,rp)*Ly
           eta_ = sqrt(u_inf/(mu_inf*x0))*y_

           if(eta_ > 20.0_rp) then
             ntot = j
             exit loop1
           endif

        enddo loop1

        allocate(eta(0:ntot), f0(1:ntot), f1(1:ntot), f2(1:ntot), stat = err)
        if(err /= 0) stop ' Error allocating blasius solution'

        !
        ! ==== compute blasius coordinate
        !
        do j = 0,ntot
                
           y_     = (j-0.5_rp)/real(n_huge,rp)*Ly
           eta(j) = sqrt(u_inf/(mu_inf*x0))*y_

        enddo

        !
        ! ==== init the solution
        !
        f0(1) = 0.0_rp
        f1(1) = 0.0_rp

        lb = 0.32_rp
        ub = 0.33_rp
        f2(1) = 0.5_rp*(lb+ub)

        toll = 1.0E-14_rp
        errb = 2*toll
        iter = 1
        do while (errb > toll .and. iter < itermax)
        
           do j = 1,ntot-1

              d_eta   =  eta(j) - eta(j-1)

              f0(j+1) = f0(j) + f1(j)         * d_eta
              f1(j+1) = f1(j) + f2(j)         * d_eta
              f2(j+1) = f2(j) - 0.5_rp*f0(j)*f2(j)*d_eta


           enddo

           if(f1(ntot) > 1.0_rp) then
             ub = 0.5_rp*(ub+lb)
             f2(1) = 0.5_rp*(lb+ub)
           elseif(f1(ntot) < 1.0_rp) then
             lb = 0.5_rp*(ub+lb)
             f2(1) = 0.5_rp*(lb+ub)
           endif
                
           errb = abs(f1(ntot) - 1.0_rp)
           iter = iter + 1

        enddo
        
        !
        ! ==== write the solution on res file
        !
        call get_unit(funit)
        filename = "DATA/"//trim(data_dir)//"/GNUPLOT/"//'blasius_exact.txt'
        open(unit = funit, file = filename, status = 'replace', access = 'sequential')

        do j = 1,ntot-1
           write(funit,*) eta(j), f1(j)
        enddo

        close(funit)

        deallocate(eta,f0,f1,f2)
        return
end subroutine exact_blasius













end module post_solutions_module
