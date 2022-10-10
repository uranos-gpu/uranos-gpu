module ic_module
use parameters_module
use storage_module
use mpi_module
use fluid_functions_module
implicit none

contains


subroutine ambient(status)
! -----------------------------------------------------------------------------------
!       
!       This subroutine zero field boundary condition
!
! -----------------------------------------------------------------------------------
        implicit none
        integer, intent(out) :: status
        
        if(rank == root) write(*,'(A)', advance = 'no') ' Initializing ambient initial condition: '
        do k = lbz, ubz
           do j = lby,uby
              do i = lbx,ubx

                 phi(i,j,k,1) = 1._rp
                   P(i,j,k)   = 1._rp 
                   U(i,j,k)   = 0._rp
                   V(i,j,k)   = 0._rp
                   W(i,j,k)   = 0._rp

              enddo
           enddo
        enddo
        if(rank == root) write(*, '(A)') 'done!'
        status = 0
return
end subroutine ambient





subroutine init_periodic_euler(dir,status)
        implicit none
        character(len=1), intent(in)  :: dir
        integer, intent(out) :: status
        
        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing periodic euler along ', dir, ': '

        select case(dir)
        case('x')
        do k = lbz, ubz      
           do j = lby,uby
              do i=lbx,ubx
                  phi(i,j,k,1) = 1._rp + 0.2_rp*sin(pi*x(i))
                  u  (i,j,k)   = Mach
                  v  (i,j,k)   = 0._rp
                  w  (i,j,k)   = 0._rp
                  p  (i,j,k)   = 1._rp
              enddo
           enddo
        enddo

        case('y')
        do k = lbz, ubz
           do j = lby,uby
              do i=lbx,ubx
                  phi(i,j,k,1) = 1._rp + 0.2_rp*sin(pi*y(j))
                  u  (i,j,k)   = 0._rp
                  v  (i,j,k)   = Mach
                  w  (i,j,k)   = 0._rp
                  p  (i,j,k)   = 1._rp
              enddo
          enddo
        enddo

        case('z')
        do k = lbz, ubz
           do j = lby,uby
              do i=lbx,ubx
                  phi(i,j,k,1) = 1._rp + 0.2_rp*sin(pi*z(k))
                  u  (i,j,k)   = 0._rp
                  v  (i,j,k)   = 0._rp
                  w  (i,j,k)   = Mach
                  p  (i,j,k)   = 1._rp
              enddo
          enddo
        enddo

        endselect
        if(rank == root) write(*, '(A)') 'done!'
        status = 0
return
end subroutine init_periodic_euler




subroutine init_shock_tube(dir, status)
        implicit none
        character(len=1), intent(in)  :: dir
        integer         , intent(out) :: status
        
        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing shock tube along ', dir, ': '

        select case(dir)
        case('x')
        do k = lbz,ubz
           do j = lby,uby        
              do i = lbx, ubx

                 if(x(i).le.0._rp)then
                     phi(i,j,k,1) = 1._rp
                       U(i,j,k)   = 0._rp
                       V(i,j,k)   = 0._rp
                       W(i,j,k)   = 0._rp
                       P(i,j,k)   = 1._rp
                 elseif(x(i).gt.0._rp) then
                     phi(i,j,k,1) = 0.125_rp
                       U(i,j,k)   = 0._rp
                       V(i,j,k)   = 0._rp
                       W(i,j,k)   = 0._rp
                       P(i,j,k)   = 0.1_rp
                 end if

              enddo
           enddo
        enddo

        case('y')
        do k=lbz,ubz
           do j=lby, uby
              do i = lbx, ubx

                 if(y(j).le.0._rp)then
                     phi(i,j,k,1) = 1._rp
                       U(i,j,k)   = 0._rp
                       V(i,j,k)   = 0._rp
                       W(i,j,k)   = 0._rp
                       P(i,j,k)   = 1._rp

                 elseif(y(j).gt.0._rp) then
                     phi(i,j,k,1) = 0.125_rp
                       U(i,j,k)   = 0._rp
                       V(i,j,k)   = 0._rp
                       W(i,j,k)   = 0._rp
                       P(i,j,k)   = 0.1_rp
                 end if

              enddo
           enddo
        enddo

        case('z')
        do k=lbz,ubz
           do j=lby, uby
              do i = lbx, ubx

                 if(z(k).le.0._rp)then
                     phi(i,j,k,1) = 1._rp
                       U(i,j,k)   = 0._rp
                       V(i,j,k)   = 0._rp
                       W(i,j,k)   = 0._rp
                       P(i,j,k)   = 1._rp

                 elseif(z(k).gt.0._rp) then
                     phi(i,j,k,1) = 0.125_rp
                       U(i,j,k)   = 0._rp
                       V(i,j,k)   = 0._rp
                       W(i,j,k)   = 0._rp
                       P(i,j,k)   = 0.1_rp
                 end if

              enddo
           enddo
        enddo

        end select
        if(rank == root) write(*, '(A)') 'done!'
        status = 0
return
end subroutine init_shock_tube


subroutine init_lax_problem(dir,status)
        implicit none
        character(len=1), intent(in)  :: dir
        integer         , intent(out) :: status
        
        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing Lax problem along ', dir,': '
        select case(dir)

        case('x')
        do k = lbz,ubz
           do j = lby,uby        
              do i = lbx, ubx
                      
                 if(x(i).le.0.0_rp)then
                     phi(i,j,k,1) = 0.445_rp
                       U(i,j,k)   = 0.698_rp
                       V(i,j,k)   = 0.0_rp
                       W(i,j,k)   = 0.0_rp
                       P(i,j,k)   = 3.528_rp
                 elseif(x(i).gt.0.0_rp) then
                     phi(i,j,k,1) = 0.5_rp
                       U(i,j,k)   = 0.0_rp
                       V(i,j,k)   = 0.0_rp
                       W(i,j,k)   = 0.0_rp
                       P(i,j,k)   = 0.571_rp
                 end if

              enddo
           enddo
        enddo

        case('y')
        do k = lbz,ubz
           do j=lby, uby
              do i = lbx, ubx

                 if(y(j).le.0.0_rp)then
                     phi(i,j,k,1) = 0.445_rp
                       U(i,j,k)   = 0.0_rp
                       V(i,j,k)   = 0.698_rp
                       W(i,j,k)   = 0.0_rp
                       P(i,j,k)   = 3.528_rp
                 elseif(y(j).gt.0.0_rp) then
                     phi(i,j,k,1) = 0.5_rp
                       U(i,j,k)   = 0.0_rp
                       V(i,j,k)   = 0.0_rp
                       W(i,j,k)   = 0.0_rp
                       P(i,j,k)   = 0.571_rp
                 end if
                      
              enddo
           enddo
        enddo

        case('z')
        do k = lbz,ubz
           do j=lby, uby
              do i = lbx, ubx

                 if(z(k).le.0.0_rp)then
                     phi(i,j,k,1) = 0.445_rp
                       U(i,j,k)   = 0.0_rp
                       V(i,j,k)   = 0.0_rp
                       W(i,j,k)   = 0.698_rp
                       P(i,j,k)   = 3.528_rp
                 elseif(z(k).gt.0.0_rp) then
                     phi(i,j,k,1) = 0.5_rp
                       U(i,j,k)   = 0.0_rp
                       V(i,j,k)   = 0.0_rp
                       W(i,j,k)   = 0.0_rp
                       P(i,j,k)   = 0.571_rp
                 end if
                      
              enddo
           enddo
        enddo

        endselect
        if(rank == root) write(*, '(A)') 'done!'
        status = 0
return
end subroutine init_lax_problem




subroutine init_linear_ode(status)
! -----------------------------------------------------------------------------------
!       
!       This subroutine initialize a zero state condition for linear ODE
!       problem y' + y = 0._rp Just to test RK scheme.
!
! -----------------------------------------------------------------------------------
        implicit none
        integer, intent(out) :: status
        
        if(rank == root) write(*,'(A)', advance = 'no') ' Initializing time linear initial condition: '
        do k = lbz, ubz
           do j = lby,uby
              do i = lbx,ubx

                 phi(i,j,k,1) = 1._rp
                   P(i,j,k)   = 1._rp 
                   U(i,j,k)   = 0._rp
                   V(i,j,k)   = 0._rp
                   W(i,j,k)   = 0._rp

              enddo
           enddo
        enddo
        if(rank == root) write(*, '(A)') 'done!'
        status = 0
return
end subroutine init_linear_ode


subroutine init_linear_advection(status)
! -----------------------------------------------------------------------------------
!       
!       This subroutine initialize a zero state condition for linear ODE
!       problem y' + y = 0._rp Just to test RK scheme.
!
! -----------------------------------------------------------------------------------
        implicit none
        integer, intent(out) :: status
        
        if(rank == root) write(*,'(A)', advance = 'no') ' Initializing time linear initial condition: '
        do k = lbz, ubz
           do j = lby,uby
              do i = lbx,ubx

                 phi(i,j,k,1) = 1._rp + 0.2_rp*sin(pi*x(i))
                   P(i,j,k)   = 1._rp 
                   U(i,j,k)   = 0._rp
                   V(i,j,k)   = 0._rp
                   W(i,j,k)   = 0._rp

              enddo
           enddo
        enddo
        if(rank == root) write(*, '(A)') 'done!'
        status = 0
return
end subroutine init_linear_advection





subroutine init_shock_wave_interaction(dir,status)
! -----------------------------------------------------------------------------------
!       Initialization of a shock-entropy wave interaction in 2D:
!       ref. "Conservative Hybrid Compact-WENO schemes for shock-turbolence
!             interaction" pag. 101
! -----------------------------------------------------------------------------------
        implicit none
        character(len=1), intent(in)  :: dir
        integer         , intent(out) :: status

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing shock-entropy wawe along ', dir, ': '
        select case(dir)
        case('x')
        do k = lbz,ubz
           do j = lby,uby
              do i=lbx,ubx

                 if(x(i).le.-4.0_rp)then
                         phi(i,j,k,1) =  3.857143_rp
                           U(i,j,k)   =  2.629369_rp
                           V(i,j,k)   =  0.0_rp
                           W(i,j,k)   =  0.0_rp
                           P(i,j,k)   = 10.333333_rp
                 elseif(x(i).gt.-4.0_rp) then
                         phi(i,j,k,1) = 1.0_rp+0.2_rp*sin(5.0_rp*(x(i)))
                           U(i,j,k)   = 0.0_rp
                           V(i,j,k)   = 0.0_rp
                           W(i,j,k)   = 0.0_rp
                           P(i,j,k)   = 1.0_rp
                 end if

              enddo
           end do
        enddo
        case('y')
        do k = lbz,ubz
           do j = lby,uby
              do i=lbx,ubx

                 if(y(j).le.-4.0_rp)then
           
                         phi(i,j,k,1) =  3.857143_rp
                           U(i,j,k)   =  0.0_rp
                           V(i,j,k)   =  2.629369_rp
                           W(i,j,k)   =  0.0_rp
                           P(i,j,k)   = 10.333333_rp
           
                 elseif(y(j).gt.-4.0_rp) then
           
                         phi(i,j,k,1) = 1.0_rp+0.2_rp*sin(5.0_rp*(y(j)))
                           U(i,j,k)   = 0.0_rp
                           V(i,j,k)   = 0.0_rp
                           W(i,j,k)   = 0.0_rp
                           P(i,j,k)   = 1.0_rp
           
                 end if

              enddo
           end do
        enddo

        case('z')
        do k = lbz,ubz
           do j = lby,uby
              do i=lbx,ubx

                 if(z(k).le.-4.0_rp)then
           
                         phi(i,j,k,1) =  3.857143_rp
                           U(i,j,k)   =  0.0_rp
                           V(i,j,k)   =  0.0_rp
                           W(i,j,k)   =  2.629369_rp
                           P(i,j,k)   = 10.333333_rp
           
                 elseif(z(k).gt.-4.0_rp) then
           
                         phi(i,j,k,1) = 1.0_rp+0.2_rp*sin(5.0_rp*(z(k)))
                           U(i,j,k)   = 0.0_rp
                           V(i,j,k)   = 0.0_rp
                           W(i,j,k)   = 0.0_rp
                           P(i,j,k)   = 1.0_rp
           
                 end if

              enddo
           end do
        enddo

        end select
        if(rank == root) write(*, '(A)') 'done!'
        status = 0
return
end subroutine init_shock_wave_interaction



subroutine init_isentropic_vortex(dir,status)
! -------------------------------------------------------------------
!       
!       This subroutine initializes an isentropic vortex.
!       Reference: 2005, Kim, 
!                  "A high-order accurate hybrid scheme using a central 
!                   flux scheme and a WENO scheme for compressible
!                   flowfield analysis"
!
! -------------------------------------------------------------------
        implicit none
        character(len=1), intent(in ) :: dir
        integer         , intent(out) :: status
        real(rp)                      :: x0, y0, alpha, beta, r2
        real(rp)                      :: exp_, x_, y_

        ! vortex geometry
        x0   = 0._rp
        y0   = 0._rp
        beta = 5._rp
        alpha = (beta**2)*(gamma0-1._rp)/(8._rp*gamma0*pi**2)

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing an isentropic vortex: '
        selectcase(dir)

        case('x')
        do k = lbz,ubz
           do j = lby,uby
              do i = lbx, ubx
                      
                 x_ = x(i) - x0
                 y_ = y(j) - y0

                 r2 = x_**2 + y_**2
                 exp_ = exp(1._rp-r2)
                 
                 phi(i,j,k,1) = (1._rp - alpha * exp_)**(1._rp    /(gamma0-1._rp))
                 P  (i,j,k)   = (1._rp - alpha * exp_)**(gamma0/(gamma0-1._rp))

                 U(i,j,k) = u_inf - beta/(2._rp*pi) * (y_ * exp_**0.5_rp)
                 V(i,j,k) =         beta/(2._rp*pi) * (x_ * exp_**0.5_rp)
                 W(i,j,k) = 0._rp

              enddo
           enddo
        enddo

        case('y')
        do k = lbz,ubz
           do j = lby,uby
              do i = lbx, ubx

                 x_ = x(i) - x0
                 y_ = y(j) - y0

                 r2 = x_**2 + y_**2
                 exp_ = exp(1._rp-r2)
                 
                 phi(i,j,k,1) = (1._rp - alpha * exp_)**(1._rp    /(gamma0-1._rp))
                 P  (i,j,k)   = (1._rp - alpha * exp_)**(gamma0/(gamma0-1._rp))

                 U(i,j,k) =         beta/(2._rp*pi) * (y_ * exp_**0.5_rp)
                 V(i,j,k) = u_inf - beta/(2._rp*pi) * (x_ * exp_**0.5_rp)
                 W(i,j,k) = 0._rp
                      
              enddo
          enddo
        enddo

        endselect
        if(rank == root) write(*, '(A)') 'done!'
        status = 0
return
end subroutine init_isentropic_vortex


subroutine init_pirozzoli_vortex(status)
        implicit none
        integer, intent(out) :: status

        ! local declarations
        real(rp)            :: x_, y_, r
        real(rp), parameter :: Mv = 0.4_rp

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing Pirozzoli vortex: '

        do k = lbz,ubz
           do j = lby,uby
              do i = lbx,ubx

                 x_ = x(i)
                 y_ = y(j)

                 r = sqrt(x_**2 + y_**2)

                 U(i,j,k) = u_inf*(1.0_rp - Mv/Mach * y_ * exp(0.5_rp*(1.0_rp-r**2)))
                 V(i,j,k) = u_inf*(      Mv/Mach * x_ * exp(0.5_rp*(1.0_rp-r**2)))
                 W(i,j,k) = 0.0_rp

                 PHI(i,j,k,1) = (1.0_rp - 0.5_rp*(gamma0-1.0_rp) * Mv**2 * exp(1.0_rp-r**2))**(1.0_rp   /(gamma0-1.0_rp))
                 P  (i,j,k)   = (1.0_rp - 0.5_rp*(gamma0-1.0_rp) * Mv**2 * exp(1.0_rp-r**2))**(gamma0/(gamma0-1.0_rp))

              enddo
           enddo
        enddo


        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return
end subroutine init_pirozzoli_vortex



subroutine init_pirozzoli_vortex_y(status)
        implicit none
        integer, intent(out) :: status

        ! local declarations
        real(rp)            :: x_, y_, r
        real(rp), parameter :: Mv = 0.4_rp

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing Pirozzoli vortex: '

        do k = lbz,ubz
           do j = lby,uby
              do i = lbx,ubx

                 x_ = x(i)
                 y_ = y(j)

                 r = sqrt(x_**2 + y_**2)

                 U(i,j,k) = u_inf*(    - Mv/Mach * y_ * exp(0.5_rp*(1.0_rp-r**2)))
                 V(i,j,k) = u_inf*(1.0_rp + Mv/Mach * x_ * exp(0.5_rp*(1.0_rp-r**2)))
                 W(i,j,k) = 0.0_rp

                 PHI(i,j,k,1) = (1.0_rp - 0.5_rp*(gamma0-1.0_rp) * Mv**2 * exp(1.0_rp-r**2))**(1.0_rp   /(gamma0-1.0_rp))
                 P  (i,j,k)   = (1.0_rp - 0.5_rp*(gamma0-1.0_rp) * Mv**2 * exp(1.0_rp-r**2))**(gamma0/(gamma0-1.0_rp))

              enddo
           enddo
        enddo


        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return
end subroutine init_pirozzoli_vortex_y






subroutine init_shock_wave(dir,status)
        implicit none
! --------------------------------------------------------------------------------
!
!       This subroutine initialize a non-stationary shock wave in a domain posizion. 
!
! --------------------------------------------------------------------------------
        character(len=1), intent(in)  :: dir
        integer         , intent(out) :: status

        ! local variables
        real(rp) :: r1, p1, T1, M1 
        real(rp) :: r2, p2, T2, M2 
        real(rp) :: vel1, vel2
        real(rp) :: x0              !< shock wave initial position

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing a shock wave : '

        M1 = Mach
        x0 = -2.0_rp
        ! --------------|-------------- !
        !     (2)       |     (1)       !   <---- M1 = MACH
        ! --------------|-------------- !
        !               ^
        !               x0                

        !
        ! ==== state 1
        !
        p1   = 1.0_rp
        r1   = 1.0_rp
        T1   = p1/r1
        vel1 = - M1 * sqrt(gamma0*T1)
        !
        ! ==== state 2
        !
        call Rankine_Hugoniot(r1,p1,T1,M1,r2,p2,T2,M2)

        vel2 = - M2 * sqrt(gamma0*T2)
       
        if(dir == 'x') then

          
          do k = lbz,ubz
             do j = lby,uby
                do i = lbx,ubx

                   if(x(i) < x0 ) then 

                           phi(i,j,k,1) = r2
                             U(i,j,k)   = vel2 - vel1
                             V(i,j,k)   = 0.0_rp
                             W(i,j,k)   = 0.0_rp
                             P(i,j,k)   = p2

                   elseif(x(i) >= x0 ) then

                           phi(i,j,k,1) = r1
                             U(i,j,k)   = vel1 - vel1
                             V(i,j,k)   = 0.0_rp
                             W(i,j,k)   = 0.0_rp
                             P(i,j,k)   = p1
                   endif

                enddo
             enddo
          enddo

        elseif(dir == 'y') then

          do k = lbz,ubz
             do j = lby,uby
                do i = lbx,ubx

                   if(y(j) < x0 ) then 

                           phi(i,j,k,1) = r2
                             U(i,j,k)   = 0.0_rp
                             V(i,j,k)   = vel2 - vel1
                             W(i,j,k)   = 0.0_rp
                             P(i,j,k)   = p2

                   elseif(y(j) >= x0 ) then

                           phi(i,j,k,1) = r1
                             U(i,j,k)   = 0.0_rp
                             V(i,j,k)   = vel1 - vel1
                             W(i,j,k)   = 0.0_rp
                             P(i,j,k)   = p1
                   endif

                enddo
             enddo
          enddo

        elseif(dir == 'z') then

          do k = lbz,ubz
             do j = lby,uby
                do i = lbx,ubx

                   if(z(k) .lt. x0 ) then 

                           phi(i,j,k,1) = r2
                             U(i,j,k)   = 0.0_rp
                             V(i,j,k)   = 0.0_rp
                             W(i,j,k)   = vel2-vel1
                             P(i,j,k)   = p2

                   elseif(z(k) .ge. x0 ) then

                           phi(i,j,k,1) = r1
                             U(i,j,k)   = 0.0_rp
                             V(i,j,k)   = 0.0_rp
                             W(i,j,k)   = vel1-vel1
                             P(i,j,k)   = p1
                   endif

                enddo
             enddo
          enddo

        endif

        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return
endsubroutine init_shock_wave




subroutine init_double_mach_reflection(status)
! --------------------------------------------------------------------------------
!
!       This subroutine initialize the double mach reflection problem for 
!       Navier-Stokes equations
!
! --------------------------------------------------------------------------------
        implicit none
        integer, intent(out) :: status

        ! local variables
        real(rp), parameter :: x0 = 1.0_rp/6.0_rp, theta = pi/6.0_rp
        real(rp)            :: r1, p1, T1, M1
        real(rp)            :: r2, p2, T2, u2, v2, M2
        real(rp)            :: vel1, vel2, tan_theta, xi, yj


        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing double mach reflection: '

        tan_theta = tan(theta)

        !
        ! ==== state 1
        !
        M1   = Mach
        p1   = 1.0_rp
        r1   = 1.4_rp
        T1   = p1/r1
        vel1 = - M1 * sqrt(gamma0*T1)
        !
        ! ==== state 2
        !
        call Rankine_Hugoniot(r1,p1,T1,M1,r2,p2,T2,M2)

        vel2 = - M2 * sqrt(gamma0*T2)
        u2   =   (vel2 - vel1)*cos(theta)
        v2   = - (vel2 - vel1)*sin(theta)
        
        do k       = lbz,ubz
           do j    = lby,uby
              do i = lbx,ubx

                 xi = x(i)
                 yj = y(j)

                 if    (xi <  x0 + yj*tan_theta) then 

                         phi(i,j,k,1) = r2
                           U(i,j,k)   = u2
                           V(i,j,k)   = v2
                           W(i,j,k)   = 0.0_rp
                           P(i,j,k)   = p2

                 elseif(xi >= x0 + yj*tan_theta) then

                         phi(i,j,k,1) = r1
                           U(i,j,k)   = 0.0_rp
                           V(i,j,k)   = 0.0_rp
                           W(i,j,k)   = 0.0_rp
                           P(i,j,k)   = p1
                 endif

              enddo
           enddo
        enddo

        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return
end subroutine init_double_mach_reflection


subroutine init_I_stokes_problem(dir,status)
        implicit none
        character(len=1), intent(in)  :: dir
        integer         , intent(out) :: status

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing I Stokes problem field: '
        
        selectcase(dir)
                
          case('x')
          do k = lbz,ubz
             do j = lby,uby
                do i = lbx,ubx

                   phi(i,j,k,1) = 1.0_rp
                     V(i,j,k)   = 0.0_rp
                     W(i,j,k)   = 0.0_rp
                     P(i,j,k)   = 1.0_rp

                   if(y(j) <= 0.0_rp) then
                     U(i,j,k)   = u_inf
                   else
                     U(i,j,k)   = 0.0_rp
                   endif

                enddo
             enddo
          enddo

        case('y')
          do k = lbz,ubz
             do j = lby,uby
                do i = lbx,ubx

                   phi(i,j,k,1) = 1.0_rp
                     U(i,j,k)   = 0.0_rp
                     W(i,j,k)   = 0.0_rp
                     P(i,j,k)   = 1.0_rp

                   if(x(i) <= 0.0_rp) then
                     V(i,j,k)   = u_inf
                   else
                     V(i,j,k)   = 0.0_rp
                   endif

                enddo
             enddo
          enddo

        case('z')
          do k = lbz,ubz
             do j = lby,uby
                do i = lbx,ubx

                   phi(i,j,k,1) = 1.0_rp
                     V(i,j,k)   = 0.0_rp
                     W(i,j,k)   = 0.0_rp
                     P(i,j,k)   = 1.0_rp

                   if(z(k) <= 0.0_rp) then
                     U(i,j,k)   = u_inf
                   else
                     U(i,j,k)   = 0.0_rp
                   endif

                enddo
             enddo
          enddo


        end select

        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return
end subroutine init_I_stokes_problem


subroutine init_II_stokes_problem(dir,status)
! --------------------------------------------------------------------------------
! 
!       This subroutine initialize the velocity profile of the II Stokes problem. 
!       The parameter of the problem are:
!
!       [1] u_ref   : reference speed (<< 1 to avoid compressibility effects)
!       [2] omega   : wall frequency
!       [3] omega_c : corrected wall frequency
!
! --------------------------------------------------------------------------------
        implicit none
        character(len=1), intent(in)  :: dir
        integer         , intent(out) :: status
        real(rp)             :: omega, omega_c

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing II Stokes problem field: '

        omega = 0.25_rp*pi
        omega_c = sqrt((omega/(2*mu_inf)))
        
        selectcase(dir)

                case('x')
                do k = lbz,ubz
                   do j = lby,uby
                      do i =lbx,ubx

                         phi(i,j,k,1) = 1._rp
                           P(i,j,k)   = 1._rp 
                           U(i,j,k)   = u_inf*exp(-omega_c*y(j)) * cos(- omega_c*y(j))
                           V(i,j,k)   = 0._rp
                           W(i,j,k)   = 0._rp

                      enddo
                   enddo
                enddo

                case('y')
                do k = lbz,ubz
                   do j = lby,uby
                      do i =lbx,ubx

                         phi(i,j,k,1) = 1._rp
                           P(i,j,k)   = 1._rp 
                           U(i,j,k)   = 0._rp
                           V(i,j,k)   = u_inf*exp(-omega_c*x(i)) * cos(- omega_c*x(i))
                           W(i,j,k)   = 0._rp

                      enddo
                   enddo
                enddo

                case('z')
                do k = lbz,ubz
                   do j = lby,uby
                      do i =lbx,ubx

                         phi(i,j,k,1) = 1._rp
                           P(i,j,k)   = 1._rp 
                           U(i,j,k)   = u_inf*exp(-omega_c*z(k)) * cos(- omega_c*z(k))
                           V(i,j,k)   = 0._rp
                           W(i,j,k)   = 0._rp

                      enddo
                   enddo
                enddo

        end select

        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return
end subroutine init_II_stokes_problem


subroutine init_steady_couette(dir,status)
        implicit none
        character(len=1), intent(in)  :: dir
        integer         , intent(out) :: status

        real(rp) :: h

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing steady Couette field: '
        selectcase(dir)

                case('x')

                h = ymax - ymin

                do k = lbz,ubz
                   do j = lby,uby
                      do i = lbx,ubx

                         phi(i,j,k,1) = 1.0_rp
                           U(i,j,k)   = u_inf * (y(j)-ymin)/h
                           V(i,j,k)   = 0._rp
                           W(i,j,k)   = 0._rp
                           P(i,j,k)   = 1.0_rp

                      enddo
                   enddo
                enddo

                case('y')

                h = xmax - xmin

                do k = lbz,ubz
                   do j = lby,uby
                      do i = lbx,ubx

                         phi(i,j,k,1) = 1.0_rp
                           U(i,j,k)   = 0._rp
                           V(i,j,k)   = u_inf * (x(i)-xmin)/h
                           W(i,j,k)   = 0._rp
                           P(i,j,k)   = 1.0_rp

                      enddo
                   enddo
                enddo

                case('z')

                h = zmax - zmin

                do k = lbz,ubz
                   do j = lby,uby
                      do i = lbx,ubx

                         phi(i,j,k,1) = 1.0_rp
                           U(i,j,k)   = u_inf * (z(k)-zmin)/h
                           V(i,j,k)   = 0._rp
                           W(i,j,k)   = 0._rp
                           P(i,j,k)   = 1.0_rp

                      enddo
                   enddo
                enddo

        endselect

        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return
end subroutine init_steady_couette





subroutine init_shock_vortex(status)
! ------------------------------------------------------------------
!
!       This subroutine initialize a shock-vortex interaction.
!       Ref: Andrey Rault & al. 
!            "Shock vortex interaction at high mach number"
!            Journal of scientific computing
!
! ------------------------------------------------------------------
        implicit none
        integer, intent(out) :: status
        real(rp) :: rho1, p1, T1, M1, u1, v1, w1 !< left - shock state
        real(rp) :: rho2, p2, T2, M2, u2, v2, w2 !< right- shock state

        ! local vars
        real(rp) :: theta, r, vm, v_th, c1, c2, alpha, beta

        ! parameters
        real(rp), parameter :: x0    = 0.5_rp   !< shock position
        real(rp), parameter :: xc    = 0.25_rp  !< vortex x center 
        real(rp), parameter :: yc    = 0.5_rp   !< vortex y center
        real(rp), parameter :: a     = 0.075_rp !< vortex inner radius
        real(rp), parameter :: b     = 0.175_rp !< vortex outer radius
        real(rp), parameter :: Mv    = 1.7_rp   !< vortex mach number of max radial speed

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing shock vortex: '
        
        M1 = Mach
        vm = Mv * sqrt(gamma0)


        ! left state
        rho1 = 1._rp
          p1 = 1._rp
          T1 = p1/rho1
          u1 = M1 * sqrt(gamma0 * T1)
          v1 = 0._rp
          w1 = 0._rp

        ! right state
        call Rankine_Hugoniot(rho1, p1, T1, M1, rho2, p2, T2, M2)

          u2 = M2 * sqrt(gamma0 * T2)
          v2 = 0._rp
          w2 = 0._rp
        
        ! temperature parameters
        alpha = (gamma0-1._rp)/gamma0 * (vm/a)**2
        beta  = (gamma0-1._rp)/gamma0 * ((vm*a)/(a**2 - b**2))**2

        c2 = T1 + beta * 2*(b**2)*log(b)
        c1 = beta * ( -b**4/(2*a**2) - 2*b**2*log(a) + 0.5_rp*a**2) - 0.5_rp*alpha * a**2 + c2
        
        do k = lbz,ubz
           do j = lby,uby
              do i = lbx, ubx

                      if(x(i) <= x0) then
                              
                              ! --- vortex field ------------------------------------------------------
                              v_th = 0._rp
                              T(i,j,k) = T1
                              r = sqrt((x(i) - xc)**2 + (y(j) - yc)**2)

                              if    (r .le. a) then
                                      v_th = vm * r/a

                                      T(i,j,k) = alpha * 0.5_rp*r**2 + c1

                              elseif(r .le. b .and. r .ge. a) then
                                      v_th = vm * a/(a**2 - b**2)*(r - b**2/r)

                                      T(i,j,k) = beta * (- 0.5_rp*b**4/r**2 - 2*b**2*log(r) + 0.5_rp*r**2) + c2

                              endif
                              ! ------------------------------------------------------------------------
                              theta  = atan2((y(j)-yc),(x(i)-xc))

                                U(i,j,k)   = u1   + v_th * (-sin(theta))
                                V(i,j,k)   = v1   + v_th * ( cos(theta))
                                W(i,j,k)   = w1
                              phi(i,j,k,1) = rho1 * (T(i,j,k)/T1) ** ((    1._rp)/(gamma0-1._rp))
                                P(i,j,k)   = p1   * (T(i,j,k)/T1) ** ((gamma0)/(gamma0-1._rp))

                      else
                              phi(i,j,k,1) = rho2
                                U(i,j,k)   = u2
                                V(i,j,k)   = v2
                                W(i,j,k)   = w2
                                P(i,j,k)   = p2
                      endif

              enddo
           enddo
        enddo
        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return
end subroutine init_shock_vortex


subroutine init_shear_layer(status)
! ----------------------------------------------------------------
!
!       This subroutine initilizes a smooth shear layer without using 
!       rondom numbers. 
!       REF: "A Validated Nonlinear Kelvin-Helmholtz Benchmark 
!             for Numerical Hydrodynamics"
! ----------------------------------------------------------------
        implicit none
        integer, intent(out) :: status
        real(rp), parameter  :: a     =  0.001_rp
        real(rp), parameter  :: sigma =  0.001_rp
        real(rp), parameter  :: y1    =  0.25_rp
        real(rp), parameter  :: y2    =  0.75_rp
        real(rp), parameter  :: AA    =  0.01_rp
        real(rp), parameter  :: p0    = 10.0_rp

        real(rp) :: htan1, htan2, exp1, exp2, yj, i_a, i_s2

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing shear layer: '

        i_a  = 1.0_rp/a
        i_s2 = (1.0_rp/sigma)**2

        do k = lbz,ubz
           do j = lby,uby

              yj = y(j)
              htan1 = tanh((yj-y1)*i_a)
              htan2 = tanh((yj-y2)*i_a)

              exp1 = exp(-(yj-y1)**2*i_s2)
              exp2 = exp(-(yj-y2)**2*i_s2)

              do i = lbx, ubx

                 phi(i,j,k,1) = 1.0_rp
                 U  (i,j,k)   = u_inf * (htan1 - htan2 - 1.0_rp)
                 V  (i,j,k)   = AA * sin(2*pi*x(i)) * (exp1 + exp2)
                 W  (i,j,k)   = 0.0_rp
                 P  (i,j,k)   = p0

              enddo
           enddo
        enddo
        

        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return
end subroutine init_shear_layer

subroutine init_piecewise_shear_layer(status)
        implicit none
        integer, intent(out) :: status
        real(rp)             :: u_n, u_s, drho

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing piecewise shear layer: '

        u_n = u_inf
        u_s = 0.5_rp*u_n
        drho = 0.2_rp
        
        do k = lbz,ubz
           do j = lby,uby
              do i = lbx,ubx
                
                 phi(i,j,k,1) = 0.5_rp*(2.0_rp - drho*tanh(100*y(j)))
                 U  (i,j,k)   = 0.5_rp*((u_n+u_s) + (u_n-u_s)*tanh(100*y(j)))
                 V  (i,j,k)   = 0.1_rp*u_inf*sin(10*pi*x(i)/lx)
                 W  (i,j,k)   = 0.0_rp
                 P  (i,j,k)   = 1.0_rp

              enddo
           enddo
        enddo

        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return
end subroutine init_piecewise_shear_layer

subroutine init_turbulent_channel(status)
! ---------------------------------------------------------
!
!       Turbulent channel initialization
!       REFERENCE: DANS.HENNINGSON, 
!                 "On turbulent spots in plane Poiseuille flow"
!                 pag 186
!
! ---------------------------------------------------------
        implicit none
        integer, intent(out) :: status

        ! local declarations
        real(rp) :: dp_dx, c_1, c_2, uy,psi_y, psi_z, exp_, xi, yj, zk
        
        if(rank == root) write(*,'(A,A,A)', advance = 'no') &
        ' Initializing turbulent channal: '

        dp_dx = -1.5_rp*2._rp*u_inf
        c_1   =   0.0_rp
        c_2   = - 0.5_rp*dp_dx
        
        do k       = lbz,ubz
           do j    = lby,uby
              do i = lbx,ubx

                 xi = x(i)
                 yj = y(j)
                 zk = z(k)

                 phi(i,j,k,1) = 1.0_rp
                 U  (i,j,k)   = 0.0_rp
                 V  (i,j,k)   = 0.0_rp
                 W  (i,j,k)   = 0.0_rp
                 P  (i,j,k)   = 1.0_rp

                 if(abs(yj) < 1.0_rp) then
                
                   uy = 0.5_rp * dp_dx * yj**2 + c_1 * yj + c_2

                   exp_ = 0.3_rp*u_inf*exp(-4*(4*xi**2 + zk**2))
                
                   psi_y = -2*yj*zk * exp_
                   psi_z = (yj**2*(8*zk**2-1.0_rp) - 8*zk**2 + 1.0_rp) * exp_
        
                   U(i,j,k) = uy
                   V(i,j,k) =   psi_z
                   W(i,j,k) = - psi_y

                 endif
              enddo
           enddo
        enddo

        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return
end subroutine init_turbulent_channel


subroutine init_nscbc_perturbation(status,dir)
        implicit none
        integer     , intent(out) :: status
        character(1), intent(in)  :: dir

        ! local
        real(rp), parameter :: x0    = 0.0_rp
        real(rp), parameter :: sigma = 0.5_rp
        real(rp)            :: x_, y_

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing perturbation: '
        
        selectcase(dir)

          case('x')
          do k = lbz,ubz
             do j = lby,uby
                do i = lbx,ubx

                   x_ = (x(i) - x0)**2/(2*sigma**2)

                   phi(i,j,k,1) = 1.0_rp + Mach * exp(-x_)
                     U(i,j,k)   = 0.0_rp + u_inf
                     V(i,j,k)   = 0.0_rp 
                     W(i,j,k)   = 0.0_rp
                     P(i,j,k)   = 1.0_rp

                enddo
             enddo
          enddo

          case('y')
          do k = lbz,ubz
             do j = lby,uby
                do i = lbx,ubx

                   y_ = (y(j) - x0)**2/(2*sigma**2)

                   phi(i,j,k,1) = 1.0_rp + Mach * exp(-y_)
                     U(i,j,k)   = 0.0_rp 
                     V(i,j,k)   = 0.0_rp + u_inf
                     W(i,j,k)   = 0.0_rp
                     P(i,j,k)   = 1.0_rp

                enddo
             enddo
          enddo

        endselect

        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return
end subroutine init_nscbc_perturbation

subroutine init_shock_bubble_2D(status)
        implicit none
        integer, intent(out) :: status

        ! local variables
        real(rp) :: r1, p1, T1, M1 
        real(rp) :: r2, p2, T2, M2 
        real(rp) :: vel1, vel2
        real(rp) :: x0              !< shock wave initial position
        real(rp), parameter :: density_ratio = 0.138_rp/1.29_rp ! Helium-air

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing shock bubble 2D: '

        ! applying Rankine-Hugoniot conditions in order to compute post-shocked field

        M1 = Mach
        x0 = -1.0_rp
        !
        ! ==== state 1
        !
        p1   = 1.0_rp
        r1   = 1.0_rp
        T1   = p1/r1
        vel1 = - M1 * sqrt(gamma0*T1)
        !
        ! ==== state 2
        !
        call Rankine_Hugoniot(r1,p1,T1,M1,r2,p2,T2,M2)

        vel2 = - M2 * sqrt(gamma0*T2)

        do k = lbz,ubz
           do j = lby,uby
              do i = lbx,ubx

                 if(x(i) .lt. x0 ) then 

                   phi(i,j,k,1) = r2
                     U(i,j,k)   = vel2 - vel1
                     V(i,j,k)   = 0.0_rp
                     W(i,j,k)   = 0.0_rp
                     P(i,j,k)   = p2

                 elseif(x(i) .ge. x0 ) then

                   if((x(i)-1.0_rp)**2 + y(j)**2 <= 0.5_rp**2) then
                     phi(i,j,k,1) = density_ratio * r1
                   else
                     phi(i,j,k,1) = r1
                   endif

                     U(i,j,k)   = 0.0_rp
                     V(i,j,k)   = 0.0_rp
                     W(i,j,k)   = 0.0_rp
                     P(i,j,k)   = p1
                 endif

              enddo
           enddo
        enddo
        
        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return
end subroutine init_shock_bubble_2D


subroutine init_ambient_supersonic(status)
        implicit none
        integer, intent(out) :: status

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing ambient supersonic: '
        
        do k = lbz,ubz
           do j = lby,uby
              do i = lbx,ubx

                 if(x(i) < -4.0_rp) then
                   U(i,j,k) = u_inf
                 else
                   U(i,j,k) = 0.0_rp
                 endif

                 phi(i,j,k,1) = 1.0_rp
                   V(i,j,k)   = 0.0_rp
                   W(i,j,k)   = 0.0_rp
                   P(i,j,k)   = 1.0_rp
                
              enddo
           enddo
        enddo

        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return
end subroutine init_ambient_supersonic

subroutine init_four_quadrants(status)
        implicit none
        integer, intent(out) :: status

        real(rp)               :: xi,yj
        real(rp), dimension(4) :: r_, u_, v_, w_, p_
        
        r_ = 1.0_rp
        u_ = 0.0_rp
        v_ = 0.0_rp
        w_ = 0.0_rp
        p_ = 1.0_rp
        if(ic == 'four_quadrants_A') then
          r_ = (/ 1.5_rp, 0.5323_rp, 0.138_rp, 0.5323_rp /)
          u_ = (/ 0.0_rp, 1.206_rp,  1.206_rp, 0.0_rp    /)
          v_ = (/ 0.0_rp, 0.0_rp,    1.206_rp, 1.206_rp  /)
          w_ = (/ 0.0_rp, 0.0_rp,    0.0_rp  , 0.0_rp    /)
          p_ = (/ 1.5_rp, 0.3_rp,    0.029_rp, 0.3_rp    /)

        elseif(ic == 'four_quadrants_B') then
          r_ = (/  1.0_rp , 2.0_rp ,  1.0_rp ,  3.0_rp  /)
          p_ = (/  1.0_rp , 1.0_rp ,  1.0_rp ,  1.0_rp  /)
          u_ = (/  0.75_rp, 0.75_rp, -0.75_rp, -0.75_rp /)
          v_ = (/ -0.5_rp , 0.5_rp ,  0.5_rp , -0.5_rp  /)
          w_ = (/ 0.0_rp  , 0.0_rp ,  0.0_rp ,  0.0_rp  /)

        else
          if(rank == root) print*, 'Four quadrants ', trim(ic), ' is not implemented.'
          call secure_stop
        endif

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing : '//"'"//trim(ic)//"' "

        do k = lbz,ubz
           do j = lby,uby
              do i = lbx,ubx

                 xi = x(i)
                 yj = y(j)
                 if    (xi > 0.0_rp .and. yj >  0.0_rp) then ! I   quadrant
                   phi(i,j,k,1) = r_(1)
                   U  (i,j,k)   = u_(1)
                   V  (i,j,k)   = v_(1)
                   W  (i,j,k)   = w_(1)
                   P  (i,j,k)   = p_(1)
                 elseif(xi < 0.0_rp .and. yj >  0.0_rp) then ! II  quadrant
                   phi(i,j,k,1) = r_(2)
                   U  (i,j,k)   = u_(2)
                   V  (i,j,k)   = v_(2)
                   W  (i,j,k)   = w_(2)
                   P  (i,j,k)   = p_(2)
                 elseif(xi <  0.0_rp .and. yj < 0.0_rp) then ! III quadrant
                   phi(i,j,k,1) = r_(3)
                   U  (i,j,k)   = u_(3)
                   V  (i,j,k)   = v_(3)
                   W  (i,j,k)   = w_(3)
                   P  (i,j,k)   = p_(3)
                 elseif(xi > 0.0_rp .and. yj < 0.0_rp) then ! IV  quadrant
                   phi(i,j,k,1) = r_(4)
                   U  (i,j,k)   = u_(4)
                   V  (i,j,k)   = v_(4)
                   W  (i,j,k)   = w_(4)
                   P  (i,j,k)   = p_(4)
                 endif

              enddo
           enddo
        enddo

        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return
end subroutine init_four_quadrants


subroutine init_shock_inflow(status)

        use math_tools_module     , only: newton_raphson
        use fluid_functions_module, only: inverse_rankine_hugoniot

        implicit none
        integer, intent(out) :: status

        ! local variables
        real(rp), parameter :: x0 = -1.0_rp                  !< shock initial position
        real(rp)            :: M1, r1, p1, T1, u1, v1, w1 !< state 1
        real(rp)            :: M2, r2, p2, T2, u2, v2, w2 !< state 2
        
        ! === check coherency
        if(trim(inflow_profile) .ne. 'shock_inflow') then
          if(rank == root) then
            print*, ' FATAL ERROR: '
            print*, ' Initial condition "shock_inflow" needs inflow_profile to be "shock_inflow" as well.'
          endif
          call secure_stop
        endif
        if(Mach>1.85_rp) then
          if(rank == root) then
          print*, 'Shock inflow with mach > 1.85_rp is not compatible with Rankine-Hugoniot'
          endif
          call secure_stop
        endif

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing a shock wave : '
        
        ! ambient condition
        r1 = 1.0_rp
        p1 = 1.0_rp
        T1 = p1/r1
        u1 = 0.0_rp
        v1 = 0.0_rp
        w1 = 0.0_rp

        ! compute post-shocked field bosed on the post shock Mach number

        call newton_raphson(inverse_rankine_hugoniot, 10000, toll_equality, Mach+1.0_rp, M1)
        call Rankine_Hugoniot(r1,p1,T1,M1,r2,p2,T2,M2)

        u2 = Mach*sqrt(gamma0*T2)
        v2 = 0.0_rp
        w2 = 0.0_rp


        do k = lbz,ubz
           do j = lby,uby
              do i = lbx,ubx

                 if     (x(i) < x0 ) then 

                   phi(i,j,k,1) = r2
                     U(i,j,k)   = u2
                     V(i,j,k)   = v2
                     W(i,j,k)   = w2
                     P(i,j,k)   = p2

                 elseif(x(i) >= x0 ) then

                   phi(i,j,k,1) = r1
                     U(i,j,k)   = u1
                     V(i,j,k)   = v1
                     W(i,j,k)   = w1
                     P(i,j,k)   = p1
                 endif

              enddo
           enddo
        enddo
        
        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return
end subroutine init_shock_inflow



subroutine init_shock_impact(status)
! --------------------------------------------------------------------------------
!
!       This subroutine initialize a steady shock wave in a domain posizion. 
!
! --------------------------------------------------------------------------------
        implicit none
        integer , intent(out) :: status

        ! local declarations
        real(rp)            :: x0
        real(rp)            :: r1, p1, T1, u1, M1
        real(rp)            :: r2, p2, T2, u2, M2

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing a steady shock wave : '

        ! check input
        if(Mach < 1.0_rp) then
          if(rank == root) then
            print*, ' Mach < 1._rp No possible steady shocks!'
          endif
          call secure_stop
        endif

        M1 = Mach
        x0 = -2.0_rp
        ! --------------|-------------- !
        !     (2)       |     (1)       !   <---- M1 = MACH
        ! --------------|-------------- !
        !               ^
        !               x0                

        !
        ! ==== state 1
        !
        p1 = 1.0_rp
        r1 = 1.0_rp
        T1 = p1/r1
        u1 = - M1 * sqrt(gamma0*T1)

        !
        ! ==== state 2
        !
        call Rankine_Hugoniot(r1,p1,T1,M1,r2,p2,T2,M2)

        u2 = - M2 * sqrt(gamma0*T2)
       
        x0  = -2.0_rp

        do k = lbz,ubz
           do j = lby,uby
              do i = lbx,ubx

                 if(x(i) <= x0 ) then 

                   phi(i,j,k,1) = r2
                     U(i,j,k)   = u2
                     V(i,j,k)   = 0.0_rp
                     W(i,j,k)   = 0.0_rp
                     P(i,j,k)   = p2

                 elseif(x(i) >= x0 ) then

                   phi(i,j,k,1) = r1
                     U(i,j,k)   = u1
                     V(i,j,k)   = 0.0_rp
                     W(i,j,k)   = 0.0_rp
                     P(i,j,k)   = p1

                 endif

              enddo
           enddo
        enddo

        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return
end subroutine init_shock_impact



subroutine init_swbli(status)
! -------------------------------------------------------------------------------------
!       
!       This subroutine initialises the SHOCK-WAVE boundary layer problem
!
! -------------------------------------------------------------------------------------
        use math_tools_module     , only: DegToRad, RadToDeg
        use fluid_functions_module, only: oblique_shock

        implicit none
        integer, intent(out) :: status

        ! parameters
        real(rp), parameter :: thetaw = 3.75_rp     !< wedge angle DEGREZ
        real(rp), parameter :: x1     = 100.0_rp    !< impingment point

        ! local variables
        real(rp) :: r1, p1, T1, u1, v1, M1 !< state 1 variables
        real(rp) :: r2, p2, T2, u2, v2, M2 !< state 2 variables
        real(rp) :: r3, p3, T3, u3, v3, M3 !< state 3 variables
        real(rp) :: beta2                  !< first shock angle
        real(rp) :: beta3                  !< secnd shock angle
        real(rp) :: theta                  !< wedge angle (in radiants)
        real(rp) :: slope1, slope2         !< slope of the first shock
        real(rp) :: shock1, shock2         !< slope of the reflected shock
        integer  :: i,j,k

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing shock wave boundary layer interaction: '
        
        !
        ! === STATE 1
        !
        M1 = Mach
        r1 = 1.0_rp
        p1 = 1.0_rp
        T1 = p1/r1
        u1 = M1 * sqrt(gamma0*T1)
        v1 = 0.0_rp

        !
        ! === STATE 2
        !
        theta = DegToRad(thetaw)
        call oblique_shock(r1,p1,T1,M1,theta,r2,p2,T2,M2,beta2)
        
        slope1 = tan(pi-beta2)

        u2 = M2 * sqrt(gamma0*T2) * cos(theta)
        v2 = M2 * sqrt(gamma0*T2) * sin(theta)

        !
        ! === STATE 3
        !
        call oblique_shock(r2,p2,T2,M2,theta,r3,p3,T3,M3,beta3)
        
        slope2 = tan(beta3 - theta)

        u3 = M3 * sqrt(gamma0*T3)
        v3 = 0.0_rp
        
        do k       = lbz,ubz
           do j    = lby,uby
              do i = lbx,ubx

                   shock1 = y(j) - ymin - slope1 * (x(i) - x1)
                   shock2 = y(j) - ymin - slope2 * (x(i) - x1)
                
                   ! REGION 1
                   if(shock1 < 0.0_rp) then
                           
                     phi(i,j,k,1) = r1
                       U(i,j,k)   = u1
                       V(i,j,k)   = v1
                       W(i,j,k)   = 0.0_rp
                       P(i,j,k)   = p1
                   
                   ! REGION 2
                   elseif(shock1 > 0.0_rp .and. shock2 < 0.0_rp) then

                     phi(i,j,k,1) =   r3
                       U(i,j,k)   =   u3
                       V(i,j,k)   =   0.0_rp
                       W(i,j,k)   =   0.0_rp
                       P(i,j,k)   =   p3

                   !! REGION 3
                   else

                     phi(i,j,k,1) =   r2
                       U(i,j,k)   =   u2
                       V(i,j,k)   = - v2
                       W(i,j,k)   =   0.0_rp
                       P(i,j,k)   =   p2

                   endif

              enddo
           enddo
        enddo

        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return
end subroutine init_swbli
        


subroutine init_LaminarBoundaryLayer(status)

        use fluid_functions_module, only: BlasiusProfile
        use FileModule

        implicit none
        integer, intent(out) :: status

        ! local declarations
        real(rp), dimension(1-GN:ny+GN) :: LaminarVel_inc
        real(rp), dimension(1-GN:ny+GN) :: LaminarVel_cmp
        real(rp), dimension(1-GN:ny+GN) :: LaminarTmp
        real(rp), dimension(1-GN:ny+GN) :: LaminarRho

        real(rp), allocatable, dimension(:) :: y_tmp
        real(rp)                            :: Tw, Rw
        type(FileType)                      :: Blasius
        integer                             :: i,j,k

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing Laminar boundary layer: '
        
        selectcase(inflow_profile)
          case('blasius')
          call BlasiusProfile(LaminarVel_inc, LaminarTmp, LaminarRho)

          case('pohlhausen')
          call PohlhausenProfile(LaminarVel_inc)
          call CompressibleCorrection(ny,.false.,Mach,Prandtl,Trat,Tw,Rw,LaminarVel_inc,LaminarVel_cmp,LaminarRho,LaminarTmp)

          case default
          print*, ' Initial condition ', trim(inflow_profile), ' is not implemented.'
          stop

        endselect

        do k       = lbz,ubz
           do j    = lby,uby
              do i = lbx,ubx

                 phi(i,j,k,1) = LaminarRho(j)
                   U(i,j,k)   = LaminarVel_inc(j)
                   V(i,j,k)   = 0.0_rp
                   W(i,j,k)   = 0.0_rp
                   P(i,j,k)   = 1.0_rp

              enddo
           enddo
        enddo
        ! 
        ! === print blasius profile
        ! 
        if(rank == root) then

          Blasius%name = 'BlasiusProfile'
          Blasius%dir  = trim(data_dir)//'/BLASIUS_PROFILE'
          call OpenNewFile(Blasius,it)
        
          allocate(y_tmp(1-GN:ny+GN))
          call compute_grid_point(y_tmp,ny,lbound(y_tmp,1),ubound(y_tmp,1),ymin,ymax,stretching_par,gridpoint_y)

          do j = 1-GN,ny+GN
             write(Blasius%unit,*) y_tmp(j), LaminarVel_inc(j), LaminarVel_cmp(j), LaminarTmp(j), LaminarRho(j)
          enddo

          deallocate(y_tmp)

        endif   

        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return
end subroutine init_LaminarBoundaryLayer




subroutine init_TurbulentBoundaryLayerMean(status)

        use fluid_functions_module, only: MuskerProfile, CompressibleCorrection, BoundaryLayerQuantities, Sutherland
        use real_to_integer_module, only: locate
        use interpolation_module  , only: polint
        use FileModule 

        implicit none
        integer, intent(out) :: status

        ! local declarations
        real(rp), dimension(:,:,:), allocatable :: tempPrim
        real(rp), dimension(:)    , allocatable :: delta_array
        real(rp), dimension(:)    , allocatable :: deltaVArray
        real(rp), dimension(:)    , allocatable :: frict_array
        real(rp), dimension(:)    , allocatable :: theta_array
        real(rp), dimension(:)    , allocatable :: y_tmp
        real(rp), dimension(0:ny)               :: UPlusInc, UplusCmp, RhoY, TmpY, velY
        real(rp), parameter                     :: toll = 1.0E-14_rp
        real(rp)                                :: tWall, rWall
        real(rp)                                :: delta, deltaV, thrat, deltaStar, theta, cf
        real(rp)                                :: dx, dy, ReTau_new, ReTau_old, err
        real(rp)                                :: ue, ReTauI, yl, errR, errU
        real(rp)                                :: ri, ui, vi_j00, vi_jm1, vi_j
        real(rp)                                :: r_, ir, u_, v_
        integer                                 :: i,j,jj,jjj, m, jg, ji
        integer                                 :: itr, itrmax = 30, ierr = 0

        real(rp) :: delta_REF, theta_REF, deltaVREF, frict_REF
        type(fileType) :: BlData
#ifdef DEBUG
        type(FileType) :: inflow, frictionFile
        real(rp)       :: u_y, mu_wall
#endif

        !
        ! === allocate temporary variables
        !
        allocate(delta_array(-GN:nx+GN), &
                 deltaVArray(-GN:nx+GN), &
                 frict_Array(-GN:nx+GN), &
                 theta_array(-GN:nx+GN), stat = ierr)
        if(ierr.ne.0) stop ' Allocation error'
        
        allocate(tempPrim(-GN:nx+GN, 1-GN:ny+GN, 3), &
                 y_tmp(0:ny), stat = ierr)
        if(ierr.ne.0) stop ' Allocation error'
        
        call compute_grid_point(y_tmp,ny,0,ny,ymin,ymax,stretching_par,gridpoint_y)
        y_tmp(0) = 0.0_rp
        

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing Turbulent boundary layer: '

        !
        ! === getting the REFERENCE inflow profile
        !
        call MuskerProfile(ny,y_tmp,ReTau,UplusInc)
        call CompressibleCorrection(ny,.true.,Mach,Prandtl,Trat,tWall,rWall,UPlusInc,UPlusCmp,RhoY,TmpY)

#ifdef DEBUG
        if(rank == root) then
          inflow%name = 'inflowProfile'
          inflow%dir  = trim(data_dir)//'/INITIAL_CONDITION'
          call OpenNewFile(inflow,it)
          do j = 0,ny
             write(inflow%unit,*) y_tmp(j), RhoY(j), uPlusCmp(j)/maxval(uPlusCmp), TmpY(j)
          enddo
          call CloseFile(inflow)
        endif
#endif

        !
        ! === computing the growth of the boundary layer
        !
        call BoundaryLayerQuantities(y_tmp,RhoY,UPlusCmp,rWall,thrat,deltaStar,cf)
        delta_REF = 1.0_rp
        deltaVREF = delta_REF/ReTau
        frict_REF = cf
        theta_REF = thrat*delta_REF

        ReTau_old = ReTau - 1.0_rp
        ReTau_new = ReTau
        dx        = Lx/real(nx,rp)

        ! get the node-0 values integrating from REF in half cell space
        itr = 0
        err = huge(1.0_rp)
        do while(err > toll .and. itr < itrmax)

           call MuskerProfile(ny,y_tmp,ReTau_new,UplusInc)
           call CompressibleCorrection(ny,.true.,Mach,Prandtl,Trat,tWall,rWall,UPlusInc,UPlusCmp,RhoY,TmpY)
           call BoundaryLayerQuantities(y_tmp,RhoY,UPlusCmp,rWall,thrat,deltaStar,cf)

           theta     = theta_REF - 0.25_rp*0.5_rp*dx*(cf+frict_REF)
           delta     = theta/thrat
           deltaV    = deltaVREF*sqrt(frict_REF/cf)
           reTau_new = delta/deltaV

           err = abs(ReTau_new - ReTau_old) 

           ReTau_old = ReTau_new
           itr = itr + 1

        enddo
        theta_array(0) = theta
        frict_array(0) = cf
        delta_array(0) = delta
        deltaVArray(0) = deltaV

        ! get the node 1 values integrating from REF in half cell space
        itr = 0
        err = huge(1.0_rp)
        do while(err > toll .and. itr < itrmax)

           call MuskerProfile(ny,y_tmp,ReTau_new,UplusInc)
           call CompressibleCorrection(ny,.true.,Mach,Prandtl,Trat,tWall,rWall,UPlusInc,UPlusCmp,RhoY,TmpY)
           call BoundaryLayerQuantities(y_tmp,RhoY,UPlusCmp,rWall,thrat,deltaStar,cf)

           theta     = theta_REF + 0.25_rp*0.5_rp*dx*(cf+frict_REF)
           delta     = theta/thrat
           deltaV    = deltaVREF*sqrt(frict_REF/cf)
           reTau_new = delta/deltaV

           err = abs(ReTau_new - ReTau_old) 

           ReTau_old = ReTau_new
           itr = itr + 1

        enddo
        theta_array(1) = theta
        frict_array(1) = cf
        delta_array(1) = delta
        deltaVArray(1) = deltaV

        !
        ! integrating from -1 to -GN
        !
        do i = -1,-GN,-1

           itr = 0
           err = huge(1.0_rp)
           do while(err > toll .and. itr < itrmax)

              call MuskerProfile(ny,y_tmp,ReTau_new,UplusInc)
              call CompressibleCorrection(ny,.true.,Mach,Prandtl,Trat,tWall,rWall,UPlusInc,UPlusCmp,RhoY,TmpY)
              call BoundaryLayerQuantities(y_tmp,RhoY,UPlusCmp,rWall,thrat,deltaStar,cf)

              theta     = theta_array(i+1) - 0.25_rp*dx*(cf+frict_array(i+1))
              delta     = theta/thrat
              deltaV    = deltaVArray(i+1)*sqrt(frict_array(i+1)/cf)
              reTau_new = delta/deltaV

              err = abs(ReTau_new - ReTau_old) 

              ReTau_old = ReTau_new
              itr = itr + 1

           enddo
           theta_array(i) = theta
           frict_array(i) = cf
           delta_array(i) = delta
           deltaVArray(i) = deltaV

        enddo
        !
        ! integrating from 2 to nx+GN
        !
        do i = 2,nx+GN

           itr = 0
           err = huge(1.0_rp)
           do while(err > toll .and. itr < itrmax)

              call MuskerProfile(ny,y_tmp,ReTau_new,UplusInc)
              call CompressibleCorrection(ny,.true.,Mach,Prandtl,Trat,tWall,rWall,UPlusInc,UPlusCmp,RhoY,TmpY)
              call BoundaryLayerQuantities(y_tmp,RhoY,UPlusCmp,rWall,thrat,deltaStar,cf)

              theta     = theta_array(i-1) + 0.25_rp*dx*(cf+frict_array(i-1))
              delta     = theta/thrat
              deltaV    = deltaVArray(i-1)*sqrt(frict_array(i-1)/cf)
              ReTau_new = delta/deltaV

              err       = abs(ReTau_new - Retau_old)
              ReTau_old = ReTau_new
              itr       = itr + 1

           enddo
           theta_array(i) = theta
           delta_array(i) = delta
           deltaVarray(i) = deltaV
           frict_array(i) = cf

        enddo

        if(rank == root) then
          bldata%name = 'blData'
          bldata%dir  = trim(data_dir)//'/INITIAL_CONDITION'
          call OpenNewFile(blData,it)
          write(blData%unit,*)   '# x    theta   delta   delta_v   cf   ReTau'
          do i = -GN,nx+GN
             write(blData%unit,*) Lx*(i-0.5_rp)/real(nx,rp), theta_array(i), &
                     delta_array(i), deltaVArray(i), frict_array(i), delta_array(i)/deltaVarray(i)
          enddo
          call CloseFile(BlData)
        endif

        !
        ! === interpolation on the mesh
        !
        do i = -GN,nx+GN
           delta  = delta_array(i)
           deltaV = deltaVArray(i)
           ReTauI = delta/deltaV

           call MuskerProfile(ny,y_tmp,ReTauI,UplusInc)
           call CompressibleCorrection(ny,.true.,Mach,Prandtl,Trat,tWall,rWall,UPlusInc,UPlusCmp,RhoY,TmpY)

           ue = maxval(UPlusCmp)
           VelY = UPlusCmp/ue

           do j = 0, ny

              yl = y_tmp(j)/delta
              call locate(y_tmp,0,ny,yl,jj)

              m = 4
              jjj = min(max(jj-(m-1)/2,1),ny+1-m)
              call polint(y_tmp(jjj),rhoY(jjj),m,yl,ri,errR)
              call polint(y_tmp(jjj),velY(jjj),m,yl,ui,errU)

              tempPrim(i,j,1) = ri
              tempPrim(i,j,2) = ri*ui
              tempPrim(i,j,3) = 0.0_rp

           enddo
        enddo

        ! === integrate rhov
        do j    = 1   , ny
           do i = 1-GN, nx+GN

              dy = y_tmp(j) - y_tmp(j-1)

              vi_j00 = (tempPrim(i,j  ,2) - tempPrim(i-1,j  ,2))/dx
              vi_jm1 = (tempPrim(i,j-1,2) - tempPrim(i-1,j-1,2))/dx
              vi_j   = 0.5_rp*(vi_j00 + vi_jm1)

              tempPrim(i,j,3) = tempPrim(i,j-1,3) - vi_j * dy

           enddo
        enddo

        ! extapolation on north face ghosts
        do    j = 1,GN
           do i = -GN,nx+GN

              jg = ny + j      !< ghost node
              ji = ny - (j-1)  !< inner node

              tempPrim(i,jg,1) =  tempPrim(i,ji,1)
              tempPrim(i,jg,2) =  tempPrim(i,ji,2)
              tempPrim(i,jg,3) =  tempPrim(i,ji,3)

           enddo
        enddo

        ! extapolation on south face ghosts
        do    j = 1,GN
           do i = -GN,nx+GN

              jg = 1 - j      !< ghost node
              ji = 1 + (j-1)  !< inner node

              tempPrim(i,jg,1) =  tempPrim(i,ji,1)
              tempPrim(i,jg,2) = -tempPrim(i,ji,2)
              tempPrim(i,jg,3) = -tempPrim(i,ji,3)

           enddo
        enddo

        !
        ! === distribute the results within the procs
        !
        do       k = lbz,ubz
           do    j = lby,uby
              do i = lbx,ubx

                 r_ = tempPrim(i,j,1)
                 ir = 1.0_rp/r_
                 u_ = tempPrim(i,j,2)*ir
                 v_ = tempPrim(i,j,3)*ir
                
                 phi(i,j,k,1) = r_
                 U  (i,j,k)   = u_inf*u_
                 V  (i,j,k)   = u_inf*v_
                 W  (i,j,k)   = 0.0_rp
                 P  (i,j,k)   = 1.0_rp

              enddo
           enddo
        enddo

#ifdef DEBUG
        ! check consistency of the friction coefficient
        frictionFile%name = 'InitFriction'//trim(str(rank))
        frictionFile%dir  = trim(data_dir)//'/INITIAL_CONDITION'
        call OpenNewFile(frictionFile,it)
        do i = lbx,ubx
           u_y     = U(i,sy,sz)/y_tmp(1)
           mu_wall = Sutherland(tWall)
           write(frictionFile%unit,*) x(i), mu_wall*mu_inf*u_y/q_inf
        enddo
        call CloseFile(frictionFile)
#endif

        deallocate(delta_array,deltaVArray,frict_Array,theta_array)
        deallocate(y_tmp,tempPrim)

        if(rank == root) write(*, '(A)') 'done!'
        status = 0

        return
end subroutine init_TurbulentBoundaryLayerMean




subroutine init_TBLNoise(status)

        use df_module
        use random_module
        use FileModule

        implicit none
        integer, intent(out) :: status

        ! local declarations
        real(rp), dimension(:,:,:,:), allocatable :: vf
        real(rp), dimension(:,:,:,:), allocatable :: uf
        real(rp), dimension(:)      , allocatable :: y_tmp
        real(rp), dimension(3)                    :: exppar1, exppar2, par1, par2
        real(rp)                                  :: um, vm , wm, Tm, rm
        real(rp)                                  :: up, vp , wp, Tp
        real(rp)                                  :: dx
        integer                                   :: err = 0, i, j, k, jg, ji, m
#ifdef DEBUG
        type(FileType)         :: reyStressFile
        real(rp), dimension(3) :: rms, rmean, rmeansq
        real(rp)               :: rey12, rey13, rey23
#endif
        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing a turbulent-consistent field: '
        
        ! allocations
        allocate(y_tmp(0:ny), stat = err)
        if(err .ne. 0) stop ' Allocation error in init_TBLNoise'

        allocate(DF%Rnd3D(3, 1-GN:nx+GN+1, 1-Df%N:ny+DF%N, 1-DF%N:nz+DF%N), stat = err)
        if(err .ne. 0) stop ' Allocation error in init_TBLNoise'

        allocate(vf(3,1-GN:nx+GN+1, 1:ny, 1:nz), stat = err)
        if(err .ne. 0) stop ' Allocation error in init_TBLNoise'

        allocate(uf(3,1-GN:nx+GN+1, 1-GN:ny+GN, 1-GN:nz+GN), stat = err)
        if(err .ne. 0) stop ' Allocation error in init_TBLNoise'
        
        ! get a temporary grid
        call compute_grid_point(y_tmp,ny,0,ny,ymin,ymax,stretching_par,gridpoint_y)
        y_tmp(0) = 0.0_rp

        ! create a 3D random field
        call DFRandomField3D(nx,ny,nz,gn,DF)

        ! filter the random data 
        call DFConvolution3D(1-GN,nx+GN+1,1,ny,1,nz,DF,vf)

        ! enforce statistics
        call DFEnforceReynoldsStresses3D(1-GN,nx+GN+1,1,ny,1,nz,vf,DF,uf)

        ! extrapolation on north face
        do j = 1,GN
           jg = ny + j      !< ghost node
           ji = ny - (j-1)  !< inner node

           uf(:,:,jg,:) = uf(:,:,ji,:)
        enddo
        ! wall on south face
        do j = 1,GN
           jg = 1 - j      !< ghost node
           ji = 1 + (j-1)  !< inner node

           uf(:,:,jg,:) = - uf(:,:,ji,:)
        enddo
        ! periodicity along z
        do k = 1,GN
           uf(:,:,:,nz+k) = uf(:,:,:,k)
           uf(:,:,:,1-k)  = uf(:,:,:,nz+1-k)
        enddo

        ! x correlation via Castro equation
        dx      = gbl_min_step(1)
        exppar1 = -0.5_rp*pi*dx/DF%xlen
        exppar2 = -pi*dx/DF%xlen
        par1    = exp(exppar1)
        par2    = sqrt(1.0_rp-exp(exppar2))

        do i=nx+GN,1,-1
           do m = 1,3
              uf(m,i,:,:) = uf(m,i+1,:,:)*par1(m) + uf(m,i,:,:)*par2(m)
           enddo
        enddo
        
        do k       = lbz,ubz
           do j    = lby,uby
              do i = lbx,ubx

                 rm = phi(i,j,k,1)
                 um = U(i,j,k)
                 vm = V(i,j,k)
                 wm = W(i,j,k)
                 Tm = 1.0_rp/rm 

                 up = u_inf*uf(1,i,j,k)
                 vp = u_inf*uf(2,i,j,k)
                 wp = u_inf*uf(3,i,j,k)
                 Tp = -0.75_rp*(gamma0-1.0_rp)/gamma0 * up * um / Tm
                
                 phi(i,j,k,1) = 1.0_rp/(Tm + Tp)
                 U(i,j,k)     = um + up
                 V(i,j,k)     = vm + vp
                 W(i,j,k)     = wm + wp

              enddo
           enddo
        enddo

        if(rank == root) write(*, '(A)') 'done!'

#ifdef DEBUG
        reyStressFile%name = 'rey_computed' 
        reyStressFile%dir  = trim(data_dir)//'/DF_REYNOLDS_STRESS'
        call OpenNewFile(ReyStressFile,it)
        do j=1,ny
         rmean   = 0._rp
         rmeansq = 0._rp
         rey12   = 0._rp
         rey13   = 0._rp
         rey23   = 0._rp
         do k=1,nz
          do i=1-GN,nx+GN
           do m=1,3
            rmean(m)   = rmean(m)   + uf(m,i,j,k)
            rmeansq(m) = rmeansq(m) + uf(m,i,j,k)**2
           enddo
           rey12 = rey12+uf(1,i,j,k)*uf(2,i,j,k)
           rey13 = rey13+uf(1,i,j,k)*uf(3,i,j,k)
           rey23 = rey23+uf(2,i,j,k)*uf(3,i,j,k)
          enddo
         enddo
         rmean   = rmean/(nz)/(nx+2*gn)
         rmeansq = rmeansq/(nz)/(nx+2*gn)
         rms     = sqrt(rmeansq-rmean*rmean)
         rey12   = rey12/(nz)/(nx+2*gn)
         rey13   = rey13/(nz)/(nx+2*gn)
         rey23   = rey23/(nz)/(nx+2*gn)
         rey12   = rey12-rmean(1)*rmean(2)
         rey13   = rey13-rmean(1)*rmean(3)
         rey23   = rey23-rmean(2)*rmean(3)
         write(reyStressFile%unit,*) &
         y_tmp(j),(rms(m)**2,m=1,3),rey12,rey13,rey23,(rmean(m),m=1,3)
        enddo
        call CloseFile(ReyStressFile)
#endif


        deallocate(y_tmp,uf,vf,DF%Rnd3D)
        status = 0

        return
end subroutine init_TBLNoise











end module ic_module
