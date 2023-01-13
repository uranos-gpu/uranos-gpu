module bc_module
! -------------------------------------------------------------
!
!       This module contains all the routine usefull to update
!       the boundary conditions of a subdomain
!
! -------------------------------------------------------------
use parameters_module
use storage_module
use mpi_module
use fluid_functions_module
use inflow_module
use mpi_comm_module
use profiling_module

contains
subroutine set_bc_conditions
        implicit none
        type(face_type), dimension(6) :: all_bound
        type(face_type)               :: bound
        integer                       :: f
        real(rp)                      :: M_Max, Twall, rfc
#ifdef TIME
        call mpi_stime(s_bcs_time)
#endif
        call StartProfRange("set_bc_conditions")


        all_bound%node = (/ sx , ex , sy , ey , sz ,  ez/)
        all_bound%norm = (/ -1 ,  1 , -1 ,  1 , -1 ,  1 /)
        all_bound%face = (/ 'W', 'E', 'S', 'N', 'B', 'F'/)
        
        !
        ! === maximum mach number computation in the whole domain
        !
        if(trim(bc(1)) == 'nscbc_inflow_relaxed') call compute_max_Mach(M_Max)
        if(trim(L1_wave) == 'poinsot_model') then
          if(bc(1) == 'nscbc_outflow' .or. &
             bc(2) == 'nscbc_outflow' .or. &
             bc(3) == 'nscbc_outflow' .or. &
             bc(4) == 'nscbc_outflow') then
             call compute_max_Mach(M_Max)
          endif
        endif

        do f = 1, nb_neigh
           if(my_neighbour(f) == MPI_PROC_NULL) then
                
             ! store everything a scalar variable
             bound%node = all_bound(f)%node
             bound%norm = all_bound(f)%norm
             bound%face = all_bound(f)%face

             selectcase(trim(bc(f)))

             case('neumann')
               call neumann(phi,bound) !yes acc

             case('adiabatic_wall')
               call adiabatic_wall(phi,bound) !yes acc

             case('nscbc_adiabatic_wall')
               if(restart_flag) then
                 call adiabatic_wall(phi,bound) !yes acc
               else
                 Twall = Trat*(1.0_rp + Prandtl**(1.0_rp/3.0_rp)*0.5_rp*(gamma0-1.0_rp)*Mach**2)
                 call nscbc_wall(phi,bound,Twall) !yes acc
               endif

             case('nscbc_isothermal_wall')
               Twall = 1.0_rp
               if(restart_flag) then
                 call isothermal_wall(phi,bound,Twall) !yes acc
               else
                 call nscbc_wall(phi,bound,Twall) !yes acc
               endif

             case('isothermal_wall'             , &
                  'dws_isothermal'              , &
                  'idws_isothermal'             , &
                  'static_dws_isothermal')

                  tWall = 1.0_rp
                  call isothermal_wall(phi,bound,Twall) !yes acc

             case('dws_adiabatic'        , &
                  'idws_adiabatic'       , &
                  'static_dws_adiabatic' , &
                  'istatic_dws_adiabatic', &
                  'adiabatic_wall_fix_temperature')

                  rfc = Prandtl**(1.0_rp/3.0_rp)
                  tWall = Trat*(1.0_rp + rfc*0.5_rp*(gamma0-1.0_rp)*Mach**2)
                  call isothermal_wall(phi,bound,Twall) !yes acc

             case('slip_wall')
               call slip_wall(phi,bound) !yes acc

             case('postshock_slipwall')
               call postshock_slipwall(phi,bound)

             case('oblique_shock')
               call oblique_shock_bc(phi,bound)

             case('flate_plate')
               call flate_plate_bc(phi,bound)

             case('travelling_shock_wave')
               call travelling_shock_wave(phi,bound)

             case('nscbc_inflow')
               if(restart_flag) then
                 call neumann(phi,bound) !yes acc
               else
                 call inflow_bc(phi,iFlow,bound)
               endif

             case('nscbc_inflow_relaxed')
               if(restart_flag) then !binary equality is avoid in restart nscbc
                 call neumann(phi,bound) !yes acc
               else
                 call nscbc_inflow_relaxed(phi,iFlow,bound,M_max)
               endif

             case('supersonic_inflow')
               if(restart_flag) then !binary equality is avoid in restart nscbc
                 call neumann(phi,bound) !yes acc
               else
                 call inflow_bc(phi,iFlow,bound) !yes acc
               endif

             case('nscbc_outflow')
                 if(Mach>1.0_rp) then
                   call supersonic_outflow(phi,bound) !yes acc
                 else
                   call outflow_nscbc_static(phi,bound)
                 endif

             case('shock_inflow')
               ! check the critical mach number
               if(Mach >= 1.0_rp) then
                 call inflow_bc(phi,iFLow,bound)
               else
                 call nscbc_inflow(phi,iFlow,bound)
               endif

             case('neumann_wall')
                call Neumann_wall(phi,bound)

             case default

                write(*,*) ' Boundary condition ', "'"//trim(bc(f))//"'", ' is not implemented.'
                stop

             endselect

           endif ! end if proc null
        enddo ! end do faces

        call EndProfRange
#ifdef TIME
        call mpi_etime(s_bcs_time,t_bcs_calls,t_bcs_time)
#endif
        return
end subroutine set_bc_conditions


subroutine mpi_bc_communications
        implicit none
#ifdef TIME
        call mpi_stime(s_mpi_time)
#endif
        call StartProfRange("mpi_bc_communications")

call  mpi_share(mpi_comm_cart,type_send_cons,type_recv_cons,my_neighbour,dims, &
                phi_bfr_send_E, phi_bfr_send_W, phi_bfr_recv_E, phi_bfr_recv_W, &
                phi_bfr_send_N, phi_bfr_send_S, phi_bfr_recv_N, phi_bfr_recv_S, &
                phi_bfr_send_B, phi_bfr_send_F, phi_bfr_recv_B, phi_bfr_recv_F, &
                phi)

        call EndProfRange
#ifdef TIME
        call mpi_etime(s_mpi_time,t_mpi_calls,t_mpi_time)
#endif
        return
end subroutine mpi_bc_communications


subroutine neumann(phi,b)
!--------------------------------------------------------------------------------
! 
!       This subroutine set Neumann omogeneous boundary condition
!       in a node.
!       INPUT: phi      < variable we want bound
!              bound    < derived data type of the boundary
!             
!--------------------------------------------------------------------------------
        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi
        type(face_type)                          , intent(in   ) :: b

        integer :: i,j,k,l
        integer :: ig, ii, jg, ji, kg,ki
        integer :: node, norm

        node = b%node
        norm = b%norm
        
        select case(b%face)
          case('E', 'W')
        
            !$acc parallel default(present)
            !$acc loop gang, vector collapse(4)
            do          l = 1,5
               do       k = lbz,ubz
                  do    j = lby,uby
                     do i = 1,GN
                        ig = node + norm*i
                        ii = node - norm*(i-1)
                        phi(ig,j,k,l) = phi(ii,j,k,l)
                     enddo
                  enddo
               enddo
            enddo
            !$acc end parallel

          case('S', 'N')
                 
            !$acc parallel default(present)
            !$acc loop gang, vector collapse(4)
            do          l = 1,5
               do       k = lbz,ubz
                  do    j = 1,GN
                     do i = lbx,ubx
                        jg = node + norm*j
                        ji = node - norm*(j-1)
                        phi(i,jg,k,l) = phi(i,ji,k,l)
                     enddo
                  enddo
               enddo
            enddo
            !$acc end parallel

          case('B', 'F')

            !$acc parallel default(present)
            !$acc loop gang, vector collapse(4)
            do          l = 1,5
               do       k = 1,GN
                  do    j = lby,uby
                     do i = lbx,ubx
                        kg = node + norm*k
                        ki = node - norm*(k-1)
                        phi(i,j,kg,l) = phi(i,j,ki,l)
                     enddo
                  enddo
               enddo
            enddo
            !$acc end parallel

        end select
        return
end subroutine neumann



subroutine adiabatic_wall(phi,b)
        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: phi
        type(face_type)                          , intent(in)    :: b

        real(rp) :: vel
        real(rp) :: r_in, u_in, v_in, w_in, p_in, T_in, irin
        real(rp) :: r_gh, u_gh, v_gh, w_gh, p_gh, T_gh
        integer  :: i,j,k
        integer  :: ig,jg,kg,ii,ji,ki
        integer  :: node, norm

        node = b%node
        norm = b%norm

        call set_wall_speed(vel,b)

        select case(b%face)

          case('E','W')
                  
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do        k = lbz,ubz
             do     j = lby,uby
                 do i = 1,GN

                    ig = node + norm*i
                    ii = node - norm*(i-1)

                    ! compute primitive variable in the internal nodes
                    r_in = phi(ii,j,k,1)
                    irin = 1._rp/r_in
                    u_in = phi(ii,j,k,2) * irin
                    v_in = phi(ii,j,k,3) * irin
                    w_in = phi(ii,j,k,4) * irin
                    p_in = (gamma0-1._rp)*(phi(ii,j,k,5)&
                          -0.5_rp*r_in*(u_in*u_in + v_in*v_in + w_in*w_in))
                    T_in = p_in*irin

                    ! compute primitives variables in the ghosts nodes
                    u_gh =         - u_in
                    v_gh = 2*vel   - v_in
                    w_gh =         - w_in
                    T_gh =           T_in  
                    p_gh = p_in
                    r_gh = p_gh/T_gh

                    phi(ig,j,k,1) = r_gh
                    phi(ig,j,k,2) = r_gh*u_gh 
                    phi(ig,j,k,3) = r_gh*v_gh 
                    phi(ig,j,k,4) = r_gh*w_gh 
                    phi(ig,j,k,5) = p_gh/(gamma0-1.0_rp) &
                                 + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)

                 enddo
             enddo
          end do
          !$acc end parallel

          case('S','N')

            !$acc parallel default(present)
            !$acc loop gang, vector collapse(3)
            do        k = lbz,ubz
               do     j = 1,GN
                   do i = lbx,ubx

                      jg = node + norm*j      !< ghost node
                      ji = node - norm*(j-1)  !< inner node
                         
                      ! compute primitive variable in the internal nodes
                      r_in = phi(i,ji,k,1)
                      irin = 1._rp/r_in
                      u_in = phi(i,ji,k,2) * irin
                      v_in = phi(i,ji,k,3) * irin
                      w_in = phi(i,ji,k,4) * irin
                      p_in = (gamma0-1._rp)*(phi(i,ji,k,5)&
                             -0.5_rp*r_in*(u_in*u_in + v_in*v_in + w_in*w_in))
                      T_in = p_in*irin

                      ! compute primitives variables in the ghosts nodes
                      u_gh = 2*vel   - u_in
                      v_gh =         - v_in
                      w_gh =         - w_in
                      T_gh =           T_in  
                      p_gh = p_in
                      r_gh = p_gh/T_gh

                      phi(i,jg,k,1) = r_gh
                      phi(i,jg,k,2) = r_gh*u_gh 
                      phi(i,jg,k,3) = r_gh*v_gh 
                      phi(i,jg,k,4) = r_gh*w_gh 
                      phi(i,jg,k,5) = p_gh/(gamma0-1.0_rp) &
                                   + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)

                   enddo
               enddo
          enddo
          !$acc end parallel

          case('B','F')

            !$acc parallel default(present)
            !$acc loop gang, vector collapse(3)
            do        k = 1,GN
               do     j = lby,uby
                   do i = lbx,ubx

                      kg = node + norm*k      !< ghost node
                      ki = node - norm*(k-1)  !< inner node

                      ! compute primitive variable in the internal nodes
                      r_in = phi(i,j,ki,1)
                      irin = 1._rp/r_in
                      u_in = phi(i,j,ki,2) * irin
                      v_in = phi(i,j,ki,3) * irin
                      w_in = phi(i,j,ki,4) * irin
                      p_in = (gamma0-1._rp)*(phi(i,j,ki,5)&
                             -0.5_rp*r_in*(u_in*u_in + v_in*v_in + w_in*w_in))
                      T_in = p_in*irin

                      ! compute primitives variables in the ghosts nodes
                      u_gh = 2*vel   - u_in
                      v_gh =         - v_in
                      w_gh =         - w_in
                      T_gh =           T_in  
                      p_gh = p_in
                      r_gh = p_gh/T_gh

                      phi(i,j,kg,1) = r_gh
                      phi(i,j,kg,2) = r_gh*u_gh 
                      phi(i,j,kg,3) = r_gh*v_gh 
                      phi(i,j,kg,4) = r_gh*w_gh 
                      phi(i,j,kg,5) = p_gh/(gamma0-1.0_rp) &
                                   + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)

                   enddo
               enddo
            enddo
          !$acc end parallel


        end select
        return
end subroutine adiabatic_wall

subroutine slip_wall(phi,b)
        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: phi
        type(face_type)                          , intent(in)    :: b

        integer :: i,j,k
        integer :: ig,jg,kg,ii,ji,ki
        integer :: node, norm

        node = b%node
        norm = b%norm

        select case(b%face)

          case('E','W')
                  
            !$acc parallel default(present)
            !$acc loop gang, vector collapse(3)
            do        k = lbz,ubz
               do     j = lby,uby
                   do i = 1,GN

                      ig = node + norm*i
                      ii = node - norm*(i-1)

                      phi(ig,j,k,1) = + phi(ii,j,k,1) 
                      phi(ig,j,k,2) = - phi(ii,j,k,2) 
                      phi(ig,j,k,3) = + phi(ii,j,k,3) 
                      phi(ig,j,k,4) = + phi(ii,j,k,4) 
                      phi(ig,j,k,5) = + phi(ii,j,k,5) 

                   enddo
               enddo
            end do
            !$acc end parallel

          case('S','N')

            !$acc parallel default(present)
            !$acc loop gang, vector collapse(3)
            do        k = lbz,ubz
               do     j = 1,GN
                   do i = lbx,ubx

                      jg = node + norm*j      !< ghost node
                      ji = node - norm*(j-1)  !< inner node

                      phi(i,jg,k,1) = + phi(i,ji,k,1) 
                      phi(i,jg,k,2) = + phi(i,ji,k,2) 
                      phi(i,jg,k,3) = - phi(i,ji,k,3) 
                      phi(i,jg,k,4) = + phi(i,ji,k,4) 
                      phi(i,jg,k,5) = + phi(i,ji,k,5) 

                   enddo
               enddo
            enddo
            !$acc end parallel

          case('B','F')

            !$acc parallel default(present)
            !$acc loop gang, vector collapse(3)
            do        k = 1,GN
               do     j = lby,uby
                   do i = lbx,ubx

                      kg = node + norm*k      !< ghost node
                      ki = node - norm*(k-1)  !< inner node

                      phi(i,j,kg,1) = + phi(i,j,ki,1) 
                      phi(i,j,kg,2) = + phi(i,j,ki,2) 
                      phi(i,j,kg,3) = + phi(i,j,ki,3) 
                      phi(i,j,kg,4) = - phi(i,j,ki,4) 
                      phi(i,j,kg,5) = + phi(i,j,ki,5) 

                   enddo
               enddo
            enddo
            !$acc end parallel

        end select
        return
end subroutine slip_wall


subroutine isothermal_wall(phi,b,Twall)
        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: phi
        type(face_type)                          , intent(in)    :: b
        real(rp)                                 , intent(in)    :: tWall

        ! local declarations
        real(rp) :: vel
        real(rp) :: r_in, u_in, v_in, w_in, p_in, T_in, irin
        real(rp) :: r_gh, u_gh, v_gh, w_gh, p_gh, T_gh
        integer  :: i,j,k,ig,ii,jg,ji,kg,ki
        integer  :: node, norm

        call set_wall_speed(vel,b)

        node = b%node
        norm = b%norm

        select case(b%face)

          case('E','W')

            !$acc parallel default(present)
            !$acc loop gang, vector collapse(3)
            do k       = lbz,ubz
               do j    = lby,uby
                  do i = 1,GN

                     ig = node + norm*i      !< ghost node
                     ii = node - norm*(i-1)  !< inner node
                        
                     ! compute primitive variable in the internal nodes
                     r_in = phi(ii,j,k,1)
                     irin = 1._rp/r_in
                     u_in = phi(ii,j,k,2) * irin
                     v_in = phi(ii,j,k,3) * irin
                     w_in = phi(ii,j,k,4) * irin
                     p_in = (gamma0-1._rp)*(phi(ii,j,k,5)&
                           -0.5_rp*r_in*(u_in*u_in + v_in*v_in + w_in*w_in))
                     T_in = p_in*irin

                     ! compute primitives variables in the ghosts nodes
                     u_gh =         - u_in
                     v_gh = 2*vel   - v_in
                     w_gh =         - w_in
                     T_gh = 2*Twall - T_in  
                     p_gh = p_in
                     r_gh = p_gh/T_gh

                     phi(ig,j,k,1) = r_gh
                     phi(ig,j,k,2) = r_gh*u_gh 
                     phi(ig,j,k,3) = r_gh*v_gh 
                     phi(ig,j,k,4) = r_gh*w_gh 
                     phi(ig,j,k,5) = p_gh/(gamma0-1.0_rp) &
                                   + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)

                  enddo
               enddo
            enddo
            !$acc end parallel

          case('S','N')
        
            !$acc parallel default(present)
            !$acc loop gang, vector collapse(3)
            do k       = lbz,ubz
               do j    = 1,GN
                  do i = lbx,ubx

                     jg = node + norm*j      !< ghost node
                     ji = node - norm*(j-1)  !< inner node
                        
                     ! compute primitive variable in the internal nodes
                     r_in = phi(i,ji,k,1)
                     irin = 1._rp/r_in
                     u_in = phi(i,ji,k,2) * irin
                     v_in = phi(i,ji,k,3) * irin
                     w_in = phi(i,ji,k,4) * irin
                     p_in = (gamma0-1._rp)*(phi(i,ji,k,5)&
                            -0.5_rp*r_in*(u_in*u_in + v_in*v_in + w_in*w_in))
                     T_in = p_in*irin

                     ! compute primitives variables in the ghosts nodes
                     u_gh = 2*vel   - u_in
                     v_gh =         - v_in
                     w_gh =         - w_in
                     T_gh = 2*Twall - T_in  
                     p_gh = p_in
                     r_gh = p_gh/T_gh

                     phi(i,jg,k,1) = r_gh
                     phi(i,jg,k,2) = r_gh*u_gh 
                     phi(i,jg,k,3) = r_gh*v_gh 
                     phi(i,jg,k,4) = r_gh*w_gh 
                     phi(i,jg,k,5) = p_gh/(gamma0-1.0_rp) &
                                   + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)

                  enddo
               enddo
            enddo
            !$acc end parallel

          case('B','F')

            !$acc parallel default(present)
            !$acc loop gang, vector collapse(3)
            do k       = 1,GN
               do j    = lby,uby
                  do i = lbx,ubx

                     kg = node + norm*k      !< ghost node
                     ki = node - norm*(k-1)  !< inner node

                     ! compute primitive variable in the internal nodes
                     r_in = phi(i,j,ki,1)
                     irin = 1._rp/r_in
                     u_in = phi(i,j,ki,2) * irin
                     v_in = phi(i,j,ki,3) * irin
                     w_in = phi(i,j,ki,4) * irin
                     p_in = (gamma0-1._rp)*(phi(i,j,ki,5)&
                            -0.5_rp*r_in*(u_in*u_in + v_in*v_in + w_in*w_in))
                     T_in = p_in*irin

                     ! compute primitives variables in the ghosts nodes
                     u_gh = 2*vel   - u_in
                     v_gh =         - v_in
                     w_gh =         - w_in
                     T_gh = 2*Twall - T_in  
                     p_gh = p_in
                     r_gh = p_gh/T_gh

                     phi(i,j,kg,1) = r_gh
                     phi(i,j,kg,2) = r_gh*u_gh 
                     phi(i,j,kg,3) = r_gh*v_gh 
                     phi(i,j,kg,4) = r_gh*w_gh 
                     phi(i,j,kg,5) = p_gh/(gamma0-1.0_rp) &
                                   + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)

                  enddo
               enddo
            enddo
            !$acc end parallel

        end select
        return
end subroutine isothermal_wall



subroutine postshock_slipwall(phi,b)
! ---------------------------------------------------------------------
!
!       This subroutine applies a postshock/slipwall boundary condition.
!       The condition is used just for the double mach reflection problem.
!
! ---------------------------------------------------------------------
        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi
        type(face_type)                          , intent(in)    :: b

        real(rp), parameter :: x0 = 1.0_rp/6.0_rp, theta = pi/6.0_rp
        real(rp)            :: r1, p1, T1, M1
        real(rp)            :: r2, p2, T2, u2, v2, w2, M2
        real(rp)            :: u_in, v_in, w_in, p_in, T_in, r_in, irin
        real(rp)            :: u_gh, v_gh, w_gh, p_gh, T_gh, r_gh
        real(rp)            :: vel1, vel2
        integer             :: jg, ji, b_node, b_norm

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
        w2   = 0.0_rp

        b_node = b%node
        b_norm = b%norm
        
        if(b%face == 'S') then
               
                !$acc parallel default(present)
                !$acc loop collapse(2)
                do k       = lbz,ubz
                   do i    = lbx,ubx
                      do j = 1,GN

                         jg = b_node + b_norm*j      !< ghost node
                         ji = b_node - b_norm*(j-1)  !< inner node

                         if(x(i) <= x0) then
                                 ! post-shoked

                                 phi(i,jg,k,1) = r2
                                 phi(i,jg,k,2) = r2 * u2
                                 phi(i,jg,k,3) = r2 * v2
                                 phi(i,jg,k,4) = r2 * w2
                                 phi(i,jg,k,5) = (p2)/(gamma0-1.0_rp)+0.5_rp*r2*(u2*u2+v2*v2+w2*w2)

                         else
                                 ! slip-wall
                                 r_in = phi(i,ji,k,1)
                                 irin = 1.0_rp/r_in
                                 u_in = phi(i,ji,k,2)*irin
                                 v_in = phi(i,ji,k,3)*irin
                                 w_in = phi(i,ji,k,4)*irin
                                 p_in = (gamma0-1._rp)*(phi(i,ji,k,5)&
                                     -0.5_rp*r_in*(u_in*u_in + v_in*v_in + w_in*w_in))
                                 T_in = p_in*irin

                                 u_gh =   u_in
                                 v_gh = - v_in
                                 w_gh = - w_in
                                 T_gh =   T_in
                                 p_gh =   p_in
                                 r_gh = p_gh/T_gh
                                 
                                 phi(i,ji,k,1) = r_gh
                                 phi(i,ji,k,2) = r_gh*u_gh 
                                 phi(i,ji,k,3) = r_gh*v_gh 
                                 phi(i,ji,k,4) = r_gh*w_gh 
                                 phi(i,ji,k,5) = p_gh/(gamma0-1.0_rp) &
                                   + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)


                         endif
                      enddo
                   enddo
                enddo
                !$acc end parallel

        else
                print*, 'Postshock/slipwall boundary non implemented for face: ', b%face
                stop
        endif
        return
end subroutine postshock_slipwall


subroutine travelling_shock_wave(phi,b)
! ---------------------------------------------------------------------
!
!       This subroutine applies a travelling shock boundary condition.
!       The condition is used just for the double mach reflection problem.
!
! ---------------------------------------------------------------------
        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi
        type(face_type)                          , intent(in)    :: b

        real(rp), parameter :: x0 = 1.0_rp/6.0_rp, theta = pi/6.0_rp
        real(rp)            :: r1, p1, T1, u1, v1, w1, M1
        real(rp)            :: r2, p2, T2, u2, v2, w2, M2
        real(rp)            :: vel1, vel2, tan_theta, xi, yjg
        integer             :: jg, b_norm, b_node

        if(b%face /= 'N') then
          if(rank == root) then
            print*, ' Travelling shock wave boundary non implemented for face: ', b%face
            stop
          endif
        endif

        tan_theta = tan(theta)

        !
        ! ==== state 1
        !
        M1   = Mach
        p1   = 1.0_rp
        r1   = 1.4_rp
        T1   = p1/r1
        vel1 = - M1 * sqrt(gamma0*T1)
        u1   = 0.0_rp
        v1   = 0.0_rp
        w1   = 0.0_rp

        !
        ! ==== state 2
        !
        call Rankine_Hugoniot(r1,p1,T1,M1,r2,p2,T2,M2)

        vel2 = - M2 * sqrt(gamma0*T2)
        u2   =   (vel2 - vel1)*cos(theta)
        v2   = - (vel2 - vel1)*sin(theta)
        w2   = 0.0_rp
        
        b_node = b%node
        b_norm = b%norm

        !$acc parallel default(present)
        !$acc loop collapse(2)
        do       k = lbz,ubz
           do    i = lbx,ubx
              do j = 1,GN

                 jg = b_node + b_norm*j      !< ghost node

                 xi  = x(i)
                 yjg = y(jg)

                 if    (xi + time*vel1*sqrt(gamma0) <= x0 + yjg*tan_theta) then 

                   phi(i,jg,k,1) = r2
                   phi(i,jg,k,2) = r2 * u2
                   phi(i,jg,k,3) = r2 * v2
                   phi(i,jg,k,4) = r2 * w2
                   phi(i,jg,k,5) = (p2)/(gamma0-1.0_rp)+0.5_rp*r2*(u2*u2+v2*v2+w2*w2)

                 elseif(xi + time*vel1*sqrt(gamma0) >= x0 + yjg*tan_theta) then

                   phi(i,jg,k,1) = r1
                   phi(i,jg,k,2) = r1 * u1
                   phi(i,jg,k,3) = r1 * v1
                   phi(i,jg,k,4) = r1 * w1
                   phi(i,jg,k,5) = (p1)/(gamma0-1.0_rp)+0.5_rp*r1*(u1*u1+v1*v1+w1*w1)

                 endif

              enddo
           enddo
        enddo
        !$acc end parallel
        
        return
end subroutine travelling_shock_wave






















subroutine outflow_nscbc_static(phi,b)
        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi
        type(face_type)                          , intent(in) :: b

        integer  :: i,j,k,l
        real(rp) :: r_1, ir1,u_1, v_1, w_1, ek1, p_1, T_1, cc1, c_1
        real(rp) :: r_2, ir2,u_2, v_2, w_2, ek2, p_2, T_2, cc2, c_2
        real(rp) :: d1k, d2k, d3k, d1kh, d2kh, d3kh

        real(rp) :: rhoh, rhohc, rhohci
        real(rp) :: r_, u_, v_, w_, p_

        real(rp), dimension(5) :: deltaw, deltav
        real(rp), dimension(5,5) :: LL, RR


        if(b%face == 'E') then
          do    k = sz,ez
             do j = sy,ey

                i = ex
                r_2 = phi(i,j,k,1)
                ir2 = 1.0_rp/r_2
                u_2 = phi(i,j,k,2)*ir2
                v_2 = phi(i,j,k,3)*ir2
                w_2 = phi(i,j,k,4)*ir2
                ek2 = 0.5_rp*(u_2*u_2 + v_2*v_2 + w_2*w_2)
                p_2 = (gamma0-1.0_rp)*(phi(i,j,k,5) - r_2*ek2)
                T_2 = p_2*ir2
                cc2 = gamma0*T_2
                c_2 = sqrt(cc2)

                i = ex-1
                r_1 = phi(i,j,k,1)
                ir1 = 1.0_rp/r_1
                u_1 = phi(i,j,k,2)*ir1
                v_1 = phi(i,j,k,3)*ir1
                w_1 = phi(i,j,k,4)*ir1
                ek1 = 0.5_rp*(u_1*u_1 + v_1*v_1 + w_1*w_1)
                p_1 = (gamma0-1.0_rp)*(phi(i,j,k,5) - r_1*ek1)
                T_1 = p_1*ir1
                cc1 = gamma0*T_1
                c_1 = sqrt(cc1)

                deltav(1) = r_2 - r_1
                deltav(2) = u_2 - u_1
                deltav(3) = v_2 - v_1
                deltav(4) = w_2 - w_1
                deltav(5) = p_2 - p_1

                ! left on nx-1
                d1k = 1.0_rp
                d2k = 0.0_rp
                d3k = 0.0_rp 
                LL(1,:) = [d1k, 0.0_rp, d3k, -d2k,-d1k/cc2]
                LL(2,:) = [d2k,-d3k, 0.0_rp,  d1k,-d2k/cc2]
                LL(3,:) = [d3k, d2k,-d1k,  0.0_rp,-d3k/cc2]
                LL(4,:) = [0.0_rp, d1k, d2k,  d3k, ir2/c_2]
                LL(5,:) = [0.0_rp,-d1k,-d2k, -d3k, ir2/c_2]

                deltaw = matmul(LL,deltav)
        
                ! right on nx
                d1kh = 0.5_rp*d1k
                d2kh = 0.5_rp*d2k
                d3kh = 0.5_rp*d3k
                
                rhoh = 0.5_rp*r_2
                rhohc = rhoh*c_2
                rhohci = rhoh/c_2

                RR(1,:) = [d1k , d2k , d3k , rhohci, rhohci]
                RR(2,:) = [0.0_rp , -d3k, d2k , d1kh  , -d1kh]
                RR(3,:) = [d3k , 0.0_rp , -d1k, d2kh  , -d2kh] 
                RR(4,:) = [-d2k, d1k , 0.0_rp , d3kh  , -d3kh]
                RR(5,:) = [0.0_rp , 0.0_rp , 0.0_rp , rhohc , rhohc]

                deltaw(5) = 0.0_rp
                if(u_2 < 0.0_rp) then
                  deltaw(1) = 0.0_rp
                  deltaw(2) = 0.0_rp
                  deltaw(3) = 0.0_rp
                endif

                deltav = matmul(RR, deltaw)

                r_ = r_2 + deltav(1)
                u_ = u_2 + deltav(2)
                v_ = v_2 + deltav(3)
                w_ = w_2 + deltav(4)
                p_ = p_2 + deltav(5)
        
                do l = 1,GN
                   phi(ex+l,j,k,1) = r_
                   phi(ex+l,j,k,2) = r_*u_
                   phi(ex+l,j,k,3) = r_*v_
                   phi(ex+l,j,k,4) = r_*w_
                   phi(ex+l,j,k,5) = p_/(gamma0-1.0_rp) + 0.5_rp*r_*(u_*u_+v_*v_+w_*w_)
                enddo


             enddo
          enddo

        elseif(b%face == 'N') then
          do k    = sz,ez
             do i = sx,ex

                j = ey
                r_2 = phi(i,j,k,1)
                ir2 = 1.0_rp/r_2
                u_2 = phi(i,j,k,2)*ir2
                v_2 = phi(i,j,k,3)*ir2
                w_2 = phi(i,j,k,4)*ir2
                ek2 = 0.5_rp*(u_2*u_2 + v_2*v_2 + w_2*w_2)
                p_2 = (gamma0-1.0_rp)*(phi(i,j,k,5) - r_2*ek2)
                T_2 = p_2*ir2
                cc2 = gamma0*T_2
                c_2 = sqrt(cc2)

                j = ey-1
                r_1 = phi(i,j,k,1)
                ir1 = 1.0_rp/r_1
                u_1 = phi(i,j,k,2)*ir1
                v_1 = phi(i,j,k,3)*ir1
                w_1 = phi(i,j,k,4)*ir1
                ek1 = 0.5_rp*(u_1*u_1 + v_1*v_1 + w_1*w_1)
                p_1 = (gamma0-1.0_rp)*(phi(i,j,k,5) - r_1*ek1)
                T_1 = p_1*ir1
                cc1 = gamma0*T_1
                c_1 = sqrt(cc1)

                deltav(1) = r_2 - r_1
                deltav(2) = u_2 - u_1
                deltav(3) = v_2 - v_1
                deltav(4) = w_2 - w_1
                deltav(5) = p_2 - p_1

                ! left on ny
                d1k = 1.0_rp
                d2k = 0.0_rp
                d3k = 0.0_rp 
                LL(1,:) = [d1k, 0.0_rp, d3k, -d2k,-d1k/cc2]
                LL(2,:) = [d2k,-d3k, 0.0_rp,  d1k,-d2k/cc2]
                LL(3,:) = [d3k, d2k,-d1k,  0.0_rp,-d3k/cc2]
                LL(4,:) = [0.0_rp, d1k, d2k,  d3k, ir2/c_2]
                LL(5,:) = [0.0_rp,-d1k,-d2k, -d3k, ir2/c_2]

                deltaw = matmul(LL,deltav)
        
                ! right on ny
                d1kh = 0.5_rp*d1k
                d2kh = 0.5_rp*d2k
                d3kh = 0.5_rp*d3k
                
                rhoh = 0.5_rp*r_2
                rhohc = rhoh*c_2
                rhohci = rhoh/c_2

                RR(1,:) = [d1k , d2k , d3k , rhohci, rhohci]
                RR(2,:) = [0.0_rp , -d3k, d2k , d1kh  , -d1kh]
                RR(3,:) = [d3k , 0.0_rp , -d1k, d2kh  , -d2kh] 
                RR(4,:) = [-d2k, d1k , 0.0_rp , d3kh  , -d3kh]
                RR(5,:) = [0.0_rp , 0.0_rp , 0.0_rp , rhohc , rhohc]

                deltaw(5) = 0.0_rp
                if(v_2 < 0.0_rp) then
                  deltaw(1) = 0.0_rp
                  deltaw(2) = 0.0_rp
                  deltaw(3) = 0.0_rp
                endif

                deltav = matmul(RR, deltaw)

                r_ = r_2 + deltav(1)
                u_ = u_2 + deltav(2)
                v_ = v_2 + deltav(3)
                w_ = w_2 + deltav(4)
                p_ = p_2 + deltav(5)
        
                do l = 1,GN
                   phi(i,ey+l,k,1) = r_
                   phi(i,ey+l,k,2) = r_*u_
                   phi(i,ey+l,k,3) = r_*v_
                   phi(i,ey+l,k,4) = r_*w_
                   phi(i,ey+l,k,5) = p_/(gamma0-1.0_rp) + 0.5_rp*r_*(u_*u_+v_*v_+w_*w_)
                enddo

             enddo
          enddo




        else
          write(*,*) ' Outflow boundary is not implemented for face', b%face
          
        endif


        return
end subroutine outflow_nscbc_static





subroutine nscbc_inflow(phi,iFlow,b)

        use math_tools_module, only: smooth_step
        use inflow_module    , only: turbulent_inflow

        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi
        type(ifl)                                , intent(inout) :: Iflow
        type(face_type)                          , intent(in)    :: b

        ! local declarations
        integer                  :: i, j,k, is, s, j0, j1
        real(rp), dimension(5)   :: q_is,L
        real(rp)                 :: rho_, u_, v_, w_, T_, c, step
        real(rp)                 :: i_dx, p_x, u_x, aik, bik, fws, irho_is, p_is
        real(rp)                 :: istep, uf, vf, wf

        ! check face
        if(b%face /= 'W') then
          if(rank == root) print*, 'Inflow is not implemented for face ', b%face
          call secure_stop
        endif

        step = 1.0_rp
        if(iFlow%SmthFlag) step = smooth_step(0.0_rp, 0.3_rp*real(itmax,rp), real(it,rp))

        if(iFlow%TurbFlag) call turbulent_inflow(iFlow)

        aik = a_rk(ik)
        bik = b_rk(ik)
        
        j0 = sy-1
        j1 = sy-1
        if(trim(inflow_profile) == 'smooth_body_inflow') then
           istep = 0.22_rp/0.032_rp
           do j = sy,ey
              if(y(j) > 0) then
                 j0 = j
                 exit
              endif
           enddo
           do j = sy,ey
              if(y(j) > istep) then
                 j1 = j 
                 exit
              endif
           enddo
        endif

        !open(unit = 10, file = "DATA/"//trim(data_dir)//'/tmp1._rptxt')
        !open(unit = 11, file = "DATA/"//trim(data_dir)//'/tmp2._rptxt')
        !do j = sy,ey
        !   write(10,*) y(j), iflow%turb(1,j,1)
        !   if(j > j1) then
        !   write(11,*) y(j), iflow%turb(1,j-j1+j0,1)
        !   else
        !   write(11,*) y(j), 0.0_rp
        !   endif
        !enddo
        !close(10)
        !stop
      
        do       k = sz,ez
           do    j = sy,ey
              do i = lbx,sx-1
                
              uf = 0.0_rp
              vf = 0.0_rp
              wf = 0.0_rp
              if(j > j1) then
                uf = iFlow%turb(1,j-j1+j0,k)
                vf = iFlow%turb(2,j-j1+j0,k)
                wf = iFlow%turb(3,j-j1+j0,k)
              endif

              rho_ = phi(i,j,k,1)
              u_   = step*(iFlow%mean(i,j,2) + uf) 
              v_   =      (iFlow%mean(i,j,3) + vf)
              w_   =      (iFlow%mean(i,j,4) + wf)
              T_   =       iFlow%mean(i,j,5)
              c    = sqrt(gamma0 * T_)

              ! downwinding FD along x-direction
              i_dx = xstep_i(i)
              
              p_x = 0.0_rp
              u_x = 0.0_rp
              do s = 0,fward_fd_order

                 is  = i+s
                 fws = fward_1(s)

                 q_is(:) = phi(is,j,k,:)
                 irho_is = 1.0_rp/q_is(1)
                 p_is    = cv_i*(q_is(5) - 0.5_rp*irho_is*(q_is(2)*q_is(2) + q_is(3)*q_is(3) + q_is(4)*q_is(4)))

                 p_x = p_x + fws * p_is
                 u_x = u_x + fws * q_is(2)*irho_is

              enddo
              p_x = p_x * i_dx
              u_x = u_x * i_dx
              
              ! charactherist waves 
              L(1) = (u_-c) * (p_x - rho_ * c * u_x)
              L(5) = L(1)
              L(2) = 0.5_rp * (gamma0 - 1.0_rp) * (L(5) + L(1))
              L(3) = 0.0_rp
              L(4) = 0.0_rp

              RHS(i,j,k,1) = - 1.0_rp/c**2 *(L(2) + 0.5_rp*(L(5) + L(1)))

              ! Evolving RHO at the bound with RK
              phi(i,j,k,1) = aik * phi_n(i,j,k,1) + bik * (phi(i,j,k,1) + dt * RHS(i,j,k,1))
                    
              ! update conservative variables on the ghost
              rho_ = phi(i,j,k,1)
              phi(i,j,k,2) = rho_ * u_
              phi(i,j,k,3) = rho_ * v_
              phi(i,j,k,4) = rho_ * w_
              phi(i,j,k,5) =(rho_ * T_)/(gamma0-1.0_rp)+0.5_rp*rho_*(u_**2+v_**2+w_**2)
              enddo
           enddo
        enddo

        return
end subroutine nscbc_inflow

subroutine nscbc_inflow_relaxed(phi,iFlow,b,M_max)
        
        use math_tools_module, only: smooth_step
        use inflow_module    , only: turbulent_inflow

        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi
        type(ifl)                                , intent(inout) :: Iflow
        type(face_type)                          , intent(in)    :: b
        real(rp)                                 , intent(in)    :: M_max

        ! local declarations
        integer                  :: i, j,k, is, s, j0, j1
        real(rp), dimension(5)   :: q_is,L, eta, d
        real(rp)                 :: rho_, irho, u_, v_, w_, T_, c, p_, ek, step
        real(rp)                 :: u0, v0, w0, T0, istep, uf, vf, wf
        real(rp)                 :: i_dx, p_x, u_x, aik, bik, fws, irho_is, p_is

        ! check face
        if(b%face /= 'W') then
          if(rank == root) print*, 'Inflow is not implemented for face ', b%face
          call secure_stop
        endif

        step = 1.0_rp
        if(iFlow%SmthFlag) step = smooth_step(0.0_rp, 0.3_rp*real(itmax,rp), real(it,rp))

        if(iFlow%TurbFlag) call turbulent_inflow(iFlow)

        aik = a_rk(ik)
        bik = b_rk(ik)

        !eta(:) = 0.5_rp
        eta(:) = 1.0_rp/(0.25_rp*0.15_rp/u_inf) ! 1/4 DF correlation time

        j0 = sy-1
        j1 = sy-1
        if(trim(inflow_profile) == 'smooth_body_inflow') then
           istep = 0.22_rp/0.032_rp
           do j = sy,ey
              if(y(j) > 0) then
                 j0 = j
                 exit
              endif
           enddo
           do j = sy,ey
              if(y(j) > istep) then
                 j1 = j 
                 exit
              endif
           enddo
        endif

        do k       = sz,ez
           do j    = sy,ey
              do i = lbx,sx-1

                 uf = 0.0_rp
                 vf = 0.0_rp
                 wf = 0.0_rp
                 if(j > j1) then
                   uf = iFlow%turb(1,j-j1+j0,k)
                   vf = iFlow%turb(2,j-j1+j0,k)
                   wf = iFlow%turb(3,j-j1+j0,k)
                 endif

                 u0 = step*(iFlow%mean(i,j,2) + uf) 
                 v0 =      (iFlow%mean(i,j,3) + vf)
                 w0 =      (iFlow%mean(i,j,4) + wf)
                 T0 =       iFlow%mean(i,j,5)

                 !
                 ! ==== access memory and save scalars
                 ! 
                 rho_ = phi(i,j,k,1)
                 irho = 1.0_rp/rho_
                 u_   = phi(i,j,k,2)*irho
                 v_   = phi(i,j,k,3)*irho
                 w_   = phi(i,j,k,4)*irho
                 p_   = (gamma0-1.0_rp)*(phi(i,j,k,5) - 0.5_rp*rho_*(u_*u_ + v_*v_ + w_*w_))
                 T_   = p_*irho
                 c    = sqrt(gamma0 * T_)
                 ! 
                 ! ==== compute downwinding finite differences
                 !
                 i_dx = xstep_i(i)
                 
                 p_x = 0.0_rp
                 u_x = 0.0_rp
                 do s = 0,fward_fd_order

                    is  = i+s
                    fws = fward_1(s)

                    q_is(:) = phi(is,j,k,:)
                    irho_is = 1.0_rp/q_is(1)
                    p_is    = cv_i*(q_is(5) - 0.5_rp*irho_is*(q_is(2)*q_is(2) + q_is(3)*q_is(3) + q_is(4)*q_is(4)))

                    p_x = p_x + fws * p_is
                    u_x = u_x + fws * q_is(2)*irho_is

                 enddo
                 p_x = p_x * i_dx
                 u_x = u_x * i_dx
                 !
                 ! ==== compute waves
                 !
                 L(1) = (u_-c) * (p_x - rho_ * c * u_x)
                 L(2) = eta(2) *                     rho_*c/Lx * (T0-T_)
                 L(3) = eta(3) *                          c/Lx * (v_-v0)
                 L(4) = eta(4) *                          c/Lx * (w_-w0)
                 L(5) = eta(5) *rho_*c**2*(1.0_rp-M_max*M_max)/Lx * (u_-u0)

                 ek   = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)

                 d(1) = 1._rp/c**2 *(L(2)+0.5_rp*(L(5)+L(1)))
                 d(2) = 0.5_rp*(L(5)+L(1))
                 d(3) = 1._rp/(2*rho_*c)*(L(5)-L(1))
                 d(4) = L(3)
                 d(5) = L(4)
                 !
                 ! ==== evolving phi via RK
                 !
                 RHS(i,j,k,1) =    d(1)
                 RHS(i,j,k,2) = u_*d(1)+rho_*d(3)
                 RHS(i,j,k,3) = v_*d(1)+rho_*d(4)
                 RHS(i,j,k,4) = w_*d(1)+rho_*d(5)
                 RHS(i,j,k,5) = ek*d(1) + d(2)/(gamma0-1.0_rp) + rho_*u_*d(3) + rho_*v_*d(4) + rho_*w_*d(5)

                 phi(i,j,k,:) = aik * phi_n(i,j,k,:) + bik * (phi(i,j,k,:) - dt * RHS(i,j,k,:))

              enddo
           enddo
        enddo

        return
end subroutine nscbc_inflow_relaxed

subroutine nscbc_outflow(phi,b,L1_wave,M_max)
! -------------------------------------------------------------------------
!
!       Characteristic outflow boundary condition (Poinsot-Lele)
!
!       INPUT:  phi     !< conservative variables
!               b       !< bound type
!               M_max   !< Maximum Mach number in the whole domain
!               L1_wave !< type of the L1 characteristic wave
!
!               L1_wave = poinsot_model   (standard model)
!               L1_wave = pirozzoli_model (optimised model)
!               L1_wave = zero_reflection (sometimes could be good)
!
! -------------------------------------------------------------------------
        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi
        type(face_type)                          , intent(in)    :: b
        real(rp)                                 , intent(in)    :: M_max
        character(*)                             , intent(in)    :: L1_wave

        ! local declarations
        real(rp), allocatable, dimension(:) :: cd_coef 
        real(rp), parameter    :: p_g = 1.0_rp
        real(rp), dimension(5) :: L, lambda, d
        integer                :: s, is, js, fdL,fdR
        integer                :: ii, ji, ig, jg
        real(rp)               :: rho_, irho, u_, v_, w_, p_, T_, c
        real(rp)               :: r_is, iris, u_is, v_is, w_is, p_is
        real(rp)               :: r_x, p_x, u_x, v_x, w_x
        real(rp)               :: r_y, p_y, u_y, v_y, w_y
        real(rp)               :: i_dx, i_dy, bw_s, aik, bik, ek

!        real(rp) :: i_dh, i2dh, ix_xi, d2xi_dx2, cs_1, cs_2
!        real(rp), dimension(3) :: vj, vel_y, velyy, div_sigma
        real(rp), dimension(5) :: rhs_vis
!        real(rp)               :: irjs, vis_y, mu_

        aik = a_rk(ik)
        bik = b_rk(ik)
        
        !
        ! === set the finite difference coefficients
        !
        selectcase(b%face)
          case('E','N','F')
            fdL = - bward_fd_order
            fdR =   0
            allocate(cd_coef(fdL:fdR))
            cd_coef = bward_1(:)

          case('W','S','B')
            fdL =   0 
            fdR =   fward_fd_order
            allocate(cd_coef(fdL:fdR))
            cd_coef = fward_1(:)
          case default
            print*, ' nscbc_outflow is not implemented for face ', b%face
            stop
        
        endselect

        selectcase(b%face)

          !
          ! ==== East and West bounds
          !
          case('E','W')
          do k = sz,ez
             do j = sy,ey
                do i = 1,1

                   ig = b%node + b%norm*i      !< ghost node
                   ii = b%node - b%norm*(i-1)  !< inner node
                  
                   ! scalarize some variables
                   rho_  = phi(ig,j,k,1)
                   irho  = 1.0_rp/rho_
                   u_    = phi(ig,j,k,2)*irho
                   v_    = phi(ig,j,k,3)*irho 
                   w_    = phi(ig,j,k,4)*irho 
                   ek    = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
                   p_    = cv_i * (phi(ig,j,k,5) - rho_*ek)
                   T_    = p_*irho
                   c     = sqrt(gamma0 * T_)

                   ! computing derivatives
                   i_dx = xstep_i(ig)

                   r_x = 0.0_rp
                   p_x = 0.0_rp
                   u_x = 0.0_rp
                   v_x = 0.0_rp
                   w_x = 0.0_rp
                   do s = fdL, fdR

                      is = ig+s

                      bw_s = cd_coef(s)

                      r_is = phi(is,j,k,1)
                      iris = 1.0_rp/r_is
                      u_is = phi(is,j,k,2)*iris
                      v_is = phi(is,j,k,3)*iris
                      w_is = phi(is,j,k,4)*iris
                      p_is = cv_i*(phi(is,j,k,5) - 0.5_rp*r_is*(u_is*u_is + v_is*v_is + w_is*w_is))

                      r_x = r_x + bw_s * r_is
                      p_x = p_x + bw_s * p_is
                      u_x = u_x + bw_s * u_is
                      v_x = v_x + bw_s * v_is
                      w_x = w_x + bw_s * w_is

                   enddo
                   r_x = r_x * i_dx
                   p_x = p_x * i_dx
                   u_x = u_x * i_dx
                   v_x = v_x * i_dx
                   w_x = w_x * i_dx
                          
                   ! Euler eigenvalues
                   lambda(1) = 0.5_rp * (u_ - c + b%norm * abs(u_ - c))
                   lambda(2) = 0.5_rp * (u_     + b%norm * abs(u_    ))
                   lambda(3) = 0.5_rp * (u_     + b%norm * abs(u_    ))
                   lambda(4) = 0.5_rp * (u_     + b%norm * abs(u_    ))
                   lambda(5) = 0.5_rp * (u_ + c + b%norm * abs(u_ + c))
          
                   ! characteristics
                   L(1) = 0.0_rp
                   if(abs(u_)/c > 1.0_rp) then
                     L(1) = lambda(1)*(p_x - rho_*c*u_x)
                   else

                     if    (L1_wave == 'poinsot_model'  ) then 
                       L(1) = L1Poinsot(Lx,p_,p_g,c,M_max)
                     elseif(L1_wave == 'pirozzoli_model') then
                       L(1) = L1Pirozzoli(u_,c,p_,p_g,i_dx)
                     elseif(L1_wave == 'zero_reflection') then
                       L(1) = 0.0_rp
                     else
                       print*, 'L1_wave', trim(L1_wave), ' is not implemented.'
                       stop
                     endif
                   endif

                   L(2) = lambda(2) * ((c**2)*r_x - p_x)
                   L(3) = lambda(3) * (v_x)
                   L(4) = lambda(4) * (w_x)
                   L(5) = lambda(5) * (p_x + rho_*c*u_x)

                   d(1) = 1._rp/c**2 *(L(2)+0.5_rp*(L(5)+L(1)))
                   d(2) = 0.5_rp*(L(5)+L(1))
                   d(3) = 1._rp/(2*rho_*c)*(L(5)-L(1))
                   d(4) = L(3)
                   d(5) = L(4)
                
                   ! viscous correction
                   rhs_vis = 0.0_rp
!                   if(viscous) then
!                     i_dh = etastep_i(j)
!                     i2dh = i_dh*i_dh
!                     vel_y = 0.0_rp
!                     velyy = 0.0_rp
!                     vis_y = 0.0_rp
!                     do s = -central_fd_order/2,central_fd_order/2
!                        cs_1 = i_dh*central_1(s)
!                        cs_2 = i2dh*central_2(s)
!
!                        irjs  = 1.0_rp/phi(ig-1,j+s,k,1)
!                        vj(1) = phi(ig-1,j+s,k,2)*irjs
!                        vj(2) = phi(ig-1,j+s,k,3)*irjs
!                        vj(3) = phi(ig-1,j+s,k,4)*irjs
!
!                        vis_y = vis_y + cs_1 * VIS(ig-1,j+s,k)
!                        vel_y = vel_y + cs_1 * vj
!                        velyy = velyy + cs_2 * vj
!                     enddo
!                     ix_xi    = 1.0_rp/(y_eta(j))
!                     d2xi_dx2 = - y_eta2(j)*ix_xi**3
!
!                     velyy = velyy*(ix_xi)**2 + vel_y * d2xi_dx2
!                     vel_y = vel_y*ix_xi
!                     vis_y = vis_y*ix_xi
!                
!                     mu_ = VIS(ig-1,j,k)
!                     div_sigma(1) = mu_*velyy(1) + vel_y(1)*vis_y
!                     div_sigma(2) = mu_*velyy(2) + vel_y(2)*vis_y
!                     div_sigma(3) = mu_*velyy(3) + vel_y(3)*vis_y
!                     rhs_vis(2)  = - mu_inf * div_sigma(1)
!                     rhs_vis(3)  = - mu_inf * div_sigma(2)
!                     rhs_vis(4)  = - mu_inf * div_sigma(3)
!                     rhs_vis(5)  = - mu_inf * (div_sigma(1)*u_ + div_sigma(2)*v_ + div_sigma(3)*w_)
!                   endif

                   RHS(ig,j,k,1) =    d(1)
                   RHS(ig,j,k,2) = u_*d(1)+rho_*d(3) + rhs_vis(2)
                   RHS(ig,j,k,3) = v_*d(1)+rho_*d(4) + rhs_vis(3)
                   RHS(ig,j,k,4) = w_*d(1)+rho_*d(5) + rhs_vis(4)
                   RHS(ig,j,k,5) = ek*d(1) + d(2)/(gamma0-1.0_rp) + rho_*u_*d(3) + rho_*v_*d(4) + rho_*w_*d(5) + rhs_vis(5)

                   ! evolving phi via RK
                   phi(ig,j,k,:) = aik * phi_n(ig,j,k,:) + bik * (phi(ig,j,k,:) - dt * RHS(ig,j,k,:))

                enddo

                do i = 2, GN
                   ig = b%node + b%norm*i              !< ghost node
                   phi(ig,j,k,:) = phi(b%node+b%norm,j,k,:)
                enddo

             enddo
          enddo
        
          !
          ! ==== South and North bounds
          !
          case('S','N')
          do k = sz,ez
             do j = 1,1
                do i = sx,ex

                   jg = b%node + b%norm*j      !< ghost node
                   ji = b%node - b%norm*(j-1)  !< inner node
                  
                   ! scalarize some variables
                   rho_  = phi(i,jg,k,1)
                   irho  = 1.0_rp/rho_
                   u_    = phi(i,jg,k,2)*irho
                   v_    = phi(i,jg,k,3)*irho 
                   w_    = phi(i,jg,k,4)*irho 
                   ek    = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
                   p_    = cv_i * (phi(i,jg,k,5) - rho_*ek)
                   T_    = p_*irho
                   c     = sqrt(gamma0 * T_)

                   ! computing derivatives
                   i_dy = ystep_i(jg)

                   r_y = 0.0_rp
                   p_y = 0.0_rp
                   u_y = 0.0_rp
                   v_y = 0.0_rp
                   w_y = 0.0_rp
                   do s = fdL, fdR

                      js = jg+s

                      bw_s = cd_coef(s)

                      r_is = phi(i,js,k,1)
                      iris = 1.0_rp/r_is
                      u_is = phi(i,js,k,2)*iris
                      v_is = phi(i,js,k,3)*iris
                      w_is = phi(i,js,k,4)*iris
                      p_is = cv_i*(phi(i,js,k,5) - 0.5_rp*r_is*(u_is*u_is + v_is*v_is + w_is*w_is))

                      r_y = r_y + bw_s * r_is
                      p_y = p_y + bw_s * p_is
                      u_y = u_y + bw_s * u_is
                      v_y = v_y + bw_s * v_is
                      w_y = w_y + bw_s * w_is

                   enddo
                   r_y = r_y * i_dy
                   p_y = p_y * i_dy
                   u_y = u_y * i_dy
                   v_y = v_y * i_dy
                   w_y = w_y * i_dy
                          
                   ! Euler eigenvalues
                   lambda(1) = 0.5_rp * (v_ - c + b%norm * abs(v_ - c))
                   lambda(2) = 0.5_rp * (v_     + b%norm * abs(v_    ))
                   lambda(3) = 0.5_rp * (v_     + b%norm * abs(v_    ))
                   lambda(4) = 0.5_rp * (v_     + b%norm * abs(v_    ))
                   lambda(5) = 0.5_rp * (v_ + c + b%norm * abs(v_ + c))
          
                   ! characteristics
                   L(1) = 0.0_rp
                   if(abs(v_)/c > 1.0_rp) then
                     L(1) = lambda(1)*(p_y - rho_*c*v_y)
                   else

                     if    (L1_wave == 'poinsot_model'  ) then 
                       L(1) = L1Poinsot(Ly,p_,p_g,c,M_max)
                     elseif(L1_wave == 'pirozzoli_model') then
                       L(1) = L1Pirozzoli(v_,c,p_,p_g,i_dy)
                     elseif(L1_wave == 'zero_reflection') then
                       L(1) = 0.0_rp
                     else
                       print*, 'L1_wave', trim(L1_wave), ' is not implemented.'
                       stop
                     endif
                   endif

                   L(2) = lambda(2) * ((c**2)*r_y - p_y)
                   L(3) = lambda(3) * (u_y)
                   L(4) = lambda(4) * (w_y)
                   L(5) = lambda(5) * (p_y + rho_*c*v_y)

                   d(1) = 1._rp/c**2 *(L(2)+0.5_rp*(L(5)+L(1)))
                   d(2) = 0.5_rp*(L(5)+L(1))
                   d(3) = L(3)
                   d(4) = 1._rp/(2*rho_*c)*(L(5)-L(1))
                   d(5) = L(4)

                   RHS(i,jg,k,1) =    d(1)
                   RHS(i,jg,k,2) = u_*d(1)+rho_*d(3)
                   RHS(i,jg,k,3) = v_*d(1)+rho_*d(4)
                   RHS(i,jg,k,4) = w_*d(1)+rho_*d(5)
                   RHS(i,jg,k,5) = ek*d(1) + d(2)/(gamma0-1.0_rp) + rho_*u_*d(3) + rho_*v_*d(4) + rho_*w_*d(5)

                   ! evolving phi via RK
                   phi(i,jg,k,:) = aik * phi_n(i,jg,k,:) + bik * (phi(i,jg,k,:) - dt * RHS(i,jg,k,:))

                enddo
             enddo

             do j    = 2,GN
                do i = sx,ex
                   jg = b%node + b%norm*j
                   phi(i,jg,k,:) = phi(i,b%node+b%norm,k,:)
                enddo
             enddo

          enddo

          case default
            print*, ' nscbc_outflow is not implemented for face ', b%face
            stop

        endselect

        deallocate(cd_coef)

        return
end subroutine nscbc_outflow

function L1Poinsot(Lx,p,p_inf,c,M_max) result(L1)
        implicit none
        real(rp), intent(in) :: Lx, p, p_inf, c, M_max
        real(rp)             :: L1

        ! local declarations
        real(rp), parameter  :: sigma = 0.27_rp
        
        L1 = sigma*(1.0_rp - M_max*M_Max) * c/Lx * (p - p_inf)

        return
end function L1Poinsot

function L1Pirozzoli(u,c,p,p_inf,i_dx) result(L1)
        implicit none
        real(rp), intent(in) :: u, c, p, p_inf, i_dx
        real(rp)             :: L1

        L1 = i_dx*(u-c)*(p_inf - p)

        return
end function L1Pirozzoli

subroutine inflow_bc(phi,iFlow,b)
! ----------------------------------------------------------------------
!
!       This subroutine provides a supersonic inflow boundary condition
!
! ----------------------------------------------------------------------
        
        use math_tools_module, only: smooth_step
        use inflow_module    , only: turbulent_inflow
        
        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi
        type(ifl)                                , intent(inout) :: iFlow
        type(face_type)                          , intent(in)    :: b

        ! local declaration
        real(rp), parameter :: alpha = 0.75_rp*(gamma0-1.0_rp)/gamma0
        real(rp) :: r_, u_, v_, w_, T_, p_, step, istep
        real(rp) :: um, vm, wm, Tm, rm
        real(rp) :: uf, vf, wf, Tf
        integer  :: i,j,k, j0, j1, l
        logical  :: logic_step = .false.

        ! check face
        if(b%face /= 'W') then
          if(rank == root) print*, 'Inflow is not implemented for face ', b%face
          call secure_stop
        endif

        step = 1.0_rp
        if(iFlow%SmthFlag) step = smooth_step(0.0_rp, 0.3_rp*real(itmax,rp), real(it,rp))
        
        if(iFlow%TurbFlag) call turbulent_inflow(iFlow)

        j0 = sy-1
        j1 = sy-1
        if(trim(inflow_profile) == 'smooth_body_inflow') then
          istep = 0.22_rp/0.032_rp
          logic_step = .true.
        endif
        if(trim(inflow_profile) == 'supersonic_ramp')    then
          istep = 0.0_rp
          logic_step = .true.
        endif

        if(logic_step) then
           !$acc parallel default(present)
           !$acc loop gang, vector
           do j = sy,ey
              if(y(j) > 0) then
                 j0 = j
                 exit
              endif
           enddo
           !$acc end parallel
           !
           !$acc parallel default(present)
           !$acc loop gang, vector
           do j = sy,ey
              if(y(j) > istep) then
                 j1 = j 
                 exit
              endif
           enddo
           !$acc end parallel
        endif
        
        !$acc parallel default(present)
        !$acc loop gang, vector collapse(2)
        do       k = sz,ez
           do    j = sy,ey
              ! Taylor hypothesis
              do l = 1,5
                 do i = lbx+1,sx-1
                    phi(i,j,k,l) = phi(i,j,k,l) - &
                    dt*c_rk(ik)*u_inf/xstep(i)*(phi(i,j,k,l) - phi(i-1,j,k,l))
                 enddo
              enddo
              do i = lbx,lbx
                
                 ! mean quantities
                 rm = iflow%mean(i,j,1)
                 um = iflow%mean(i,j,2)
                 vm = iflow%mean(i,j,3)
                 wm = iflow%mean(i,j,4)
                 Tm = iflow%mean(i,j,5)

                 ! turbulent quantities
                 uf = 0.0_rp
                 vf = 0.0_rp
                 wf = 0.0_rp
                 if(j > j1) then
                   uf = u_inf*iFlow%turb(1,j-j1+j0,k)
                   vf = u_inf*iFlow%turb(2,j-j1+j0,k)
                   wf = u_inf*iFlow%turb(3,j-j1+j0,k)
                 endif
                 Tf = -alpha * uf * um / Tm
                
                 ! global quantities
                 u_ = step*(um + uf)
                 v_ =       vm + vf
                 w_ =       wm + wf
                 T_ =       Tm + Tf
                 if(Mach > 1.0_rp) then
                   p_ = 1.0_rp
                   r_ = p_/T_
                 else
                   r_ = phi(sx,j,k,1)
                   p_ = r_*T_
                 endif
        
                 ! assign to conservative variables
                 phi(i,j,k,1) = r_
                 phi(i,j,k,2) = r_ * u_
                 phi(i,j,k,3) = r_ * v_
                 phi(i,j,k,4) = r_ * w_
                 phi(i,j,k,5) = p_/(gamma0-1.0_rp) + 0.5_rp*r_*(u_*u_ + v_*v_ + w_*w_)

              enddo
           enddo
        enddo
        !$acc end parallel

        return
end subroutine inflow_bc





subroutine compute_max_Mach(M_Max)
        implicit none
        real(rp), intent(out) :: M_max
        real(rp)              :: u_, v_, w_, T_, my_M_Max

        my_M_Max = 0.0_rp
        do k = sz,ez
           do j = sy,ey
              do i = sx,ex
                
                 u_ = U(i,j,k)
                 v_ = V(i,j,k)
                 w_ = W(i,j,k)
                 T_ = T(i,j,k)

                 my_M_max = max(((u_*u_ + v_*v_ + w_*w_)/(gamma0*T_)), my_M_max)

              enddo
           enddo
        enddo
        my_M_max = sqrt(my_M_max)

        call MPI_allreduce(my_M_max, M_max, 1, MPI_RP, MPI_MAX, mpi_comm_cart, mpi_err)

        return
end subroutine compute_max_Mach

subroutine set_wall_speed(vel,b)
        use math_tools_module, only: smooth_step

        implicit none
        type(face_type), intent(in)  :: b
        real(rp)       , intent(out) :: vel
        real(rp)                     :: step

        vel = 0.0_rp

        step = smooth_step(0.0_rp, 0.3_rp*real(itmax,rp), real(it,rp))

        ! II stokes wall speed
        if    (ic == 'II_stokes_problem_x') then
          vel = u_inf*cos(0.25_rp*pi*time)
        elseif(ic == 'II_stokes_problem_y') then
          vel = u_inf*cos(0.25_rp*pi*time)
        elseif(ic == 'II_stokes_problem_z') then
          vel = u_inf*cos(0.25_rp*pi*time)

        ! coutte wall speed
        elseif(ic == 'couette_x' .and. b%face == 'N') then
          vel = step*u_inf
        elseif(ic == 'couette_y' .and. b%face == 'E') then
          vel = step*u_inf
        elseif(ic == 'couette_z' .and. b%face == 'F') then
          vel = step*u_inf

        ! stationary couette flow
        elseif(ic == 'steady_couette_x' .and. b%face == 'N') then
          vel = u_inf
        elseif(ic == 'steady_couette_y' .and. b%face == 'E') then
          vel = u_inf
        elseif(ic == 'steady_couette_z' .and. b%face == 'F') then
          vel = u_inf

        endif

        return
end subroutine set_wall_speed















subroutine oblique_shock_bc(phi,b)
! ----------------------------------------------------------------------
!       
!       This subroutine provides the correct upper boundary condition
!       for the swbli problem
!
! ----------------------------------------------------------------------
        
        use parameters_module
        use math_tools_module     , only: DegToRad
        use fluid_functions_module, only: oblique_shock

        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: phi
        type(face_type)                          , intent(in)    :: b

        ! parameters
        real(rp) :: thetaw   !< wedge angle DEGREZ (3.75_rp)
        real(rp) :: x1       !< impingement pont (DeGrez 100)

        ! local declaration
        real(rp), parameter :: bfr = 0.2_rp
        real(rp) :: r1, p1, T1, u1, v1, M1 !< left - shock state
        real(rp) :: r2, p2, T2, u2, v2, M2 !< right- shock state
        real(rp) :: beta                   !< shock angle
        real(rp) :: theta, slope1, x0, shock1
        integer  :: i,j,k,jg

        if(b%face /= 'N') then
          print*, 'Oblique_shock is not implemented for face ', b%face
          stop
        endif

        thetaw = wedge_angle    !< wedge angle DEGREZ (3.75_rp)
        x1     = impin_point    !< impingement pont (DeGrez 100)

        !
        ! STATE 1
        !
        M1 = Mach
        r1 = 1.0_rp
        p1 = 1.0_rp
        T1 = p1/r1
        u1 = M1 * sqrt(gamma0*T1)
        v1 = 0.0_rp
        
        !
        ! STATE 2
        !
        theta = DegToRad(thetaw)
        call oblique_shock(r1,p1,T1,M1,theta,r2,p2,T2,M2,beta)

        u2 = M2 * sqrt(gamma0*T2) * cos(theta)
        v2 = M2 * sqrt(gamma0*T2) * sin(theta)

        ! angular coefficient
        slope1 = tan(pi-beta)
        x0     = x1 + ymax/slope1      !< shock initial position
        
        !
        ! forcing the ny point to respect R-H around the shock
        !
        do k       = lbz,ubz
           do j    = ny+1,uby
              do i = lbx,ubx

                   shock1 = y(j) - slope1 * (x(i) - x1)

                   if(    shock1 < 0.0_rp .and. x(i) > x0 - bfr) then
                           
                     phi(i,j,k,1) = r1
                     phi(i,j,k,2) = r1*u1
                     phi(i,j,k,3) = 0.0_rp
                     phi(i,j,k,4) = 0.0_rp
                     phi(i,j,k,5) = p1/(gamma0-1.0_rp) + 0.5_rp*r1*u1*u1

                   elseif(shock1 > 0.0_rp .and. x(i) < x0 + bfr) then

                     phi(i,j,k,1) =  r2
                     phi(i,j,k,2) =  r2*u2
                     phi(i,j,k,3) = -r2*v2
                     phi(i,j,k,4) = 0.0_rp
                     phi(i,j,k,5) = p2/(gamma0-1.0_rp) + 0.5_rp*r2*(u2*u2+v2*v2)

                   endif

              enddo
           enddo
        enddo
        !
        ! extrapoling on the ghosts elsewhere
        !
        do       k = lbz,ubz
           do    j = 1,GN
              do i = lbx,ubx

                 jg = b%node + b%norm*j      !< ghost node
                 shock1 = y(jg) - slope1 * (x(i) - x1)

                 if(    shock1 < 0.0_rp .and. x(i) < x0 - bfr) then
                   phi(i,jg,k,:) = phi(i,b%node,k,:)
                 elseif(shock1 > 0.0_rp .and. x(i) > x0 + bfr) then
                   phi(i,jg,k,:) = phi(i,b%node,k,:)
                 endif
                 
              enddo
           enddo
        enddo



        return
end subroutine oblique_shock_bc


subroutine flate_plate_bc(phi,b)
        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: phi
        type(face_type)                          , intent(in)    :: b

        ! local declarations
        real(rp), parameter :: x0 = 0.0_rp
        integer             :: i,j,k,jg, ji

        if(b%face /= 'S') then
          print*, 'Flate plate is not implemented for face ', b%face
          stop
        endif

        do  k      = lbz,ubz
           do j    = 1  ,GN
              do i = lbx,ubx
                      
                 jg = b%node + b%norm*j      !< ghost node
                 ji = b%node - b%norm*(j-1)  !< inner node

                 if(x(i) < x0) then !< symmetry
                   phi(i,jg,k,1) =    phi(i,ji,k,1) 
                   phi(i,jg,k,2) =    phi(i,ji,k,2) 
                   phi(i,jg,k,3) =  - phi(i,ji,k,3) 
                   phi(i,jg,k,4) =    phi(i,ji,k,4) 
                   phi(i,jg,k,5) =    phi(i,ji,k,5) 

                 else               !< adiabatic no-slip wall
                   phi(i,jg,k,1) =   phi(i,ji,k,1) 
                   phi(i,jg,k,2) = - phi(i,ji,k,2) 
                   phi(i,jg,k,3) = - phi(i,ji,k,3) 
                   phi(i,jg,k,4) = - phi(i,ji,k,4) 
                   phi(i,jg,k,5) =   phi(i,ji,k,5) 

                 endif

              enddo
           enddo
        enddo

        return
end subroutine flate_plate_bc



subroutine nscbc_wall(phi,b,Twall)
        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: phi
        type(face_type)                          , intent(in)    :: b
        real(rp)                                 , intent(in)    :: Twall

        real(rp) :: L1, d1
        real(rp) :: aik, bik
        real(rp) :: r_gh, irgh, u_gh, v_gh, w_gh, ekgh, p_gh, T_gh, c_gh
        real(rp) :: r_in, irin, u_in, v_in, w_in, ekin, p_in, T_in
        real(rp) :: r_js, irjs, u_js, v_js, w_js, p_js
        real(rp) :: cs, i_st, p_y, v_y
        integer  :: i,j,k,ji, jg, s, js, fdL, fdR
        integer  :: node, norm

        aik = a_rk(ik)
        bik = b_rk(ik)

        node = b%node
        norm = b%norm
       
        selectcase(b%face)
        case('N')
        
        !$acc parallel default(present)
        !$acc loop gang, vector collapse(2)
        do        k = lbz,ubz
           do     i = lbx,ubx
              !$acc loop seq
               do j = 1,GN

                 jg = node + norm*j      !< ghost node
                 ji = node - norm*(j-1)  !< inner node

                 ! compute scalar on internal nodes
                 r_in = phi(i,ji,k,1)
                 irin = 1.0_rp/r_in
                 u_in = phi(i,ji,k,2) * irin
                 v_in = phi(i,ji,k,3) * irin
                 w_in = phi(i,ji,k,4) * irin
                 ekin = 0.5_rp*(u_in*u_in + v_in*v_in + w_in*w_in)
                 p_in = (gamma0-1._rp)*(phi(i,ji,k,5)-r_in*ekin)
                 T_in = p_in*irin

                 ! compute scalar on ghost nodes 
                 r_gh = phi(i,jg,k,1)
                 irgh = 1.0_rp/r_gh
                 u_gh = phi(i,jg,k,2)*irgh
                 v_gh = phi(i,jg,k,3)*irgh 
                 w_gh = phi(i,jg,k,4)*irgh
                 ekgh = 0.5_rp*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)
                 p_gh = cv_i * (phi(i,jg,k,5) - r_gh*ekgh)
                 T_gh = p_gh*irgh
                 c_gh = sqrt(gamma0 * T_gh)

                 ! computing derivatives 
                 i_st = ystep_i(jg)

                 p_y = 0.0_rp
                 v_y = 0.0_rp
                 !$acc loop seq
                 do s = - bward_fd_order,0

                    js = jg+s

                    cs = bward_1(s)

                    r_js = phi(i,js,k,1)
                    irjs = 1.0_rp/r_js
                    u_js = phi(i,js,k,2)*irjs
                    v_js = phi(i,js,k,3)*irjs
                    w_js = phi(i,js,k,4)*irjs
                    p_js = cv_i*(phi(i,js,k,5) - 0.5_rp*r_js*(u_js*u_js + v_js*v_js + w_js*w_js))

                    p_y = p_y + cs * p_js
                    v_y = v_y + cs * v_js

                 enddo
                 p_y = p_y * i_st
                 v_y = v_y * i_st
                        
                 !!! evolving RHO via RK
                 L1            = (v_gh + norm*c_gh) * (p_y + norm*r_gh*c_gh*v_y)
                 d1            = 1._rp/c_gh**2 *L1
                 RHS(i,jg,k,1) = d1
                 phi(i,jg,k,1) = aik * phi_n(i,jg,k,1) + bik * (phi(i,jg,k,1) - dt * RHS(i,jg,k,1))
        
                 ! update ghosts
                 r_gh = phi(i,jg,k,1)
                 irgh = 1.0_rp/r_gh
                 u_gh =         - u_in
                 v_gh =         - v_in
                 w_gh =         - w_in
                 T_gh = 2*twall - T_in  
                 p_gh = r_gh*T_gh
                
                 phi(i,jg,k,1) = r_gh
                 phi(i,jg,k,2) = r_gh*u_gh 
                 phi(i,jg,k,3) = r_gh*v_gh 
                 phi(i,jg,k,4) = r_gh*w_gh 
                 phi(i,jg,k,5) = p_gh/(gamma0-1.0_rp) + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)

              enddo
           enddo
        enddo
        !$acc end parallel

        case('S')

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(2)
        do       k = lbz,ubz
           do    i = lbx,ubx
              !$acc loop seq
              do j = 1,GN

                 jg = node + norm*j      !< ghost node
                 ji = node - norm*(j-1)  !< inner node

                 ! compute scalar on internal nodes
                 r_in = phi(i,ji,k,1)
                 irin = 1.0_rp/r_in
                 u_in = phi(i,ji,k,2) * irin
                 v_in = phi(i,ji,k,3) * irin
                 w_in = phi(i,ji,k,4) * irin
                 ekin = 0.5_rp*(u_in*u_in + v_in*v_in + w_in*w_in)
                 p_in = (gamma0-1._rp)*(phi(i,ji,k,5)-r_in*ekin)
                 T_in = p_in*irin

                 ! compute scalar on ghost nodes 
                 r_gh = phi(i,jg,k,1)
                 irgh = 1.0_rp/r_gh
                 u_gh = phi(i,jg,k,2)*irgh
                 v_gh = phi(i,jg,k,3)*irgh 
                 w_gh = phi(i,jg,k,4)*irgh
                 ekgh = 0.5_rp*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)
                 p_gh = cv_i * (phi(i,jg,k,5) - r_gh*ekgh)
                 T_gh = p_gh*irgh
                 c_gh = sqrt(gamma0 * T_gh)

                 ! computing derivatives 
                 i_st = ystep_i(jg)

                 p_y = 0.0_rp
                 v_y = 0.0_rp
                 !$acc loop seq
                 do s = 0,fward_fd_order

                    js = jg+s

                    cs = fward_1(s)

                    r_js = phi(i,js,k,1)
                    irjs = 1.0_rp/r_js
                    u_js = phi(i,js,k,2)*irjs
                    v_js = phi(i,js,k,3)*irjs
                    w_js = phi(i,js,k,4)*irjs
                    p_js = cv_i*(phi(i,js,k,5) - 0.5_rp*r_js*(u_js*u_js + v_js*v_js + w_js*w_js))

                    p_y = p_y + cs * p_js
                    v_y = v_y + cs * v_js

                 enddo
                 p_y = p_y * i_st
                 v_y = v_y * i_st
                        
                 !!! evolving RHO via RK
                 L1            = (v_gh + norm*c_gh) * (p_y + norm*r_gh*c_gh*v_y)
                 d1            = 1._rp/c_gh**2 *L1
                 RHS(i,jg,k,1) = d1
                 phi(i,jg,k,1) = aik * phi_n(i,jg,k,1) + bik * (phi(i,jg,k,1) - dt * RHS(i,jg,k,1))
        
                 ! update ghosts
                 r_gh = phi(i,jg,k,1)
                 irgh = 1.0_rp/r_gh
                 u_gh =         - u_in
                 v_gh =         - v_in
                 w_gh =         - w_in
                 T_gh = 2*twall - T_in  
                 p_gh = r_gh*T_gh
                
                 phi(i,jg,k,1) = r_gh
                 phi(i,jg,k,2) = r_gh*u_gh 
                 phi(i,jg,k,3) = r_gh*v_gh 
                 phi(i,jg,k,4) = r_gh*w_gh 
                 phi(i,jg,k,5) = p_gh/(gamma0-1.0_rp) + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)

              enddo
           enddo
        enddo
        !$acc end parallel

        case default 
          print*, 'nscbc_wall is not implemented for face ', trim(b%face)
          stop

        endselect

        return
end subroutine nscbc_wall








subroutine supersonic_outflow(phi,b)
        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: phi
        type(face_type)                          , intent(in)    :: b

        integer :: i,j,k,ig, jg
        integer :: bnode, bnorm

        bnode = b%node
        bnorm = b%norm
        
        selectcase(b%face)
          case('E')

          !$acc parallel default(present)
          !$acc loop gang, vector collapse(2)
          do       k = sz,ez
             do    j = sy,ey
                !$acc loop seq
                do i = 1,GN

                   ig = bnode + bnorm*i  !< ghost node

                   phi(ig,j,k,:) = phi(bnode,j,k,:)

                enddo
             enddo
          enddo
          !$acc end parallel
        case('N')

          !$acc parallel default(present)
          !$acc loop gang, vector collapse(2)
          do       k = sz,ez
             do    i = sx,ex
                !$acc loop seq
                do j = 1,GN

                   jg = bnode + bnorm*j  !< ghost node

                   phi(i,jg,k,:) = phi(i,bnode,k,:)

                enddo
             enddo
          enddo
          !$acc end parallel

        case default
          print*, ' supersonic outflow is not implemented for face ', trim(b%face)
          stop

        endselect

        return
end subroutine supersonic_outflow





subroutine Neumann_wall(phi,b)

        use math_tools_module, only: DegToRad

        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: phi
        type(face_type)                          , intent(in)    :: b

        real(rp), dimension(5) :: c
        real(rp)               :: x0, l1, l2
        real(rp), parameter    :: r1     = 52.80_rp
        real(rp), parameter    :: r2     = 32.36_rp
        real(rp), parameter    :: L0     = 150.0_rp
        real(rp), parameter    :: alpha1 = 10
        real(rp), parameter    :: alpha2 = 22
        integer                :: i,j,k, jg, ji
        
        ! geometrical locations
        l1 = r1*cos(DegToRad(alpha1))
        l2 = r2*cos(DegToRad(alpha2))
        x0 = (l1+l2)/L0

        c = (/1.0_rp, 1.0_rp, -1.0_rp, 1.0_rp, 1.0_rp/)
        if(viscous) c = (/1.0_rp, -1.0_rp, -1.0_rp, -1.0_rp, 1.0_rp/)

        selectcase(b%face)
          case('N')

          do       k = lbz,ubz
             do    j = 1,GN
                do i = lbx,ubx

                   jg = b%node + b%norm*j
                   ji = b%node - b%norm*(j-1)

                   if(x(i) < x0) then
                     phi(i,jg,k,:) = phi(i,ji,k,:) 

                   else
                     phi(i,jg,k,:) =  c(:) * phi(i,ji,k,:) 
                   endif

                enddo
             enddo
          enddo

        case default
          print*, ' supersonic outflow is not implemented for face ', trim(b%face)
          stop

        endselect

        return
end subroutine Neumann_wall







subroutine subsonic_inflow_bernardini1(phi,iFlow,b)

        use math_tools_module, only: smooth_step
        use inflow_module    , only: turbulent_inflow 

        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi
        type(ifl)                                , intent(inout) :: Iflow
        type(face_type)                          , intent(in)    :: b

        ! local declarations
        real(rp) :: step
        integer  :: i,j,k,ig

        real(rp) :: r_l, irl, rul, rvl, rwl, rel
        real(rp) :: u_l, v_l, w_l, ekl, T_l, p_l, ccl, c_l

        real(rp) :: r_r, irr, rur, rvr, rwr, rer
        real(rp) :: u_r, v_r, w_r, ekr, T_r, p_r, ccr, c_r

        real(rp) :: rmach1, fmach, deltaw5, deltav1, deltav2, deltav5
        real(rp) :: r_, u_, v_, w_, p_, gmach

        step = 1.0_rp
        if(iFlow%SmthFlag) step = smooth_step(0.0_rp, 0.3_rp*real(itmax,rp), real(it,rp))

        if(iFlow%TurbFlag) call turbulent_inflow(iFlow)

        if(b%face == 'W') then

          do       k = sz,ez
             do    j = sy,ey

               i = 1
               r_l = phi(i,j,k,1)
               irl = 1.0_rp/r_l
               rul = phi(i,j,k,2)
               rvl = phi(i,j,k,3)
               rwl = phi(i,j,k,4)
               rel = phi(i,j,k,5)

               u_l = rul*irl
               v_l = rvl*irl
               w_l = rwl*irl
               ekl = 0.5_rp*(u_l*u_l + v_l*v_l + w_l*w_l)
               p_l = (gamma0-1.0_rp)*(rel - r_l*ekl)
               T_l = p_l*irl
               ccl = gamma0*t_l
               c_l = sqrt(ccl)
               rMach1 = u_l/c_l

               i = 2
               r_r = phi(i,j,k,1)
               irr = 1.0_rp/r_r
               rur = phi(i,j,k,2)
               rvr = phi(i,j,k,3)
               rwr = phi(i,j,k,4)
               rer = phi(i,j,k,5)

               u_r = rur*irr
               v_r = rvr*irr
               w_r = rwr*irr
               ekr = 0.5_rp*(u_r*u_r + v_r*v_r + w_r*w_r)
               p_r = (gamma0-1.0_rp)*(rer - r_r*ekr)
               T_r = p_r*irr
               ccr = gamma0*t_r
               c_r = sqrt(ccr)

               deltaw5 = (p_r - p_l)*irl/c_l - (u_r - u_l)
               !fmach = rmach1*rmach1+1.0_rp
               !fmach = fmach/(rmach1+1.0_rp)**2
               !deltav1 = r_l/c_l*deltaw5*fmach
               !deltav5 = r_l*c_l*deltaw5*fmach
               !fmach = rmach1
               !fmach = fmach/(rmach1+1.0_rp)**2
               !deltav2 = -4*deltaw5*fmach
               fmach = 1.0_rp-gamma0*rmach1+(gamma0-1.0_rp)*rmach1**2
               fmach = fmach/(1.0_rp+gamma0*rmach1+(gamma0-1.0_rp)*rmach1**2)
               gmach = 1.0_rp-rmach1 - fmach*(rmach1+1.0_rp)
               gmach = gmach/rmach1
               deltav1 = (0.5_rp*r_l/c_l*(1.0_rp+gmach+fmach))*deltaw5
               deltav2 = -0.5_rp*deltaw5*(1.0_rp-fmach)
               deltav5 = (0.5_rp*r_l/c_l*(1.0_rp+fmach))*deltaw5

               r_ = r_l - deltav1
               u_ = u_l - deltav2
               v_ = 0.0_rp
               w_ = 0.0_rp
               p_ = p_l - deltav5

               do i = 1,GN

                  ig = b%node + b%norm*i

                  phi(ig,j,k,1) = r_
                  phi(ig,j,k,2) = r_*u_
                  phi(ig,j,k,3) = r_*v_
                  phi(ig,j,k,4) = r_*w_
                  phi(ig,j,k,5) = p_/(gamma0-1.0_rp)+ 0.5_rp*r_*(u_*u_ + v_*v_ + w_*w_)
               enddo

             enddo
        enddo

        else
          if(rank == root) print*, 'Inflow is not implemented for face ', b%face
          call secure_stop
        endif

        return
end subroutine subsonic_inflow_bernardini1



subroutine subsonic_inflow_bernardini2(phi,b)
        use inflow_module

        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi
        type(face_type)                          , intent(in) :: b

        integer :: i,j,k,l
        real(rp) :: r_1, ir1,u_1, v_1, w_1, ek1, p_1, T_1, cc1, c_1
        real(rp) :: r_2, ir2,u_2, v_2, w_2, ek2, p_2, T_2, cc2, c_2
        real(rp) :: d1k, d2k, d3k, d1kh, d2kh, d3kh

        real(rp) :: rhoh, rhohc, rhohci
        real(rp) :: r_, u_, v_, w_, p_, T_

        real(rp), dimension(5) :: deltaw, deltav
        real(rp), dimension(5,5) :: LL, RR


        if(b%face == 'W') then
          do    k = sz,ez
             do j = sy,ey

                i = sx+1
                r_2 = phi(i,j,k,1)
                ir2 = 1.0_rp/r_2
                u_2 = phi(i,j,k,2)*ir2
                v_2 = phi(i,j,k,3)*ir2
                w_2 = phi(i,j,k,4)*ir2
                ek2 = 0.5_rp*(u_2*u_2 + v_2*v_2 + w_2*w_2)
                p_2 = (gamma0-1.0_rp)*(phi(i,j,k,5) - r_2*ek2)
                T_2 = p_2*ir2
                cc2 = gamma0*T_2
                c_2 = sqrt(cc2)

                i = sx
                r_1 = phi(i,j,k,1)
                ir1 = 1.0_rp/r_1
                u_1 = phi(i,j,k,2)*ir1
                v_1 = phi(i,j,k,3)*ir1
                w_1 = phi(i,j,k,4)*ir1
                ek1 = 0.5_rp*(u_1*u_1 + v_1*v_1 + w_1*w_1)
                p_1 = (gamma0-1.0_rp)*(phi(i,j,k,5) - r_1*ek1)
                T_1 = p_1*ir1
                cc1 = gamma0*T_1
                c_1 = sqrt(cc1)

                deltav(1) = r_2 - r_1
                deltav(2) = u_2 - u_1
                deltav(3) = v_2 - v_1
                deltav(4) = w_2 - w_1
                deltav(5) = p_2 - p_1

                ! left on 1
                d1k = 1.0_rp
                d2k = 0.0_rp
                d3k = 0.0_rp 
                LL(1,:) = [d1k, 0.0_rp, d3k, -d2k,-d1k/cc1]
                LL(2,:) = [d2k,-d3k, 0.0_rp,  d1k,-d2k/cc1]
                LL(3,:) = [d3k, d2k,-d1k,  0.0_rp,-d3k/cc1]
                LL(4,:) = [0.0_rp, d1k, d2k,  d3k, ir1/c_1]
                LL(5,:) = [0.0_rp,-d1k,-d2k, -d3k, ir1/c_1]

                deltaw = matmul(LL,deltav)
        
                ! right on nx
                d1kh = 0.5_rp*d1k
                d2kh = 0.5_rp*d2k
                d3kh = 0.5_rp*d3k
                
                rhoh = 0.5_rp*r_1
                rhohc = rhoh*c_1
                rhohci = rhoh/c_1

                RR(1,:) = [d1k , d2k , d3k , rhohci, rhohci]
                RR(2,:) = [0.0_rp , -d3k, d2k , d1kh  , -d1kh]
                RR(3,:) = [d3k , 0.0_rp , -d1k, d2kh  , -d2kh] 
                RR(4,:) = [-d2k, d1k , 0.0_rp , d3kh  , -d3kh]
                RR(5,:) = [0.0_rp , 0.0_rp , 0.0_rp , rhohc , rhohc]
        
                deltaw(1) = 0.0_rp
                deltaw(2) = 0.0_rp
                deltaw(3) = 0.0_rp
                deltaw(4) = 0.0_rp

                deltav = matmul(RR, deltaw)
                
                u_ = iflow%mean(0,j,2)
                v_ = iflow%mean(0,j,3)
                w_ = iflow%mean(0,j,4)
                T_ = iflow%mean(0,j,5)
                !p_ = p_1 - deltav(5)
                !r_ = p_/T_
                p_ = 1.0_rp
                r_ = 1.0_rp

                do l = 1,GN
                   phi(sx-l,j,k,1) = r_
                   phi(sx-l,j,k,2) = r_*u_
                   phi(sx-l,j,k,3) = r_*v_
                   phi(sx-l,j,k,4) = r_*w_
                   phi(sx-l,j,k,5) = p_/(gamma0-1.0_rp) + 0.5_rp*r_*(u_*u_+v_*v_+w_*w_)
                enddo


             enddo
          enddo

        else
          write(*,*) ' Outflow boundary is not implemented for face', b%face
          
        endif


        return
end subroutine subsonic_inflow_bernardini2



































end module bc_module
