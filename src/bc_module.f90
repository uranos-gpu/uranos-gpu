module bc_module
use parameters_module
use storage_module
use mpi_module
use fluid_functions_module
use mpi_comm_module
use profiling_module
use thermodynamics_module
use mesh_module, only: deta_dx, deta_dy

contains
subroutine set_bc_conditions
        implicit none
        type(face_type), dimension(6) :: all_bound
        type(face_type)               :: bound
        integer                       :: f
        real(rp)                      :: Twall

        call StartProfRange("set_bc_conditions")

        all_bound%node = (/ sx , ex , sy , ey , sz ,  ez/)
        all_bound%norm = (/ -1 ,  1 , -1 ,  1 , -1 ,  1 /)
        all_bound%face = (/ 'W', 'E', 'S', 'N', 'B', 'F'/)

        do f = 1, nb_neigh
           if(my_neighbour(f) == MPI_PROC_NULL) then

             bound%node = all_bound(f)%node
             bound%norm = all_bound(f)%norm
             bound%face = all_bound(f)%face

             selectcase(trim(bc(f)))

             case('neumann')
               call neumann(phi,bound)

             case('adiabatic_wall')
               call adiabatic_wall(phi,bound)

             case('isothermal_wall')

                  tWall = 1.0_rp
                  call isothermal_wall(phi,bound,Twall)

             case('diabatic_wall')

                  tWall = Trat*T_recovery
                  call isothermal_wall(phi,bound,Twall)

             case('nscbc_outflow', 'nscbc_inflow_generalized', 'extrapolation')
               extrapolation_bc(f) = .true.
               call extrapolation(phi,bound)

             case default

                write(*,*) ' Boundary condition ', "'"//trim(bc(f))//"'", ' is not implemented.'
                stop

             endselect

           endif
        enddo

        call EndProfRange

        return
end subroutine set_bc_conditions

subroutine mpi_bc_communications
        implicit none

        call StartProfRange("mpi_bc_communications")

call  mpi_share(mpi_comm_cart,type_send_cons,type_recv_cons,my_neighbour,dims, &
                phi_bfr_send_E, phi_bfr_send_W, phi_bfr_recv_E, phi_bfr_recv_W, &
                phi_bfr_send_N, phi_bfr_send_S, phi_bfr_recv_N, phi_bfr_recv_S, &
                phi_bfr_send_B, phi_bfr_send_F, phi_bfr_recv_B, phi_bfr_recv_F, &
                phi)

        call EndProfRange

        return
end subroutine mpi_bc_communications

subroutine neumann(phi,b)
        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi
        type(face_type)                          , intent(in   ) :: b

        integer :: i,j,k,l
        integer :: node, norm

        node = b%node
        norm = b%norm

        select case(b%face)
          case('W','E')

            !$acc parallel default(present)
            !$acc loop gang, vector collapse(4)
            do          l = 1,5
               do       k = lbz,ubz
                  do    j = lby,uby
                     do i = 1,GN
                        phi(node+norm*i,j,k,l) = phi(node-norm*(i-1),j,k,l)
                     enddo
                  enddo
               enddo
            enddo
            !$acc end parallel

          case('S','N')

            !$acc parallel default(present)
            !$acc loop gang, vector collapse(4)
            do          l = 1,5
               do       k = lbz,ubz
                  do    j = 1,GN
                     do i = lbx,ubx
                        phi(i,node+norm*j,k,l) = phi(i,node-norm*(j-1),k,l)
                     enddo
                  enddo
               enddo
            enddo
            !$acc end parallel

          case('B','F')

            !$acc parallel default(present)
            !$acc loop gang, vector collapse(4)
            do          l = 1,5
               do       k = 1,GN
                  do    j = lby,uby
                     do i = lbx,ubx
                        phi(i,j,node+norm*k,l) = phi(i,j,node-norm*(k-1),l)
                     enddo
                  enddo
               enddo
            enddo
            !$acc end parallel

        end select

end subroutine neumann

subroutine adiabatic_wall(phi,b)
        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: phi
        type(face_type)                          , intent(in)    :: b

        real(rp) :: gm1, igm1
        real(rp) :: r_in, u_in, v_in, w_in, p_in, T_in, irin, e_in
        real(rp) :: r_gh, u_gh, v_gh, w_gh, p_gh, T_gh
        integer  :: i,j,k
        integer  :: ig,jg,kg,ii,ji,ki
        integer  :: node, norm

        node = b%node
        norm = b%norm

        gm1  = gamma_Tref - 1.0_rp
        igm1 = 1.0_rp/gm1

        select case(b%face)

          case('E','W')

            if (thermally_perfect_gas > 0) then
               !$acc parallel default(present)
               !$acc loop gang, vector collapse(3)
               do        k = lbz,ubz
                  do     j = lby,uby
                     do  i = 1,GN

                        ig = node + norm*i
                        ii = node - norm*(i-1)

                        r_in = phi(ii,j,k,1)
                        irin = 1._rp/r_in
                        u_in = phi(ii,j,k,2) * irin
                        v_in = phi(ii,j,k,3) * irin
                        w_in = phi(ii,j,k,4) * irin
                        e_in = phi(ii,j,k,5) * irin - 0.5_rp*(u_in*u_in+v_in*v_in+w_in*w_in)
                        T_in = inverse_temperature(e_in)
                        p_in = r_in * T_in

                        u_gh =         - u_in
                        v_gh = - v_in
                        w_gh =         - w_in
                        T_gh =           T_in
                        p_gh = p_in
                        r_gh = p_gh/T_gh

                        phi(ig,j,k,1) = r_gh
                        phi(ig,j,k,2) = r_gh*u_gh
                        phi(ig,j,k,3) = r_gh*v_gh
                        phi(ig,j,k,4) = r_gh*w_gh
                        phi(ig,j,k,5) = r_gh*enthalpy_t(T_gh) - p_gh &
                                      + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)
                     enddo
                  enddo
               enddo
               !$acc end parallel
            else
               !$acc parallel default(present)
               !$acc loop gang, vector collapse(3)
               do        k = lbz,ubz
                  do     j = lby,uby
                     do  i = 1,GN

                        ig = node + norm*i
                        ii = node - norm*(i-1)

                        r_in = phi(ii,j,k,1)
                        irin = 1._rp/r_in
                        u_in = phi(ii,j,k,2) * irin
                        v_in = phi(ii,j,k,3) * irin
                        w_in = phi(ii,j,k,4) * irin
                        p_in = gm1*(phi(ii,j,k,5) &
                             - 0.5_rp*r_in*(u_in*u_in + v_in*v_in + w_in*w_in))
                        T_in = p_in*irin

                        u_gh =         - u_in
                        v_gh = - v_in
                        w_gh =         - w_in
                        T_gh =           T_in
                        p_gh = p_in
                        r_gh = p_gh/T_gh

                        phi(ig,j,k,1) = r_gh
                        phi(ig,j,k,2) = r_gh*u_gh
                        phi(ig,j,k,3) = r_gh*v_gh
                        phi(ig,j,k,4) = r_gh*w_gh
                        phi(ig,j,k,5) = p_gh*igm1 &
                                      + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)
                     enddo
                  enddo
               enddo
               !$acc end parallel
            endif

          case('S','N')

            if (thermally_perfect_gas > 0) then
               !$acc parallel default(present)
               !$acc loop gang, vector collapse(3)
               do        k = lbz,ubz
                  do     j = 1,GN
                     do  i = lbx,ubx

                        jg = node + norm*j
                        ji = node - norm*(j-1)

                        r_in = phi(i,ji,k,1)
                        irin = 1._rp/r_in
                        u_in = phi(i,ji,k,2) * irin
                        v_in = phi(i,ji,k,3) * irin
                        w_in = phi(i,ji,k,4) * irin
                        e_in = phi(i,ji,k,5) * irin - 0.5_rp*(u_in*u_in+v_in*v_in+w_in*w_in)
                        T_in = inverse_temperature(e_in)
                        p_in = r_in * T_in

                        u_gh = - u_in
                        v_gh =         - v_in
                        w_gh =         - w_in
                        T_gh =           T_in
                        p_gh = p_in
                        r_gh = p_gh/T_gh

                        phi(i,jg,k,1) = r_gh
                        phi(i,jg,k,2) = r_gh*u_gh
                        phi(i,jg,k,3) = r_gh*v_gh
                        phi(i,jg,k,4) = r_gh*w_gh
                        phi(i,jg,k,5) = r_gh*enthalpy_t(T_gh) - p_gh &
                                      + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)
                     enddo
                  enddo
               enddo
               !$acc end parallel
            else
               !$acc parallel default(present)
               !$acc loop gang, vector collapse(3)
               do        k = lbz,ubz
                  do     j = 1,GN
                     do  i = lbx,ubx

                        jg = node + norm*j
                        ji = node - norm*(j-1)

                        r_in = phi(i,ji,k,1)
                        irin = 1._rp/r_in
                        u_in = phi(i,ji,k,2) * irin
                        v_in = phi(i,ji,k,3) * irin
                        w_in = phi(i,ji,k,4) * irin
                        p_in = gm1*(phi(i,ji,k,5) &
                             - 0.5_rp*r_in*(u_in*u_in + v_in*v_in + w_in*w_in))
                        T_in = p_in*irin

                        u_gh = - u_in
                        v_gh =         - v_in
                        w_gh =         - w_in
                        T_gh =           T_in
                        p_gh = p_in
                        r_gh = p_gh/T_gh

                        phi(i,jg,k,1) = r_gh
                        phi(i,jg,k,2) = r_gh*u_gh
                        phi(i,jg,k,3) = r_gh*v_gh
                        phi(i,jg,k,4) = r_gh*w_gh
                        phi(i,jg,k,5) = p_gh*igm1 &
                                      + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)
                     enddo
                  enddo
               enddo
               !$acc end parallel
            endif

          case('B','F')

            if (thermally_perfect_gas > 0) then
               !$acc parallel default(present)
               !$acc loop gang, vector collapse(3)
               do        k = 1,GN
                  do     j = lby,uby
                     do  i = lbx,ubx

                        kg = node + norm*k
                        ki = node - norm*(k-1)

                        r_in = phi(i,j,ki,1)
                        irin = 1._rp/r_in
                        u_in = phi(i,j,ki,2) * irin
                        v_in = phi(i,j,ki,3) * irin
                        w_in = phi(i,j,ki,4) * irin
                        e_in = phi(i,j,ki,5) * irin - 0.5_rp*(u_in*u_in+v_in*v_in+w_in*w_in)
                        T_in = inverse_temperature(e_in)
                        p_in = r_in * T_in

                        u_gh = - u_in
                        v_gh =         - v_in
                        w_gh =         - w_in
                        T_gh =           T_in
                        p_gh = p_in
                        r_gh = p_gh/T_gh

                        phi(i,j,kg,1) = r_gh
                        phi(i,j,kg,2) = r_gh*u_gh
                        phi(i,j,kg,3) = r_gh*v_gh
                        phi(i,j,kg,4) = r_gh*w_gh
                        phi(i,j,kg,5) = r_gh*enthalpy_t(T_gh) - p_gh &
                                      + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)
                     enddo
                  enddo
               enddo
               !$acc end parallel
            else
               !$acc parallel default(present)
               !$acc loop gang, vector collapse(3)
               do        k = 1,GN
                  do     j = lby,uby
                     do  i = lbx,ubx

                        kg = node + norm*k
                        ki = node - norm*(k-1)

                        r_in = phi(i,j,ki,1)
                        irin = 1._rp/r_in
                        u_in = phi(i,j,ki,2) * irin
                        v_in = phi(i,j,ki,3) * irin
                        w_in = phi(i,j,ki,4) * irin
                        p_in = gm1*(phi(i,j,ki,5) &
                             - 0.5_rp*r_in*(u_in*u_in + v_in*v_in + w_in*w_in))
                        T_in = p_in*irin

                        u_gh = - u_in
                        v_gh =         - v_in
                        w_gh =         - w_in
                        T_gh =           T_in
                        p_gh = p_in
                        r_gh = p_gh/T_gh

                        phi(i,j,kg,1) = r_gh
                        phi(i,j,kg,2) = r_gh*u_gh
                        phi(i,j,kg,3) = r_gh*v_gh
                        phi(i,j,kg,4) = r_gh*w_gh
                        phi(i,j,kg,5) = p_gh*igm1 &
                                      + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)
                     enddo
                  enddo
               enddo
               !$acc end parallel
            endif

        end select
        return
end subroutine adiabatic_wall

subroutine isothermal_wall(phi,b,Twall)
        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: phi
        type(face_type)                          , intent(in)    :: b
        real(rp)                                 , intent(in)    :: Twall

        real(rp) :: gm1, igm1
        real(rp) :: r_in, u_in, v_in, w_in, p_in, T_in, irin, e_in
        real(rp) :: r_gh, u_gh, v_gh, w_gh, p_gh, T_gh
        integer  :: i,j,k
        integer  :: ig,jg,kg,ii,ji,ki
        integer  :: node, norm

        real(rp), parameter :: T_floor = 0.1_rp

        node = b%node
        norm = b%norm

        gm1  = gamma_Tref - 1.0_rp
        igm1 = 1.0_rp/gm1

        select case(b%face)
          case('W','E')

            if (thermally_perfect_gas > 0) then
               !$acc parallel default(present)
               !$acc loop gang, vector collapse(3)
               do        k = lbz,ubz
                  do     j = lby,uby
                     do  i = 1,GN

                        ig = node + norm*i
                        ii = node - norm*(i-1)

                        r_in = phi(ii,j,k,1)
                        irin = 1._rp/r_in
                        u_in = phi(ii,j,k,2) * irin
                        v_in = phi(ii,j,k,3) * irin
                        w_in = phi(ii,j,k,4) * irin
                        e_in = phi(ii,j,k,5) * irin - 0.5_rp*(u_in*u_in+v_in*v_in+w_in*w_in)
                        T_in = inverse_temperature(e_in)
                        p_in = r_in * T_in

                        u_gh =         - u_in
                        v_gh = - v_in
                        w_gh =         - w_in
                        T_gh = max(2*Twall - T_in, T_floor*Twall)
                        p_gh = p_in
                        r_gh = p_gh/T_gh

                        phi(ig,j,k,1) = r_gh
                        phi(ig,j,k,2) = r_gh*u_gh
                        phi(ig,j,k,3) = r_gh*v_gh
                        phi(ig,j,k,4) = r_gh*w_gh
                        phi(ig,j,k,5) = r_gh*enthalpy_t(T_gh) - p_gh &
                                      + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)
                     enddo
                  enddo
               enddo
               !$acc end parallel
            else
               !$acc parallel default(present)
               !$acc loop gang, vector collapse(3)
               do        k = lbz,ubz
                  do     j = lby,uby
                     do  i = 1,GN

                        ig = node + norm*i
                        ii = node - norm*(i-1)

                        r_in = phi(ii,j,k,1)
                        irin = 1._rp/r_in
                        u_in = phi(ii,j,k,2) * irin
                        v_in = phi(ii,j,k,3) * irin
                        w_in = phi(ii,j,k,4) * irin
                        p_in = gm1*(phi(ii,j,k,5) &
                             - 0.5_rp*r_in*(u_in*u_in + v_in*v_in + w_in*w_in))
                        T_in = p_in*irin

                        u_gh =         - u_in
                        v_gh = - v_in
                        w_gh =         - w_in
                        T_gh = max(2*Twall - T_in, T_floor*Twall)
                        p_gh = p_in
                        r_gh = p_gh/T_gh

                        phi(ig,j,k,1) = r_gh
                        phi(ig,j,k,2) = r_gh*u_gh
                        phi(ig,j,k,3) = r_gh*v_gh
                        phi(ig,j,k,4) = r_gh*w_gh
                        phi(ig,j,k,5) = p_gh*igm1 &
                                      + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)
                     enddo
                  enddo
               enddo
               !$acc end parallel
            endif

          case('S','N')

            if (thermally_perfect_gas > 0) then
               !$acc parallel default(present)
               !$acc loop gang, vector collapse(3)
               do        k = lbz,ubz
                  do     j = 1,GN
                     do  i = lbx,ubx

                        jg = node + norm*j
                        ji = node - norm*(j-1)

                        r_in = phi(i,ji,k,1)
                        irin = 1._rp/r_in
                        u_in = phi(i,ji,k,2) * irin
                        v_in = phi(i,ji,k,3) * irin
                        w_in = phi(i,ji,k,4) * irin
                        e_in = phi(i,ji,k,5) * irin - 0.5_rp*(u_in*u_in+v_in*v_in+w_in*w_in)
                        T_in = inverse_temperature(e_in)
                        p_in = r_in * T_in

                        u_gh = - u_in
                        v_gh =         - v_in
                        w_gh =         - w_in
                        T_gh = max(2*Twall - T_in, T_floor*Twall)
                        p_gh = p_in
                        r_gh = p_gh/T_gh

                        phi(i,jg,k,1) = r_gh
                        phi(i,jg,k,2) = r_gh*u_gh
                        phi(i,jg,k,3) = r_gh*v_gh
                        phi(i,jg,k,4) = r_gh*w_gh
                        phi(i,jg,k,5) = r_gh*enthalpy_t(T_gh) - p_gh &
                                      + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)
                     enddo
                  enddo
               enddo
               !$acc end parallel
            else
               !$acc parallel default(present)
               !$acc loop gang, vector collapse(3)
               do        k = lbz,ubz
                  do     j = 1,GN
                     do  i = lbx,ubx

                        jg = node + norm*j
                        ji = node - norm*(j-1)

                        r_in = phi(i,ji,k,1)
                        irin = 1._rp/r_in
                        u_in = phi(i,ji,k,2) * irin
                        v_in = phi(i,ji,k,3) * irin
                        w_in = phi(i,ji,k,4) * irin
                        p_in = gm1*(phi(i,ji,k,5) &
                             - 0.5_rp*r_in*(u_in*u_in + v_in*v_in + w_in*w_in))
                        T_in = p_in*irin

                        u_gh = - u_in
                        v_gh =         - v_in
                        w_gh =         - w_in
                        T_gh = max(2*Twall - T_in, T_floor*Twall)
                        p_gh = p_in
                        r_gh = p_gh/T_gh

                        phi(i,jg,k,1) = r_gh
                        phi(i,jg,k,2) = r_gh*u_gh
                        phi(i,jg,k,3) = r_gh*v_gh
                        phi(i,jg,k,4) = r_gh*w_gh
                        phi(i,jg,k,5) = p_gh*igm1 &
                                      + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)
                     enddo
                  enddo
               enddo
               !$acc end parallel
            endif

          case('B','F')

            if (thermally_perfect_gas > 0) then
               !$acc parallel default(present)
               !$acc loop gang, vector collapse(3)
               do        k = 1,GN
                  do     j = lby,uby
                     do  i = lbx,ubx

                        kg = node + norm*k
                        ki = node - norm*(k-1)

                        r_in = phi(i,j,ki,1)
                        irin = 1._rp/r_in
                        u_in = phi(i,j,ki,2) * irin
                        v_in = phi(i,j,ki,3) * irin
                        w_in = phi(i,j,ki,4) * irin
                        e_in = phi(i,j,ki,5) * irin - 0.5_rp*(u_in*u_in+v_in*v_in+w_in*w_in)
                        T_in = inverse_temperature(e_in)
                        p_in = r_in * T_in

                        u_gh = - u_in
                        v_gh =         - v_in
                        w_gh =         - w_in
                        T_gh = max(2*Twall - T_in, T_floor*Twall)
                        p_gh = p_in
                        r_gh = p_gh/T_gh

                        phi(i,j,kg,1) = r_gh
                        phi(i,j,kg,2) = r_gh*u_gh
                        phi(i,j,kg,3) = r_gh*v_gh
                        phi(i,j,kg,4) = r_gh*w_gh
                        phi(i,j,kg,5) = r_gh*enthalpy_t(T_gh) - p_gh &
                                      + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)
                     enddo
                  enddo
               enddo
               !$acc end parallel
            else
               !$acc parallel default(present)
               !$acc loop gang, vector collapse(3)
               do        k = 1,GN
                  do     j = lby,uby
                     do  i = lbx,ubx

                        kg = node + norm*k
                        ki = node - norm*(k-1)

                        r_in = phi(i,j,ki,1)
                        irin = 1._rp/r_in
                        u_in = phi(i,j,ki,2) * irin
                        v_in = phi(i,j,ki,3) * irin
                        w_in = phi(i,j,ki,4) * irin
                        p_in = gm1*(phi(i,j,ki,5) &
                             - 0.5_rp*r_in*(u_in*u_in + v_in*v_in + w_in*w_in))
                        T_in = p_in*irin

                        u_gh = - u_in
                        v_gh =         - v_in
                        w_gh =         - w_in
                        T_gh = max(2*Twall - T_in, T_floor*Twall)
                        p_gh = p_in
                        r_gh = p_gh/T_gh

                        phi(i,j,kg,1) = r_gh
                        phi(i,j,kg,2) = r_gh*u_gh
                        phi(i,j,kg,3) = r_gh*v_gh
                        phi(i,j,kg,4) = r_gh*w_gh
                        phi(i,j,kg,5) = p_gh*igm1 &
                                      + 0.5_rp*r_gh*(u_gh*u_gh + v_gh*v_gh + w_gh*w_gh)
                     enddo
                  enddo
               enddo
               !$acc end parallel
            endif

        end select

end subroutine isothermal_wall

subroutine extrapolation(phi,b)
        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: phi
        type(face_type)                          , intent(in)    :: b

        integer :: i,j,k,l
        integer :: node, norm

        node = b%node
        norm = b%norm

        select case(b%face)

          case('E','W')

            !$acc parallel default(present)
            !$acc loop gang, vector collapse(4)
            do          l = 1,5
               do       k = lbz,ubz
                  do    j = lby,uby
                     do i = 1,GN
                        phi(node+norm*i,j,k,l) = phi(node,j,k,l)
                     enddo
                  enddo
               enddo
            enddo
            !$acc end parallel

          case('S','N')

            !$acc parallel default(present)
            !$acc loop gang, vector collapse(4)
            do          l = 1,5
               do       k = lbz,ubz
                  do    j = 1,GN
                     do i = lbx,ubx
                        phi(i,node+norm*j,k,l) = phi(i,node,k,l)
                     enddo
                  enddo
               enddo
            enddo
            !$acc end parallel

          case('B','F')

            !$acc parallel default(present)
            !$acc loop gang, vector collapse(4)
            do          l = 1,5
               do       k = 1,GN
                  do    j = lby,uby
                     do i = lbx,ubx
                        phi(i,j,node+norm*k,l) = phi(i,j,node,l)
                     enddo
                  enddo
               enddo
            enddo
            !$acc end parallel

          case default
            print*, ' extrapolation is not implemented for face ', trim(b%face)
            stop

        end select

        return
end subroutine extrapolation

end module bc_module
