module post_bc_module
use parameters_module
use storage_module
use post_storage_module
use bc_module

implicit none
contains
subroutine apply_pst_bc_conditions
        implicit none
        type(face_type), dimension(6) :: all_bound
        type(face_type)               :: bound
        real(rp)                      :: tWall, rfc
        integer                       :: f

        all_bound%node = (/ sx , ex , sy , ey , sz ,  ez/)
        all_bound%norm = (/ -1 ,  1 , -1 ,  1 , -1 ,  1 /)
        all_bound%face = (/ 'W', 'E', 'S', 'N', 'B', 'F'/)

        do f = 1,6
                
           ! store everything a scalar variable
           bound%node = all_bound(f)%node
           bound%norm = all_bound(f)%norm
           bound%face = all_bound(f)%face

           if    (trim(bc(f)) == 'periodic') then
             call periodic(phi,bound)

           elseif(trim(bc(f)) == 'neumann' .or. trim(bc(f)) == 'nscbc_outflow') then
             call neumann(phi,bound)

           elseif(trim(bc(f)) == 'adiabatic_wall' .or. trim(bc(f)) == 'nscbc_adiabatic_wall') then
             call adiabatic_wall(phi,bound)

           elseif(trim(bc(f)) == 'isothermal_wall'.or. &
                  trim(bc(f)) == 'nscbc_isothermal_wall') then
             tWall = 1.0_rp
             call isothermal_wall(phi,bound,tWall)

           elseif(trim(bc(f)) == 'slip_wall') then
             call slip_wall(phi,bound)

           elseif(trim(bc(f)) == 'postshock_slipwall') then
             call postshock_slipwall(phi,bound)

           elseif(trim(bc(f)) == 'oblique_shock') then
             call oblique_shock_bc(phi,bound)

           elseif(trim(bc(f)) == 'flate_plate') then
             call flate_plate_bc(phi,bound)

           elseif(trim(bc(f)) == 'travelling_shock_wave') then
             call travelling_shock_wave(phi,bound)

           elseif(trim(bc(f)) == 'subsonic_inflow'      .or. &
                  trim(bc(f)) == 'nscbc_inflow'         .or. &
                  trim(bc(f)) == 'nscbc_inflow_relaxed' .or. &
                  trim(bc(f)) == 'supersonic_inflow') then
             call neumann(phi,bound)

           elseif(trim(bc(f)) == 'shock_inflow') then
               call neumann(phi,bound)

           elseif(trim(bc(f)) == 'neumann_wall') then
               call neumann_wall(phi,bound)

           elseif(trim(bc(f)) == 'isothermal_wall'       .or. &
                  trim(bc(f)) == 'dws_isothermal'        .or. &
                  trim(bc(f)) == 'idws_isothermal'       .or. &
                  trim(bc(f)) == 'static_dws_isothermal') then
               tWall = 1.0_rp
               call isothermal_wall(phi,bound,tWall)

           elseif(trim(bc(f)) == 'dws_adiabatic'        .or. &
                  trim(bc(f)) == 'idws_adiabatic'       .or. &
                  trim(bc(f)) == 'static_dws_adiabatic' .or. &
                  trim(bc(f)) == 'istatic_dws_adiabatic') then

               rfc = Prandtl**(1.0_rp/3.0_rp)
               tWall = Trat*(1.0_rp + rfc*0.5_rp*(gamma0-1.0_rp)*Mach**2)
               call isothermal_wall(phi,bound,tWall)

           else
             write(*,*) ' Boundary condition ', "'"//trim(bc(f))//"'", ' is not implemented.'
             stop

           endif
        enddo

        return
end subroutine apply_pst_bc_conditions


subroutine periodic(phi,bound)
        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi
        type(face_type)                          , intent(in)    :: bound
        integer                                                  :: i,j,k

        select case(bound%face)
                case('E')

                  do i = 1,GN
                     phi(sx-i,:,:,:) = phi(ex+1-i,:,:,:)
                     phi(ex+i,:,:,:) = phi(sx-1+i,:,:,:)
                  enddo

                case('S')

                  do j=1,GN
                     phi(:,sy-j,:,:) = phi(:,ey+1-j,:,:)
                     phi(:,ey+j,:,:) = phi(:,sy-1+j,:,:)
                  end do

                case('B')

                  do k=1,GN
                     phi(:,:,sz-k,:) = phi(:,:,ez+1-k,:)
                     phi(:,:,ez+k,:) = phi(:,:,sz-1+k,:)
                  end do

        end select
        return
end subroutine periodic





end module post_bc_module
