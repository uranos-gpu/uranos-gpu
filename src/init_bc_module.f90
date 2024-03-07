module init_bc_module
use parameters_module
use storage_module
use ic_module
use random_module
implicit none

public init_boundary_conditions

contains
subroutine init_boundary_conditions
        
        use bc_module      , only: set_bc_conditions, mpi_bc_communications
        use rhs_module     , only: update_all, update_all_ghosts, mpi_wait_conservatives

        implicit none
        integer :: err = 0
        integer :: i,j,k
        real(rp) :: r_, u_, v_, w_, p_

        if(len_trim(restart_file).eq.0) then
                
                ! ambient
                if    (ic == 'ambient'        .or. &
                       ic == 'hTurb'          .or. &
                       ic == 'KolmogorovFlow' .or. &
                       ic == 'blades'         .or. &
                       ic == 'smooth_body')  then
                        call ambient(err)

                ! Periodic Euler
                elseif(ic == 'periodic_euler_x') then
                        call init_periodic_euler('x',err)
                elseif(ic == 'periodic_euler_y') then
                        call init_periodic_euler('y',err)
                elseif(ic == 'periodic_euler_z') then
                        call init_periodic_euler('z',err)
                
                ! shock tube
                elseif(ic == 'shock_tube_x') then
                        call init_shock_tube('x',err)
                elseif(ic == 'shock_tube_y') then
                        call init_shock_tube('y',err)
                elseif(ic == 'shock_tube_z') then
                        call init_shock_tube('z',err)

                ! Lax shock tube
                elseif(ic == 'lax_problem_x') then
                        call init_lax_problem('x',err)
                elseif(ic == 'lax_problem_y') then
                        call init_lax_problem('y',err)
                elseif(ic == 'lax_problem_z') then
                        call init_lax_problem('z',err)

                ! Shock-wave interaction
                elseif(ic == 'shock_wave_interaction_x') then
                        call init_shock_wave_interaction('x',err)
                elseif(ic == 'shock_wave_interaction_y') then
                        call init_shock_wave_interaction('y',err)
                elseif(ic == 'shock_wave_interaction_z') then
                        call init_shock_wave_interaction('z',err)

                ! isentropic vortex
                elseif(ic == 'isentropic_vortex_x') then
                        call init_isentropic_vortex('x',err)
                elseif(ic == 'isentropic_vortex_y') then
                        call init_isentropic_vortex('y',err)
                elseif(ic == 'pirozzoli_vortex') then
                        call init_pirozzoli_vortex(err)
                elseif(ic == 'pirozzoli_vortex_y') then
                        call init_pirozzoli_vortex_y(err)

                ! Shock-wave
                elseif(ic == 'shock_wave_x') then
                        call init_shock_wave('x',err)
                elseif(ic == 'shock_wave_y') then
                        call init_shock_wave('y',err)
                elseif(ic == 'shock_wave_z') then
                        call init_shock_wave('z',err)
                elseif(ic == 'shock_inflow') then
                        call init_shock_inflow(err)
                elseif(ic == 'shock_impact') then
                        call init_shock_impact(err)

                ! double mach reflection
                elseif(ic == 'double_mach_reflection') then
                        call init_double_mach_reflection(err)

                elseif(ic == 'four_quadrants_A' .or.  &
                       ic == 'four_quadrants_B'     ) then
                        call init_four_quadrants(err)

                ! couette and Poiseuille flow
                elseif(ic == 'couette_x' .or. &
                       ic == 'couette_y' .or. &
                       ic == 'couette_z' .or. &
                       ic == 'poiseuille_x' .or. &
                       ic == 'poiseuille_y' .or. &
                       ic == 'poiseuille_z' .or. &
                       ic == 'inflow_poiseuille') then
                        call ambient(err)
        
                ! steady couette initialization
                elseif(ic == 'steady_couette_x') then
                        call init_steady_couette('x',err)
                elseif(ic == 'steady_couette_y') then
                        call init_steady_couette('y',err)
                elseif(ic == 'steady_couette_z') then
                        call init_steady_couette('z',err)

                ! I stokes problems
                elseif(ic == 'I_stokes_problem_x') then
                        call init_I_stokes_problem('x',err)
                elseif(ic == 'I_stokes_problem_y') then
                        call init_I_stokes_problem('y',err)
                elseif(ic == 'I_stokes_problem_z') then
                        call init_I_stokes_problem('z',err)
        
                ! II stokes problem
                elseif(ic == 'II_stokes_problem_x') then
                        call init_II_stokes_problem('x',err)
                elseif(ic == 'II_stokes_problem_y') then
                        call init_II_stokes_problem('y',err)
                elseif(ic == 'II_stokes_problem_z') then
                        call init_II_stokes_problem('z',err)

                elseif(ic == 'shock_vortex') then
                        call init_shock_vortex(err)

                elseif(ic == 'turbulent_channel') then
                        !call init_turbulent_channel(err)
                        call init_turbulent_chennel_random(err)

                elseif(ic == 'ambient_supersonic') then
                        call init_ambient_supersonic(err)
                
                ! shear layer
                elseif(ic == 'shear_layer') then
                        call init_shear_layer(err)
                elseif(ic == 'piecewise_shear_layer') then
                        call init_piecewise_shear_layer(err)

                elseif(ic == 'shock_bubble_2D') then
                        call init_shock_bubble_2D(err)

                elseif(ic == 'nscbc_perturbation_x') then
                        call init_nscbc_perturbation(err,'x')
                elseif(ic == 'nscbc_perturbation_y') then
                        call init_nscbc_perturbation(err,'y')

                elseif(ic == 'blasius') then
                        call ambient(err)

                elseif(ic == 'supersonic_intake') then
                        call ambient(err)

                elseif(ic == 'swbli') then
                        call init_swbli(err)
                        !call init_LaminarBoundaryLayer(err)

                elseif(ic == 'turbulent_BL') then
                        call init_TurbulentBoundaryLayerMean(err)
                        call init_TBLNoise(err)

                ! linear ODE y'+y = 0 initialization
                elseif(ic == 'linear_ode') then
                        call init_linear_ode(err)

                ! linear advection u_t + u_x = 0 initialization
                elseif(ic == 'linear_advection') then
                        call init_linear_advection(err)

                ! rectangular cylinder (barc) initialization
                elseif(ic == 'barc') then
                        call ambient(err)

                else
                        if(rank == root) print*, ' Initial boundary condition ', "'"//trim(ic)//"'", ' is not implemented.'
                        call secure_stop
                endif
                
                ! initialize conservative variable
                do       k = lbz,ubz
                   do    j = lby,uby
                      do i = lbx,ubx
                           
                           r_ = phi(i,j,k,1)
                           u_ = U(i,j,k)
                           v_ = V(i,j,k)
                           w_ = W(i,j,k)
                           p_ = P(i,j,k)

                           phi(i,j,k,2) = r_ * u_
                           phi(i,j,k,3) = r_ * v_
                           phi(i,j,k,4) = r_ * w_
                           phi(i,j,k,5) = p_/(gamma0-1._rp)+0.5_rp*r_*(u_*u_+v_*v_+w_*w_)

                           T(i,j,k) = p_/r_

                      enddo
                   enddo
               enddo

               if(viscous) then
                 do       k = lbz,ubz
                    do    j = lby,uby
                       do i = lbx,ubx

                          VIS(i,j,k) = laminar_viscosity(T(i,j,k),tref,vis_flag)
                          LMD(i,j,k) = k_inf * VIS(i,j,k)

                       enddo
                    enddo
                 enddo
               endif

        ! ============= init from file ============== !
        else 
          restart_flag = .true. 

          if(trim(restart_file) == 'restart') call Getrestart(restart_file)

          ! read conservative field from a file
          call mpi_read(restart_file,err)
        
          !$acc data copy(phi,U,V,W,P,T,VIS,LMD) &
          !$acc copyin(phi_bfr_send_E, phi_bfr_send_W, phi_bfr_recv_E, phi_bfr_recv_W) &
          !$acc copyin(phi_bfr_send_N, phi_bfr_send_S, phi_bfr_recv_N, phi_bfr_recv_S) &
          !$acc copyin(phi_bfr_send_B, phi_bfr_send_F, phi_bfr_recv_B, phi_bfr_recv_F) &
          !$acc copyin(bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W) &
          !$acc copyin(bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S) &
          !$acc copyin(bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F) &
          !$acc copyin(my_neighbour,bc)

          ! shared mpi boundary
          call mpi_bc_communications
          ! apply istantaneous boundary conditions
          call set_bc_conditions

          ! update all the other fields
          call update_all

          ! compute time from restart
          time_from_restart = time

          restart_flag = .false.

          !$acc end data

          if(restartStat) call init_previous_stats(restart_file)

        endif

        if(err /= 0) then
          if(rank == root) print*, ' Initialization error. Program exits with', err
          call secure_stop
        endif

        return
end subroutine init_boundary_conditions



subroutine GetRestart(restart_file)

        use FileModule

        implicit none
        character(dl), intent(out) :: restart_file

        character(dl) :: temp_file
        integer       :: f_unit = 10, err = 0

        ! save ls in a temporary file
        temp_file = 'DATA/'//trim(data_dir)//'/contents'//trim(str(rank))//'.txt'
        call execute_command_line('ls DATA/'//trim(data_dir)//'/BINARY/ > '//temp_file)

        open(unit = f_unit,file = temp_file, action = "read", access = "sequential", iostat = err)
        if(err .ne. 0) stop 'error in GetRestart'

        do
         read(f_unit,FMT='(a)',iostat=err) restart_file
         if (err/=0) EXIT
        end do

        restart_file = 'DATA/'//trim(data_dir)//'/BINARY/'//trim(restart_file)

        call execute_command_line('rm -f '//trim(temp_file))
        close(f_unit)

        call MPI_BARRIER(mpi_comm_world,err)

        return
end subroutine GetRestart



subroutine mpi_read(filename,err)
        implicit none
        character(*), intent(in)  :: filename
        integer     , intent(out) :: err

        ! local declarations
        integer                             :: fh
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(kind = MPI_OFFSET_KIND)     :: initial_displacement
        integer, parameter                  :: array_rank=3
        integer, dimension(array_rank)      :: shape_array, shape_sub_array, start_coord
        integer, dimension(array_rank)      :: shape_view_array, shape_sub_view_array, start_view_coord
        integer                             :: type_sub_array = 0, type_sub_view_array = 0

        CALL MPI_FILE_OPEN(mpi_comm_cart, filename, &
                 MPI_MODE_RDONLY                  , &
                 MPI_INFO_NULL, fh, mpi_err)

        if(mpi_err .ne. 0) then
                if(rank == root) write(*,'(1x,A,A,A)') ' The file ', trim(filename), ' does not exists.'
          call secure_stop
        endif

        if(rank == root) then
          write(*,'(1x,A,A)') ' Reading initial boundary condition from ', "'"//trim(filename)//"'"
        endif

        ! subarrays shape
        shape_array     = SHAPE(phi(:,:,:,1))
        shape_sub_array = SHAPE(phi(sx:ex, sy:ey, sz:ez,1))
        start_coord     = (/GN, GN, GN/)

        call mpi_type_create_subarray(array_rank, shape_array, shape_sub_array, start_coord, &
                                      MPI_ORDER_FORTRAN, MPI_RP              , &
                                      type_sub_array, mpi_err)

        call mpi_type_commit(type_sub_array, mpi_err)

        ! view array shape
        shape_view_array     = (/nx, ny, nz/)
        shape_sub_view_array = SHAPE(phi(sx:ex, sy:ey, sz:ez,1))
        start_view_coord     =  (/sx-1 , sy-1 , sz-1/)

        call mpi_type_create_subarray(array_rank, shape_view_array, shape_sub_view_array, start_view_coord, &
                                      MPI_ORDER_FORTRAN, MPI_RP                             , &
                                      type_sub_view_array, mpi_err)
        call mpi_type_commit(type_sub_view_array, mpi_err)

        CALL MPI_FILE_READ_ALL(fh, nx , 1, MPI_INTEGER, status, mpi_err)
        CALL MPI_FILE_READ_ALL(fh, ny , 1, MPI_INTEGER, status, mpi_err)
        CALL MPI_FILE_READ_ALL(fh, nz , 1, MPI_INTEGER, status, mpi_err)

        CALL MPI_FILE_READ_ALL(fh, it   , 1, MPI_INTEGER, status, mpi_err)
        CALL MPI_FILE_READ_ALL(fh, time , 1, MPI_RP, status, mpi_err)
        CALL MPI_FILE_READ_ALL(fh, dt   , 1, MPI_RP, status, mpi_err)

        initial_displacement = 4*ip + 2*rp
        CALL MPI_FILE_SET_VIEW(fh, initial_displacement, MPI_RP, &
             type_sub_view_array, "native", MPI_INFO_NULL, mpi_err)

        CALL MPI_FILE_READ_ALL(fh, phi(:,:,:,1), 1, type_sub_array, status, mpi_err)
        CALL MPI_FILE_READ_ALL(fh, phi(:,:,:,2), 1, type_sub_array, status, mpi_err)
        CALL MPI_FILE_READ_ALL(fh, phi(:,:,:,3), 1, type_sub_array, status, mpi_err)
        CALL MPI_FILE_READ_ALL(fh, phi(:,:,:,4), 1, type_sub_array, status, mpi_err)
        CALL MPI_FILE_READ_ALL(fh, phi(:,:,:,5), 1, type_sub_array, status, mpi_err)

        CALL MPI_FILE_CLOSE(fh, mpi_err)
        CALL check_mpi('mpi_read', mpi_err)

        ! check for tmax
        if(time >= tmax) then
          if(rank == root) write(*,'(1x,A)') ' The restarting time is greater than total simulation time.'
          call secure_stop
        endif
        
        err = 0
        return
end subroutine mpi_read



subroutine init_previous_stats(restart_file)
        implicit none
        character(*), intent(in) :: restart_file
        character(dl) :: iteration, restartStatFile, restartStatWmles

        integer, dimension(4) :: lenght
        integer :: lenghtFile, sum_len
        
        !
        ! getting the statitics file corresponding to the restart binary
        !
        write(iteration,'(i7.7)') it
        lenghtFile = len_trim(restart_file)
        lenght(1) = 4   !.bin
        lenght(2) = 7   ! iteration
        lenght(3) = len_trim(output_file_name)
        lenght(4) = 7   ! BINARY
        sum_len = sum(lenght)
        
        restartStatFile = restart_file(1:lenghtFile-sum_len)//'BINARY_STAT/stats'//&
                          trim(iteration)//'.bin'

        selectcase(ic)
        case('turbulent_channel','poiseuille_x') 
          call mpi_read_stat1D(restartStatFile,sy,ey,ny,nvAve1D,vmean1D)

        case default
          call mpi_read_stat2D(restartStatFile,sx,ex,sy,ey,nx,ny,vmean2D)

          if(wmles) then
           restartStatWmles = restart_file(1:lenghtFile-sum_len)//&
                              'BINARY_WMLES_STAT/stats'//trim(iteration)//'.bin'
           call mpi_read_stat1D(restartStatWmles,sx,ex,nx,nVWMLESData,vmean1D_wmles)
          endif

        endselect
        
        itStat = itStat - 1 


        return
end subroutine init_previous_stats


subroutine mpi_read_stat1D(filename,sid,eid,n,nStats,vmean)
        implicit none
        real(rp), dimension(:,:), allocatable, intent(inout) :: vmean
        character(*)                         , intent(in)    :: filename
        integer                              , intent(in)    :: sid, eid, n, nStats

        ! local declarations
        integer, parameter             :: array_rank = 1
        integer, parameter             :: intprec = MPI_INTEGER

        integer, dimension(array_rank) :: shape_array
        integer, dimension(array_rank) :: shape_sub_array
        integer, dimension(array_rank) :: start_coord

        integer, dimension(array_rank) :: shape_view_array
        integer, dimension(array_rank) :: shape_sub_view_array
        integer, dimension(array_rank) :: start_view_coord

        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(kind = MPI_OFFSET_KIND)     :: initial_displacement

        integer :: type_sub_array, type_sub_view_array
        integer :: fh, err = 0, l

        !
        ! === open statistic file
        !
        CALL MPI_FILE_OPEN(mpi_comm_cart, filename, &
                 MPI_MODE_RDONLY                  , &
                 MPI_INFO_NULL, fh, err)

        if(err.ne.0) then
          if(rank == root) write(*,'(1x,A)') ' No previous statistics found!'
          CALL MPI_FILE_CLOSE(fh, err)
          return
        endif

        if(rank == root) write(*,'(1x,A,A)') &
          ' Reading statistics from ', "'"//trim(filename)//"'"

        !
        ! === subarrays
        !
        shape_array     = SHAPE(vmean(:,1))
        shape_sub_array = SHAPE(vmean(sid:eid,1))
        start_coord     = (/GN/)

        call mpi_type_create_subarray(array_rank, shape_array, shape_sub_array, &
                                      start_coord, MPI_ORDER_FORTRAN, MPI_RP , &
                                      type_sub_array, err)
        call mpi_type_commit(type_sub_array, err)
        !
        ! === view array
        !
        shape_view_array     = (/n/)
        shape_sub_view_array = SHAPE(vmean(sid:eid,1))
        start_view_coord     =  (/sid-1/)

        call mpi_type_create_subarray(array_rank, shape_view_array, &
                                      shape_sub_view_array, start_view_coord, &
                                      MPI_ORDER_FORTRAN, MPI_RP, &
                                      type_sub_view_array, err)
        call mpi_type_commit(type_sub_view_array, err)
        !
        ! === reading the file
        !
        CALL MPI_FILE_READ_ALL(fh, itStat , 1, intprec, status, err)
        initial_displacement = 1*ip 

        CALL MPI_FILE_SET_VIEW(fh, initial_displacement, MPI_RP, &
             type_sub_view_array, "native", MPI_INFO_NULL, err)

        do l = 1,nStats
           CALL MPI_FILE_READ_ALL(fh, vmean(:,l), 1, type_sub_array, status, err)
        enddo

        CALL MPI_FILE_CLOSE(fh, err)


        return
end subroutine mpi_read_stat1D



subroutine mpi_read_stat2D(filename,sx,ex,sy,ey,nx,ny,vmean)
        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: vmean
        character(*)                           , intent(in)    :: filename
        integer                                , intent(in)    :: sx, ex, sy,ey, nx, ny

        ! local declarations
        integer, parameter             :: array_rank = 2
        integer, parameter             :: intprec = MPI_INTEGER

        integer, dimension(array_rank) :: shape_array
        integer, dimension(array_rank) :: shape_sub_array
        integer, dimension(array_rank) :: start_coord

        integer, dimension(array_rank) :: shape_view_array
        integer, dimension(array_rank) :: shape_sub_view_array
        integer, dimension(array_rank) :: start_view_coord

        integer, dimension(MPI_STATUS_SIZE) :: status
        integer(kind = MPI_OFFSET_KIND)     :: initial_displacement

        integer :: type_sub_array, type_sub_view_array
        integer :: fh, err = 0, l
        
        !
        ! === open statistic file
        !
        CALL MPI_FILE_OPEN(mpi_comm_cart, filename, &
                 MPI_MODE_RDONLY                  , &
                 MPI_INFO_NULL, fh, err)

        if(err.ne.0) then
          if(rank == root) write(*,'(1x,A)') ' No previous statistics found! '
          CALL MPI_FILE_CLOSE(fh, err)
          return
        endif

        if(rank == root) write(*,'(1x,A,A)') &
          ' Reading statistics from ', "'"//trim(filename)//"'"

        !
        ! === subarrays
        !
        shape_array     = SHAPE(vmean(:,:,1))
        shape_sub_array = SHAPE(vmean(sx:ex,sy:ey,1))
        start_coord     = (/GN,GN/)

        call mpi_type_create_subarray(array_rank, shape_array, shape_sub_array, &
                                      start_coord, MPI_ORDER_FORTRAN, MPI_RP  , &
                                      type_sub_array, err)
        call mpi_type_commit(type_sub_array, err)
        !
        ! === view array
        !
        shape_view_array     = (/nx, ny/)
        shape_sub_view_array = SHAPE(vmean(sx:ex, sy:ey,1))
        start_view_coord     =  (/sx-1, sy-1/)

        call mpi_type_create_subarray(array_rank, shape_view_array, &
                                      shape_sub_view_array, start_view_coord, &
                                      MPI_ORDER_FORTRAN, MPI_RP, &
                                      type_sub_view_array, err)
        call mpi_type_commit(type_sub_view_array, err)
        !
        ! === reading the file
        !
        CALL MPI_FILE_READ_ALL(fh, itStat , 1, intprec, status, err)
        initial_displacement = 1*ip 

        CALL MPI_FILE_SET_VIEW(fh, initial_displacement, MPI_RP, &
             type_sub_view_array, "native", MPI_INFO_NULL, err)

        do l = 1,nvAve2D
           CALL MPI_FILE_READ_ALL(fh, vmean(:,:,l), 1, type_sub_array, status, err)
        enddo

        CALL MPI_FILE_CLOSE(fh, err)


        return
end subroutine mpi_read_stat2D









end module init_bc_module
