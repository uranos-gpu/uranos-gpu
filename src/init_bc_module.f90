module init_bc_module
use parameters_module
use storage_module
use ic_module
use random_module
use thermodynamics_module
use mesh_module, only: deta_dx, deta_dy
implicit none

public init_boundary_conditions

contains
subroutine init_boundary_conditions

        use bc_module      , only: set_bc_conditions, mpi_bc_communications
        use rhs_module           , only: mpi_wait_conservatives
        use update_variable_module, only: update_all

        implicit none
        character(dl) :: restart_file_stat
        integer :: err = 0, exit_code
        integer :: i,j,k,l
        real(rp) :: r_, u_, v_, w_, p_, cp_

        if(len_trim(restart_file).eq.0) then

                if    (ic == 'ambient' .or. &
                       ic == 'KolmogorovFlow') then
                        call ambient(err)

                elseif(ic == 'TGV') then
                        call TGV(err)

                elseif(ic == 'periodic_euler_x') then
                        call init_periodic_euler('x',err)
                elseif(ic == 'periodic_euler_y') then
                        call init_periodic_euler('y',err)
                elseif(ic == 'periodic_euler_z') then
                        call init_periodic_euler('z',err)

                elseif(ic == 'shock_tube_x') then
                        call init_shock_tube('x',err)
                elseif(ic == 'shock_tube_y') then
                        call init_shock_tube('y',err)
                elseif(ic == 'shock_tube_z') then
                        call init_shock_tube('z',err)

                elseif(ic == 'isentropic_vortex_x') then
                        call init_isentropic_vortex('x',err)
                elseif(ic == 'isentropic_vortex_y') then
                        call init_isentropic_vortex('y',err)

                elseif(ic == 'turbulent_channel') then
                        call init_turbulent_channel(err)

                elseif(ic == 'turbulent_BL') then

                        call init_TurbulentBoundaryLayerMean(err)
                        call init_TBLNoise_fast(err)

                else
                        if(rank == root) print*, ' Initial boundary condition ', "'"//trim(ic)//"'", ' is not implemented.'
                        call secure_stop
                endif

                if(thermally_perfect_gas > 0) then
                   do       k = lbz,ubz
                      do    j = lby,uby
                         do i = lbx,ubx

                              r_       = phi(i,j,k,1)
                              u_       = U(i,j,k)
                              v_       = V(i,j,k)
                              w_       = W(i,j,k)
                              p_       = P(i,j,k)

                              T(i,j,k)             = p_/r_

                              cp_                  = heat_capacity_t(T(i,j,k))
                              ENTHALPY(i,j,k)      = enthalpy_t(T(i,j,k))
                              GAMMA_T(i,j,k)       = cp_ / (cp_ - 1.0_rp)

                              phi(i,j,k,2) = r_ * u_
                              phi(i,j,k,3) = r_ * v_
                              phi(i,j,k,4) = r_ * w_
                              phi(i,j,k,5) = r_ * ENTHALPY(i,j,k) - p_ +0.5_rp*r_*(u_*u_+v_*v_+w_*w_)

                         enddo
                      enddo
                   enddo
                else
                   do       k = lbz,ubz
                      do    j = lby,uby
                         do i = lbx,ubx

                              r_       = phi(i,j,k,1)
                              u_       = U(i,j,k)
                              v_       = V(i,j,k)
                              w_       = W(i,j,k)
                              p_       = P(i,j,k)

                              T(i,j,k)             = p_/r_

                              ENTHALPY(i,j,k) = cp0 * T(i,j,k)
                              GAMMA_T(i,j,k)       = gamma0

                              phi(i,j,k,2) = r_ * u_
                              phi(i,j,k,3) = r_ * v_
                              phi(i,j,k,4) = r_ * w_
                              phi(i,j,k,5) = r_ * ENTHALPY(i,j,k) - p_ +0.5_rp*r_*(u_*u_+v_*v_+w_*w_)

                         enddo
                      enddo
                   enddo
               endif

               if(viscous) then
                 do       k = lbz,ubz
                    do    j = lby,uby
                       do i = lbx,ubx

                          VIS(i,j,k) = laminar_viscosity(T(i,j,k),tref,vis_flag)
                          if(thermally_perfect_gas > 0) then
                             LMD(i,j,k) = cp_local(T(i,j,k)) * iPrandtl * VIS(i,j,k)
                          else
                             LMD(i,j,k) = k_inf * VIS(i,j,k)
                          endif

                       enddo
                    enddo
                 enddo
               endif

        else
          restart_flag = .true.

          if(trim(restart_file) == 'restart') call GetRestart(restart_file)

          call mpi_read(restart_file,err)

          if(dims == 2) then
            do l = 1,5
               do k = lbz,ubz
                  if(k /= sz) phi(:,:,k,l) = phi(:,:,sz,l)
               enddo
            enddo
          endif

          !$acc data copy(phi,U,V,W,P,T,VIS,LMD) &
          !$acc copy(GAMMA_T, ENTHALPY) &
          !$acc copyin(phi_bfr_send_E, phi_bfr_send_W, phi_bfr_recv_E, phi_bfr_recv_W) &
          !$acc copyin(phi_bfr_send_N, phi_bfr_send_S, phi_bfr_recv_N, phi_bfr_recv_S) &
          !$acc copyin(phi_bfr_send_B, phi_bfr_send_F, phi_bfr_recv_B, phi_bfr_recv_F) &
          !$acc copyin(bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W) &
          !$acc copyin(bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S) &
          !$acc copyin(bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F) &
          !$acc copyin(my_neighbour,bc,x,y,z) &
          !$acc copyin(deta_dx,deta_dy)

          call mpi_bc_communications
          call set_bc_conditions

          call update_all

          time_from_restart = time

          restart_flag = .false.

          !$acc end data

          if(restartStat) then
            call GetRestartStat(restart_file_stat,exit_code)
            call init_previous_stats(restart_file_stat,exit_code)
          endif

        endif

        if(err /= 0) then
          if(rank == root) print*, ' Initialization error. Program exits with', err
          call secure_stop
        endif

        return
end subroutine init_boundary_conditions

subroutine GetRestart(restart_file)
    use FileModule
    use mpi

    implicit none

    character(dl), intent(out) :: restart_file

    character(dl) :: temp_file, line
    integer       :: f_unit   = 10
    integer       :: err      = 0
    integer       :: ierr_mpi = 0
    integer       :: cmdstat  = 0
    logical       :: have_file

    restart_file = ''
    have_file    = .false.

    if (rank == 0) then
        temp_file = 'DATA/'//trim(data_dir)//'/contents_tmp.txt'

        call execute_command_line( &
            'rm -f DATA/'//trim(data_dir)//'/BINARY/*.lock 2> /dev/null', &
            cmdstat = cmdstat )

        call execute_command_line( &
            'ls -1v DATA/'//trim(data_dir)//'/BINARY 2> /dev/null | grep "\.bin$" > '//trim(temp_file), &
            cmdstat = cmdstat )

        if (cmdstat == 0) then
            open(unit=f_unit, file=temp_file, action="read", status="old", &
                 access="sequential", iostat=err)

            if (err == 0) then
                do
                    read(f_unit,'(a)',iostat=err) line
                    if (err /= 0) exit
                    if (len_trim(line) > 0) then
                        restart_file = trim(line)
                        have_file = .true.
                    end if
                end do
                close(f_unit)
            end if
        end if

        call execute_command_line('rm -f '//trim(temp_file), cmdstat = cmdstat)

        if (.not. have_file) restart_file = ''
    end if

    call MPI_BCAST(restart_file, len(restart_file), MPI_CHARACTER, &
                   0, mpi_comm_world, ierr_mpi)

    call MPI_BARRIER(mpi_comm_world, ierr_mpi)

    if (len_trim(restart_file) > 0) then
        restart_file = 'DATA/'//trim(data_dir)//'/BINARY/'//trim(restart_file)
    end if

end subroutine GetRestart

subroutine GetRestartStat(restart_file_stat, exit_code)
    use mpi
    implicit none
    character(dl), intent(out) :: restart_file_stat
    integer      , intent(out) :: exit_code

    character(dl)        :: stat_dir, tmp_file, line
    integer, allocatable :: k_n(:), k_it(:)
    integer              :: funit, ios, n_s, it_s, nl, k, l, ierr
    logical              :: exists

    restart_file_stat = ''
    exit_code         = 1

    if(rank == root) then
      stat_dir = 'DATA/'//trim(data_dir)//'/BINARY_STAT'
      inquire(file=trim(stat_dir)//'/stats_index.txt', exist=exists)

      if(exists) then
        open(newunit=funit, file=trim(stat_dir)//'/stats_index.txt', action='read', status='old')
        nl = 0
        do
           read(funit,*,iostat=ios) n_s, it_s
           if(ios /= 0) exit
           nl = nl + 1
        enddo
        allocate(k_n(nl), k_it(nl))
        rewind(funit)
        k = 0
        do l = 1,nl
           read(funit,*) n_s, it_s
           do while(k > 0)
              if(k_it(k) < it_s) exit
              k = k - 1
           enddo
           k = k + 1
           k_n(k) = n_s; k_it(k) = it_s
        enddo
        close(funit)
        do while(k > 0)
           if(k_it(k) <= it) exit
           k = k - 1
        enddo
        if(k > 0) then
          write(restart_file_stat,'(A,I0.7,A)') 'stats', k_n(k), '.bin'
          exit_code = 0
        else
          exit_code = 2
        endif
        deallocate(k_n, k_it)
      else
        tmp_file = 'DATA/'//trim(data_dir)//'/contents_stat.txt'
        call execute_command_line('ls -1v '//trim(stat_dir)//' 2>/dev/null | grep -E "^stats[0-9]+\.bin$" > ' &
                                  //trim(tmp_file), exitstat=ios)
        open(newunit=funit, file=trim(tmp_file), action='read', status='old', iostat=ios)
        if(ios == 0) then
          do
             read(funit,'(A)',iostat=ios) line
             if(ios /= 0) exit
             if(len_trim(line) > 0) restart_file_stat = trim(line)
          enddo
          close(funit, status='delete')
        endif
        if(len_trim(restart_file_stat) > 0) then
          exit_code = 0
          write(*,'(1x,3A)') ' WARNING: no stats_index.txt, statistics restarted from ', trim(restart_file_stat), &
                             ' without checking it against the restart iteration'
        endif
      endif
    endif

    call MPI_BCAST(exit_code, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(restart_file_stat, len(restart_file_stat), MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)

    if(exit_code == 2) then
      if(rank == root) write(*,'(1x,A,I0,A)') ' ERROR: no statistics sample at or before the restart iteration ', &
                                              it, '. Run with restartStat = .false.'
      call secure_stop
    endif

    return
end subroutine GetRestartStat

subroutine mpi_read(filename,err)
        implicit none
        character(*), intent(in)  :: filename
        integer     , intent(out) :: err

        integer                             :: fh, ne
        integer                             :: int_size, real_size
        integer                             :: open_err_max, alr_err
        integer(kind = MPI_OFFSET_KIND)     :: initial_displacement
        integer, parameter                  :: array_rank=3
        integer, dimension(array_rank)      :: shape_array, shape_sub_array, start_coord
        integer, dimension(array_rank)      :: shape_view_array, shape_sub_view_array, start_view_coord
        integer                             :: type_sub_array, type_sub_view_array
        integer                             :: nx_file, ny_file, nz_file

        integer, parameter                  :: NIO = 14
        integer, dimension(NIO)             :: errs
        character(len=28), dimension(NIO)   :: lbls

        CALL MPI_FILE_OPEN(mpi_comm_cart, filename, &
                 MPI_MODE_RDONLY, MPI_INFO_NULL, fh, mpi_err)

        call MPI_Allreduce(mpi_err, open_err_max, 1, MPI_INTEGER, MPI_MAX, mpi_comm_cart, alr_err)

        if(open_err_max .ne. 0) then
                if(rank == root) write(*,'(1x,A,A,A)') ' The file ', trim(filename), ' does not exists.'
          call secure_stop
        endif

        if(rank == root) then
          write(*,'(1x,A,A)') ' Reading initial boundary condition from ', "'"//trim(filename)//"'"
        endif

        ne = 0

        shape_array     = SHAPE(phi(:,:,:,1))
        shape_sub_array = SHAPE(phi(sx:ex, sy:ey, sz:ez,1))
        start_coord     = (/GN, GN, GN/)

        ne=ne+1; lbls(ne)='type_create_subarray loc'
        call mpi_type_create_subarray(array_rank, shape_array, shape_sub_array, start_coord, &
                                      MPI_ORDER_FORTRAN, MPI_RP, type_sub_array, errs(ne))
        ne=ne+1; lbls(ne)='type_commit loc'
        call mpi_type_commit(type_sub_array, errs(ne))

        shape_view_array     = (/nx, ny, nz/)
        shape_sub_view_array = SHAPE(phi(sx:ex, sy:ey, sz:ez,1))
        start_view_coord     =  (/sx-1 , sy-1 , sz-1/)

        ne=ne+1; lbls(ne)='type_create_subarray view'
        call mpi_type_create_subarray(array_rank, shape_view_array, shape_sub_view_array, start_view_coord, &
                                      MPI_ORDER_FORTRAN, MPI_RP, type_sub_view_array, errs(ne))
        ne=ne+1; lbls(ne)='type_commit view'
        call mpi_type_commit(type_sub_view_array, errs(ne))

        ne=ne+1; lbls(ne)='read_all nx'
        CALL MPI_FILE_READ_ALL(fh, nx_file, 1, MPI_INTEGER, MPI_STATUS_IGNORE, errs(ne))
        ne=ne+1; lbls(ne)='read_all ny'
        CALL MPI_FILE_READ_ALL(fh, ny_file, 1, MPI_INTEGER, MPI_STATUS_IGNORE, errs(ne))
        ne=ne+1; lbls(ne)='read_all nz'
        CALL MPI_FILE_READ_ALL(fh, nz_file, 1, MPI_INTEGER, MPI_STATUS_IGNORE, errs(ne))

        call check_mpi_vec('mpi_read: header', errs(1:ne), lbls(1:ne))

        if(nx_file /= nx .or. ny_file /= ny .or. nz_file /= nz) then
          if(rank == root) then
            write(*,'(1x,A)')       ' Restart file grid does not match the input deck:'
            write(*,'(1x,A,3I8)')   '   file  (nx,ny,nz) = ', nx_file, ny_file, nz_file
            write(*,'(1x,A,3I8)')   '   input (nx,ny,nz) = ', nx     , ny     , nz
          endif
          call secure_stop
        endif

        ne = 0

        ne=ne+1; lbls(ne)='read_all it'
        CALL MPI_FILE_READ_ALL(fh, it,   1, MPI_INTEGER, MPI_STATUS_IGNORE, errs(ne))
        ne=ne+1; lbls(ne)='read_all time'
        CALL MPI_FILE_READ_ALL(fh, time, 1, MPI_RP,      MPI_STATUS_IGNORE, errs(ne))
        ne=ne+1; lbls(ne)='read_all dt'
        CALL MPI_FILE_READ_ALL(fh, dt,   1, MPI_RP,      MPI_STATUS_IGNORE, errs(ne))

        ne=ne+1; lbls(ne)='type_size MPI_INTEGER'
        call MPI_TYPE_SIZE(MPI_INTEGER, int_size,  errs(ne))
        ne=ne+1; lbls(ne)='type_size MPI_RP'
        call MPI_TYPE_SIZE(MPI_RP,      real_size, errs(ne))

        initial_displacement = 4*int_size + 2*real_size
        ne=ne+1; lbls(ne)='file_set_view'
        CALL MPI_FILE_SET_VIEW(fh, initial_displacement, MPI_RP, type_sub_view_array, "native", MPI_INFO_NULL, errs(ne))

        ne=ne+1; lbls(ne)='read_all phi(1)'
        CALL MPI_FILE_READ_ALL(fh, phi(:,:,:,1), 1, type_sub_array, MPI_STATUS_IGNORE, errs(ne))
        ne=ne+1; lbls(ne)='read_all phi(2)'
        CALL MPI_FILE_READ_ALL(fh, phi(:,:,:,2), 1, type_sub_array, MPI_STATUS_IGNORE, errs(ne))
        ne=ne+1; lbls(ne)='read_all phi(3)'
        CALL MPI_FILE_READ_ALL(fh, phi(:,:,:,3), 1, type_sub_array, MPI_STATUS_IGNORE, errs(ne))
        ne=ne+1; lbls(ne)='read_all phi(4)'
        CALL MPI_FILE_READ_ALL(fh, phi(:,:,:,4), 1, type_sub_array, MPI_STATUS_IGNORE, errs(ne))
        ne=ne+1; lbls(ne)='read_all phi(5)'
        CALL MPI_FILE_READ_ALL(fh, phi(:,:,:,5), 1, type_sub_array, MPI_STATUS_IGNORE, errs(ne))

        ne=ne+1; lbls(ne)='file_close'
        CALL MPI_FILE_CLOSE(fh, errs(ne))

        ne=ne+1; lbls(ne)='type_free loc'
        call MPI_TYPE_FREE(type_sub_array,      errs(ne))
        ne=ne+1; lbls(ne)='type_free view'
        call MPI_TYPE_FREE(type_sub_view_array, errs(ne))

        call check_mpi_vec('mpi_read: body', errs(1:ne), lbls(1:ne))

        if(time >= tmax) then
          if(rank == root) write(*,'(1x,A)') ' The restarting time is greater than total simulation time.'
          call secure_stop
        endif

        err = 0
        return
end subroutine mpi_read

subroutine init_previous_stats(restart_file_stat,code)

        use FileModule

        implicit none
        character(*), intent(in) :: restart_file_stat
        integer     , intent(in) :: code

        character(dl) :: statFile

        if(code == 0) then

           statFile = 'DATA/'//trim(data_dir)//'/BINARY_STAT/'//trim(restart_file_stat)

           selectcase(stat_kind)
           case('TFT')
             call mpi_read_stat1D(statFile,sy,ey,ny,nvAve1D,vmean1D)

           case('FFT')
             call mpi_read_stat2D(statFile,sx,ex,sy,ey,nx,ny,vmean2D)

           case('FTT')
             call mpi_read_stat2D(statFile,sx,ex,sy,ey,nx,ny,vmean2D)

           endselect

        else
          if(rank == root) write(*,'(1x,A)') ' No previous statistics found!'
          return

        endif

        return
end subroutine init_previous_stats

subroutine mpi_read_stat1D(filename,sid,eid,n,nStats,vmean)
        implicit none
        real(rp), dimension(:,:), allocatable, intent(inout) :: vmean
        character(*)                         , intent(in)    :: filename
        integer                              , intent(in)    :: sid, eid, n, nStats

        integer, parameter             :: array_rank = 1

        integer, dimension(array_rank) :: shape_array
        integer, dimension(array_rank) :: shape_sub_array
        integer, dimension(array_rank) :: start_coord

        integer, dimension(array_rank) :: shape_view_array
        integer, dimension(array_rank) :: shape_sub_view_array
        integer, dimension(array_rank) :: start_view_coord

        integer(kind = MPI_OFFSET_KIND)     :: initial_displacement

        integer :: type_sub_array, type_sub_view_array
        integer :: fh, err, l, int_size, real_size, ne, nOnDisk, nRead
        integer :: open_err_max, alr_err
        integer(kind = MPI_OFFSET_KIND)     :: fsize

        integer, parameter                  :: NIO = 13
        integer, dimension(NIO)             :: errs
        character(len=28), dimension(NIO)   :: lbls

        err = 0
        CALL MPI_FILE_OPEN(mpi_comm_cart, filename, &
                 MPI_MODE_RDONLY, MPI_INFO_NULL, fh, err)

        call MPI_Allreduce(err, open_err_max, 1, MPI_INTEGER, MPI_MAX, mpi_comm_cart, alr_err)

        if(open_err_max.ne.0) then
          if(rank == root) write(*,'(1x,A)') ' Error in reading previous statistics!'
          return
        endif

        if(rank == root) write(*,'(1x,A,A)') &
          ' Reading statistics from ', "'"//trim(filename)//"'"

        ne = 0

        shape_array     = SHAPE(vmean(:,1))
        shape_sub_array = SHAPE(vmean(sid:eid,1))
        start_coord     = (/GN/)

        ne=ne+1; lbls(ne)='type_create_subarray loc'
        call mpi_type_create_subarray(array_rank, shape_array, shape_sub_array, &
                                      start_coord, MPI_ORDER_FORTRAN, MPI_RP, type_sub_array, errs(ne))
        ne=ne+1; lbls(ne)='type_commit loc'
        call mpi_type_commit(type_sub_array, errs(ne))

        shape_view_array     = (/n/)
        shape_sub_view_array = SHAPE(vmean(sid:eid,1))
        start_view_coord     =  (/sid-1/)

        ne=ne+1; lbls(ne)='type_create_subarray view'
        call mpi_type_create_subarray(array_rank, shape_view_array, shape_sub_view_array, start_view_coord, &
                                      MPI_ORDER_FORTRAN, MPI_RP, type_sub_view_array, errs(ne))
        ne=ne+1; lbls(ne)='type_commit view'
        call mpi_type_commit(type_sub_view_array, errs(ne))

        ne=ne+1; lbls(ne)='read_all itStat'
        CALL MPI_FILE_READ_ALL(fh, itStat, 1, MPI_INTEGER, MPI_STATUS_IGNORE, errs(ne))

        ne=ne+1; lbls(ne)='type_size MPI_INTEGER'
        call MPI_TYPE_SIZE(MPI_INTEGER, int_size, errs(ne))
        ne=ne+1; lbls(ne)='type_size MPI_RP'
        call MPI_TYPE_SIZE(MPI_RP, real_size, errs(ne))

        ne=ne+1; lbls(ne)='file_get_size'
        call MPI_FILE_GET_SIZE(fh, fsize, errs(ne))
        nOnDisk = int((fsize - int(int_size,MPI_OFFSET_KIND)) /                       &
                      (int(n,MPI_OFFSET_KIND)*int(real_size,MPI_OFFSET_KIND)))
        if((nOnDisk /= nStats .and. .not.(nStats == nvAve1D .and. nOnDisk == 27)) .or. &
           fsize /= int(int_size,MPI_OFFSET_KIND) +                                   &
                    int(n,MPI_OFFSET_KIND)*int(nOnDisk,MPI_OFFSET_KIND)*int(real_size,MPI_OFFSET_KIND)) then
          if(rank == root) write(*,'(1x,4A)') ' ERROR: ', "'"//trim(filename)//"'", &
               ' is not a 1D statistics file of this case.', ' Run with restartStat = .false.'
          call secure_stop
        endif
        nRead = min(nOnDisk, nStats)
        if(nRead < nStats .and. rank == root)                                         &
          write(*,'(1x,A,I0,A,I0,A)') ' WARNING: statistics file holds ', nOnDisk,    &
               ' slots, this build expects ', nStats, ' -- the new slots start from zero and their'// &
               ' means are biased low (divided by all samples).'

        initial_displacement = int_size
        ne=ne+1; lbls(ne)='file_set_view'
        CALL MPI_FILE_SET_VIEW(fh, initial_displacement, MPI_RP, type_sub_view_array, "native", MPI_INFO_NULL, errs(ne))

        ne=ne+1; lbls(ne)='read_all slabs'
        errs(ne) = 0
        do l = 1,nRead
           CALL MPI_FILE_READ_ALL(fh, vmean(:,l), 1, type_sub_array, MPI_STATUS_IGNORE, err)
           errs(ne) = ior(errs(ne), abs(err))
        enddo

        ne=ne+1; lbls(ne)='file_close'
        CALL MPI_FILE_CLOSE(fh, errs(ne))

        ne=ne+1; lbls(ne)='type_free loc'
        call MPI_TYPE_FREE(type_sub_array,      errs(ne))
        ne=ne+1; lbls(ne)='type_free view'
        call MPI_TYPE_FREE(type_sub_view_array, errs(ne))

        call check_mpi_vec('mpi_read_stat1D', errs(1:ne), lbls(1:ne))

        return
end subroutine mpi_read_stat1D

subroutine mpi_read_stat2D(filename,sx,ex,sy,ey,nx,ny,vmean)
        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: vmean
        character(*)                           , intent(in)    :: filename
        integer                                , intent(in)    :: sx, ex, sy,ey, nx, ny

        integer, parameter             :: array_rank = 2

        integer, dimension(array_rank) :: shape_array
        integer, dimension(array_rank) :: shape_sub_array
        integer, dimension(array_rank) :: start_coord

        integer, dimension(array_rank) :: shape_view_array
        integer, dimension(array_rank) :: shape_sub_view_array
        integer, dimension(array_rank) :: start_view_coord

        integer(kind = MPI_OFFSET_KIND)     :: initial_displacement

        integer :: type_sub_array, type_sub_view_array
        integer :: fh, err, l, int_size, ne
        integer :: open_err_max, alr_err
        integer :: real_size, nOnDisk, nRead
        integer(kind = MPI_OFFSET_KIND)     :: fsize

        integer, parameter                  :: NIO = 14
        integer, dimension(NIO)             :: errs
        character(len=28), dimension(NIO)   :: lbls

        err = 0
        CALL MPI_FILE_OPEN(mpi_comm_cart, filename, &
                 MPI_MODE_RDONLY, MPI_INFO_NULL, fh, err)

        call MPI_Allreduce(err, open_err_max, 1, MPI_INTEGER, MPI_MAX, mpi_comm_cart, alr_err)

        if(open_err_max.ne.0) then
          if(rank == root) write(*,'(1x,A)') ' No previous statistics found! '
          return
        endif

        if(rank == root) write(*,'(1x,A,A)') &
          ' Reading statistics from ', "'"//trim(filename)//"'"

        ne = 0

        shape_array     = SHAPE(vmean(:,:,1))
        shape_sub_array = SHAPE(vmean(sx:ex,sy:ey,1))
        start_coord     = (/GN,GN/)

        ne=ne+1; lbls(ne)='type_create_subarray loc'
        call mpi_type_create_subarray(array_rank, shape_array, shape_sub_array, &
                                      start_coord, MPI_ORDER_FORTRAN, MPI_RP, type_sub_array, errs(ne))
        ne=ne+1; lbls(ne)='type_commit loc'
        call mpi_type_commit(type_sub_array, errs(ne))

        shape_view_array     = (/nx, ny/)
        shape_sub_view_array = SHAPE(vmean(sx:ex, sy:ey,1))
        start_view_coord     =  (/sx-1, sy-1/)

        ne=ne+1; lbls(ne)='type_create_subarray view'
        call mpi_type_create_subarray(array_rank, shape_view_array, shape_sub_view_array, start_view_coord, &
                                      MPI_ORDER_FORTRAN, MPI_RP, type_sub_view_array, errs(ne))
        ne=ne+1; lbls(ne)='type_commit view'
        call mpi_type_commit(type_sub_view_array, errs(ne))

        ne=ne+1; lbls(ne)='read_all itStat'
        CALL MPI_FILE_READ_ALL(fh, itStat, 1, MPI_INTEGER, MPI_STATUS_IGNORE, errs(ne))

        ne=ne+1; lbls(ne)='type_size MPI_INTEGER'
        call MPI_TYPE_SIZE(MPI_INTEGER, int_size, errs(ne))
        ne=ne+1; lbls(ne)='type_size MPI_RP'
        call MPI_TYPE_SIZE(MPI_RP, real_size, errs(ne))

        ne=ne+1; lbls(ne)='file_get_size'
        call MPI_FILE_GET_SIZE(fh, fsize, errs(ne))
        nOnDisk = int((fsize - int(int_size,MPI_OFFSET_KIND)) /                       &
                      (int(nx,MPI_OFFSET_KIND)*int(ny,MPI_OFFSET_KIND)*int(real_size,MPI_OFFSET_KIND)))
        if((nOnDisk /= nvAve2D .and. nOnDisk /= 94) .or. fsize /= int(int_size,MPI_OFFSET_KIND) + &
           int(nOnDisk,MPI_OFFSET_KIND)*int(nx,MPI_OFFSET_KIND)*int(ny,MPI_OFFSET_KIND)*int(real_size,MPI_OFFSET_KIND)) then
          if(rank == root) write(*,'(1x,4A)') ' ERROR: ', "'"//trim(filename)//"'", &
               ' is not a complete 2D statistics file of this case (crashed write or other grid).', &
               ' Remove it from BINARY_STAT/stats_index.txt to use the previous sample, or run with restartStat = .false.'
          call secure_stop
        endif
        nRead   = min(nOnDisk, nvAve2D)
        if(nRead < nvAve2D .and. rank == root)                                        &
          write(*,'(1x,A,I0,A,I0,A)') ' WARNING: statistics file holds ', nOnDisk,    &
               ' slots, this build expects ', nvAve2D, ' -- the new slots start from zero and their'// &
               ' means are biased low (divided by all samples).'

        initial_displacement = int_size
        ne=ne+1; lbls(ne)='file_set_view'
        CALL MPI_FILE_SET_VIEW(fh, initial_displacement, MPI_RP, type_sub_view_array, "native", MPI_INFO_NULL, errs(ne))

        ne=ne+1; lbls(ne)='read_all planes'
        errs(ne) = 0
        do l = 1,nRead
           CALL MPI_FILE_READ_ALL(fh, vmean(:,:,l), 1, type_sub_array, MPI_STATUS_IGNORE, err)
           errs(ne) = ior(errs(ne), abs(err))
        enddo

        ne=ne+1; lbls(ne)='file_close'
        CALL MPI_FILE_CLOSE(fh, errs(ne))

        ne=ne+1; lbls(ne)='type_free loc'
        call MPI_TYPE_FREE(type_sub_array,      errs(ne))
        ne=ne+1; lbls(ne)='type_free view'
        call MPI_TYPE_FREE(type_sub_view_array, errs(ne))

        call check_mpi_vec('mpi_read_stat2D', errs(1:ne), lbls(1:ne))

        return
end subroutine mpi_read_stat2D

end module init_bc_module
