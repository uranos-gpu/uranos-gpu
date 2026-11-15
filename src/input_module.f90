module input_module
use parameters_module
use mpi_module

implicit none
interface open_unit
  module procedure open_unit
endinterface

character(dl) :: input_file

contains

subroutine parse_command_line_args
        implicit none

        if(command_argument_count() < 1) then
          if(rank == root) then
            print*, ' ERROR: missing input data file argument.'
            print*, '        usage: Uranos.exe <input.dat> [restart.bin]'
          endif
          call secure_stop
        endif

        call get_command_argument(1, input_file)
        call get_command_argument(2, restart_file)

        if(rank == root) then
          print*, ' Reading input data from file ', "'"//trim(input_file)//"'", '.'
          if(len_trim(restart_file) > 0) then
            print*, ' Restart file set to ', "'"//trim(restart_file)//"'", '.'
          endif
        endif

        return
end subroutine parse_command_line_args

subroutine read_input_namelist
        implicit none
        integer              :: iunit, io_err
        character(len = 256) :: io_msg

        open(newunit = iunit, file = input_file, &
             action = 'read', status = 'old', &
             iostat = io_err, iomsg = io_msg)
        if(io_err /= 0) then
          if(rank == root) then
            print*, ' ERROR: cannot open input file "'//trim(input_file)//'".'
            print*, '        iomsg = ', trim(io_msg)
          endif
          call secure_stop
        endif

        read(iunit, nml = input_list, iostat = io_err, iomsg = io_msg)
        if(io_err /= 0) then
          if(rank == root) then
            print*, ' ERROR: failed to read /input_list/ from "'//trim(input_file)//'".'
            print*, '        iomsg = ', trim(io_msg)
          endif
          close(iunit)
          call secure_stop
        endif

        close(iunit)

        !$acc update device(Tref)

        if(.not. StFlg) StOut = 0

        return
end subroutine read_input_namelist

subroutine set_reference_quantities
        use thermodynamics_module, only: gamma_Tref, cp_Tref, recovery_temperature
        implicit none

        mu_inf     = sqrt(gamma_Tref)*Mach/Reynolds
        u_inf      = sqrt(gamma_Tref)*Mach
        q_inf      = 0.5_rp*gamma_Tref*Mach**2
        k_inf      = cp_Tref/Prandtl
        iPrandtl   = 1.0_rp/Prandtl
        T_recovery = recovery_temperature()

        return
end subroutine set_reference_quantities

subroutine resolve_scheme_flags
        implicit none
        real(rp) :: kp_max, CFL_max

        if(logical_CFL) then
          if(CFL <= 0.0_rp) then
            if(rank == root) print*, ' ERROR: CFL must be > 0 (got ', CFL, ').'
            call secure_stop
          endif
          if(viscous .and. Fourier <= 0.0_rp) then
            if(rank == root) print*, ' ERROR: Fourier must be > 0 (got ', Fourier, ').'
            call secure_stop
          endif
          if(rank == root) then
            selectcase(fd_order)
              case(2)     ; kp_max = 1.0000_rp
              case(4)     ; kp_max = 1.3722_rp
              case default; kp_max = 1.5860_rp
            endselect
            CFL_max = sqrt(3.0_rp)/kp_max
            if(CFL > CFL_max) then
              print*, ' WARNING: CFL = ', CFL, ' exceeds the SSP-RK3 imaginary-axis limit'
              print*, '          for fd_order = ', fd_order, ' : CFL_max = ', CFL_max
            endif
            if(viscous .and. Fourier > 0.5_rp) &
              print*, ' WARNING: Fourier = ', Fourier, ' is close to the SSP-RK3 real-axis limit.'
          endif
        endif

        hybrid_weno = trim(scheme) == 'hybrid_weno5'  .or. &
                      trim(scheme) == 'hybrid_wenoEP' .or. &
                      trim(scheme) == 'hybrid_tenoEP' .or. &
                      trim(scheme) == 'hybrid_tenoaEP'

        if(hybrid_weno) then
          if(weno_order /= 3 .and. weno_order /= 5 .and. weno_order /= 7) then
            if(rank == root) print*, ' ERROR: weno_order must be 3, 5 or 7 (got ', weno_order, ').'
            call secure_stop
          endif
          select case(trim(scheme))
            case('hybrid_tenoEP')
              if(weno_order /= 5 .and. weno_order /= 7) then
                if(rank == root) print*, ' ERROR: scheme "hybrid_tenoEP" requires weno_order = 5 or 7 (got ', weno_order, ').'
                call secure_stop
              endif
            case('hybrid_tenoaEP')
              if(weno_order /= 5) then
                if(rank == root) print*, ' ERROR: scheme "hybrid_tenoaEP" requires weno_order = 5 (got ', weno_order, ').'
                call secure_stop
              endif
          end select
        endif

        weno_num = (weno_order+1)/2

        select case(trim(ep_type))
          case('JP'  ); ep_id = 1
          case('KEEP'); ep_id = 2
          case('AEC' )
            if(thermally_perfect_gas > 0) then
              ep_id = 4
            else
              ep_id = 3
            endif
          case('PEP' )
            if(thermally_perfect_gas > 0) then
              ep_id = 5
            else
              if(rank == root) print*, ' ERROR: ep_type "PEP" is restricted to a thermally-perfect'// &
                                       ' gas (thermally_perfect_gas = 1). Its calorically-perfect'// &
                                       ' path is unvalidated; use JP, KEEP or AEC for a CP gas.'
              call secure_stop
            endif
          case default
            if(rank == root) print*, ' ERROR: ep_type "'//trim(ep_type)//'" is not implemented.'// &
                                     ' Valid: JP, KEEP, AEC, PEP.'
            call secure_stop
        end select

        return
end subroutine resolve_scheme_flags

subroutine open_unit(iunit,filename,status,access,position)

        implicit none
        integer                     , intent(in) :: iunit
        character(len = *)          , intent(in) :: filename
        character(len = *), optional, intent(in) :: status
        character(len = *), optional, intent(in) :: access
        character(len = *), optional, intent(in) :: position
        integer                                  :: ios = 0

        if    (present(status)  ) then
          open(unit = iunit, file = filename, status = status, iostat = ios)

        elseif(present(access)  ) then
          open(unit = iunit, file = filename, access = access,iostat = ios)

        elseif(present(position)) then
          open(unit = iunit, file = filename, position = position,iostat = ios)

        elseif(present(status).and.present(access)  ) then
          open(unit = iunit, file = filename, status = status, access = access,iostat = ios)

        elseif(present(status).and.present(position)) then
          open(unit = iunit, file = filename, status = status, position = position,iostat = ios)

        elseif(present(access).and.present(position)) then
          open(unit = iunit, file = filename, access = access, position = position,iostat = ios)

        elseif(present(access).and.present(status).and.present(position)) then
          open(unit = iunit, file = filename, status = status, access = access, position = position,iostat = ios)

        else
          open(unit = iunit, file = filename)
        endif

        if(ios /= 0) then

          write(*,'(a)') ''
          write(*,'(a)') ' FATAL ERROR!'
          write(*,'(a)') ' Could not open file "' //trim(filename)// '"'

          stop

        endif

        return
end subroutine open_unit

subroutine make_saving_directories
        use FileModule

        implicit none
        logical :: dir_exist

        if(rank == root) then

          if(data_dir.eq.'temp') then
            inquire(file = 'DATA/'//trim(data_dir), exist= dir_exist)
            if(dir_exist) then
              call execute_command_line('rm -r DATA/'//trim(data_dir))
              write(*,*) " Previous DATA/"//trim(data_dir)//" has been removed."
            endif
          endif

          call mkdir('DATA/'//trim(data_dir))
          call mkdir('DATA/'//trim(data_dir)//'/BINARY')
          if(StOut > 0) call mkdir('DATA/'//trim(data_dir)//'/BINARY_STAT')
          if(turb_inflow) call mkdir('DATA/'//trim(data_dir)//'/BINARY_DF')

        endif
        return
end subroutine make_saving_directories

end module
