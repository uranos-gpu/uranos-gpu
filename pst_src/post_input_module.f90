module post_input_module
use post_storage_module
use post_shell_tools_module
use FileModule
implicit none


contains
subroutine read_input_file
! -----------------------------------------------------------------
!       
!       This subroutine reads the same input file of the Uranos code
!
! -----------------------------------------------------------------
        implicit none
        integer       :: io_err = 0
        character(dl) :: input_file

        call get_command_argument(1,input_file)
        call get_command_argument(2,restart_file)

        if(len_trim(input_file) == 0) then
                
            print*, ' Input data file missed!. Program will be stopped.' 
            stop

        else
          
            write(*,'(A,A,A)') ' Reading input data from file ', "'"//trim(input_file)//"'", '.'

        endif

        ! Read input file
        open(unit = 1, file = input_file, IOSTAT = io_err)
        if(io_err .ne. 0) then
                
          print*, ' Error in reading file ', trim(input_file), '.'
          stop

        endif

        read(1, nml = input_list)
        close(1)

        ! adimensional grups
        mu_inf = sqrt(gamma0)*Mach/Reynolds
         u_inf = sqrt(gamma0)*Mach
         q_inf = 0.5_rp*gamma0*Mach**2

        k_inf = gamma0/((gamma0-1._rp) * Prandtl)


        ! Set inflow flag
        inflow = .false.
        if(trim(bc(1)) == 'supersonic_inflow'     .or. &
           trim(bc(1)) == 'subsonic_inflow'       .or. &
           trim(bc(1)) == 'nscbc_inflow'          .or. &
           trim(bc(1)) == 'nscbc_inflow_relaxed'  .or. &
           trim(bc(1)) == 'shock_inflow')         then

           inflow = .true.

        endif

        ! set hybrid_weno flag
        if(scheme == 'hybrid_wenoEP') hybrid_weno = .true.
        if(scheme == 'hybrid_tenoEP') hybrid_weno = .true.

        return
end subroutine read_input_file

subroutine read_input_dir
! --------------------------------------------------------------------
!       
!       This subroutine read the contents of the simulation directory
!
! --------------------------------------------------------------------
        implicit none
        
        main_dir%name = 'DATA/'//trim(data_dir)//'/BINARY'
        
        call CheckDirExists(main_dir)
            write(*,'(A,A,A)') ' Reading input data from directory ', &
                               "'"//trim(main_dir%name)//"'", '.'

        call GetDirContent(main_dir)
        file_end = main_dir%Files

        ! get the file for restart
        file_ini = 1
        if(len_trim(restart_file) .ne. 0) then

          do f = 1, size(main_dir%content)

             ifile = trim(main_dir%name)//'/'//trim(main_dir%content(f))
             
             if(trim(ifile) == trim(restart_file)) then
                     file_ini = f
                     exit
             endif

          enddo

          if(file_ini == 1) then
            print*, "Restart file ", trim(restart_file), ' does not exists'
            stop
          endif

        endif

        return
end subroutine read_input_dir


subroutine read_binary_file(filename)
        implicit none
        character(len=*), intent(in) :: filename

        ! local declaration
        integer :: iostatus, f_unit = 10

        open(  unit   = f_unit,           &
             & file   = filename,         & 
             & STATUS = 'old',            &
             & ACCESS = 'stream',         &
             & FORM   = 'unformatted',    &
             & IOSTAT = iostatus) 

        if(iostatus.ne.0) then
                write(*,'(3A)') ' File ', "'"//trim(filename)//"'", ' does not exist.'
                STOP
        endif
        
        ! update time old
        time_old = time_new

        read(f_unit) nx
        read(f_unit) ny
        read(f_unit) nz
        read(f_unit) it
        read(f_unit) time_new
        read(f_unit) dt

        read(f_unit) phi(sx:ex,sy:ey,sz:ez,1)
        read(f_unit) phi(sx:ex,sy:ey,sz:ez,2)
        read(f_unit) phi(sx:ex,sy:ey,sz:ez,3)
        read(f_unit) phi(sx:ex,sy:ey,sz:ez,4)
        read(f_unit) phi(sx:ex,sy:ey,sz:ez,5)

        close(f_unit)

        ! computing file dt
        file_dt = time_new - time_old   !< file time step

        ! computing actual time
        time    = time + file_dt        !< current time (exactly equals to time_new)
        
        ! computing time from the restart file
        if(f == file_ini) time_restart = time
        time_from_restart = time - time_restart
        
return
end subroutine read_binary_file


subroutine make_saving_directories
! ----------------------------------------------------------------------
!
!       This subroutine create the saving directories
!
! ----------------------------------------------------------------------
        implicit none
        logical :: dir_exist

        ! --- create VTK subdirectory --- !
        if(VTK) then

                inquire(file = 'DATA/'//trim(data_dir)//'/VTK', exist= dir_exist)
                if(dir_exist.eqv..false.) then
                        call execute_command_line('mkdir -p DATA/'//trim(data_dir)//'/VTK')
                        write(*,*) " Directory DATA/"//trim(data_dir)//"/VTK has been created."
                endif

        endif

        ! --- create NPY subdirectory --- !
        if(NPY) then

                inquire(file = 'DATA/'//trim(data_dir)//'/NPY', exist= dir_exist)
                if(dir_exist.eqv..false.) then
                        call execute_command_line('mkdir -p DATA/'//trim(data_dir)//'/NPY')
                        write(*,*) " Directory DATA/"//trim(data_dir)//"/NPY has been created."
                endif

        endif

        ! --- create VTK2D subdirectory --- !
        if(VTK2D) then

                inquire(file = 'DATA/'//trim(data_dir)//'/VTK2D', exist= dir_exist)
                if(dir_exist.eqv..false.) then
                        call execute_command_line('mkdir -p DATA/'//trim(data_dir)//'/VTK2D')
                        write(*,*) " Directory DATA/"//trim(data_dir)//"/VTK2D has been created."
                endif

        endif
        
        ! --- create GNUPLOT subdirectory --- !
        if(GNUPLOT) then

                inquire(file = 'DATA/'//trim(data_dir)//'/GNUPLOT', exist= dir_exist)
                if(dir_exist.eqv..false.) then
                        call execute_command_line('mkdir -p DATA/'//trim(data_dir)//'/GNUPLOT')
                        write(*,*) " Directory DATA/"//trim(data_dir)//"/GNUPLOT has been created."
                endif

        endif

        ! --- create GRID subdirectory --- !
        inquire(file = 'DATA/'//trim(data_dir)//'/GRID', exist= dir_exist)
        if(dir_exist.eqv..false.) then
                call execute_command_line('mkdir -p DATA/'//trim(data_dir)//'/GRID')
                write(*,*) " Directory DATA/"//trim(data_dir)//"/GRID has been created."
        endif

        return
end subroutine make_saving_directories







end module post_input_module
