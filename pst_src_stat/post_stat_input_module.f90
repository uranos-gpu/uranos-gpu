module post_stat_input_module
use post_stat_storage_module
use FileModule
implicit none

private
public ReadInputFile, ReadInputDirectoryContent, ReadBinaryStatFile, &
       ReadBinaryWmlesStatFile, GetRestartFile

contains
subroutine ReadInputFile
        implicit none
        integer       :: io_err = 0
        character(dl) :: input_file

        call get_command_argument(1,input_file)
        call get_command_argument(2,restart_file)

        if(len_trim(input_file) == 0) then
            stop ' Input data file missed!. Program will be stopped.' 
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

        wmles = .false.
        if(trim(bc(3)) == 'dws_isothermal'          .or. &
           trim(bc(3)) == 'dws_adiabatic'           .or. &
           trim(bc(3)) == 'idws_adiabatic'          .or. &
           trim(bc(3)) == 'static_dws_adiabatic'    .or. &
           trim(bc(3)) == 'istatic_dws_adiabatic'   .or. &
           trim(bc(3)) == 'vorticity_dws_adiabatic' .or. &
           trim(bc(3)) == 'aws_adiabatic'           .or. &
           trim(bc(3)) == 'aws_isothermal') then
           wmles = .true.
        endif

        return
end subroutine ReadInputFile


subroutine ReadInputDirectoryContent(bindir,dir)

        implicit none
        character(*) , intent(in)  :: bindir
        type(DirType), intent(out) :: dir

        dir%name = 'DATA/'//trim(data_dir)//'/'//trim(bindir)

        call CheckDirExists(dir)

        write(*,'(A,A)') ' Reading input data from directory ', &
                         "'"//trim(dir%name)//"'"

        call GetDirContent(dir)

        return
end subroutine ReadInputDirectoryContent


subroutine GetRestartFile(dir,file_ini,file_end)
        implicit none
        type(DirType), intent(in)  :: dir
        integer      , intent(out) :: file_ini, file_end

        character(500) :: ifile
        integer        :: i

        file_end = dir%Files

        ! get the restart file
        file_ini = 1
        if(len_trim(restart_file) .ne. 0) then

          if(trim(restart_file) == 'restart') then
            file_ini = file_end
          else

            do i = 1, dir%Files

               ifile = trim(dir%name)//'/'//trim(dir%content(i))
               
               if(trim(ifile) == trim(restart_file)) then
                 file_ini = i
                 exit
               endif
            enddo

          endif
        endif
        
        if(len_trim(restart_file) .ne. 0) then
          if(file_ini == 1) then
            print*, "Restart file ", trim(restart_file), ' does not exists'      
            stop
          endif
        endif

        return
end subroutine GetRestartFile


subroutine ReadBinaryStatFile(filename)
        implicit none
        character(*), intent(in) :: filename

        ! local declaration
        integer :: l, err = 0, funit = 10

        open(unit = funit, file = filename, status = 'old', access = 'stream', &
             form = 'unformatted', iostat=err)
        if(err .ne. 0) then
          write(*,'(3A)') ' File ', "'"//trim(filename)//"'", ' does not exist.'
          stop
        endif

        selectcase(pcharAll)

          case('FFT','FTT')

          read(funit) itStat
          do l = 1,nvAve2D
             read(funit) vmean2D(sx:ex,sy:ey,l)
          enddo

          case('TFT')

          read(funit) itStat
          do l = 1,nvAve1D
             read(funit) vmean1D(sy:ey,l)
          enddo

        endselect
        
        close(funit)

        return
end subroutine ReadBinaryStatFile


subroutine ReadBinaryWmlesStatFile(filename)
        implicit none
        character(*), intent(in) :: filename

        ! local declarations
        integer :: l, err = 0, funit = 10, dummy

        open(unit = funit, file = filename, status = 'old', access = 'stream', &
             form = 'unformatted', iostat = err)

        selectcase(ic)

          case('turbulent_BL','swbli')

          read(funit) dummy
          do l = 1, nVWMLESData
             read(funit) vmean1D_wmles(sx:ex,l)
          enddo

        end select

        close(funit)

        return
end subroutine ReadBinaryWmlesStatFile





end module post_stat_input_module
