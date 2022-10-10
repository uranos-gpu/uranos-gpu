module post_shell_tools_module
use parameters_module
implicit none


contains
subroutine get_what_is_inside_a_directory(dir_name,dir_contents,files)
! ----------------------------------------------------------------------
!       
!       This subroutine gets the name of a directory and its content.
!
! ----------------------------------------------------------------------
        implicit none
        character(dl)                           , intent(in ) :: dir_name
        character(dl), allocatable, dimension(:), intent(out) :: dir_contents
        integer                                 , intent(out) :: files

        ! local declaration
        character(dl) :: temp_file
        integer       :: f_unit = 10, line, fake_var, io_err
        
        ! save ls in a temporary file
        temp_file = 'contents.txt'
        call execute_command_line('ls '//trim(dir_name)//' > '//temp_file)

        ! --- open temp_file and read it line by line
        open(unit   = f_unit,       &
             file   = temp_file,    &
             action = "read",       &
             access = "sequential", &
             iostat = io_err)

        if(io_err .ne. 0) then
                write(*,*)' Access denied to ', temp_file
                stop
        endif

        ! how many files are there?
        files = 0
        do
                read(f_unit,FMT='(a)',iostat=io_err) fake_var
                if (io_err/=0) EXIT
                files = files+1
        end do
        rewind(f_unit)

        ! --- save the dir content in a variable
        allocate(dir_contents(files))

        do line = 1,files
           read(f_unit,'(a)') dir_contents(line)
        enddo

        call execute_command_line('rm -f '//trim(temp_file))
        close(f_unit)

        return
end subroutine get_what_is_inside_a_directory







end module post_shell_tools_module
