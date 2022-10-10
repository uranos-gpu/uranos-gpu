module FileModule
implicit none

public
type :: FileType
  character(500) :: name
  character(500) :: dir
  character(10)  :: iter
  integer        :: unit
  integer        :: err = 0
  logical        :: isopen = .false.
endtype FileType

type DirType
  character(500)                            :: name
  character(500), allocatable, dimension(:) :: content
  integer                                   :: files
end type DirType

type(FileType) :: file
type(DirType)  :: directory

contains


subroutine GetFileUnit(file)
! ------------------------------------------------
!       This function get the unit of a file
! ------------------------------------------------

        implicit none
        type(FileType), intent(inout) :: file
        integer :: i
        
        file%unit = 10
        
        do i = 10, 99
        
          !if(i /= 5 .and. i /= 6 .and. i /= 9) then
        
            inquire(unit = i, opened = file%isopen, iostat = file%err)
        
            if(file%err == 0) then
              if(.not. file%isopen) then

                file%unit = i
                return

              end if
            end if
          !end if
        end do
        
        return
end subroutine GetFileUnit


subroutine OpenNewFile(file,it)
! -----------------------------------------------------------------
!       This subroutine open a new file
! -----------------------------------------------------------------

        implicit none
        type(FileType), intent(inout) :: file
        integer       , intent(in)    :: it
        logical                       :: dir_exist

        file%err = 0
        !
        ! === Set file directory
        !
        inquire(file = 'DATA/'//trim(file%dir), exist= dir_exist)
        if(.not. dir_exist) then

          call execute_command_line('mkdir -p DATA/'//trim(file%dir))

        endif
        !
        ! === Set file name
        !
        write(file%iter,'(i7.7)') it

        file%name = 'DATA/'//trim(file%dir)//'/'//trim(file%name)//trim(file%iter)//'.txt'

        !
        ! === Open the file
        !
        call GetFileUnit(file)
        open(file = file%name, unit = file%unit, status = 'unknown', iostat = file%err)

        if(file%err .ne. 0) then
          print*, ' FATAL ERROR! File ', trim(file%name), ' has not been opened'
          stop
        endif

        return
end subroutine OpenNewFile


subroutine OpenVTKFile(file,it)
        implicit none
        type(FileType), intent(inout) :: file
        integer       , intent(in)    :: it
        logical                       :: dir_exist


        file%err = 0
        !
        ! === Set file directory
        !
        inquire(file = 'DATA/'//trim(file%dir), exist= dir_exist)
        if(.not. dir_exist) then

          call execute_command_line('mkdir -p DATA/'//trim(file%dir))

        endif
        !
        ! === Set file name
        !
        write(file%iter,'(i7.7)') it

        file%name = 'DATA/'//trim(file%dir)//'/'//trim(file%name)//trim(file%iter)//'.vtk'

        !
        ! === Open the file
        !
        call GetFileUnit(file)

        open(file   = file%name         , &
             unit   = file%unit         , &
             form   = 'unformatted'     , &
             action = 'write'           , &
             access = 'stream'          , &
             iostat = file%err          , &
             convert = 'big_endian')

        if(file%err .ne. 0) then
          print*, ' FATAL ERROR! File ', trim(file%name), ' has not been opened'
          stop
        endif


        return
end subroutine OpenVTKFile





subroutine AppendToFile(file,it)
        implicit none
        type(FileType), intent(inout) :: file
        integer       , intent(in)    :: it
        logical                       :: dir_exist
        logical                       :: file_exist

        file%err = 0
        !
        ! === Set file directory
        !
        inquire(file = 'DATA/'//trim(file%dir), exist= dir_exist)
        if(.not. dir_exist) then

          call execute_command_line('mkdir -p DATA/'//trim(file%dir))

        endif
        !
        ! === Set the file name
        !
        file%name = 'DATA/'//trim(file%dir)//'/'//trim(file%name)//'.txt'
        !
        ! === Open the file
        !
        call GetFileUnit(file)
        inquire(file = file%name, exist= file_exist)
        if(it == 0) then

          if(file_exist) call execute_command_line('rm -f '//trim(file%name))
          open(file = file%name, unit = file%unit, status = 'new', iostat = file%err)

        else
          open(file = file%name, unit = file%unit, status = 'unknown', position='append',iostat = file%err)

        endif

        if(file%err .ne. 0) then
          print*, ' FATAL ERROR! File ', trim(file%name), ' has not been opened'
          stop
        endif

        return
end subroutine AppendToFile


subroutine CloseFile(file)
        implicit none
        type(FileType), intent(in) :: file

        flush(file%unit)
        close(file%unit)

        return
end subroutine CloseFile




subroutine CheckDirExists(directory)
        implicit none
        type(DirType), intent(in) :: directory
        logical                   :: dir_exist

        inquire(file = trim(directory%name), exist= dir_exist)
        if(.not.dir_exist) then

          write(*,'(A)')     ' FATAL ERROR! '
          write(*,'(A,A,A)') ' The directory ',trim(directory%name), 'does not exist'
          write(*,'(A)')     ' The program will be stopped'
          stop

        endif
        return
end subroutine CheckDirExists




subroutine GetDirContent(directory)
        implicit none
        type(DirType), intent(inout) :: directory

        ! local
        character(500) :: temp_file
        integer        :: f_unit = 10, io_err = 0
        integer        :: line, fake_var


        ! save ls in a temporary file
        temp_file = 'contents.txt'
        call execute_command_line('ls '//trim(directory%name)//' > '//temp_file)

        ! --- open temp_file and read it line by line
        open(unit   = f_unit,       &
             file   = temp_file,    &
             action = "read",       &
             access = "sequential", &
             iostat = io_err)
        if(io_err .ne. 0) stop ' FATAL ERROR! GetDirContent Fails'

        ! how many files are there?
        directory%files = 0
        do
                read(f_unit,FMT='(a)',iostat=io_err) fake_var
                if (io_err/=0) EXIT
                directory%files = directory%files+1
        end do
        rewind(f_unit)

        ! --- save the dir content in a variable
        allocate(directory%content(directory%files))

        do line = 1,directory%files
           read(f_unit,'(a)') directory%content(line)
        enddo

        call execute_command_line('rm -f '//trim(temp_file))
        close(f_unit)


        return
end subroutine GetDirContent




subroutine mkdir(directory)
        implicit none
        character(*), intent(in) :: directory

        ! local declarations
        logical :: dir_exist

        inquire(file = directory, exist= dir_exist)
        if(.not.dir_exist) then
          call execute_command_line('mkdir -p '//trim(directory))
          write(*,*) " Directory "//trim(directory)//" has been created."
        endif

        return
end subroutine mkdir
















subroutine GnuplotHeader(file,sx,ex,sy,ey,sz,ez,time,ic)

        implicit none
        type(FileType)              , intent(in) :: file
        integer                     , intent(in) :: sx, ex, sy, ey, sz, ez
        real(8)                     , intent(in) :: time
        character(*)                , intent(in) :: ic

        write(file%unit,'(A)'      ) " # ================================= #"
        write(file%unit,'(A,A)'    ) " # test     = ", trim(ic)
        write(file%unit,'(A,e18.6)') " # time     = ", time
        write(file%unit,'(A,I3)'   ) " # x_points = ", ex-sx+1
        write(file%unit,'(A,I3)'   ) " # y_points = ", ey-sy+1
        write(file%unit,'(A,I3)'   ) " # z_points = ", ez-sz+1
        write(file%unit,'(A)'      ) " # ================================= #"


        return
end subroutine GnuplotHeader


function str(num) result (string)
        integer, intent(in) :: num
        character(20)       :: string
        
        write(string,*) num

        string = adjustl(string)
        
        return
end function str





end module FileModule

