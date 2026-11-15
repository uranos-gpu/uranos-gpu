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

        implicit none
        type(FileType), intent(inout) :: file
        integer :: i

        file%unit = 10

        do i = 10, 99

            inquire(unit = i, opened = file%isopen, iostat = file%err)

            if(file%err == 0) then
              if(.not. file%isopen) then

                file%unit = i
                return

              end if
            end if
        end do

        return
end subroutine GetFileUnit

subroutine OpenNewFile(file,it)

        implicit none
        type(FileType), intent(inout) :: file
        integer       , intent(in)    :: it
        logical                       :: dir_exist

        file%err = 0
        inquire(file = 'DATA/'//trim(file%dir), exist= dir_exist)
        if(.not. dir_exist) then

          call execute_command_line('mkdir -p DATA/'//trim(file%dir))

        endif
        write(file%iter,'(i0.7)') it

        file%name = 'DATA/'//trim(file%dir)//'/'//trim(file%name)//trim(file%iter)//'.txt'

        call GetFileUnit(file)
        open(file = file%name, unit = file%unit, status = 'unknown', iostat = file%err)

        if(file%err .ne. 0) then
          print*, ' FATAL ERROR! File ', trim(file%name), ' has not been opened'
          stop
        endif

        return
end subroutine OpenNewFile

subroutine AppendToFile(file,it)
        implicit none
        type(FileType), intent(inout) :: file
        integer       , intent(in)    :: it
        logical                       :: dir_exist
        logical                       :: file_exist

        file%err = 0
        inquire(file = 'DATA/'//trim(file%dir), exist= dir_exist)
        if(.not. dir_exist) then

          call execute_command_line('mkdir -p DATA/'//trim(file%dir))

        endif
        file%name = 'DATA/'//trim(file%dir)//'/'//trim(file%name)//'.txt'
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

subroutine mkdir(directory)
        implicit none
        character(*), intent(in) :: directory

        logical :: dir_exist

        inquire(file = directory, exist= dir_exist)
        if(.not.dir_exist) then
          call execute_command_line('mkdir -p '//trim(directory))
          write(*,*) " Directory "//trim(directory)//" has been created."
        endif

        return
end subroutine mkdir

subroutine prepare_row_file(fname, hdr, it_now)
        use mpi
        implicit none
        character(*), intent(in) :: fname, hdr
        integer     , intent(in) :: it_now
        character(len=len(hdr)+1):: line
        character(len=12)        :: c_it
        integer(kind=MPI_OFFSET_KIND) :: fsize, keep, reclen
        integer                  :: unit, ios, nrow, lo, hi, md, it_r, err
        logical                  :: exists

        reclen = len(hdr) + 1
        inquire(file=trim(fname), exist=exists, size=fsize)
        if(it_now > 0 .and. exists) then
          line = ''
          if(fsize >= reclen) then
            open(newunit=unit, file=trim(fname), access='stream', form='unformatted', &
                 status='old', action='read')
            read(unit, pos=1) line
            close(unit)
          endif
          if(line(1:len(hdr)) == hdr .and. line(len(hdr)+1:) == achar(10)) then
            nrow = int(fsize/reclen) - 1
            open(newunit=unit, file=trim(fname), access='stream', form='unformatted', &
                 status='old', action='read')
            lo = 1; hi = nrow + 1
            do while(lo < hi)
               md = (lo + hi)/2
               read(unit, pos=int(md,MPI_OFFSET_KIND)*reclen + 1) c_it
               read(c_it,*,iostat=ios) it_r
               if(ios == 0 .and. it_r <= it_now) then
                 lo = md + 1
               else
                 hi = md
               endif
            enddo
            close(unit)
            keep = int(lo,MPI_OFFSET_KIND)*reclen
            if(keep < fsize) call cut_file(fname, keep, fsize, it_now)
            return
          endif
        endif
        if(exists) call move_aside(fname)

        open(newunit=unit, file=trim(fname), status='replace', iostat=ios)
        if(ios /= 0) then
          print*, ' ERROR: cannot create ', trim(fname)
          call MPI_Abort(MPI_COMM_WORLD, 1, err)
        endif
        write(unit,'(A)') hdr
        close(unit)

        return
end subroutine prepare_row_file

subroutine prepare_record_file(fname, hdr, reclen, it_now, fend)
        use mpi
        implicit none
        character(*)                 , intent(in)  :: fname
        integer, dimension(:)        , intent(in)  :: hdr
        integer(kind=MPI_OFFSET_KIND), intent(in)  :: reclen
        integer                      , intent(in)  :: it_now
        integer(kind=MPI_OFFSET_KIND), intent(out) :: fend
        integer, dimension(size(hdr))  :: old
        integer(kind=MPI_OFFSET_KIND)  :: fsize, hlen, nrec, r
        integer(8)                     :: it_r, it_p
        integer                        :: unit, ios, err
        logical                        :: exists

        hlen = 4*size(hdr)
        inquire(file=trim(fname), exist=exists, size=fsize)
        if(it_now > 0 .and. exists) then
          old = -1
          open(newunit=unit, file=trim(fname), access='stream', form='unformatted', &
               status='old', action='read')
          if(fsize >= hlen) read(unit, pos=1) old
          if(all(old == hdr)) then
            nrec = (fsize - hlen)/reclen
            fend = hlen
            it_p = -1
            do r = 1,nrec
               read(unit, pos=fend + 1) it_r
               if(it_r <= it_p .or. it_r > it_now) exit
               it_p = it_r
               fend = fend + reclen
            enddo
            close(unit)
            if(fend < fsize) call cut_file(fname, fend, fsize, it_now)
            return
          endif
          close(unit)
        endif
        if(exists) call move_aside(fname)

        open(newunit=unit, file=trim(fname), access='stream', form='unformatted', &
             status='replace', iostat=ios)
        if(ios /= 0) then
          print*, ' ERROR: cannot create ', trim(fname)
          call MPI_Abort(MPI_COMM_WORLD, 1, err)
        endif
        write(unit) hdr
        close(unit)
        fend = hlen

        return
end subroutine prepare_record_file

subroutine cut_file(fname, keep, fsize, it_now)
        use mpi
        implicit none
        character(*)                 , intent(in) :: fname
        integer(kind=MPI_OFFSET_KIND), intent(in) :: keep, fsize
        integer                      , intent(in) :: it_now
        character(len=1000) :: cmd
        integer             :: ios, err

        write(cmd,'(A,I0,A)') 'truncate -s ', keep, ' '//trim(fname)
        call execute_command_line(trim(cmd), exitstat=ios)
        if(ios /= 0) then
          print*, ' ERROR: cannot cut ', trim(fname), ' back to it =', it_now
          call MPI_Abort(MPI_COMM_WORLD, 1, err)
        endif
        print*, ' ', trim(fname), ': dropped ', fsize - keep, ' bytes beyond the restart iteration'

        return
end subroutine cut_file

subroutine move_aside(fname)
        use mpi
        implicit none
        character(*), intent(in) :: fname
        character(len=1000) :: cmd
        integer             :: n, ios, err
        logical             :: exists

        n = 0
        cmd = trim(fname)//'.old'
        inquire(file=trim(cmd), exist=exists)
        do while(exists)
           n = n + 1
           write(cmd,'(A,I0)') trim(fname)//'.old', n
           inquire(file=trim(cmd), exist=exists)
        enddo
        call execute_command_line('mv '//trim(fname)//' '//trim(cmd), exitstat=ios)
        if(ios /= 0) then
          print*, ' ERROR: cannot move aside ', trim(fname)
          call MPI_Abort(MPI_COMM_WORLD, 1, err)
        endif
        print*, ' WARNING: ', trim(fname), ' is not continued by this run, moved to ', trim(cmd)

        return
end subroutine move_aside

function str(num) result (string)
        integer, intent(in) :: num
        character(20)       :: string

        write(string,*) num

        string = adjustl(string)

        return
end function str

end module FileModule
