module vtk_utils_module
use parameters_module, only: rp, data_dir, output_file_name
use FileModule

interface WriteVTKHeader
  module procedure WriteVTKHeader2D, WriteVTKHeader3D
end interface

interface ScalarVTK
  module procedure ScalarVTK_BINARY2D, ScalarVTK_BINARY2Dcodim, ScalarVTK_BINARY3D, ScalarVTK_BINARY3Dcodim, &
                   Int1scalarVTK_BINARY3D
end interface

interface VectorVTK
  module procedure vectorVTK_BINARY2D, vectorVTK_BINARY3D
end interface

private
public WriteVTKHeader, ScalarVTK, VectorVTK

contains

subroutine WriteVTKHeader2D(file,x,y,lo,up)
        implicit none
        type(FileType)                     , intent(in) :: file
        real(rp), dimension(:), allocatable, intent(in) :: x,y
        integer , dimension(2)             , intent(in) :: lo, up

        character(len=1), parameter :: newline = achar(10)
        integer                     :: funit, err = 0, nx, ny
        character(len=100)          :: s_buffer
        character(len=30)           :: buf2
        character(len=256)          :: GRID_header

        funit = file%unit
        nx = up(1) - lo(1) + 1
        ny = up(2) - lo(2) + 1

        write(unit = funit, iostat = err) '# vtk DataFile Version 3.0' // newline
        write(unit = funit, iostat = err) 'test file' // newline
        write(unit = funit, iostat = err) 'BINARY' // newline
        write(unit = funit, iostat = err) newline

        ! === grid header
        write(s_buffer, fmt = '(A)'     , iostat = err) 'DATASET RECTILINEAR_GRID'
        write(unit = funit, iostat = err) trim(s_buffer) // newline
        write(s_buffer, FMT = '(A,3I12)', iostat = err) 'DIMENSIONS', nx, ny, 1
        write(unit = funit, iostat = err) trim(s_buffer) // newline
        
        ! === x coordinate
        write(buf2,'(i8)') nx
        GRID_header="X_COORDINATES "//trim(adjustl(buf2))//" double"//newline
        write(unit = funit) trim(GRID_header),real(x(lo(1):up(1)),rp),newline

        ! === y coordinate
        write(buf2,'(i8)') ny
        GRID_header="Y_COORDINATES "//trim(adjustl(buf2))//" double"//newline
        write(unit = funit) trim(GRID_header),real(y(lo(2):up(2)),rp),newline

        ! ==== z coordinate
        write(buf2,'(i8)') 1
        GRID_header="Z_COORDINATES "//trim(adjustl(buf2))//" double"//newline
        write(unit = funit) trim(GRID_header),0.0_rp,newline

        ! ==== writing the data
        write(unit = funit, iostat = err) newline
        write(s_buffer, FMT = '(A,I12)', iostat = err) 'POINT_DATA', nx*ny
        write(unit = funit, iostat = err) trim(s_buffer) // newline


        return
end subroutine WriteVTKHeader2D




subroutine WriteVTKHeader3D(file,x,y,z,id,lo,up)
        implicit none
        type(FileType)                     , intent(in) :: file
        real(rp), dimension(:), allocatable, intent(in) :: x,y,z
        integer , dimension(3)             , intent(in) :: id, lo, up

        character(len=1), parameter :: newline = achar(10)
        integer                     :: funit, err = 0
        character(len=100)          :: s_buffer
        character(len=30)           :: buf2
        character(len=256)          :: GRID_header

        funit = file%unit

        write(unit = funit, iostat = err) '# vtk DataFile Version 3.0' // newline
        write(unit = funit, iostat = err) 'test file' // newline
        write(unit = funit, iostat = err) 'BINARY' // newline
        write(unit = funit, iostat = err) newline

        ! === grid header
        write(s_buffer, fmt = '(A)'     , iostat = err) 'DATASET RECTILINEAR_GRID'
        write(unit = funit, iostat = err) trim(s_buffer) // newline
        write(s_buffer, FMT = '(A,3I12)', iostat = err) 'DIMENSIONS', id(1), id(2), id(3)
        write(unit = funit, iostat = err) trim(s_buffer) // newline
        
        ! === x coordinate
        write(buf2,'(i8)') id(1)
        GRID_header="X_COORDINATES "//trim(adjustl(buf2))//" double"//newline
        write(unit = funit) trim(GRID_header),real(x(lo(1):up(1)),rp),newline

        ! === y coordinate
        write(buf2,'(i8)') id(2)
        GRID_header="Y_COORDINATES "//trim(adjustl(buf2))//" double"//newline
        write(unit = funit) trim(GRID_header),real(y(lo(2):up(2)),rp),newline

        ! ==== z coordinate
        write(buf2,'(i8)') id(3)
        GRID_header="Z_COORDINATES "//trim(adjustl(buf2))//" double"//newline
        if(id(3)==1) then
          write(unit = funit) trim(GRID_header),0.0_rp,newline
        else
          write(unit = funit) trim(GRID_header),real(z(lo(3):up(3)),rp),newline
        endif

        ! ==== writing the data
        write(unit = funit, iostat = err) newline
        write(s_buffer, FMT = '(A,I12)', iostat = err) 'POINT_DATA', id(1)*id(2)*id(3)
        write(unit = funit, iostat = err) trim(s_buffer) // newline


        return
end subroutine WriteVTKHeader3D



subroutine scalarVTK_BINARY2D(file, lo, up, var_name, var)
        implicit none
        integer , dimension(2)               , intent(in) :: lo, up
        real(rp), allocatable, dimension(:,:), intent(in) :: var
        type(FileType)                       , intent(in) :: file
        character(len=*)                     , intent(in) :: var_name
        
        ! local declarations
        character(len=1) , parameter  :: newline = char(10)
        character(len=100)            :: s_buffer
        integer                       :: iostatus
        
        write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS '//trim(var_name)//' double', 1
        write(unit = file%unit, iostat = iostatus) trim(s_buffer) // newline
        write(unit = file%unit, iostat = iostatus) 'LOOKUP_TABLE default' // newline
        
        write(file%unit) var(lo(1):up(1),lo(2):up(2))

        return
end subroutine scalarVTK_BINARY2D


subroutine scalarVTK_BINARY2Dcodim(file, lo, up, var_name, var, codim)
        implicit none
        integer , dimension(2)                 , intent(in) :: lo, up
        real(rp), allocatable, dimension(:,:,:), intent(in) :: var
        type(FileType)                         , intent(in) :: file
        character(len=*)                       , intent(in) :: var_name
        integer                                , intent(in) :: codim
        
        ! local declarations
        character(len=1) , parameter  :: newline = char(10)
        character(len=100)            :: s_buffer
        integer                       :: iostatus
        
        write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS '//trim(var_name)//' double', 1
        write(unit = file%unit, iostat = iostatus) trim(s_buffer) // newline
        write(unit = file%unit, iostat = iostatus) 'LOOKUP_TABLE default' // newline
        
        write(file%unit) var(lo(1):up(1),lo(2):up(2),codim)

        return
end subroutine scalarVTK_BINARY2Dcodim


subroutine scalarVTK_BINARY3D(file, lo, up, var_name, var)
        implicit none
        integer , dimension(3)                 , intent(in) :: lo, up
        real(rp), allocatable, dimension(:,:,:), intent(in) :: var
        type(FileType)                         , intent(in) :: file
        character(len=*)                       , intent(in) :: var_name
        
        ! local declarations
        character(len=1) , parameter  :: newline = char(10)
        character(len=100)            :: s_buffer
        integer                       :: iostatus
        
        write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS '//trim(var_name)//' double', 1
        write(unit = file%unit, iostat = iostatus) trim(s_buffer) // newline
        write(unit = file%unit, iostat = iostatus) 'LOOKUP_TABLE default' // newline
        
        write(file%unit) var(lo(1):up(1),lo(2):up(2),lo(3):up(3))

        return
end subroutine scalarVTK_BINARY3D


subroutine scalarVTK_BINARY3Dcodim(file, lo, up, var_name, var, codim)
        implicit none
        integer , dimension(3)                   , intent(in) :: lo, up
        real(rp), allocatable, dimension(:,:,:,:), intent(in) :: var
        type(FileType)                           , intent(in) :: file
        character(len=*)                         , intent(in) :: var_name
        integer                                  , intent(in) :: codim
        
        ! local declarations
        character(len=1) , parameter  :: newline = char(10)
        character(len=100)            :: s_buffer
        integer                       :: iostatus
        
        write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS '//trim(var_name)//' double', 1
        write(unit = file%unit, iostat = iostatus) trim(s_buffer) // newline
        write(unit = file%unit, iostat = iostatus) 'LOOKUP_TABLE default' // newline
        
        write(file%unit) var(lo(1):up(1),lo(2):up(2),lo(3):up(3),codim)

        return
end subroutine scalarVTK_BINARY3Dcodim


subroutine Int1scalarVTK_BINARY3D(file, lo, up, var_name, var)
        implicit none
        integer   , dimension(3)                 , intent(in) :: lo, up
        integer(1), allocatable, dimension(:,:,:), intent(in) :: var
        type(FileType)                           , intent(in) :: file
        character(len=*)                         , intent(in) :: var_name
        
        ! local declarations
        character(len=1) , parameter  :: newline = char(10)
        character(len=100)            :: s_buffer
        integer                       :: iostatus
        
        write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'SCALARS '//trim(var_name)//' double', 1
        write(unit = file%unit, iostat = iostatus) trim(s_buffer) // newline
        write(unit = file%unit, iostat = iostatus) 'LOOKUP_TABLE default' // newline
        
        write(file%unit) real(var(lo(1):up(1),lo(2):up(2),lo(3):up(3)),rp)

        return
end subroutine Int1scalarVTK_BINARY3D


subroutine vectorVTK_BINARY2D(file,lo,up,vec_name,vec_x,vec_y,vec_z)
        implicit none
        integer, dimension(2)                , intent(in) :: lo,up
        real(rp), allocatable, dimension(:,:), intent(in) :: vec_x, vec_y, vec_z
        type(FileType)                       , intent(in) :: file
        character(len=*)                     , intent(in) :: vec_name

        ! local declarations
        character(len=1) , parameter  :: newline = char(10)
        character(len=100)            :: s_buffer
        integer                       :: iostatus, i,j
        
        write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'VECTORS '//trim(vec_name)//' double'
        write(unit = file%unit, iostat = iostatus) trim(s_buffer) // newline
       
        do j    = lo(2),up(2)
           do i = lo(1),up(1)

              write(file%unit) vec_x(i,j), vec_y(i,j), vec_z(i,j)

           enddo
        enddo

        return
end subroutine vectorVTK_BINARY2D


subroutine vectorVTK_BINARY3D(file, lo,up,vec_name, vec_x, vec_y, vec_z)
        implicit none
        integer, dimension(3)                  , intent(in) :: lo,up
        real(rp), allocatable, dimension(:,:,:), intent(in) :: vec_x, vec_y, vec_z
        type(FileType)                         , intent(in) :: file
        character(len=*)                       , intent(in) :: vec_name

        ! local declarations
        character(len=1) , parameter  :: newline = char(10)
        character(len=100)            :: s_buffer
        integer                       :: iostatus, i,j,k
        
        write(s_buffer, FMT = '(A,I12)', iostat = iostatus) 'VECTORS '//trim(vec_name)//' double'
        write(unit = file%unit, iostat = iostatus) trim(s_buffer) // newline
       
        do k       = lo(3),up(3)
           do j    = lo(2),up(2)
              do i = lo(1),up(1)

                 write(file%unit) vec_x(i,j,k), vec_y(i,j,k), vec_z(i,j,k)

              enddo
           enddo
        enddo

        return
end subroutine vectorVTK_BINARY3D












end module vtk_utils_module
