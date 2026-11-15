module output_restart
use parameters_module
use mpi_module
use storage_module
use profiling_module, only: StartProfRange, EndProfRange
implicit none

private
public mpi_write, mpi_write_stat1D, mpi_write_stat2D

contains
subroutine mpi_write
        implicit none

        character(dl)                       :: filename, iteration
        integer                             :: fh
        integer                             :: int_size, real_size, header_count
        integer                             :: open_err_max, alr_err, ne
        integer(kind = MPI_OFFSET_KIND)     :: initial_displacement, header_offset
        integer, parameter                  :: array_rank=3
        integer, dimension(array_rank)      :: shape_array, shape_sub_array, start_coord
        integer, dimension(array_rank)      :: shape_view_array, shape_sub_view_array, start_view_coord
        integer                             :: type_sub_array, type_sub_view_array

        integer, parameter                  :: NIO = 22
        integer, dimension(NIO)             :: errs
        character(len=28), dimension(NIO)   :: lbls

        call StartProfRange('mpi_write')

        write(iteration, '(i0.7)') it
        filename = "DATA/"//trim(data_dir)//"/BINARY/"//trim(output_file_name)//trim(iteration)//'.bin'

        call StartProfRange('mpi_write_open')

        CALL MPI_FILE_OPEN(mpi_comm_cart, filename, &
                 MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, mpi_err)

        call MPI_Allreduce(mpi_err, open_err_max, 1, MPI_INTEGER, MPI_MAX, mpi_comm_cart, alr_err)

        call EndProfRange

        if(open_err_max .ne. 0) then
          if(rank == root) write(*,'(A,A)') 'MPI_IO error in writing file ', trim(filename)
          call secure_stop
        endif

        ne = 0

        call StartProfRange('mpi_write_truncate')
        ne=ne+1; lbls(ne)='file_set_size'
        call MPI_FILE_SET_SIZE(fh, 0_MPI_OFFSET_KIND, errs(ne))
        call EndProfRange

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

        ne=ne+1; lbls(ne)='type_size MPI_INTEGER'
        call MPI_TYPE_SIZE(MPI_INTEGER, int_size, errs(ne))
        ne=ne+1; lbls(ne)='type_size MPI_RP'
        call MPI_TYPE_SIZE(MPI_RP, real_size, errs(ne))

        header_count = 0
        if(rank == root) header_count = 1

        header_offset = 0
        ne=ne+1; lbls(ne)='write_at_all nx'
        call MPI_FILE_WRITE_AT_ALL(fh, header_offset, nx, header_count, MPI_INTEGER, MPI_STATUS_IGNORE, errs(ne))
        header_offset = header_offset + int_size
        ne=ne+1; lbls(ne)='write_at_all ny'
        call MPI_FILE_WRITE_AT_ALL(fh, header_offset, ny, header_count, MPI_INTEGER, MPI_STATUS_IGNORE, errs(ne))
        header_offset = header_offset + int_size
        ne=ne+1; lbls(ne)='write_at_all nz'
        call MPI_FILE_WRITE_AT_ALL(fh, header_offset, nz, header_count, MPI_INTEGER, MPI_STATUS_IGNORE, errs(ne))
        header_offset = header_offset + int_size
        ne=ne+1; lbls(ne)='write_at_all it'
        call MPI_FILE_WRITE_AT_ALL(fh, header_offset, it, header_count, MPI_INTEGER, MPI_STATUS_IGNORE, errs(ne))
        header_offset = header_offset + int_size
        ne=ne+1; lbls(ne)='write_at_all time'
        call MPI_FILE_WRITE_AT_ALL(fh, header_offset, time, header_count, MPI_RP, MPI_STATUS_IGNORE, errs(ne))
        header_offset = header_offset + real_size
        ne=ne+1; lbls(ne)='write_at_all dt'
        call MPI_FILE_WRITE_AT_ALL(fh, header_offset, dt, header_count, MPI_RP, MPI_STATUS_IGNORE, errs(ne))

        initial_displacement = 4*int_size + 2*real_size
        call StartProfRange('mpi_write_write')
        ne=ne+1; lbls(ne)='file_set_view'
        CALL MPI_FILE_SET_VIEW(fh, initial_displacement, MPI_RP, type_sub_view_array, "native", MPI_INFO_NULL, errs(ne))

        ne=ne+1; lbls(ne)='write_all phi(1)'
        CALL MPI_FILE_WRITE_ALL(fh, phi(:,:,:,1), 1, type_sub_array, MPI_STATUS_IGNORE, errs(ne))
        ne=ne+1; lbls(ne)='write_all phi(2)'
        CALL MPI_FILE_WRITE_ALL(fh, phi(:,:,:,2), 1, type_sub_array, MPI_STATUS_IGNORE, errs(ne))
        ne=ne+1; lbls(ne)='write_all phi(3)'
        CALL MPI_FILE_WRITE_ALL(fh, phi(:,:,:,3), 1, type_sub_array, MPI_STATUS_IGNORE, errs(ne))
        ne=ne+1; lbls(ne)='write_all phi(4)'
        CALL MPI_FILE_WRITE_ALL(fh, phi(:,:,:,4), 1, type_sub_array, MPI_STATUS_IGNORE, errs(ne))
        ne=ne+1; lbls(ne)='write_all phi(5)'
        CALL MPI_FILE_WRITE_ALL(fh, phi(:,:,:,5), 1, type_sub_array, MPI_STATUS_IGNORE, errs(ne))

        call EndProfRange
        call StartProfRange('mpi_write_close')
        ne=ne+1; lbls(ne)='file_close'
        CALL MPI_FILE_CLOSE(fh, errs(ne))
        call EndProfRange

        ne=ne+1; lbls(ne)='type_free loc'
        call MPI_TYPE_FREE(type_sub_array,      errs(ne))
        ne=ne+1; lbls(ne)='type_free view'
        call MPI_TYPE_FREE(type_sub_view_array, errs(ne))

        call check_mpi_vec('mpi_write', errs(1:ne), lbls(1:ne))

        call EndProfRange

        return
end subroutine mpi_write

subroutine mpi_write_stat1D(sid,eid,n,nstats,dirname,vmean, writer_comm, is_writer)
        use FileModule

        implicit none
        integer                              , intent(in) :: sid, eid, n, nstats
        character(*)                         , intent(in) :: dirname
        real(rp), dimension(:,:), allocatable, intent(in) :: vmean
        integer                              , intent(in) :: writer_comm
        logical                              , intent(in) :: is_writer

        integer, parameter :: array_rank=1

        integer, dimension(array_rank)      :: shape_array
        integer, dimension(array_rank)      :: shape_sub_array
        integer, dimension(array_rank)      :: start_coord

        integer, dimension(array_rank)      :: shape_view_array
        integer, dimension(array_rank)      :: shape_sub_view_array
        integer, dimension(array_rank)      :: start_view_coord

        integer(kind = MPI_OFFSET_KIND)     :: initial_displacement

        character(dl) :: iteration, filename, saving_dir
        integer       :: fh, err, l, int_size, header_count, ne
        integer       :: open_err_max, alr_err
        integer       :: type_sub_array, type_sub_view_array

        integer, parameter                :: NIO = 12
        integer, dimension(NIO)           :: errs
        character(len=28), dimension(NIO) :: lbls

        call StartProfRange('mpi_write_stat1D')

        err = 0
        write(iteration, '(i0.7)') itStat
        saving_dir = "DATA/"//trim(data_dir)//"/"//trim(dirname)
        filename = trim(saving_dir)//"/stats"//trim(iteration)//'.bin'

        errs = 0
        ne   = 0

        if(is_writer) then

        call StartProfRange('mpi_write_stat1D_open')

        CALL MPI_FILE_OPEN(writer_comm, filename, &
                 MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, err)

        call MPI_Allreduce(err, open_err_max, 1, MPI_INTEGER, MPI_MAX, writer_comm, alr_err)

        call EndProfRange

        if(open_err_max.ne.0) then
          if(rank == root) write(*,'(A,A)') 'MPI_IO error in writing file ', trim(filename)
          call secure_stop
        endif

        call StartProfRange('mpi_write_stat1D_truncate')
        ne=ne+1; lbls(ne)='file_set_size'
        call MPI_FILE_SET_SIZE(fh, 0_MPI_OFFSET_KIND, errs(ne))
        call EndProfRange

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
        start_view_coord     = (/sid-1/)

        ne=ne+1; lbls(ne)='type_create_subarray view'
        call mpi_type_create_subarray(array_rank, shape_view_array, shape_sub_view_array, start_view_coord, &
                                      MPI_ORDER_FORTRAN, MPI_RP, type_sub_view_array, errs(ne))
        ne=ne+1; lbls(ne)='type_commit view'
        call mpi_type_commit(type_sub_view_array, errs(ne))

        ne=ne+1; lbls(ne)='type_size MPI_INTEGER'
        call MPI_TYPE_SIZE(MPI_INTEGER, int_size, errs(ne))

        header_count = 0
        if(rank == root) header_count = 1
        ne=ne+1; lbls(ne)='write_at_all itStat'
        call MPI_FILE_WRITE_AT_ALL(fh, 0_MPI_OFFSET_KIND, itStat, header_count, MPI_INTEGER, MPI_STATUS_IGNORE, errs(ne))

        initial_displacement = int_size
        call StartProfRange('mpi_write_stat1D_write')
        ne=ne+1; lbls(ne)='file_set_view'
        CALL MPI_FILE_SET_VIEW(fh, initial_displacement, MPI_RP, type_sub_view_array, "native", MPI_INFO_NULL, errs(ne))

        ne=ne+1; lbls(ne)='write_all slabs'
        errs(ne) = 0
        do l = 1,nstats
           CALL MPI_FILE_WRITE_ALL(fh,vmean(:,l),1,type_sub_array,MPI_STATUS_IGNORE,err)
           errs(ne) = ior(errs(ne), abs(err))
        enddo

        call EndProfRange
        call StartProfRange('mpi_write_stat1D_close')
        ne=ne+1; lbls(ne)='file_close'
        CALL MPI_FILE_CLOSE(fh, errs(ne))
        call EndProfRange

        ne=ne+1; lbls(ne)='type_free loc'
        call MPI_TYPE_FREE(type_sub_array,      errs(ne))
        ne=ne+1; lbls(ne)='type_free view'
        call MPI_TYPE_FREE(type_sub_view_array, errs(ne))

        endif

        call check_mpi_vec('mpi_write_stat1D', errs(1:ne), lbls(1:ne))

        call EndProfRange

        return
end subroutine mpi_write_stat1D

subroutine mpi_write_stat2D(sx,ex,sy,ey,nx,ny,vmean, writer_comm, is_writer)
        implicit none
        integer                                , intent(in) :: sx,ex,sy,ey
        integer                                , intent(in) :: nx,ny
        real(rp), dimension(:,:,:), allocatable, intent(in) :: vmean
        integer                                , intent(in) :: writer_comm
        logical                                , intent(in) :: is_writer

        integer, parameter :: array_rank=2

        integer, dimension(array_rank)      :: shape_array
        integer, dimension(array_rank)      :: shape_sub_array
        integer, dimension(array_rank)      :: start_coord

        integer, dimension(array_rank)      :: shape_view_array
        integer, dimension(array_rank)      :: shape_sub_view_array
        integer, dimension(array_rank)      :: start_view_coord

        integer(kind = MPI_OFFSET_KIND)     :: initial_displacement

        character(dl) :: iteration, filename
        integer       :: fh, err, l, int_size, header_count, ne
        integer       :: open_err_max, alr_err
        integer       :: type_sub_array, type_sub_view_array

        integer, parameter                :: NIO = 12
        integer, dimension(NIO)           :: errs
        character(len=28), dimension(NIO) :: lbls

        call StartProfRange('mpi_write_stat2D')

        err  = 0
        errs = 0
        ne   = 0

        if(is_writer) then
        write(iteration, '(i0.7)') itStat
        filename = "DATA/"//trim(data_dir)//"/BINARY_STAT/stats"&
                //trim(iteration)//'.bin'

        call StartProfRange('mpi_write_stat2D_open')

        CALL MPI_FILE_OPEN(writer_comm, filename, &
                 MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, err)

        call MPI_Allreduce(err, open_err_max, 1, MPI_INTEGER, MPI_MAX, writer_comm, alr_err)

        call EndProfRange

        if(open_err_max.ne.0) then
          if(rank == root) write(*,'(A,A)') 'MPI_IO error in writing file ', trim(filename)
          call secure_stop
        endif

        call StartProfRange('mpi_write_stat2D_truncate')

        ne=ne+1; lbls(ne)='file_set_size'
        call MPI_FILE_SET_SIZE(fh, 0_MPI_OFFSET_KIND, errs(ne))

        call EndProfRange

        shape_array     = SHAPE(vmean(:,:,1))
        shape_sub_array = SHAPE(vmean(sx:ex,sy:ey,1))
        start_coord     = (/GN,GN/)

        ne=ne+1; lbls(ne)='type_create_subarray loc'
        call mpi_type_create_subarray(array_rank, shape_array, shape_sub_array, &
                                      start_coord, MPI_ORDER_FORTRAN, MPI_RP, type_sub_array, errs(ne))
        ne=ne+1; lbls(ne)='type_commit loc'
        call mpi_type_commit(type_sub_array, errs(ne))

        shape_view_array     = (/nx,ny/)
        shape_sub_view_array = SHAPE(vmean(sx:ex,sy:ey,1))
        start_view_coord     = (/sx-1,sy-1/)

        ne=ne+1; lbls(ne)='type_create_subarray view'
        call mpi_type_create_subarray(array_rank, shape_view_array, shape_sub_view_array, start_view_coord, &
                                      MPI_ORDER_FORTRAN, MPI_RP, type_sub_view_array, errs(ne))
        ne=ne+1; lbls(ne)='type_commit view'
        call mpi_type_commit(type_sub_view_array, errs(ne))

        ne=ne+1; lbls(ne)='type_size MPI_INTEGER'
        call MPI_TYPE_SIZE(MPI_INTEGER, int_size, errs(ne))

        header_count = 0
        if(rank == root) header_count = 1
        ne=ne+1; lbls(ne)='write_at_all itStat'
        call MPI_FILE_WRITE_AT_ALL(fh, 0_MPI_OFFSET_KIND, itStat, header_count, MPI_INTEGER, MPI_STATUS_IGNORE, errs(ne))

        initial_displacement = int_size

        call StartProfRange('mpi_write_stat2D_write')
        ne=ne+1; lbls(ne)='file_set_view'
        CALL MPI_FILE_SET_VIEW(fh, initial_displacement, MPI_RP, type_sub_view_array, "native", MPI_INFO_NULL, errs(ne))

        ne=ne+1; lbls(ne)='write_all planes'
        errs(ne) = 0
        do l = 1,nvAve2D
           CALL MPI_FILE_WRITE_ALL(fh,vmean(:,:,l),1,type_sub_array,MPI_STATUS_IGNORE,err)
           errs(ne) = ior(errs(ne), abs(err))
        enddo
        call EndProfRange

        call StartProfRange('mpi_write_stat2D_close')
        ne=ne+1; lbls(ne)='file_close'
        CALL MPI_FILE_CLOSE(fh, errs(ne))
        call EndProfRange

        ne=ne+1; lbls(ne)='type_free loc'
        call MPI_TYPE_FREE(type_sub_array,      errs(ne))
        ne=ne+1; lbls(ne)='type_free view'
        call MPI_TYPE_FREE(type_sub_view_array, errs(ne))

        endif

        call check_mpi_vec('mpi_write_stat2D', errs(1:ne), lbls(1:ne))

        call EndProfRange

        return
end subroutine mpi_write_stat2D

end module output_restart
