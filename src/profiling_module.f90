module profiling_module
  use nvtx
#ifdef PROFCPU
  use mpi
  use parameters_module, only: data_dir
  use FileModule,        only: GetFileUnit, FileType
#endif
  implicit none
  private

  public :: StartProfRange, EndProfRange, FinalizeProfiling

#if defined(PROFCPU) && defined(NVTX)
#  error "profiling_module: PROFCPU and NVTX are mutually exclusive"
#endif

#ifdef PROFCPU
  integer, parameter :: PROF_NAME_LEN   = 64
  integer, parameter :: PROF_MAX_RANGES = 512
  integer, parameter :: PROF_MAX_DEPTH  = 64

  type :: prof_range_t
     character(len=PROF_NAME_LEN) :: name     = ' '
     integer(kind=8)              :: ncalls   = 0_8
     real(kind=8)                 :: tot_time = 0.0_8
     real(kind=8)                 :: min_time = huge(1.0_8)
     real(kind=8)                 :: max_time = 0.0_8
  end type prof_range_t

  type(prof_range_t), save :: registry(PROF_MAX_RANGES)
  integer,            save :: n_registered = 0

  integer,      save :: stack_depth       = 0
  integer,      save :: stack_id(PROF_MAX_DEPTH)
  real(kind=8), save :: stack_t0(PROF_MAX_DEPTH)

  integer(kind=8), save :: n_end_without_start = 0_8
  integer(kind=8), save :: n_stack_overflow    = 0_8
  integer(kind=8), save :: n_registry_overflow = 0_8
#endif

contains

  subroutine StartProfRange(char)
    implicit none
    character(*), intent(in) :: char
#ifdef PROFCPU
    integer :: id
#endif

    if (len_trim(char) == 0) return

#ifdef NVTX
    call nvtxStartRange(char)
#endif

#ifdef PROFCPU
    id = prof_register_or_find(char)
    if (id <= 0) return
    if (stack_depth >= PROF_MAX_DEPTH) then
       n_stack_overflow = n_stack_overflow + 1_8
       return
    end if
    stack_depth           = stack_depth + 1
    stack_id(stack_depth) = id
    stack_t0(stack_depth) = MPI_Wtime()
#endif
  end subroutine StartProfRange

  subroutine EndProfRange
    implicit none
#ifdef PROFCPU
    integer      :: id
    real(kind=8) :: dt
#endif

#ifdef NVTX
    call nvtxEndRange
#endif

#ifdef PROFCPU
    if (stack_depth <= 0) then
       n_end_without_start = n_end_without_start + 1_8
       return
    end if
    id = stack_id(stack_depth)
    dt = MPI_Wtime() - stack_t0(stack_depth)
    stack_depth = stack_depth - 1

    registry(id)%ncalls   = registry(id)%ncalls + 1_8
    registry(id)%tot_time = registry(id)%tot_time + dt
    if (dt < registry(id)%min_time) registry(id)%min_time = dt
    if (dt > registry(id)%max_time) registry(id)%max_time = dt
#endif
  end subroutine EndProfRange

  subroutine FinalizeProfiling
    implicit none
#ifdef PROFCPU
    integer, parameter :: master = 0
    integer :: ierr, myrank, nranks, i, j, n_glob, io_err
    real(kind=8)    :: loc_tot(PROF_MAX_RANGES), loc_min(PROF_MAX_RANGES), loc_max(PROF_MAX_RANGES)
    integer(kind=8) :: loc_ncalls(PROF_MAX_RANGES)
    real(kind=8)    :: g_sum(PROF_MAX_RANGES),  g_tmax(PROF_MAX_RANGES)
    real(kind=8)    :: g_cmax(PROF_MAX_RANGES), g_cmin(PROF_MAX_RANGES)
    integer(kind=8) :: g_ncalls(PROF_MAX_RANGES)
    integer         :: order(PROF_MAX_RANGES)
    logical         :: mpi_is_init, dir_exist
    real(kind=8)    :: mean_ms
    type(FileType)  :: prof_file
    character(len=1024) :: log_dir, log_path
    character(len=PROF_NAME_LEN) :: g_names(PROF_MAX_RANGES)
    character(len=PROF_NAME_LEN), allocatable :: my_names(:), all_names(:)
    integer, allocatable :: counts(:), displs(:)

    call MPI_Initialized(mpi_is_init, ierr)
    if (.not. mpi_is_init) return

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)

    allocate(counts(nranks), displs(nranks))
    call MPI_Gather(n_registered, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, &
                    master, MPI_COMM_WORLD, ierr)
    allocate(my_names(max(n_registered,1)))
    do i = 1, n_registered
       my_names(i) = registry(i)%name
    end do
    if (myrank == master) then
       displs(1) = 0
       do i = 2, nranks
          displs(i) = displs(i-1) + counts(i-1)*PROF_NAME_LEN
       end do
       allocate(all_names(max(sum(counts),1)))
    else
       allocate(all_names(1))
    end if
    call MPI_Gatherv(my_names, n_registered*PROF_NAME_LEN, MPI_CHARACTER, &
                     all_names, counts*PROF_NAME_LEN, displs, MPI_CHARACTER, &
                     master, MPI_COMM_WORLD, ierr)
    n_glob = 0
    if (myrank == master) then
       do i = 1, sum(counts)
          if (any(g_names(1:n_glob) == all_names(i))) cycle
          if (n_glob >= PROF_MAX_RANGES) then
             n_registry_overflow = n_registry_overflow + 1_8
             cycle
          end if
          n_glob = n_glob + 1
          g_names(n_glob) = all_names(i)
       end do
    end if
    call MPI_Bcast(n_glob, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if (n_glob <= 0) return
    call MPI_Bcast(g_names, n_glob*PROF_NAME_LEN, MPI_CHARACTER, master, MPI_COMM_WORLD, ierr)

    loc_tot    = 0.0_8
    loc_min    = huge(1.0_8)
    loc_max    = 0.0_8
    loc_ncalls = 0_8
    do i = 1, n_registered
       do j = 1, n_glob
          if (g_names(j) /= registry(i)%name) cycle
          loc_tot(j)    = registry(i)%tot_time
          loc_min(j)    = registry(i)%min_time
          loc_max(j)    = registry(i)%max_time
          loc_ncalls(j) = registry(i)%ncalls
          exit
       end do
    end do

    call MPI_Reduce(loc_tot,    g_sum,    n_glob, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, master, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(loc_tot,    g_tmax,   n_glob, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, master, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(loc_max,    g_cmax,   n_glob, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, master, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(loc_min,    g_cmin,   n_glob, MPI_DOUBLE_PRECISION, &
                    MPI_MIN, master, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(loc_ncalls, g_ncalls, n_glob, MPI_INTEGER8,         &
                    MPI_SUM, master, MPI_COMM_WORLD, ierr)

    if (myrank == master) then
       do i = 1, n_glob; order(i) = i; end do
       call prof_sort_desc(g_tmax, order, n_glob)

       log_dir  = 'DATA/'//trim(data_dir)
       log_path = trim(log_dir)//'/profiling.log'
       inquire(file=trim(log_dir), exist=dir_exist)
       if (.not. dir_exist) call execute_command_line('mkdir -p '//trim(log_dir))

       call GetFileUnit(prof_file)
       open(file=trim(log_path), unit=prof_file%unit, status='replace', &
            action='write', iostat=io_err)
       if (io_err /= 0) then
          prof_file%unit = 6
          write(6,'(A,A,A)') 'WARN: could not open ', trim(log_path), &
                             ' -- profiling summary follows on stdout'
       end if

       write(prof_file%unit,'(A)') repeat('=', 108)
       write(prof_file%unit,'(A)') &
         'URANOS CPU profiling summary (ordered by slowest-rank total)'
       write(prof_file%unit,'(A,I0)')        '  MPI ranks          : ', nranks
       write(prof_file%unit,'(A,I0,A,I0,A)') '  Ranges registered  : ', &
                                             n_glob, '  (max ', PROF_MAX_RANGES, ')'
       if (n_end_without_start > 0_8 .or. n_stack_overflow > 0_8 &
           .or. n_registry_overflow > 0_8) then
          write(prof_file%unit,'(A,I0,A,I0,A,I0)')                      &
            '  WARN unbalanced End=', n_end_without_start,              &
            '  stack-overflow=',      n_stack_overflow,                 &
            '  registry-overflow=',   n_registry_overflow
       end if
       write(prof_file%unit,'(A)') repeat('-', 108)
       write(prof_file%unit,'(A6,1X,A32,1X,A14,1X,A12,1X,A12,1X,A10,1X,A10,1X,A10)') &
          'rank','range','calls_total','tot_max[s]','tot_sum[s]',       &
          'mean[ms]','min[ms]','max[ms]'
       write(prof_file%unit,'(A)') repeat('-', 108)
       do j = 1, n_glob
          i = order(j)
          if (g_ncalls(i) <= 0_8) cycle
          mean_ms = 1.0e3_8 * g_sum(i) / real(g_ncalls(i), kind=8)
          write(prof_file%unit,                                              &
             '(I6,1X,A32,1X,I14,1X,F12.4,1X,F12.4,1X,F10.3,1X,F10.3,1X,F10.3)') &
               j,                                               &
               g_names(i),                                      &
               g_ncalls(i),                                     &
               g_tmax(i),                                       &
               g_sum(i),                                        &
               mean_ms,                                         &
               1.0e3_8 * g_cmin(i),                             &
               1.0e3_8 * g_cmax(i)
       end do
       write(prof_file%unit,'(A)') repeat('=', 108)

       if (io_err == 0) then
          flush(prof_file%unit)
          close(prof_file%unit)
          write(6,'(A,A)') ' Profiling summary written to ', trim(log_path)
       end if
    end if

    n_registered        = 0
    stack_depth         = 0
    n_end_without_start = 0_8
    n_stack_overflow    = 0_8
    n_registry_overflow = 0_8
#endif
  end subroutine FinalizeProfiling

#ifdef PROFCPU
  function prof_register_or_find(char) result(id)
    implicit none
    character(*), intent(in) :: char
    integer :: id, i
    do i = 1, n_registered
       if (trim(registry(i)%name) == trim(char)) then
          id = i
          return
       end if
    end do
    if (n_registered >= PROF_MAX_RANGES) then
       n_registry_overflow = n_registry_overflow + 1_8
       id = 0
       return
    end if
    n_registered = n_registered + 1
    registry(n_registered)%name     = char
    registry(n_registered)%ncalls   = 0_8
    registry(n_registered)%tot_time = 0.0_8
    registry(n_registered)%min_time = huge(1.0_8)
    registry(n_registered)%max_time = 0.0_8
    id = n_registered
  end function prof_register_or_find

  subroutine prof_sort_desc(key, order, n)
    implicit none
    real(kind=8), intent(in)    :: key(:)
    integer,      intent(inout) :: order(:)
    integer,      intent(in)    :: n
    integer :: i, j, tmp
    do i = 2, n
       tmp = order(i)
       j = i - 1
       do while (j >= 1)
          if (key(order(j)) >= key(tmp)) exit
          order(j+1) = order(j)
          j = j - 1
       end do
       order(j+1) = tmp
    end do
  end subroutine prof_sort_desc
#endif

end module profiling_module
