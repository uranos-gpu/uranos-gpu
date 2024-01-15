module mpi_module
use parameters_module
use mpi
#ifdef _OPENACC
use openacc
#endif
use allocate_module
implicit none

! cpu variables
integer            :: rank              !< process rank
integer            :: nprocs            !< number of precesses
integer            :: mpi_err  = 0      !< mpi handling error
integer, parameter :: root     = 0      !< root process
logical            :: mpi_flag = .false.!< mpi_handle

! only openMP required
!$ integer            :: required, provided

! topology variables
integer                   :: mpi_comm_cart     !< cartesian communicator
integer, parameter        :: ndims = 3         !< number of dimension
integer, parameter        :: ndims_cons = 4    !< number of dimension of conservative variables

integer :: comm1Dx , comm1Dy , comm1Dz
integer :: comm2Dxy, comm2Dxz, comm2Dyz

integer, dimension(ndims) :: cart_coord        !< coordinates of the local domain
logical, dimension(ndims) :: periods           !< Topology periodicity

! neighbours variables
integer, parameter           :: nb_neigh = 6
integer, parameter, private  :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
integer, dimension(nb_neigh) :: my_neighbour

! number of ghost nodes
integer, parameter :: GN = 4

! integer inner subdomain positions
integer :: sx, ex       !< start x, end x
integer :: sy, ey       !< start y, end y
integer :: sz, ez       !< start z, end z

! derived data type
integer :: lbx, ubx     !< lower bound and upper bound along x (sx-GN, ex+GN)
integer :: lby, uby     !< lower bound and upper bound along y (sy-GN, ey+GN)
integer :: lbz, ubz     !< lower bound and upper bound along z (sz-GN, ez+GN)

! derived data type
integer, dimension(nb_neigh) :: type_send_cons, type_recv_cons  !< mpi derived data for conservative variable
integer, dimension(nb_neigh) :: type_send_flag, type_recv_flag  !< mpi derived data for weno flag
integer, dimension(nb_neigh) :: type_send_prim, type_recv_prim  !< mpi derived data for primitive variable

! buffers for integer flags communications
integer(1), allocatable, dimension(:,:,:) :: int_bfr_send_E, int_bfr_send_W, int_bfr_recv_E, int_bfr_recv_W
integer(1), allocatable, dimension(:,:,:) :: int_bfr_send_N, int_bfr_send_S, int_bfr_recv_N, int_bfr_recv_S
integer(1), allocatable, dimension(:,:,:) :: int_bfr_send_B, int_bfr_send_F, int_bfr_recv_B, int_bfr_recv_F

! buffer for conservative variable communications
real(rp), allocatable, dimension(:,:,:,:) :: phi_bfr_send_E, phi_bfr_send_W, phi_bfr_recv_E, phi_bfr_recv_W
real(rp), allocatable, dimension(:,:,:,:) :: phi_bfr_send_N, phi_bfr_send_S, phi_bfr_recv_N, phi_bfr_recv_S
real(rp), allocatable, dimension(:,:,:,:) :: phi_bfr_send_B, phi_bfr_send_F, phi_bfr_recv_B, phi_bfr_recv_F

! buffer for 3D variable
real(rp), allocatable, dimension(:,:,:) :: bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W
real(rp), allocatable, dimension(:,:,:) :: bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S
real(rp), allocatable, dimension(:,:,:) :: bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F

! buffer for 3D int1 variables
integer(1), allocatable, dimension(:,:,:) :: i13D_bfr_send_E, i13D_bfr_send_W, i13D_bfr_recv_E, i13D_bfr_recv_W
integer(1), allocatable, dimension(:,:,:) :: i13D_bfr_send_N, i13D_bfr_send_S, i13D_bfr_recv_N, i13D_bfr_recv_S
integer(1), allocatable, dimension(:,:,:) :: i13D_bfr_send_B, i13D_bfr_send_F, i13D_bfr_recv_B, i13D_bfr_recv_F

! statuses tables
integer, allocatable, dimension(:) :: req_array_xx, req_array_yz

integer :: num_dev, dev_id

contains
subroutine init_mpi
! ------------------------------------------------------
!       MPI initialization
! ------------------------------------------------------
        implicit none
#ifdef _OPENACC
        integer :: local_comm, local_rank
#endif

#ifdef OPENMP
        required = MPI_THREAD_FUNNELED
        call MPI_Init_thread(required, provided, mpi_err)
#else
        call MPI_Init(mpi_err)
#endif
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, mpi_err)
        call MPI_Comm_size(MPI_COMM_WORLD, nprocs, mpi_err)

        mpi_flag = .true.

#ifdef _OPENACC
        call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, &
                                 0, MPI_INFO_NULL, local_comm,mpi_err)
        call MPI_Comm_rank(local_comm, local_rank, mpi_err)
#ifdef NVIDIA
        call acc_init(acc_device_nvidia)
        num_dev = acc_get_num_devices(acc_device_nvidia)
        dev_id = mod(local_rank,num_dev)
        call acc_set_device_num(dev_id,acc_device_nvidia)
#endif
#ifdef AMD
        call acc_init(acc_device_host)
        dev_id = local_rank
        call acc_set_device_num(dev_id,acc_device_host)
#endif
        print*, 'Rank ', local_rank, 'is associated to GPU', dev_id
#endif

        call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
        
        return
end subroutine init_mpi


subroutine init_cartesian_topology
! ------------------------------------------------------
!
!       Creation of the Cartesian topology and distribution of the
!       procs over the sub-domains.
!       
! ------------------------------------------------------
        implicit none
        integer :: i
        logical, parameter :: reorganisation = .false.


        ! check for periodic boundary conditions ------------------------------------
        do i = 1,5,2
           if(trim(bc(i)) == 'periodic') then
             if(trim(bc(i+1)) .ne. trim(bc(i))) then
               if(rank == root) print*, ' Periodic boundary condition must be in pair.'
               call secure_stop
             endif
           endif
        enddo

        ! setting periodic directions
        periods(:) = .false.
        if(bc(W) == bc(E) .and. bc(W) == 'periodic') then
                periods(1) = .true.
        endif
        if(bc(S) == bc(N) .and. bc(S) == 'periodic') then
                periods(2) = .true.
        endif
        if(bc(B) == bc(F) .and. bc(B) == 'periodic') then
                periods(3) = .true.
        endif

        ! ================== Find the number of processes on each dimension ========================== !
        if(rank == root) then
                write(*,*) '--------------------------------------------------- ' 
                write(*,'(1x,A,i4,A)') ' Executing URANOS with', nprocs, ' procs.'
        endif

        if    (sum(cart_dims) /= 0) then
                
                ! check for 2D problems
                if(dims == 2 .and. cart_dims(3) /= 1) then
                   if(rank == root) write(*,*) ' ERROR: 2-dimensional problems accept just one proc along-z.'
                   call MPI_FINALIZE(mpi_err)
                   stop
                endif

                ! check for precessors number compatibility
                if(cart_dims(1)*cart_dims(2)*cart_dims(3) /= nprocs) then
                   if(rank == root) write(*,*) ' ERROR: The number of processors specified is not compatible.'
                   call MPI_FINALIZE(mpi_err)
                   stop
                endif

                if(rank == root) print*, ' The processors has distributed by USERs specification as follow:'

        elseif(sum(cart_dims) == 0) then
                
                if(dims == 2) cart_dims(3) = 1

                if(rank == root) print*, ' The processors will be AUTOMATICALLY distribuited as follow:'
                
        endif
        
        call MPI_dims_create(nprocs,ndims,cart_dims,mpi_err)

        ! Creation of the Cartesian Communicator over the cartesian processes --------
        mpi_comm_cart = MPI_COMM_NULL
        CALL MPI_cart_create(MPI_COMM_WORLD, ndims, cart_dims, periods, &
                reorganisation, mpi_comm_cart, mpi_err)

        call check_mpi('init_cartesian_topology',mpi_err)
        
        return
end subroutine init_cartesian_topology



subroutine init_subcartesian_topology
        implicit none
        logical, dimension(3) :: TFF, FTF, FFT, TTF, TFT, FTT
        logical, parameter    :: T = .true., F = .false.
        integer               :: err = 0

        TFF(1) = T
        TFF(2) = F
        TFF(3) = F

        FTF(1) = F
        FTF(2) = T
        FTF(3) = F

        FFT(1) = F
        FFT(2) = F
        FFT(3) = T

        ! 1D sub-communicators
        comm1Dx = MPI_COMM_NULL
        comm1Dy = MPI_COMM_NULL
        comm1Dz = MPI_COMM_NULL
       
        call MPI_CART_SUB(mpi_comm_cart, TFF, comm1Dx, err)
        call MPI_CART_SUB(mpi_comm_cart, FTF, comm1Dy, err)
        call MPI_CART_SUB(mpi_comm_cart, FFT, comm1Dz, err)

        TTF(1) = T
        TTF(2) = T
        TTF(3) = F

        TFT(1) = T
        TFT(2) = F
        TFT(3) = T

        FTT(1) = F
        FTT(2) = T
        FTT(3) = T

        ! 2D sub-communicators
        comm2Dxy = MPI_COMM_NULL
        comm2Dxz = MPI_COMM_NULL
        comm2Dyz = MPI_COMM_NULL

        call MPI_CART_SUB(mpi_comm_cart, TTF, comm2Dxy, err)
        call MPI_CART_SUB(mpi_comm_cart, TFT, comm2Dxz, err)
        call MPI_CART_SUB(mpi_comm_cart, FTT, comm2Dyz, err)

        return
end subroutine init_subcartesian_topology


subroutine init_subdomain_boundaries
! ------------------------------------------------------
!       Computation of the local grid boundary indexes
! ------------------------------------------------------
        implicit none
        integer :: min_nodes_lcl, min_nodes_gbl
        integer :: min_ghost_lcl, min_ghost_gbl

        call MPI_cart_coords(mpi_comm_cart,rank,ndims,cart_coord,mpi_err)
        
        ! x-axis inner limits
        sx = (cart_coord(1)*nx)/cart_dims(1)+1
        ex = ((cart_coord(1)+1)*nx)/cart_dims(1)

        ! y-axis inner limits
        sy = (cart_coord(2)*ny)/cart_dims(2)+1
        ey = ((cart_coord(2)+1)*ny)/cart_dims(2)

        ! z-axis inner limits
        if(dims == 2) nz = 1  !< ensure the maximum speed up for 2D problems

        sz = (cart_coord(3)*nz)/cart_dims(3)+1
        ez = ((cart_coord(3)+1)*nz)/cart_dims(3)

        ! computing lower and upper bounds
        lbx = sx - GN 
        ubx = ex + GN 

        lby = sy - GN 
        uby = ey + GN 

        lbz = sz - GN 
        ubz = ez + GN 

        !! check node distribution
        !if(rank == root) print*, 'rank ','sx ','ex ','sy ','ey ','sz ','ez'
        !write(*,'(7i4)')  rank, sx, ex, sy, ey, sz, ez

        !
        ! ==== Compute the proc with maximum ghost/node ratio
        !
        if    (dims == 2) then
          min_nodes_lcl = (ex-sx+1)*(ey-sy+1)
          min_ghost_lcl = (ex-sx+1)*2*GN + (ey-sy)*2*GN
        elseif(dims == 3) then
          min_nodes_lcl = (ex-sx+1)*(ey-sy+1)*(ez-sz+1)
          min_ghost_lcl = (ex-sx+1)*(ey-sy+1)*2*GN + (ex-sx+1)*(ez-sz+1)*2*GN + (ez-sz+1)*(ey-sy+1)*2*GN
        endif
        call MPI_allreduce(min_nodes_lcl, min_nodes_gbl, 1, MPI_INT, MPI_MIN, mpi_comm_cart, mpi_err)
        call MPI_allreduce(min_ghost_lcl, min_ghost_gbl, 1, MPI_INT, MPI_MIN, mpi_comm_cart, mpi_err)

        if(rank == root) then
                write(*,*)
                write(*,'(1x,A,i4)')   ' Problem dimensions  : ', dims
                write(*,*)
                write(*,'(1x,A,i4)')   ' Procs in x-direction: ', cart_dims(1)
                write(*,'(1x,A,i4)')   ' Procs in y-direction: ', cart_dims(2)
                write(*,'(1x,A,i4)')   ' Procs in z-direction: ', cart_dims(3)
                write(*,*)
                write(*,'(1x,A,i4)')   ' Point in x-direction: ', nx
                write(*,'(1x,A,i4)')   ' Point in y-direction: ', ny
                write(*,'(1x,A,i4)')   ' Point in z-direction: ', nz
                write(*,*)
                write(*,*)   ' Min local nodes num : ', min_nodes_gbl
                write(*,*)   ' Min local ghost num : ', min_ghost_gbl
                write(*,*)   ' Ghost/Nodes ratio % : ', real(min_ghost_gbl,rp)/real(min_nodes_gbl,rp)*100
                write(*,*)
                if(mpi_opt_level == 1) write(*,*) ' MPI communications handled with SENDRECV'
                if(mpi_opt_level == 2) write(*,*) ' MPI communications handled with ISEND/IRECV + BFR'
                if(mpi_opt_level == 3) write(*,*) ' MPI communications handled with ISEND/IRECV + DDT'
                if(cuda_aware)         write(*,*) ' CUDA-AWARE MPI enabled'
                write(*,*)
                write(*,'(1x,A,A)')    ' Advection scheme    : ', trim(scheme)
                write(*,'(1x,A,A)')    ' Diffusion scheme    : ', trim(diffusion_scheme)
                write(*,'(1x,A,i4)')   ' FD accurary order   : ', fd_order
                write(*,'(1x,A,i4)')   ' WENO accurary order : ', weno_order
                if(viscous) then
                write(*,'(1x,A,e18.6)') ' Mach number         : ', Mach
                write(*,'(1x,A,e18.6)') ' Reynolds number     : ', Reynolds
                endif
                if(les) then
                write(*,*)
                write(*,'(1x,A,A)')     ' LES model            : ', trim(sgs_model)
                endif
                write(*,*) '--------------------------------------------------- ' 
        endif



        
        if    (dims == 2) then

          if(ex-sx < 2 .or. ey-sy < 2) then
            if(rank == root) then
                    print*, ' ERROR: the number of points per single proc is too small.'
            endif
            call secure_stop
            stop
          endif

        elseif(dims == 3) then

          if(ex-sx < 2 .or. ey-sy < 2 .or. ez-sz < 2) then
            if(rank == root) then
                    print*, ' ERROR: the number of points per single proc is too small.'
            endif
            call secure_stop
            stop
          endif

        endif

        return
end subroutine init_subdomain_boundaries

subroutine init_mpi_neighbours
! ----------------------------------------------------
!       Find the neighbours of all the procs
! ----------------------------------------------------
        implicit none
        integer, parameter, dimension(3) :: dirs = (/0,1,2/)
        integer, parameter               :: step = 1

        ! initialize neighbours
        my_neighbour(:) = MPI_PROC_NULL

        ! get my estern and western neighbours
        call MPI_cart_shift(mpi_comm_cart, dirs(1), step, my_neighbour(W), my_neighbour(E), mpi_err)

        ! get my southern and northern neighbours
        call MPI_cart_shift(mpi_comm_cart, dirs(2), step, my_neighbour(S), my_neighbour(N), mpi_err)

        ! get my behind and front neighbours
        call MPI_cart_shift(mpi_comm_cart, dirs(3), step, my_neighbour(B), my_neighbour(F), mpi_err)
        
        !if(rank == root) print*, 'rank', 'W', ' E', ' S', ' N', ' B', ' F'
        !print*, rank, my_neighbour(W), my_neighbour(E), my_neighbour(S), my_neighbour(N), my_neighbour(B), my_neighbour(F) 

        return
end subroutine init_mpi_neighbours

subroutine init_mpi_exchange_data
! ------------------------------------------------
!       initialize mpi derived data type
! ------------------------------------------------
        implicit none
        
        ! ==== init mpi derived data data for conservative variables
        call init_mpi_buffers

        call init_mpi_type_field(ndims_cons,5,MPI_RP, type_send_cons, type_recv_cons)
        
        if    (dims == 2) then
          allocate(req_array_xx(4), req_array_yz(4))
        elseif(dims == 3) then
          allocate(req_array_xx(4), req_array_yz(8))
        endif

        ! ==== init mpi derived data type for hybrid weno flag
        if(hybrid_weno) then
          call init_mpi_type_field(ndims,1,MPI_INTEGER1, type_send_flag, type_recv_flag)
        endif

        ! ==== init mpi derived data type for primitive variables in LES
        call init_mpi_type_field(ndims,1,MPI_RP,type_send_prim, type_recv_prim)

        return
end subroutine init_mpi_exchange_data


subroutine init_mpi_type_field(field_dims,field_co_dims,old_field_type, type_field_send, type_field_recv)
! ---------------------------------------------------------------------------------
!
!       This subroutine create the send/recv derived data type for a general 3D field
!       with general old_type.
!
!       input : field_dims               !< field dimension
!               field_co_dims            !< field co-dimension (only for 4D fields)
!               old_field_type           !< old type of the field
!       output: type_field_send          !< handle array for send types
!               type_field_recv          !< handle array for recv types
!
! ---------------------------------------------------------------------------------
        implicit none
        integer                     , intent(in)    :: field_dims
        integer                     , intent(in)    :: field_co_dims
        integer                     , intent(in)    :: old_field_type
        integer, dimension(nb_neigh), intent(inout) :: type_field_send, type_field_recv

        ! 3D fields size and subsizes
        integer, dimension(field_dims) :: old_field_size   !< original array sizes along dirs
        integer, dimension(field_dims) :: size_sub_field_x !< subsize along East   and West
        integer, dimension(field_dims) :: size_sub_field_y !< subsize along South  and North
        integer, dimension(field_dims) :: size_sub_field_z !< subsize along Behind and Front
        
        ! starting position of 3D fields
        integer, dimension(field_dims) :: start_field_send_W, start_field_send_E
        integer, dimension(field_dims) :: start_field_recv_W, start_field_recv_E
        integer, dimension(field_dims) :: start_field_send_S, start_field_send_N
        integer, dimension(field_dims) :: start_field_recv_S, start_field_recv_N
        integer, dimension(field_dims) :: start_field_send_B, start_field_send_F
        integer, dimension(field_dims) :: start_field_recv_B, start_field_recv_F
        
        if(field_dims == 3) then

                old_field_size   = (/gn+ex-sx+1+gn, gn+ey-sy+1+gn, gn+ez-sz+1+gn/)
                size_sub_field_x = (/GN           , gn+ey-sy+1+gn, gn+ez-sz+1+gn/)
                size_sub_field_y = (/gn+ex-sx+1+gn, GN           , gn+ez-sz+1+gn/)
                size_sub_field_z = (/gn+ex-sx+1+gn, gn+ey-sy+1+gn, GN           /)

                ! starting positions
                start_field_send_W = (/GN        , 0, 0/)
                start_field_recv_E = (/ex-sx+1+GN, 0, 0/)
                start_field_send_E = (/ex-sx+1   , 0, 0/)
                start_field_recv_W = (/0         , 0, 0/)

                ! starting positions
                start_field_send_S = (/0, GN        , 0/)
                start_field_recv_N = (/0, ey-sy+1+GN, 0/)
                start_field_send_N = (/0, ey-sy+1   , 0/)
                start_field_recv_S = (/0, 0         , 0/)

                ! starting positions
                start_field_send_B = (/0, 0, GN        /)
                start_field_recv_F = (/0, 0, ez-sz+1+GN/)
                start_field_send_F = (/0, 0, ez-sz+1   /)
                start_field_recv_B = (/0, 0, 0         /)

        elseif(field_dims == 4) then

                ! size of the orginal structure
                old_field_size   = (/gn+ex-sx+1+gn, gn+ey-sy+1+gn, gn+ez-sz+1+gn, field_co_dims/)
                size_sub_field_x = (/GN           , gn+ey-sy+1+gn, gn+ez-sz+1+gn, field_co_dims/)
                size_sub_field_y = (/gn+ex-sx+1+gn, GN           , gn+ez-sz+1+gn, field_co_dims/)
                size_sub_field_z = (/gn+ex-sx+1+gn, gn+ey-sy+1+gn, GN           , field_co_dims/)

                ! starting positions
                start_field_send_W = (/GN        , 0, 0, 0/)
                start_field_recv_E = (/ex-sx+1+GN, 0, 0, 0/)
                start_field_send_E = (/ex-sx+1   , 0, 0, 0/)
                start_field_recv_W = (/0         , 0, 0, 0/)

                ! starting positions
                start_field_send_S = (/0, GN        , 0, 0/)
                start_field_recv_N = (/0, ey-sy+1+GN, 0, 0/)
                start_field_send_N = (/0, ey-sy+1   , 0, 0/)
                start_field_recv_S = (/0, 0         , 0, 0/)

                ! starting positions
                start_field_send_B = (/0, 0, GN        , 0/)
                start_field_recv_F = (/0, 0, ez-sz+1+GN, 0/)
                start_field_send_F = (/0, 0, ez-sz+1   , 0/)
                start_field_recv_B = (/0, 0, 0         , 0/)

        endif

        ! creation of subarray for East-West communications
        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_x,start_field_send_E,old_field_type,type_field_send(E))
        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_x,start_field_send_W,old_field_type,type_field_send(W))

        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_x,start_field_recv_E,old_field_type,type_field_recv(E))
        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_x,start_field_recv_W,old_field_type,type_field_recv(W))

        ! creation of subarray for South-North communications
        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_y,start_field_send_S,old_field_type,type_field_send(S))
        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_y,start_field_send_N,old_field_type,type_field_send(N))

        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_y,start_field_recv_S,old_field_type,type_field_recv(S))
        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_y,start_field_recv_N,old_field_type,type_field_recv(N))

        ! creation of subarray for Backward-Farward communications
        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_z,start_field_send_B,old_field_type,type_field_send(B))
        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_z,start_field_send_F,old_field_type,type_field_send(F))

        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_z,start_field_recv_B,old_field_type,type_field_recv(B))
        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_z,start_field_recv_F,old_field_type,type_field_recv(F))

        return
end subroutine init_mpi_type_field


subroutine MPI_create_subarray(dims,old_size,subsize,start_position,old_type,new_type)
        implicit none
        integer              , intent(in) :: dims
        integer, dimension(:), intent(in) :: old_size
        integer, dimension(:), intent(in) :: subsize
        integer, dimension(:), intent(in) :: start_position
        integer              , intent(in) :: old_type
        integer                           :: new_type
        
        ! create the subarray
        call MPI_type_create_subarray(dims,old_size, subsize, start_position , &
                                      MPI_ORDER_FORTRAN, old_type, new_type  , &
                                      mpi_err)
        ! commit the subarray
        call MPI_type_commit(new_type,mpi_err)

        ! mpi check error
        call check_mpi('create subarray', mpi_err)

        return
end subroutine MPI_create_subarray






subroutine mpi_wait_procs(req_array)

        implicit none
        integer, allocatable, dimension(:), intent(inout) :: req_array
        integer                                           :: err = 0

        call MPI_waitall(size(req_array), req_array, MPI_STATUSES_IGNORE, err)

        if(err /= 0) then
          print*, ' MPI DEADLOCK ... The program will be stopped!'
          stop
        endif

        return
end subroutine mpi_wait_procs




subroutine init_mpi_buffers
        implicit none
        
        !4D
        call AllocateReal(phi_bfr_send_W,sx,sx+GN-1, lby,uby, lbz,ubz, 1,5)
        call AllocateReal(phi_bfr_send_E,ex-GN+1,ex, lby,uby, lbz,ubz, 1,5)
        call AllocateReal(phi_bfr_recv_W,sx-GN,sx-1, lby,uby, lbz,ubz, 1,5)
        call AllocateReal(phi_bfr_recv_E,ex+1,ex+GN, lby,uby, lbz,ubz, 1,5)

        call AllocateReal(phi_bfr_send_S,lbx,ubx, sy,sy+GN-1, lbz,ubz, 1,5)
        call AllocateReal(phi_bfr_send_N,lbx,ubx, ey-GN+1,ey, lbz,ubz, 1,5)
        call AllocateReal(phi_bfr_recv_S,lbx,ubx, sy-GN,sy-1, lbz,ubz, 1,5)
        call AllocateReal(phi_bfr_recv_N,lbx,ubx, ey+1,ey+GN, lbz,ubz, 1,5)

        call AllocateReal(phi_bfr_send_B,lbx,ubx, lby,uby, sz,sz+GN-1, 1,5)
        call AllocateReal(phi_bfr_send_F,lbx,ubx, lby,uby, ez-GN+1,ez, 1,5)
        call AllocateReal(phi_bfr_recv_B,lbx,ubx, lby,uby, sz-GN,sz-1, 1,5)
        call AllocateReal(phi_bfr_recv_F,lbx,ubx, lby,uby, ez+1,ez+GN, 1,5)

        !3D
        call AllocateReal(bfr_send_W,sx,sx+GN-1, lby,uby, lbz,ubz)
        call AllocateReal(bfr_send_E,ex-GN+1,ex, lby,uby, lbz,ubz)
        call AllocateReal(bfr_recv_W,sx-GN,sx-1, lby,uby, lbz,ubz)
        call AllocateReal(bfr_recv_E,ex+1,ex+GN, lby,uby, lbz,ubz)

        call AllocateReal(bfr_send_S,lbx,ubx, sy,sy+GN-1, lbz,ubz)
        call AllocateReal(bfr_send_N,lbx,ubx, ey-GN+1,ey, lbz,ubz)
        call AllocateReal(bfr_recv_S,lbx,ubx, sy-GN,sy-1, lbz,ubz)
        call AllocateReal(bfr_recv_N,lbx,ubx, ey+1,ey+GN, lbz,ubz)

        call AllocateReal(bfr_send_B,lbx,ubx, lby,uby, sz,sz+GN-1)
        call AllocateReal(bfr_send_F,lbx,ubx, lby,uby, ez-GN+1,ez)
        call AllocateReal(bfr_recv_B,lbx,ubx, lby,uby, sz-GN,sz-1)
        call AllocateReal(bfr_recv_F,lbx,ubx, lby,uby, ez+1,ez+GN)

        ! 3D int1
        call AllocateInteger(i13D_bfr_send_W,sx,sx+GN-1, lby,uby, lbz,ubz)
        call AllocateInteger(i13D_bfr_send_E,ex-GN+1,ex, lby,uby, lbz,ubz)
        call AllocateInteger(i13D_bfr_recv_W,sx-GN,sx-1, lby,uby, lbz,ubz)
        call AllocateInteger(i13D_bfr_recv_E,ex+1,ex+GN, lby,uby, lbz,ubz)

        call AllocateInteger(i13D_bfr_send_S,lbx,ubx, sy,sy+GN-1, lbz,ubz)
        call AllocateInteger(i13D_bfr_send_N,lbx,ubx, ey-GN+1,ey, lbz,ubz)
        call AllocateInteger(i13D_bfr_recv_S,lbx,ubx, sy-GN,sy-1, lbz,ubz)
        call AllocateInteger(i13D_bfr_recv_N,lbx,ubx, ey+1,ey+GN, lbz,ubz)

        call AllocateInteger(i13D_bfr_send_B,lbx,ubx, lby,uby, sz,sz+GN-1)
        call AllocateInteger(i13D_bfr_send_F,lbx,ubx, lby,uby, ez-GN+1,ez)
        call AllocateInteger(i13D_bfr_recv_B,lbx,ubx, lby,uby, sz-GN,sz-1)
        call AllocateInteger(i13D_bfr_recv_F,lbx,ubx, lby,uby, ez+1,ez+GN)
                 
        return
end subroutine init_mpi_buffers







subroutine mpi_safe_allreduce(local,glbal,mpi_operation,mpi_comm,err)
! ------------------------------------------------------------------------
!       
!       This subroutine implements a binary safe all reduce operation
!       splitting it in three mpi different operations.
!
! ------------------------------------------------------------------------
        implicit none
        real(rp), intent(in)  :: local
        integer , intent(in)  :: mpi_operation
        integer , intent(in)  :: mpi_comm
        real(rp), intent(out) :: glbal
        integer , intent(out) :: err

        real(rp), dimension(nprocs) :: glb_a
        
        ! gather local to global array
        call mpi_gather(local,1,MPI_RP, &
                        glb_a,1,MPI_RP, &
                        root, mpi_comm,err)
        
        ! perform the computation just on root
        if(rank == root) then

          selectcase(mpi_operation)

                case(mpi_min)
                    glbal = minval(glb_a)

                case(mpi_max)
                    glbal = maxval(glb_a)

                case(mpi_sum)
                    glbal = sum(glb_a)

          endselect

        endif

        ! broadcast the computation
        call mpi_bcast(glbal,1,MPI_RP,root,mpi_comm,err)

        return
end subroutine mpi_safe_allreduce







subroutine mpi_communicate4D_int1(tp_send,tp_recv,var)
        implicit none
        integer(1), allocatable, dimension(:,:,:,:) , intent(inout) :: var
        integer                , dimension(nb_neigh), intent(in)    :: tp_send
        integer                , dimension(nb_neigh), intent(in)    :: tp_recv

        ! local declaration
        integer, parameter                             :: n_request1 = 4, n_request2 = 8, tag = 0
        integer, dimension(n_request1)                 :: request1
        integer, dimension(n_request2)                 :: request2
        integer, dimension(MPI_STATUS_SIZE,n_request1) :: tab_stat1
        integer, dimension(MPI_STATUS_SIZE,n_request2) :: tab_stat2

        ! east - west communications
        call MPI_ISEND(var,1,tp_send(E),my_neighbour(E),tag+1,mpi_comm_cart,request1(1),mpi_err)
        call MPI_IRECV(var,1,tp_recv(W),my_neighbour(W),tag+1,mpi_comm_cart,request1(2),mpi_err)

        call MPI_ISEND(var,1,tp_send(W),my_neighbour(W),tag+2,mpi_comm_cart,request1(3),mpi_err)
        call MPI_IRECV(var,1,tp_recv(E),my_neighbour(E),tag+2,mpi_comm_cart,request1(4),mpi_err)

        ! wait till communications end along x cause corners
        call MPI_waitall(n_request1, request1, tab_stat1, mpi_err)

        ! nord-south communications
        call MPI_ISEND(var,1,tp_send(N),my_neighbour(N),tag+3,mpi_comm_cart,request2(1),mpi_err)
        call MPI_IRECV(var,1,tp_recv(S),my_neighbour(S),tag+3,mpi_comm_cart,request2(2),mpi_err)

        call MPI_ISEND(var,1,tp_send(S),my_neighbour(S),tag+4,mpi_comm_cart,request2(3),mpi_err)
        call MPI_IRECV(var,1,tp_recv(N),my_neighbour(N),tag+4,mpi_comm_cart,request2(4),mpi_err)

        ! backward-forwvar communications
        call MPI_ISEND(var,1,tp_send(F),my_neighbour(F),tag+5,mpi_comm_cart,request2(5),mpi_err)
        call MPI_IRECV(var,1,tp_recv(B),my_neighbour(B),tag+5,mpi_comm_cart,request2(6),mpi_err)

        call MPI_ISEND(var,1,tp_send(B),my_neighbour(B),tag+6,mpi_comm_cart,request2(7),mpi_err)
        call MPI_IRECV(var,1,tp_recv(F),my_neighbour(F),tag+6,mpi_comm_cart,request2(8),mpi_err)
         
        ! wait till communications end
        call MPI_waitall(n_request2, request2, tab_stat2, mpi_err)

        return
end subroutine mpi_communicate4D_int1



subroutine init_mpi_integer_buffers
        implicit none
        integer :: err = 0
        integer, dimension(6) :: sS, eS
        integer, dimension(6) :: sR, eR

        ! 
        ! === indeces
        !
        sS = (/sx-GN, ex+1 , sy-GN, ey+1 , sz-GN, ez+1 /)             
        eS = (/sx-1 , ex+GN, sy-1 , ey+GN, sz-1 , ez+GN/)             

        sR = (/sx     , ex-GN+1, sy     , ey-GN+1, sz     , ez-GN+1/)             
        eR = (/sx+GN-1, ex     , sy+GN-1, ey     , sz+GN-1, ez     /)


        allocate(int_bfr_send_W(sS(W):eS(W), lby:uby    , lbz:ubz    ), &
                 int_bfr_send_E(sS(E):eS(E), lby:uby    , lbz:ubz    ), &
                 int_bfr_send_S(lbx:ubx    , sS(S):eS(S), lbz:ubz    ), &
                 int_bfr_send_N(lbx:ubx    , sS(N):eS(N), lbz:ubz    ), &
                 int_bfr_send_B(lbx:ubx    , lby:uby    , sS(B):eS(B)), &
                 int_bfr_send_F(lbx:ubx    , lby:uby    , sS(F):eS(F)), &

                 int_bfr_recv_W(sR(W):eR(W), lby:uby    , lbz:ubz    ), &
                 int_bfr_recv_E(sR(E):eR(E), lby:uby    , lbz:ubz    ), &
                 int_bfr_recv_S(lbx:ubx    , sR(S):eR(S), lbz:ubz    ), &
                 int_bfr_recv_N(lbx:ubx    , sR(N):eR(N), lbz:ubz    ), &
                 int_bfr_recv_B(lbx:ubx    , lby:uby    , sR(B):eR(B)), &
                 int_bfr_recv_F(lbx:ubx    , lby:uby    , sR(F):eR(F)), stat = err)

        if(err .ne. 0) then
          if(rank == root) print*, ' Allocation error in init_mpi_integer_buffers'
          call secure_stop
        endif

        return
end subroutine init_mpi_integer_buffers


subroutine communicate_int_flag(my_int_flag)
        implicit none
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: my_int_flag
        
        ! local declarations
        integer, parameter    :: tag = 0, dtype = mpi_integer1
        integer, dimension(3) :: bfr_dims
        integer, dimension(6) :: sS, eS
        integer, dimension(6) :: sR, eR
        integer               :: bfr_size, err = 0
        
        !
        ! === buffer dimensions
        !
        bfr_dims = shape(int_bfr_send_E)
        bfr_size = bfr_dims(1)*bfr_dims(2)*bfr_dims(3)

        ! 
        ! === indeces
        !
        sS = (/sx-GN, ex+1 , sy-GN, ey+1 , sz-GN, ez+1 /)             
        eS = (/sx-1 , ex+GN, sy-1 , ey+GN, sz-1 , ez+GN/)             

        sR = (/sx     , ex-GN+1, sy     , ey-GN+1, sz     , ez-GN+1/)             
        eR = (/sx+GN-1, ex     , sy+GN-1, ey     , sz+GN-1, ez     /)


        !
        ! === E/W communications
        !
        int_bfr_send_W = my_int_flag(sS(W):eS(W),:,:)
        call MPI_ISEND(int_bfr_send_W,bfr_size,dtype,my_neighbour(W),tag+1,mpi_comm_cart,req_array_xx(1), err)
        call MPI_IRECV(int_bfr_recv_E,bfr_size,dtype,my_neighbour(E),tag+1,mpi_comm_cart,req_array_xx(2), err)

        int_bfr_send_E = my_int_flag(sS(E):eS(E),:,:)
        call MPI_ISEND(int_bfr_send_E,bfr_size,dtype,my_neighbour(E),tag+2,mpi_comm_cart,req_array_xx(3), err)
        call MPI_IRECV(int_bfr_recv_W,bfr_size,dtype,my_neighbour(W),tag+2,mpi_comm_cart,req_array_xx(4), err)

        ! wait till communications end (corner needed)
        call mpi_wait_procs(req_array_xx)

        !
        ! === N/S communications
        !
        int_bfr_send_S = my_int_flag(:,sS(S):eS(S),:)
        call MPI_ISEND(int_bfr_send_S,bfr_size,dtype,my_neighbour(S),tag+3,mpi_comm_cart,req_array_yz(1), err)
        call MPI_IRECV(int_bfr_recv_N,bfr_size,dtype,my_neighbour(N),tag+3,mpi_comm_cart,req_array_yz(2), err)

        int_bfr_send_N = my_int_flag(:,sS(N):eS(N),:)
        call MPI_ISEND(int_bfr_send_N,bfr_size,dtype,my_neighbour(N),tag+4,mpi_comm_cart,req_array_yz(3), err)
        call MPI_IRECV(int_bfr_recv_S,bfr_size,dtype,my_neighbour(S),tag+4,mpi_comm_cart,req_array_yz(4), err)

        !
        ! === B/F communications
        !
        if(dims == 3) then
          int_bfr_send_B = my_int_flag(:,:,sS(B):eS(B))
          call MPI_ISEND(int_bfr_send_B,bfr_size,dtype,my_neighbour(B),tag+5,mpi_comm_cart,req_array_yz(5), err)
          call MPI_IRECV(int_bfr_recv_F,bfr_size,dtype,my_neighbour(F),tag+5,mpi_comm_cart,req_array_yz(6), err)

          int_bfr_send_F = my_int_flag(:,:,sS(F):eS(F))
          call MPI_ISEND(int_bfr_send_F,bfr_size,dtype,my_neighbour(F),tag+6,mpi_comm_cart,req_array_yz(7), err)
          call MPI_IRECV(int_bfr_recv_B,bfr_size,dtype,my_neighbour(B),tag+6,mpi_comm_cart,req_array_yz(8), err)
        endif

        ! wait till communications end
        call mpi_wait_procs(req_array_yz)
        
        if(my_neighbour(W) > 0) then
        my_int_flag(sR(W):eR(W),:,:) = min(my_int_flag(sR(W):eR(W),:,:), int_bfr_recv_W(sR(W):eR(W),:,:))
        endif
        if(my_neighbour(E) > 0) then
        my_int_flag(sR(E):eR(E),:,:) = min(my_int_flag(sR(E):eR(E),:,:), int_bfr_recv_E(sR(E):eR(E),:,:))
        endif
        
        if(my_neighbour(S) > 0) then
        my_int_flag(:,sR(S):eR(S),:) = min(my_int_flag(:,sR(S):eR(S),:), int_bfr_recv_S(:,sR(S):eR(S),:))
        endif
        if(my_neighbour(N) > 0) then
        my_int_flag(:,sR(N):eR(N),:) = min(my_int_flag(:,sR(N):eR(N),:), int_bfr_recv_N(:,sR(N):eR(N),:))
        endif

        if(my_neighbour(B) > 0) then
        my_int_flag(:,:,sR(B):eR(B)) = min(my_int_flag(:,:,sR(B):eR(B)), int_bfr_recv_B(:,:,sR(B):eR(B)))
        endif
        if(my_neighbour(F) > 0) then
        my_int_flag(:,:,sR(F):eR(F)) = min(my_int_flag(:,:,sR(F):eR(F)), int_bfr_recv_F(:,:,sR(F):eR(F)))
        endif

        return
end subroutine communicate_int_flag


subroutine mpi_flag_splitting(lflag,old_comm,old_rank,new_comm,new_rank,new_nprocs,err)
! ------------------------------------------------------------------------------------
!       This finds whichs processes are flagged inside an old mpi
!       communicator and builds a new mpi communicator for them. 
! ------------------------------------------------------------------------------------
        implicit none
        logical, intent(in ) :: lflag      !< logical flag
        integer, intent(in ) :: old_comm   !< old mpi communicator
        integer, intent(in ) :: old_rank   !< old mpi rank
        integer, intent(out) :: new_comm   !< new mpi communicator
        integer, intent(out) :: new_rank   !< new mpi rank
        integer, intent(out) :: new_nprocs !< new communicator n of procs
        integer, intent(out) :: err        !< error

        ! local declarations
        integer :: color = MPI_UNDEFINED   !< by default the color is undefined
        
        ! determine the color in respect of the flag
        if(lflag) color = 1
        call mpi_comm_split(old_comm, color, old_rank, new_comm, err)
       
        ! determine the rank and the procs with color = 1
        if(lflag) then
          call mpi_comm_rank(new_comm, new_rank  , err)
          call mpi_comm_size(new_comm, new_nprocs, err)
          !print*, old_rank, new_rank, new_nprocs
        endif

        return
end subroutine mpi_flag_splitting





subroutine mpi_plot_scl(lcl_x,lcl_f,n,my_str,my_end,my_id,comm,nprocs,root,outfile,it)

        use parameters_module, only: rp
        use FileModule  
        use mpi

        implicit none
        real(rp), dimension(:), allocatable, intent(in)    :: lcl_x
        real(rp), dimension(:), allocatable, intent(in)    :: lcl_f
        integer                            , intent(in)    :: n
        integer                            , intent(in)    :: my_str, my_end, my_id
        integer                            , intent(in)    :: comm, nprocs, root,it
        type(FileType)                     , intent(inout) :: outfile

        ! local declarations
        real(rp), allocatable, dimension(:) :: gbl_f
        real(rp), allocatable, dimension(:) :: gbl_x
        real(rp), allocatable, dimension(:) :: buff
        real(rp), allocatable, dimension(:) :: bufx
        integer                             :: i, err = 0
        integer                             :: tag = 400, pid, sr, er
        
        if(my_id == root) then

          allocate(gbl_f(1:n), gbl_x(1:n), stat = err)
          if(err .ne. 0) stop ' Allocation error in mpi_plot'
          
          gbl_f(my_str:my_end) = lcl_f(my_str:my_end)
          gbl_x(my_str:my_end) = lcl_x(my_str:my_end)

          do pid = 1, nprocs-1

             ! recieve start and end indices
             call MPI_Recv(sr, 1, mpi_int, pid, tag  , comm, mpi_status_ignore, err)
             call MPI_Recv(er, 1, mpi_int, pid, tag+1, comm, mpi_status_ignore, err)

             ! recieve the piece of data
             allocate(buff(sr:er), bufx(sr:er))

             call MPI_Recv(buff, size(buff,1), MPI_RP, pid, tag+2, comm, mpi_status_ignore, err)
             call MPI_Recv(bufx, size(bufx,1), MPI_RP, pid, tag+3, comm, mpi_status_ignore, err)

             gbl_f(sr:er) = buff
             gbl_x(sr:er) = bufx

             deallocate(buff,bufx)

          enddo
                
          ! write to the file
          call OpenNewFile(outFile,it)
          do i = 1,n
             write(outfile%unit,*) gbl_x(i), gbl_f(i)
          enddo
          call CloseFile(outFile)

          deallocate(gbl_x,gbl_f)

        else
        
          ! send start and end indices
          call MPI_Send(my_str, 1, mpi_int, root, tag  , comm, err)
          call MPI_Send(my_end, 1, mpi_int, root, tag+1, comm, err)
          
          allocate(buff(my_str:my_end), bufx(my_str:my_end))

          bufx = lcl_x(my_str:my_end)
          buff = lcl_f(my_str:my_end)

          call MPI_Send(buff, size(buff,1), MPI_RP, root, tag+2, comm, err)
          call MPI_Send(bufx, size(bufx,1), MPI_RP, root, tag+3, comm, err)

          deallocate(buff,bufx)

        endif

        return
end subroutine mpi_plot_scl








subroutine mpi_plot_vec(lcl_x,lcl_f,n,my_str,my_end,my_id,comm,nprocs,root,outfile,it)

        use parameters_module, only: rp
        use FileModule  
        use mpi

        implicit none
        real(rp), dimension(:)  , allocatable, intent(in)    :: lcl_x
        real(rp), dimension(:,:), allocatable, intent(in)    :: lcl_f
        integer                              , intent(in)    :: n
        integer                              , intent(in)    :: my_str, my_end, my_id
        integer                              , intent(in)    :: comm, nprocs, root,it
        type(FileType)                       , intent(inout) :: outFile

        ! local declarations
        real(rp), allocatable, dimension(:,:) :: gbl_f
        real(rp), allocatable, dimension(:)   :: gbl_x
        real(rp), allocatable, dimension(:,:) :: buff
        real(rp), allocatable, dimension(:)   :: bufx
        integer                               :: i, err = 0
        integer                               :: tag = 400, pid, sr, er, lb,ub
        
        lb = lbound(lcl_f,2)
        ub = ubound(lcl_f,2)

        if(my_id == root) then

          allocate(gbl_f(1:n,lb:ub), gbl_x(1:n), stat = err)
          if(err .ne. 0) stop ' Allocation error in mpi_plot'
          
          gbl_f(my_str:my_end,lb:ub) = lcl_f(my_str:my_end,lb:ub)
          gbl_x(my_str:my_end)       = lcl_x(my_str:my_end)

          do pid = 1, nprocs-1

             ! recieve start and end indices
             call MPI_Recv(sr, 1, mpi_int, pid, tag  , comm, mpi_status_ignore, err)
             call MPI_Recv(er, 1, mpi_int, pid, tag+1, comm, mpi_status_ignore, err)

             ! recieve the piece of data
             allocate(buff(sr:er,lb:ub), bufx(sr:er))

             call MPI_Recv(buff, size(buff,1)*size(buff,2), MPI_RP, pid, tag+2, comm, mpi_status_ignore, err)
             call MPI_Recv(bufx, size(bufx,1)             , MPI_RP, pid, tag+3, comm, mpi_status_ignore, err)

             gbl_f(sr:er,lb:ub) = buff
             gbl_x(sr:er)       = bufx

             deallocate(buff,bufx)

          enddo
                
          ! write to the file
          call OpenNewFile(outFile,it)
          do i = 1,n
             write(outfile%unit,'(100e18.9)') gbl_x(i), gbl_f(i,lb:ub)
          enddo
          call CloseFile(outFile)

          deallocate(gbl_x,gbl_f)

        else
        
          ! send start and end indices
          call MPI_Send(my_str, 1, mpi_int, root, tag  , comm, err)
          call MPI_Send(my_end, 1, mpi_int, root, tag+1, comm, err)
          
          allocate(buff(my_str:my_end, lb:ub), bufx(my_str:my_end))

          bufx = lcl_x(my_str:my_end)
          buff = lcl_f(my_str:my_end,lb:ub)

          call MPI_Send(buff, size(buff,1)*size(buff,2), MPI_RP, root, tag+2, comm, err)
          call MPI_Send(bufx, size(bufx,1)             , MPI_RP, root, tag+3, comm, err)

          deallocate(buff,bufx)

        endif

        return
end subroutine mpi_plot_vec



























subroutine end_mpi()
! ------------------------------------------------------
!       MPI finalization
! ------------------------------------------------------
        implicit none
        integer :: i
       
        !4D
        call DeallocateReal(phi_bfr_send_E)
        call DeallocateReal(phi_bfr_send_W)
        call DeallocateReal(phi_bfr_recv_E)
        call DeallocateReal(phi_bfr_recv_W)
        call DeallocateReal(phi_bfr_send_N)
        call DeallocateReal(phi_bfr_send_S)
        call DeallocateReal(phi_bfr_recv_N)
        call DeallocateReal(phi_bfr_recv_S)
        call DeallocateReal(phi_bfr_send_B)
        call DeallocateReal(phi_bfr_send_F)
        call DeallocateReal(phi_bfr_recv_B)
        call DeallocateReal(phi_bfr_recv_F)
        
        !3D
        call DeallocateReal(bfr_send_E)
        call DeallocateReal(bfr_send_W)
        call DeallocateReal(bfr_recv_E)
        call DeallocateReal(bfr_recv_W)
        call DeallocateReal(bfr_send_N)
        call DeallocateReal(bfr_send_S)
        call DeallocateReal(bfr_recv_N)
        call DeallocateReal(bfr_recv_S)
        call DeallocateReal(bfr_send_B)
        call DeallocateReal(bfr_send_F)
        call DeallocateReal(bfr_recv_B)
        call DeallocateReal(bfr_recv_F)

        !3D int1
        call DeallocateInteger(i13D_bfr_send_E)
        call DeallocateInteger(i13D_bfr_send_W)
        call DeallocateInteger(i13D_bfr_recv_E)
        call DeallocateInteger(i13D_bfr_recv_W)
        call DeallocateInteger(i13D_bfr_send_N)
        call DeallocateInteger(i13D_bfr_send_S)
        call DeallocateInteger(i13D_bfr_recv_N)
        call DeallocateInteger(i13D_bfr_recv_S)
        call DeallocateInteger(i13D_bfr_send_B)
        call DeallocateInteger(i13D_bfr_send_F)
        call DeallocateInteger(i13D_bfr_recv_B)
        call DeallocateInteger(i13D_bfr_recv_F)


        do i = 1, nb_neigh
           call MPI_Type_free(type_send_cons(i),mpi_err)
           call MPI_Type_free(type_recv_cons(i),mpi_err)
        enddo

        if(hybrid_weno) then
          do i = 1, nb_neigh
             call MPI_Type_free(type_send_flag(i),mpi_err)
             call MPI_Type_free(type_recv_flag(i),mpi_err)
          enddo
        endif

        ! free the new communicator
        call MPI_COMM_FREE(mpi_comm_cart,mpi_err)
        
        ! close MPI..
        call MPI_Finalize(mpi_err)
        
        return
end subroutine end_mpi


subroutine secure_stop
! ------------------------------------------------------
!       Stop the program safetly
! ------------------------------------------------------
        implicit none
        
        if(rank == root) then
                print*, ' The program will be stopped.'
                write(*,*)
        endif

        call end_mpi

        stop

        return
end subroutine secure_stop

subroutine check_mpi(sub_prog_name,check)
! ------------------------------------------------------
!       Check for MPI errors
! ------------------------------------------------------
        integer     , intent(in) :: check
        character(*), intent(in) :: sub_prog_name

        if(check .ne. 0) then
          if(rank == root) write(*,*) ' MPI problem in "', sub_prog_name, '"'
          call secure_stop
        end if
        return
end subroutine check_mpi



end module mpi_module
