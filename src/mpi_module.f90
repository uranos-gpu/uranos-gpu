module mpi_module
use parameters_module
use mpi
#ifdef _OPENACC
use openacc
#endif
use allocate_module
implicit none

integer            :: rank
integer            :: nprocs
integer            :: mpi_err  = 0
integer, parameter :: root     = 0
logical            :: mpi_flag = .false.

!$ integer            :: required, provided

integer                   :: mpi_comm_cart
integer                   :: mpi_group_cart
integer                   :: mpi_group_E, mpi_group_W, mpi_group_S, mpi_group_N, mpi_group_B, mpi_group_F
integer, parameter        :: ndims = 3
integer, parameter        :: ndims_cons = 4

integer :: comm1Dx , comm1Dy , comm1Dz
integer :: comm2Dxy, comm2Dxz, comm2Dyz

integer, dimension(ndims) :: cart_coord
logical, dimension(ndims) :: periods

integer, parameter           :: nb_neigh = 6
integer, parameter, private  :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
integer, dimension(nb_neigh) :: my_neighbour

integer, parameter :: GN = 4

integer :: sx, ex
integer :: sy, ey
integer :: sz, ez

integer :: lbx, ubx
integer :: lby, uby
integer :: lbz, ubz

integer, dimension(nb_neigh) :: type_send_cons, type_recv_cons
integer, dimension(nb_neigh) :: type_send_flag, type_recv_flag
integer, dimension(nb_neigh) :: type_send_prim, type_recv_prim

real(rp), allocatable, dimension(:,:,:,:) :: phi_bfr_send_E, phi_bfr_send_W, phi_bfr_recv_E, phi_bfr_recv_W
real(rp), allocatable, dimension(:,:,:,:) :: phi_bfr_send_N, phi_bfr_send_S, phi_bfr_recv_N, phi_bfr_recv_S
real(rp), allocatable, dimension(:,:,:,:) :: phi_bfr_send_B, phi_bfr_send_F, phi_bfr_recv_B, phi_bfr_recv_F

integer :: win_phi_bfr_E, win_phi_bfr_W
integer :: win_phi_bfr_N, win_phi_bfr_S
integer :: win_phi_bfr_B, win_phi_bfr_F

real(rp), allocatable, dimension(:,:,:) :: bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W
real(rp), allocatable, dimension(:,:,:) :: bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S
real(rp), allocatable, dimension(:,:,:) :: bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F

integer :: win_bfr_E, win_bfr_W
integer :: win_bfr_N, win_bfr_S
integer :: win_bfr_B, win_bfr_F

integer(1), allocatable, dimension(:,:,:) :: i13D_bfr_send_E, i13D_bfr_send_W, i13D_bfr_recv_E, i13D_bfr_recv_W
integer(1), allocatable, dimension(:,:,:) :: i13D_bfr_send_N, i13D_bfr_send_S, i13D_bfr_recv_N, i13D_bfr_recv_S
integer(1), allocatable, dimension(:,:,:) :: i13D_bfr_send_B, i13D_bfr_send_F, i13D_bfr_recv_B, i13D_bfr_recv_F

integer :: win_i13D_bfr_E, win_i13D_bfr_W
integer :: win_i13D_bfr_N, win_i13D_bfr_S
integer :: win_i13D_bfr_B, win_i13D_bfr_F

integer, allocatable, dimension(:) :: req_array_xx, req_array_yz

integer :: num_dev, dev_id

contains
subroutine init_mpi
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
        print*, 'Rank ', local_rank, 'is associated to GPU', dev_id
#endif

        call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)

        return
end subroutine init_mpi

subroutine init_cartesian_topology
        implicit none
        integer :: i, comm_tiled
        logical, parameter :: reorganisation = .false.

        do i = 1,5,2
           if(trim(bc(i)) == 'periodic') then
             if(trim(bc(i+1)) .ne. trim(bc(i))) then
               if(rank == root) print*, ' Periodic boundary condition must be in pair.'
               call secure_stop
             endif
           endif
        enddo

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

        if(rank == root) then
                write(*,*) '--------------------------------------------------- '
                write(*,'(1x,A,i4,A)') ' Executing URANOS with', nprocs, ' procs.'
        endif

        if(dims == 2) then
                if(cart_dims(3) == 0) cart_dims(3) = 1
                if(cart_dims(3) /= 1) then
                   if(rank == root) write(*,*) ' ERROR: 2-dimensional problems accept just one proc along-z.'
                   call MPI_FINALIZE(mpi_err)
                   stop
                endif
        endif

        if(all(cart_dims /= 0)) then

                if(cart_dims(1)*cart_dims(2)*cart_dims(3) /= nprocs) then
                   if(rank == root) write(*,*) ' ERROR: The number of processors specified is not compatible.'
                   call MPI_FINALIZE(mpi_err)
                   stop
                endif

                if(rank == root) print*, ' The processors has distributed by USERs specification as follow:'

        else

                call min_halo_surface_dims(nprocs, nx, ny, nz, cart_dims)

                if(rank == root) print*, ' The processors are distributed to minimise the halo surface as follow:'

        endif

        mpi_comm_cart = MPI_COMM_NULL
        call node_tile_comm(cart_dims, comm_tiled)

        CALL MPI_cart_create(comm_tiled, ndims, cart_dims, periods, &
                reorganisation, mpi_comm_cart, mpi_err)

        call MPI_Comm_rank(mpi_comm_cart, rank, mpi_err)

        call check_mpi('init_cartesian_topology',mpi_err)

        return
end subroutine init_cartesian_topology

subroutine min_halo_surface_dims(np, nxg, nyg, nzg, dims)
        implicit none
        integer, intent(in)    :: np, nxg, nyg, nzg
        integer, intent(inout) :: dims(3)
        integer                :: px, py, pz, best(3)
        real(8)                :: cost, best_cost

        best      = 0
        best_cost = huge(1.0_8)

        do px = 1, np
           if(mod(np,px) /= 0) cycle
           if(dims(1) /= 0 .and. px /= dims(1)) cycle
           do py = 1, np/px
              if(mod(np/px,py) /= 0) cycle
              if(dims(2) /= 0 .and. py /= dims(2)) cycle
              pz = np/(px*py)
              if(dims(3) /= 0 .and. pz /= dims(3)) cycle

              if(px > 1 .and. nxg/px < GN) cycle
              if(py > 1 .and. nyg/py < GN) cycle
              if(pz > 1 .and. nzg/pz < GN) cycle

              cost = real(nyg,8)*real(nzg,8)/real(py*pz,8) &
                   + real(nxg,8)*real(nzg,8)/real(px*pz,8) &
                   + real(nxg,8)*real(nyg,8)/real(px*py,8)

              if(cost < best_cost) then
                 best_cost = cost
                 best      = (/px, py, pz/)
              endif
           enddo
        enddo

        if(best(1) == 0) then
           if(rank == root) write(*,*) ' ERROR: no processor grid is compatible with cart_dims and the stencil width.'
           call MPI_FINALIZE(mpi_err)
           stop
        endif

        dims = best

        return
end subroutine min_halo_surface_dims

subroutine node_tile_comm(pdims, newcomm)
        implicit none
        integer, intent(in)  :: pdims(3)
        integer, intent(out) :: newcomm
        integer :: lcomm, nloc, lrank, node, loc, crank
        integer :: gx, gy, gz, bx, by, bz, nnx, nny, nnz
        integer :: I_, J_, K_, ii, jj, kk, ci, cj, ck
        integer :: tile(3)
        real(8) :: cost, best_cost
        logical :: ok, ok_all

        newcomm = MPI_COMM_WORLD
        tile    = 1

        call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, lcomm, mpi_err)
        call MPI_Comm_size(lcomm, nloc , mpi_err)
        call MPI_Comm_rank(lcomm, lrank, mpi_err)
        call MPI_Comm_free(lcomm, mpi_err)

        ok = (mod(nprocs,nloc) == 0) .and. (lrank == mod(rank,nloc))
        call MPI_Allreduce(ok, ok_all, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, mpi_err)
        if(.not. ok_all .or. nloc == 1) return

        best_cost = huge(1.0_8)
        do gx = 1, nloc
           if(mod(nloc,gx) /= 0 .or. mod(pdims(1),gx) /= 0) cycle
           do gy = 1, nloc/gx
              if(mod(nloc/gx,gy) /= 0 .or. mod(pdims(2),gy) /= 0) cycle
              gz = nloc/(gx*gy)
              if(mod(pdims(3),gz) /= 0) cycle
              bx = (nx/pdims(1))*gx
              by = (ny/pdims(2))*gy
              bz = (nz/pdims(3))*gz
              cost = 0.0_8
              if(pdims(1)/gx > 1) cost = cost + real(by,8)*real(bz,8)
              if(pdims(2)/gy > 1) cost = cost + real(bx,8)*real(bz,8)
              if(pdims(3)/gz > 1) cost = cost + real(bx,8)*real(by,8)
              if(cost < best_cost) then
                 best_cost = cost
                 tile      = (/gx, gy, gz/)
              endif
           enddo
        enddo
        if(best_cost == huge(1.0_8)) return

        gx = tile(1); gy = tile(2); gz = tile(3)
        nnx = pdims(1)/gx; nny = pdims(2)/gy; nnz = pdims(3)/gz
        node = rank/nloc
        loc  = mod(rank,nloc)
        I_ = node/(nny*nnz); J_ = mod(node/nnz, nny); K_ = mod(node, nnz)
        ii = loc/(gy*gz)   ; jj = mod(loc/gz, gy)   ; kk = mod(loc, gz)
        ci = I_*gx + ii; cj = J_*gy + jj; ck = K_*gz + kk
        crank = (ci*pdims(2) + cj)*pdims(3) + ck

        call MPI_Comm_split(MPI_COMM_WORLD, 0, crank, newcomm, mpi_err)

        if(rank == root) write(*,'(1x,A,i3,A,3i3)') ' Ranks per node: ', nloc, '   node tile (x,y,z): ', tile

        return
end subroutine node_tile_comm

subroutine init_subcartesian_topology
        implicit none
        logical, dimension(3) :: TFF, FTF, FFT, TTF, TFT, FTT
        logical, parameter    :: T = .true., F = .false.
        integer               :: err

        TFF(1) = T
        TFF(2) = F
        TFF(3) = F

        FTF(1) = F
        FTF(2) = T
        FTF(3) = F

        FFT(1) = F
        FFT(2) = F
        FFT(3) = T

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

        comm2Dxy = MPI_COMM_NULL
        comm2Dxz = MPI_COMM_NULL
        comm2Dyz = MPI_COMM_NULL

        call MPI_CART_SUB(mpi_comm_cart, TTF, comm2Dxy, err)
        call MPI_CART_SUB(mpi_comm_cart, TFT, comm2Dxz, err)
        call MPI_CART_SUB(mpi_comm_cart, FTT, comm2Dyz, err)

        return
end subroutine init_subcartesian_topology

subroutine init_subdomain_boundaries
        implicit none

        call MPI_cart_coords(mpi_comm_cart,rank,ndims,cart_coord,mpi_err)

        sx = (cart_coord(1)*nx)/cart_dims(1)+1
        ex = ((cart_coord(1)+1)*nx)/cart_dims(1)

        sy = (cart_coord(2)*ny)/cart_dims(2)+1
        ey = ((cart_coord(2)+1)*ny)/cart_dims(2)

        if(dims == 2) nz = 1

        sz = (cart_coord(3)*nz)/cart_dims(3)+1
        ez = ((cart_coord(3)+1)*nz)/cart_dims(3)

        lbx = sx - GN
        ubx = ex + GN

        lby = sy - GN
        uby = ey + GN

        lbz = sz - GN
        ubz = ez + GN

        if(mpi_opt_level == 4) then
                if (.not.(cuda_aware)) then
                        if(rank == root) then
                                print*, ' ERROR: ONE-SIDED communication requires CUDA-AWARE'
                        endif
                endif
                if(mod(ex-sx+1,cart_dims(1)) /= 0) then
                        if(rank == root) then
                                print*, ' ERROR: the number of points is not divisible by the number of procs in x-direction.'
                        endif
                        stop
                endif
                if(mod(ey-sy+1,cart_dims(2)) /= 0) then
                        if(rank == root) then
                                print*, ' ERROR: the number of points is not divisible by the number of procs in y-direction.'
                        endif
                        stop
                endif
                if(mod(ez-sz+1,cart_dims(3)) /= 0) then
                        if(rank == root) then
                                print*, ' ERROR: the number of points is not divisible by the number of procs in z-direction.'
                        endif
                        stop
                endif
        endif

        if(rank == root) then
                write(*,*)
                write(*,'(1x,A,i4)')   ' Problem dimensions  : ', dims
                write(*,*)
                if (thermally_perfect_gas == 1) then
                write(*,'(1x,A,i4)') ' The gas is modeled as THERMALLY PERFECT.'
                else
                write(*,'(1x,A,i4)') ' The gas is modeled as CALORICALLY PERFECT.'
                endif
                write(*,'(1x,A,i4)')   ' Procs in x-direction: ', cart_dims(1)
                write(*,'(1x,A,i4)')   ' Procs in y-direction: ', cart_dims(2)
                write(*,'(1x,A,i4)')   ' Procs in z-direction: ', cart_dims(3)
                write(*,*)
                write(*,'(1x,A,i4)')   ' Point in x-direction: ', nx
                write(*,'(1x,A,i4)')   ' Point in y-direction: ', ny
                write(*,'(1x,A,i4)')   ' Point in z-direction: ', nz
                write(*,*)
                if(len_trim(mesh_file) > 0) then
                write(*,'(1x,A,A)')    ' Grid topology       :  CURVILINEAR, mesh file = ', trim(mesh_file)
                else
                write(*,'(1x,A)')      ' Grid topology       :  CARTESIAN (separable)'
                endif
                write(*,*)
                if(mpi_opt_level == 1) write(*,*) ' MPI communications handled with SENDRECV'
                if(mpi_opt_level == 2) write(*,*) ' MPI communications handled with ISEND/IRECV + BFR'
                if(mpi_opt_level == 3) write(*,*) ' MPI communications handled with ISEND/IRECV + DDT'
                if(mpi_opt_level == 4) write(*,*) ' MPI communications handled with ONE-SIDED'
                if(cuda_aware)         write(*,*) ' CUDA-AWARE MPI enabled'
                write(*,*)
                write(*,'(1x,A,A)')    ' Advection scheme    : ', trim(scheme)
                write(*,'(1x,A,A)')    ' ep_type             : ', trim(ep_type)
                write(*,'(1x,A,A)')    ' Diffusion scheme    : ', trim(diffusion_scheme)
                write(*,'(1x,A,i4)')   ' FD accurary order   : ', fd_order
                write(*,'(1x,A,i4)')   ' WENO accurary order : ', weno_order
                if(viscous) then
                write(*,'(1x,A,e18.6)') ' Mach number         : ', Mach
                write(*,'(1x,A,e18.6)') ' Reynolds number     : ', Reynolds
                endif
                if(logical_CFL) then
                write(*,'(1x,A,e18.6)') ' CFL number          : ', CFL
                if(viscous) &
                write(*,'(1x,A,e18.6)') ' Fourier number      : ', Fourier
                endif
                if(les) then
                write(*,*)
                write(*,'(1x,A,A)')     ' LES model            : ', trim(sgs_model)
                if(Yoshizawa) then
                write(*,'(1x,A)')       ' Yoshizawa            :  Active'
                endif
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
        implicit none
        integer, parameter, dimension(3) :: dirs = (/0,1,2/)
        integer, parameter               :: step = 1

        my_neighbour(:) = MPI_PROC_NULL

        call MPI_cart_shift(mpi_comm_cart, dirs(1), step, my_neighbour(W), my_neighbour(E), mpi_err)

        call MPI_cart_shift(mpi_comm_cart, dirs(2), step, my_neighbour(S), my_neighbour(N), mpi_err)

        call MPI_cart_shift(mpi_comm_cart, dirs(3), step, my_neighbour(B), my_neighbour(F), mpi_err)

        return
end subroutine init_mpi_neighbours

subroutine init_mpi_exchange_data
        implicit none

        call init_mpi_buffers

        call init_mpi_type_field(ndims_cons,5,MPI_RP, type_send_cons, type_recv_cons)

        if    (dims == 2) then
          allocate(req_array_xx(4), req_array_yz(4))
        elseif(dims == 3) then
          allocate(req_array_xx(4), req_array_yz(8))
        endif

        if(hybrid_weno) then
          call init_mpi_type_field(ndims,1,MPI_INTEGER1, type_send_flag, type_recv_flag)
        endif

        call init_mpi_type_field(ndims,1,MPI_RP,type_send_prim, type_recv_prim)

        return
end subroutine init_mpi_exchange_data

subroutine init_mpi_type_field(field_dims,field_co_dims,old_field_type, type_field_send, type_field_recv)
        implicit none
        integer                     , intent(in)    :: field_dims
        integer                     , intent(in)    :: field_co_dims
        integer                     , intent(in)    :: old_field_type
        integer, dimension(nb_neigh), intent(inout) :: type_field_send, type_field_recv

        integer, dimension(field_dims) :: old_field_size
        integer, dimension(field_dims) :: size_sub_field_x
        integer, dimension(field_dims) :: size_sub_field_y
        integer, dimension(field_dims) :: size_sub_field_z

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

                start_field_send_W = (/GN        , 0, 0/)
                start_field_recv_E = (/ex-sx+1+GN, 0, 0/)
                start_field_send_E = (/ex-sx+1   , 0, 0/)
                start_field_recv_W = (/0         , 0, 0/)

                start_field_send_S = (/0, GN        , 0/)
                start_field_recv_N = (/0, ey-sy+1+GN, 0/)
                start_field_send_N = (/0, ey-sy+1   , 0/)
                start_field_recv_S = (/0, 0         , 0/)

                start_field_send_B = (/0, 0, GN        /)
                start_field_recv_F = (/0, 0, ez-sz+1+GN/)
                start_field_send_F = (/0, 0, ez-sz+1   /)
                start_field_recv_B = (/0, 0, 0         /)

        elseif(field_dims == 4) then

                old_field_size   = (/gn+ex-sx+1+gn, gn+ey-sy+1+gn, gn+ez-sz+1+gn, field_co_dims/)
                size_sub_field_x = (/GN           , gn+ey-sy+1+gn, gn+ez-sz+1+gn, field_co_dims/)
                size_sub_field_y = (/gn+ex-sx+1+gn, GN           , gn+ez-sz+1+gn, field_co_dims/)
                size_sub_field_z = (/gn+ex-sx+1+gn, gn+ey-sy+1+gn, GN           , field_co_dims/)

                start_field_send_W = (/GN        , 0, 0, 0/)
                start_field_recv_E = (/ex-sx+1+GN, 0, 0, 0/)
                start_field_send_E = (/ex-sx+1   , 0, 0, 0/)
                start_field_recv_W = (/0         , 0, 0, 0/)

                start_field_send_S = (/0, GN        , 0, 0/)
                start_field_recv_N = (/0, ey-sy+1+GN, 0, 0/)
                start_field_send_N = (/0, ey-sy+1   , 0, 0/)
                start_field_recv_S = (/0, 0         , 0, 0/)

                start_field_send_B = (/0, 0, GN        , 0/)
                start_field_recv_F = (/0, 0, ez-sz+1+GN, 0/)
                start_field_send_F = (/0, 0, ez-sz+1   , 0/)
                start_field_recv_B = (/0, 0, 0         , 0/)

        endif

        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_x,start_field_send_E,old_field_type,type_field_send(E))
        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_x,start_field_send_W,old_field_type,type_field_send(W))

        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_x,start_field_recv_E,old_field_type,type_field_recv(E))
        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_x,start_field_recv_W,old_field_type,type_field_recv(W))

        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_y,start_field_send_S,old_field_type,type_field_send(S))
        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_y,start_field_send_N,old_field_type,type_field_send(N))

        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_y,start_field_recv_S,old_field_type,type_field_recv(S))
        call MPI_create_subarray(field_dims,old_field_size,size_sub_field_y,start_field_recv_N,old_field_type,type_field_recv(N))

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

        call MPI_type_create_subarray(dims,old_size, subsize, start_position , &
                                      MPI_ORDER_FORTRAN, old_type, new_type  , &
                                      mpi_err)
        call MPI_type_commit(new_type,mpi_err)

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

subroutine init_mpi_windows
        implicit none
        integer(kind=MPI_ADDRESS_KIND) :: win_size
        integer, parameter             :: group_size = 1
        integer, dimension(1)          :: ranks

        call MPI_COMM_GROUP(mpi_comm_cart, mpi_group_cart, mpi_err)

        if (my_neighbour(E) >= 0) then
                ranks(1) = my_neighbour(E)
                call MPI_GROUP_INCL(mpi_group_cart, group_size, ranks, mpi_group_E, mpi_err)
        endif

        if (my_neighbour(W) >= 0) then
                ranks(1) = my_neighbour(W)
                call MPI_GROUP_INCL(mpi_group_cart, group_size, ranks, mpi_group_W, mpi_err)
        endif

        if (my_neighbour(N) >= 0) then
                ranks(1) = my_neighbour(N)
                call MPI_GROUP_INCL(mpi_group_cart, group_size, ranks, mpi_group_N, mpi_err)
        endif

        if (my_neighbour(S) >= 0) then
                ranks(1) = my_neighbour(S)
                call MPI_GROUP_INCL(mpi_group_cart, group_size, ranks, mpi_group_S, mpi_err)
        endif

        if (my_neighbour(B) >= 0) then
                ranks(1) = my_neighbour(B)
                call MPI_GROUP_INCL(mpi_group_cart, group_size, ranks, mpi_group_B, mpi_err)
        endif

        if (my_neighbour(F) >= 0) then
                ranks(1) = my_neighbour(F)
                call MPI_GROUP_INCL(mpi_group_cart, group_size, ranks, mpi_group_F, mpi_err)
        endif

        if (mpi_opt_level == 4) then
                !$acc host_data use_device(phi_bfr_recv_W) &
                !$acc           use_device(phi_bfr_recv_E) &
                !$acc           use_device(phi_bfr_recv_S) &
                !$acc           use_device(phi_bfr_recv_N) &
                !$acc           use_device(phi_bfr_recv_B) &
                !$acc           use_device(phi_bfr_recv_F)
                win_size = product(shape(phi_bfr_recv_W)) * rp
                call MPI_WIN_CREATE(phi_bfr_recv_E, win_size, 1, MPI_INFO_NULL, mpi_comm_cart, win_phi_bfr_E, mpi_err)
                call MPI_WIN_CREATE(phi_bfr_recv_W, win_size, 1, MPI_INFO_NULL, mpi_comm_cart, win_phi_bfr_W, mpi_err)

                win_size = product(shape(phi_bfr_recv_S)) * rp
                call MPI_WIN_CREATE(phi_bfr_recv_S, win_size, 1, MPI_INFO_NULL, mpi_comm_cart, win_phi_bfr_S, mpi_err)
                call MPI_WIN_CREATE(phi_bfr_recv_N, win_size, 1, MPI_INFO_NULL, mpi_comm_cart, win_phi_bfr_N, mpi_err)

                win_size = product(shape(phi_bfr_recv_B)) * rp
                call MPI_WIN_CREATE(phi_bfr_recv_B, win_size, 1, MPI_INFO_NULL, mpi_comm_cart, win_phi_bfr_B, mpi_err)
                call MPI_WIN_CREATE(phi_bfr_recv_F, win_size, 1, MPI_INFO_NULL, mpi_comm_cart, win_phi_bfr_F, mpi_err)

                call MPI_Win_fence(0, win_phi_bfr_E, mpi_err)
                call MPI_Win_fence(0, win_phi_bfr_W, mpi_err)
                call MPI_Win_fence(0, win_phi_bfr_S, mpi_err)
                call MPI_Win_fence(0, win_phi_bfr_N, mpi_err)
                call MPI_Win_fence(0, win_phi_bfr_B, mpi_err)
                call MPI_Win_fence(0, win_phi_bfr_F, mpi_err)
                !$acc end host_data

                !$acc host_data use_device(bfr_recv_W) &
                !$acc           use_device(bfr_recv_E) &
                !$acc           use_device(bfr_recv_S) &
                !$acc           use_device(bfr_recv_N) &
                !$acc           use_device(bfr_recv_B) &
                !$acc           use_device(bfr_recv_F)

                win_size = product(shape(bfr_recv_W)) * rp
                call MPI_WIN_CREATE(bfr_recv_E, win_size, rp, MPI_INFO_NULL, mpi_comm_cart, win_bfr_E, mpi_err)
                call MPI_WIN_CREATE(bfr_recv_W, win_size, rp, MPI_INFO_NULL, mpi_comm_cart, win_bfr_W, mpi_err)

                win_size = product(shape(bfr_recv_S)) * rp
                call MPI_WIN_CREATE(bfr_recv_S, win_size, rp, MPI_INFO_NULL, mpi_comm_cart, win_bfr_S, mpi_err)
                call MPI_WIN_CREATE(bfr_recv_N, win_size, rp, MPI_INFO_NULL, mpi_comm_cart, win_bfr_N, mpi_err)

                win_size = product(shape(bfr_recv_B)) * rp
                call MPI_WIN_CREATE(bfr_recv_B, win_size, rp, MPI_INFO_NULL, mpi_comm_cart, win_bfr_B, mpi_err)
                call MPI_WIN_CREATE(bfr_recv_F, win_size, rp, MPI_INFO_NULL, mpi_comm_cart, win_bfr_F, mpi_err)

                call MPI_Win_fence(0, win_bfr_E, mpi_err)
                call MPI_Win_fence(0, win_bfr_W, mpi_err)
                call MPI_Win_fence(0, win_bfr_S, mpi_err)
                call MPI_Win_fence(0, win_bfr_N, mpi_err)
                call MPI_Win_fence(0, win_bfr_B, mpi_err)
                call MPI_Win_fence(0, win_bfr_F, mpi_err)
                !$acc end host_data

                !$acc host_data use_device(i13D_bfr_recv_W) &
                !$acc           use_device(i13D_bfr_recv_E) &
                !$acc           use_device(i13D_bfr_recv_S) &
                !$acc           use_device(i13D_bfr_recv_N) &
                !$acc           use_device(i13D_bfr_recv_B) &
                !$acc           use_device(i13D_bfr_recv_F)

                win_size = product(shape(i13D_bfr_recv_W)) * 1
                call MPI_WIN_CREATE(i13D_bfr_recv_E, win_size, 1, MPI_INFO_NULL, mpi_comm_cart, win_i13D_bfr_E, mpi_err)
                call MPI_WIN_CREATE(i13D_bfr_recv_W, win_size, 1, MPI_INFO_NULL, mpi_comm_cart, win_i13D_bfr_W, mpi_err)

                win_size = product(shape(i13D_bfr_recv_S)) * 1
                call MPI_WIN_CREATE(i13D_bfr_recv_S, win_size, 1, MPI_INFO_NULL, mpi_comm_cart, win_i13D_bfr_S, mpi_err)
                call MPI_WIN_CREATE(i13D_bfr_recv_N, win_size, 1, MPI_INFO_NULL, mpi_comm_cart, win_i13D_bfr_N, mpi_err)

                win_size = product(shape(i13D_bfr_recv_B)) * 1
                call MPI_WIN_CREATE(i13D_bfr_recv_B, win_size, 1, MPI_INFO_NULL, mpi_comm_cart, win_i13D_bfr_B, mpi_err)
                call MPI_WIN_CREATE(i13D_bfr_recv_F, win_size, 1, MPI_INFO_NULL, mpi_comm_cart, win_i13D_bfr_F, mpi_err)

                call MPI_Win_fence(0, win_i13D_bfr_E, mpi_err)
                call MPI_Win_fence(0, win_i13D_bfr_W, mpi_err)
                call MPI_Win_fence(0, win_i13D_bfr_S, mpi_err)
                call MPI_Win_fence(0, win_i13D_bfr_N, mpi_err)
                call MPI_Win_fence(0, win_i13D_bfr_B, mpi_err)
                call MPI_Win_fence(0, win_i13D_bfr_F, mpi_err)
                !$acc end host_data

        endif

        return
end subroutine init_mpi_windows

subroutine mpi_safe_allreduce(local,glbal,mpi_operation,mpi_comm,err)
        implicit none
        real(rp), intent(in)  :: local
        integer , intent(in)  :: mpi_operation
        integer , intent(in)  :: mpi_comm
        real(rp), intent(out) :: glbal
        integer , intent(out) :: err

        real(rp), dimension(nprocs) :: glb_a

        call mpi_gather(local,1,MPI_RP, &
                        glb_a,1,MPI_RP, &
                        root, mpi_comm,err)

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

        call mpi_bcast(glbal,1,MPI_RP,root,mpi_comm,err)

        return
end subroutine mpi_safe_allreduce

subroutine mpi_flag_splitting(lflag,old_comm,old_rank,new_comm,new_rank,new_nprocs,err)
        implicit none
        logical, intent(in ) :: lflag
        integer, intent(in ) :: old_comm
        integer, intent(in ) :: old_rank
        integer, intent(out) :: new_comm
        integer, intent(out) :: new_rank
        integer, intent(out) :: new_nprocs
        integer, intent(out) :: err

        integer :: color = MPI_UNDEFINED

        if(lflag) color = 1
        call mpi_comm_split(old_comm, color, old_rank, new_comm, err)

        if(lflag) then
          call mpi_comm_rank(new_comm, new_rank  , err)
          call mpi_comm_size(new_comm, new_nprocs, err)
        endif

        return
end subroutine mpi_flag_splitting

subroutine end_mpi()
        implicit none
        integer :: i

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

        call MPI_COMM_FREE(mpi_comm_cart,mpi_err)

        call MPI_Finalize(mpi_err)

        return
end subroutine end_mpi

subroutine secure_stop
        implicit none

        if(rank == root) then
                print*, ' The program will be stopped.'
                write(*,*)
        endif

        flush(6)

        call MPI_Abort(MPI_COMM_WORLD, 1, mpi_err)

        stop

        return
end subroutine secure_stop

subroutine check_mpi(sub_prog_name,check)
        integer     , intent(in) :: check
        character(*), intent(in) :: sub_prog_name

        if(check .ne. 0) then
          if(rank == root) write(*,*) ' MPI problem in "', sub_prog_name, '"'
          call secure_stop
        end if
        return
end subroutine check_mpi

subroutine check_mpi_vec(routine, errs, labels)
        character(*), intent(in)               :: routine
        integer, dimension(:), intent(in)      :: errs
        character(*), dimension(:), intent(in) :: labels

        integer :: i, first_fail, local_max, global_max, alr_err

        first_fail = 0
        local_max  = 0
        do i = 1, size(errs)
          if(errs(i) /= 0) then
            if(first_fail == 0) first_fail = i
            local_max = max(local_max, abs(errs(i)))
          endif
        enddo

        call MPI_Allreduce(local_max, global_max, 1, MPI_INTEGER, MPI_MAX, mpi_comm_cart, alr_err)

        if(global_max /= 0) then
          if(first_fail > 0) then
            write(*,'(1x,A,I0,A,A,A,A,A,I0)') '[rank ', rank, '] MPI error in "', &
              trim(routine), '" at step "', trim(labels(first_fail)), '", code=', errs(first_fail)
          endif
          call secure_stop
        endif

        return
end subroutine check_mpi_vec

end module mpi_module
