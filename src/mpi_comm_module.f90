module mpi_comm_module
use parameters_module, only: rp, MPI_RP, mpi_opt_level, cuda_aware
use mpi_module, only: sx,ex,sy,ey,sz,ez,lbx,ubx,lby,uby,lbz,ubz,GN
use mpi

implicit none
private
public mpi_share, mpi_share_int1


interface mpi_share
  module procedure mpi_share4D, mpi_share3D
endinterface


interface mpi_share_int1
  module procedure mpi_shareInt3D
endinterface


contains

subroutine mpi_share4D(mpi_comm_cart,type_send_cons,type_recv_cons,my_neighbour,dims, &
                phi_bfr_send_E, phi_bfr_send_W, phi_bfr_recv_E, phi_bfr_recv_W, &
                phi_bfr_send_N, phi_bfr_send_S, phi_bfr_recv_N, phi_bfr_recv_S, &
                phi_bfr_send_B, phi_bfr_send_F, phi_bfr_recv_B, phi_bfr_recv_F, &
                phi)

        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi_bfr_send_E
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi_bfr_send_W
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi_bfr_recv_E
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi_bfr_recv_W
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi_bfr_send_N
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi_bfr_send_S
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi_bfr_recv_N
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi_bfr_recv_S
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi_bfr_send_B
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi_bfr_send_F
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi_bfr_recv_B
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: phi_bfr_recv_F
        integer              , dimension(6)      , intent(in)    :: type_send_cons
        integer              , dimension(6)      , intent(in)    :: type_recv_cons
        integer              , dimension(6)      , intent(in)    :: my_neighbour
        integer                                  , intent(in)    :: dims
        integer                                  , intent(in)    :: mpi_comm_cart
        
        selectcase(mpi_opt_level)
        case(1) !SENDRECV
        
          if(cuda_aware) then
          call mpi_communicate4D_real8_sendrecv_cuda_aware(mpi_comm_cart,&
                my_neighbour,dims,&
                phi_bfr_send_E, phi_bfr_send_W, phi_bfr_recv_E, phi_bfr_recv_W, &
                phi_bfr_send_N, phi_bfr_send_S, phi_bfr_recv_N, phi_bfr_recv_S, &
                phi_bfr_send_B, phi_bfr_send_F, phi_bfr_recv_B, phi_bfr_recv_F, &
                phi)

          else
          call mpi_communicate4D_real8_sendrecv(mpi_comm_cart,my_neighbour,dims,&
                phi_bfr_send_E, phi_bfr_send_W, phi_bfr_recv_E, phi_bfr_recv_W, &
                phi_bfr_send_N, phi_bfr_send_S, phi_bfr_recv_N, phi_bfr_recv_S, &
                phi_bfr_send_B, phi_bfr_send_F, phi_bfr_recv_B, phi_bfr_recv_F, &
                phi)
          endif
        
        case(2) !ISEND/IRECV + BUFFERS
        
          if(cuda_aware) then
          call mpi_communicate4D_real8_bfr_cuda_aware(mpi_comm_cart,&
                my_neighbour,dims    , &
                phi_bfr_send_E, phi_bfr_send_W, phi_bfr_recv_E, phi_bfr_recv_W, &
                phi_bfr_send_N, phi_bfr_send_S, phi_bfr_recv_N, phi_bfr_recv_S, &
                phi_bfr_send_B, phi_bfr_send_F, phi_bfr_recv_B, phi_bfr_recv_F, &
                phi)
          else
          call mpi_communicate4D_real8_bfr(mpi_comm_cart,my_neighbour,dims    , &
                phi_bfr_send_E, phi_bfr_send_W, phi_bfr_recv_E, phi_bfr_recv_W, &
                phi_bfr_send_N, phi_bfr_send_S, phi_bfr_recv_N, phi_bfr_recv_S, &
                phi_bfr_send_B, phi_bfr_send_F, phi_bfr_recv_B, phi_bfr_recv_F, &
                phi)
          endif

        case(3)
          ! ISEND/IRECV + DERIVED DATA TYPES
        
          if(cuda_aware) then
          call mpi_communicate4D_real8_cuda_aware(mpi_comm_cart,&
                     type_send_cons,type_recv_cons,my_neighbour,dims,phi)
          else
          call mpi_communicate4D_real8(mpi_comm_cart,&
                     type_send_cons,type_recv_cons,my_neighbour,dims,phi)
          endif
        case default

          print*, 'mpi_opt_level equal to ', mpi_opt_level, ' is not implemented'
          stop

        endselect

        return
end subroutine mpi_share4D



subroutine mpi_share3D(mpi_comm_cart,type_send_prim,type_recv_prim,my_neighbour,dims, &
                bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                var)

        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: var
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: bfr_send_E
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: bfr_send_W
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: bfr_recv_E
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: bfr_recv_W
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: bfr_send_N
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: bfr_send_S
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: bfr_recv_N
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: bfr_recv_S
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: bfr_send_B
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: bfr_send_F
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: bfr_recv_B
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: bfr_recv_F
        integer              , dimension(6)    , intent(in)    :: type_send_prim
        integer              , dimension(6)    , intent(in)    :: type_recv_prim
        integer              , dimension(6)    , intent(in)    :: my_neighbour
        integer                                , intent(in)    :: dims
        integer                                , intent(in)    :: mpi_comm_cart
        
        selectcase(mpi_opt_level)
        case(1) !SENDRECV

          if(cuda_aware) then
          call mpi_communicate3D_real8_sendrecv_cuda_aware(mpi_comm_cart,&
                my_neighbour,dims,&
                bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                var)

          else
          call mpi_communicate3D_real8_sendrecv(mpi_comm_cart,my_neighbour,dims,&
                bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                var)
          endif
        
        case(2) !ISEND/IRECV + BUFFERS
        
          if(cuda_aware) then
          call mpi_communicate3D_real8_bfr_cuda_aware(mpi_comm_cart, &
                my_neighbour,dims, &
                bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                var)
          else
          call mpi_communicate3D_real8_bfr(mpi_comm_cart,my_neighbour,dims, &
                bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                var)
          endif

        case(3)
          ! ISEND/IRECV + DERIVED DATA TYPES
          if(cuda_aware) then      
          call mpi_communicate3D_real8_cuda_aware(mpi_comm_cart,&
                     type_send_prim,type_recv_prim,my_neighbour,dims,var)
          else
          call mpi_communicate3D_real8(mpi_comm_cart,&
                     type_send_prim,type_recv_prim,my_neighbour,dims,var)
          endif
        case default

          print*, 'mpi_opt_level equal to ', mpi_opt_level, ' is not implemented'
          stop

        endselect

        return
end subroutine mpi_share3D




subroutine mpi_shareInt3D(mpi_comm_cart,type_send_prim,type_recv_prim,my_neighbour,dims, &
                bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                var)

        implicit none
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: var
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: bfr_send_E
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: bfr_send_W
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: bfr_recv_E
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: bfr_recv_W
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: bfr_send_N
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: bfr_send_S
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: bfr_recv_N
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: bfr_recv_S
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: bfr_send_B
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: bfr_send_F
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: bfr_recv_B
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: bfr_recv_F
        integer                , dimension(6)    , intent(in)    :: type_send_prim
        integer                , dimension(6)    , intent(in)    :: type_recv_prim
        integer                , dimension(6)    , intent(in)    :: my_neighbour
        integer                                  , intent(in)    :: dims
        integer                                  , intent(in)    :: mpi_comm_cart
        
        selectcase(mpi_opt_level)
        case(1) !SENDRECV
        if(cuda_aware) then
        call mpi_communicateInt3D_sendrecv_cuda_aware(mpi_comm_cart,&
                        my_neighbour,dims,&
                        bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                        bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                        bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                        var)
        else
        call mpi_communicateInt3D_sendrecv(mpi_comm_cart,my_neighbour,dims,&
                        bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                        bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                        bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                        var)
        endif

        case(2) !ISEND/IRECV + BUFFERS
        if(cuda_aware) then
        call mpi_communicateInt3D_isend_irecv_bfr_cuda_aware(mpi_comm_cart,&
                        my_neighbour,dims,&
                        bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                        bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                        bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                        var)
        else
        call mpi_communicateInt3D_isend_irecv_bfr(mpi_comm_cart,my_neighbour,dims,&
                        bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                        bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                        bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                        var)
        endif

        case(3) ! ISEND/IRECV + DERIVED DATA TYPES
         if(cuda_aware) then
         call mpi_communicateInt3D_isend_irecv_ddt_cuda_aware(mpi_comm_cart,&
                    type_send_prim,type_recv_prim,my_neighbour,dims,var)
         else
         call mpi_communicateInt3D_isend_irecv_ddt(mpi_comm_cart,&
                    type_send_prim,type_recv_prim,my_neighbour,dims,var)
         endif
        case default

          print*, 'mpi_opt_level equal to ', mpi_opt_level, ' is not implemented'
          stop

        endselect

        return
end subroutine mpi_shareInt3D





























































subroutine mpi_communicate3D_real8(comm,tp_send,tp_recv,neigh,dims,var)

        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: var     
        integer              , dimension(6)    , intent(in)    :: tp_send
        integer              , dimension(6)    , intent(in)    :: tp_recv
        integer              , dimension(6)    , intent(in)    :: neigh
        integer                                , intent(in)    :: dims
        integer                                , intent(in)    :: comm

        ! local declarations
        integer, dimension(4)       :: req_xx
        integer, dimension(2**dims) :: req_yz

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0
        
        !$acc update host(var)
        !
        ! === E/W communications
        !
        call MPI_IRECV(var,1,tp_recv(W),neigh(W),tag+1,comm,req_xx(1),err)
        call MPI_IRECV(var,1,tp_recv(E),neigh(E),tag+2,comm,req_xx(2),err)

        call MPI_ISEND(var,1,tp_send(E),neigh(E),tag+1,comm,req_xx(3),err)
        call MPI_ISEND(var,1,tp_send(W),neigh(W),tag+2,comm,req_xx(4),err)

        ! wait till communications end (corner needed)
        call MPI_waitall(size(req_xx), req_xx, MPI_STATUSES_IGNORE, err)

        !
        ! === N/S communications
        !
        call MPI_IRECV(var,1,tp_recv(S),neigh(S),tag+3,comm,req_yz(1),err)
        call MPI_IRECV(var,1,tp_recv(N),neigh(N),tag+4,comm,req_yz(2),err)

        call MPI_ISEND(var,1,tp_send(N),neigh(N),tag+3,comm,req_yz(3),err)
        call MPI_ISEND(var,1,tp_send(S),neigh(S),tag+4,comm,req_yz(4),err)

        !
        ! === B/F communications
        !
        if(dims == 3) then
        call MPI_IRECV(var,1,tp_recv(B),neigh(B),tag+5,comm,req_yz(5),err)
        call MPI_IRECV(var,1,tp_recv(F),neigh(F),tag+6,comm,req_yz(6),err)

        call MPI_ISEND(var,1,tp_send(F),neigh(F),tag+5,comm,req_yz(7),err)
        call MPI_ISEND(var,1,tp_send(B),neigh(B),tag+6,comm,req_yz(8),err)
        endif

        call MPI_waitall(size(req_yz), req_yz, MPI_STATUSES_IGNORE, err)

        !$acc update device(var)

        return
end subroutine mpi_communicate3D_real8






subroutine mpi_communicate3D_real8_cuda_aware(comm,tp_send,tp_recv,neigh,dims,var)

        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: var     
        integer              , dimension(6)    , intent(in)    :: tp_send
        integer              , dimension(6)    , intent(in)    :: tp_recv
        integer              , dimension(6)    , intent(in)    :: neigh
        integer                                , intent(in)    :: dims
        integer                                , intent(in)    :: comm

        ! local declarations
        integer, dimension(4)       :: req_xx
        integer, dimension(2**dims) :: req_yz

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0
       
        !$acc host_data use_device(var)
        !
        ! === E/W communications
        !
        call MPI_IRECV(var,1,tp_recv(W),neigh(W),tag+1,comm,req_xx(1),err)
        call MPI_IRECV(var,1,tp_recv(E),neigh(E),tag+2,comm,req_xx(2),err)

        call MPI_ISEND(var,1,tp_send(E),neigh(E),tag+1,comm,req_xx(3),err)
        call MPI_ISEND(var,1,tp_send(W),neigh(W),tag+2,comm,req_xx(4),err)

        ! wait till communications end (corner needed)
        call MPI_waitall(size(req_xx), req_xx, MPI_STATUSES_IGNORE, err)

        !
        ! === N/S communications
        !
        call MPI_IRECV(var,1,tp_recv(S),neigh(S),tag+3,comm,req_yz(1),err)
        call MPI_IRECV(var,1,tp_recv(N),neigh(N),tag+4,comm,req_yz(2),err)

        call MPI_ISEND(var,1,tp_send(N),neigh(N),tag+3,comm,req_yz(3),err)
        call MPI_ISEND(var,1,tp_send(S),neigh(S),tag+4,comm,req_yz(4),err)

        !
        ! === B/F communications
        !
        if(dims == 3) then
        call MPI_IRECV(var,1,tp_recv(B),neigh(B),tag+5,comm,req_yz(5),err)
        call MPI_IRECV(var,1,tp_recv(F),neigh(F),tag+6,comm,req_yz(6),err)

        call MPI_ISEND(var,1,tp_send(F),neigh(F),tag+5,comm,req_yz(7),err)
        call MPI_ISEND(var,1,tp_send(B),neigh(B),tag+6,comm,req_yz(8),err)
        endif

        call MPI_waitall(size(req_yz), req_yz, MPI_STATUSES_IGNORE, err)

        !$acc end host_data

        return
end subroutine mpi_communicate3D_real8_cuda_aware




































subroutine mpi_communicate4D_real8(comm,tp_send,tp_recv,neigh,dims,var)

        implicit none
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: var     
        integer              , dimension(6)       , intent(in)    :: tp_send
        integer              , dimension(6)       , intent(in)    :: tp_recv
        integer              , dimension(6)       , intent(in)    :: neigh
        integer                                   , intent(in)    :: dims
        integer                                   , intent(in)    :: comm

        ! local 
        integer, dimension(4)       :: req_xx
        integer, dimension(2**dims) :: req_yz

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0

        !$acc update host(var)
        !
        ! === E/W communications
        !
        call MPI_IRECV(var,1,tp_recv(W),neigh(W),tag+1,comm,req_xx(1),err)
        call MPI_IRECV(var,1,tp_recv(E),neigh(E),tag+2,comm,req_xx(2),err)

        call MPI_ISEND(var,1,tp_send(E),neigh(E),tag+1,comm,req_xx(3),err)
        call MPI_ISEND(var,1,tp_send(W),neigh(W),tag+2,comm,req_xx(4),err)

        ! wait till communications end (corner needed)
        call MPI_waitall(size(req_xx), req_xx, MPI_STATUSES_IGNORE, err)

        !
        ! === N/S communications
        !
        call MPI_IRECV(var,1,tp_recv(S),neigh(S),tag+3,comm,req_yz(1),err)
        call MPI_IRECV(var,1,tp_recv(N),neigh(N),tag+4,comm,req_yz(2),err)

        call MPI_ISEND(var,1,tp_send(N),neigh(N),tag+3,comm,req_yz(3),err)
        call MPI_ISEND(var,1,tp_send(S),neigh(S),tag+4,comm,req_yz(4),err)

        !
        ! === B/F communications
        !
        if(dims == 3) then
        call MPI_IRECV(var,1,tp_recv(B),neigh(B),tag+5,comm,req_yz(5),err)
        call MPI_IRECV(var,1,tp_recv(F),neigh(F),tag+6,comm,req_yz(6),err)

        call MPI_ISEND(var,1,tp_send(F),neigh(F),tag+5,comm,req_yz(7),err)
        call MPI_ISEND(var,1,tp_send(B),neigh(B),tag+6,comm,req_yz(8),err)
        endif

        call MPI_waitall(size(req_yz), req_yz, MPI_STATUSES_IGNORE, err)

        !$acc update device(var)

        return
end subroutine mpi_communicate4D_real8



subroutine mpi_communicate4D_real8_cuda_aware(comm,tp_send,tp_recv,neigh,dims,var)

        implicit none
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: var     
        integer              , dimension(6)       , intent(in)    :: tp_send
        integer              , dimension(6)       , intent(in)    :: tp_recv
        integer              , dimension(6)       , intent(in)    :: neigh
        integer                                   , intent(in)    :: dims
        integer                                   , intent(in)    :: comm

        ! local 
        integer, dimension(4)       :: req_xx
        integer, dimension(2**dims) :: req_yz

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0

        !$acc host_data use_device(var)
        !
        ! === E/W communications
        !
        call MPI_IRECV(var,1,tp_recv(W),neigh(W),tag+1,comm,req_xx(1),err)
        call MPI_IRECV(var,1,tp_recv(E),neigh(E),tag+2,comm,req_xx(2),err)

        call MPI_ISEND(var,1,tp_send(E),neigh(E),tag+1,comm,req_xx(3),err)
        call MPI_ISEND(var,1,tp_send(W),neigh(W),tag+2,comm,req_xx(4),err)

        ! wait till communications end (corner needed)
        call MPI_waitall(size(req_xx), req_xx, MPI_STATUSES_IGNORE, err)

        !
        ! === N/S communications
        !
        call MPI_IRECV(var,1,tp_recv(S),neigh(S),tag+3,comm,req_yz(1),err)
        call MPI_IRECV(var,1,tp_recv(N),neigh(N),tag+4,comm,req_yz(2),err)

        call MPI_ISEND(var,1,tp_send(N),neigh(N),tag+3,comm,req_yz(3),err)
        call MPI_ISEND(var,1,tp_send(S),neigh(S),tag+4,comm,req_yz(4),err)

        !
        ! === B/F communications
        !
        if(dims == 3) then
        call MPI_IRECV(var,1,tp_recv(B),neigh(B),tag+5,comm,req_yz(5),err)
        call MPI_IRECV(var,1,tp_recv(F),neigh(F),tag+6,comm,req_yz(6),err)

        call MPI_ISEND(var,1,tp_send(F),neigh(F),tag+5,comm,req_yz(7),err)
        call MPI_ISEND(var,1,tp_send(B),neigh(B),tag+6,comm,req_yz(8),err)
        endif

        call MPI_waitall(size(req_yz), req_yz, MPI_STATUSES_IGNORE, err)

        !$acc end host_data

        return
end subroutine mpi_communicate4D_real8_cuda_aware




















subroutine mpi_communicate4D_real8_bfr(comm,neigh,dims,&
                                       bfr_s_E, bfr_s_W, bfr_r_E, bfr_r_W, &
                                       bfr_s_N, bfr_s_S, bfr_r_N, bfr_r_S, &
                                       bfr_s_B, bfr_s_F, bfr_r_B, bfr_r_F, &
                                       var)

        implicit none
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: var     
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_E
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_W
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_E
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_W
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_N
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_S
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_N
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_S
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_B
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_F
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_B
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_F
        integer              , dimension(6)       , intent(in)    :: neigh
        integer                                   , intent(in)    :: dims
        integer                                   , intent(in)    :: comm

        ! local declaration
        integer, parameter    :: dtype = MPI_RP
        integer, dimension(4) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer, dimension(4) :: req_xx
        integer, dimension(4) :: req_yy
        integer, dimension(4) :: req_zz
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: i,j,k,l,err = 0, ii, jj, kk

        !
        ! === buffer dimensions
        !
        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)*bfr_dims_x(4)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)*bfr_dims_y(4)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)*bfr_dims_z(4)


        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN    
                    ii = ex-GN+i  !ex+1-GN,ex
                    bfr_s_E(ii,j,k,l) = var(ii,j,k,l)
                    ii = sx-1+i   !sx,sx+GN-1
                    bfr_s_W(ii,j,k,l) = var(ii,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN    
                 do i = lbx,ubx
                    jj = ey-GN+j !ey+1-GN,ey
                    bfr_s_N(i,jj,k,l) = var(i,jj,k,l)
                    jj = sy-1+j  !sy,sy+GN-1
                    bfr_s_S(i,jj,k,l) = var(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN 
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = ez-GN+k  !ez+1-GN,ez
                    bfr_s_F(i,j,kk,l) = var(i,j,kk,l)
                    kk = sz-1+k   !sz,sz+GN-1
                    bfr_s_B(i,j,kk,l) = var(i,j,kk,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc wait(1)
        !$acc update host(bfr_s_E) async(4)
        !$acc update host(bfr_s_W) async(5)

        !$acc wait(2)
        !$acc update host(bfr_s_N) async(6)
        !$acc update host(bfr_s_S) async(7)

        !$acc wait(3)
        !$acc update host(bfr_s_F) async(8)
        !$acc update host(bfr_s_B) async(9)

        !
        ! === east - west communications
        !
        call MPI_IRECV(bfr_r_W,bfr_size_x,dtype,neigh(W),tag+1,comm,req_xx(2),err)
        call MPI_IRECV(bfr_r_E,bfr_size_x,dtype,neigh(E),tag+2,comm,req_xx(4),err)

        !$acc wait(4)
        call MPI_ISEND(bfr_s_E,bfr_size_x,dtype,neigh(E),tag+1,comm,req_xx(1),err)
        !$acc wait(5)
        call MPI_ISEND(bfr_s_W,bfr_size_x,dtype,neigh(W),tag+2,comm,req_xx(3),err)

        ! wait till communications end (corner needed)
        call MPI_waitall(size(req_xx), req_xx, MPI_STATUSES_IGNORE, err)
        !$acc update device(bfr_r_W,bfr_r_E) async(10)

        !
        ! ==== nord-south communications
        !
        call MPI_IRECV(bfr_r_S,bfr_size_y,dtype,neigh(S),tag+3,comm,req_yy(2),err)
        call MPI_IRECV(bfr_r_N,bfr_size_y,dtype,neigh(N),tag+4,comm,req_yy(4),err)

        !$acc wait(6)
        call MPI_ISEND(bfr_s_N,bfr_size_y,dtype,neigh(N),tag+3,comm,req_yy(1),err)
        !$acc wait(7)
        call MPI_ISEND(bfr_s_S,bfr_size_y,dtype,neigh(S),tag+4,comm,req_yy(3),err)

        call MPI_waitall(size(req_yy), req_yy, MPI_STATUSES_IGNORE, err)
        !$acc update device(bfr_r_N,bfr_r_S) async(11)

        !
        ! ==== backward-forwvar communications
        !
        call MPI_IRECV(bfr_r_B,bfr_size_z,dtype,neigh(B),tag+5,comm,req_zz(2),err)
        call MPI_IRECV(bfr_r_F,bfr_size_z,dtype,neigh(F),tag+6,comm,req_zz(4),err)

        !$acc wait(8)
        call MPI_ISEND(bfr_s_F,bfr_size_z,dtype,neigh(F),tag+5,comm,req_zz(1),err)
        !$acc wait(9)
        call MPI_ISEND(bfr_s_B,bfr_size_z,dtype,neigh(B),tag+6,comm,req_zz(3),err)

        call MPI_waitall(size(req_zz), req_zz, MPI_STATUSES_IGNORE, err)
        !$acc update device(bfr_r_F,bfr_r_B) async(12)

        !$acc wait(10)
        if(neigh(W)>=0) then
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN 
                    ii = sx-GN+i-1  !sx-GN,sx-1
                    var(ii,j,k,l) = bfr_r_W(ii,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(neigh(E)>=0) then
        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN 
                    ii = ex+i       !ex+1,ex+GN
                    var(ii,j,k,l) = bfr_r_E(ii,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait(11)
        if(neigh(S)>=0) then
        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN 
                 do i = lbx,ubx
                    jj = sy-GN+j-1    !sy-GN,sy-1
                    var(i,jj,k,l) = bfr_r_S(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(neigh(N)>=0) then
        !$acc parallel default(present) async(4)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN 
                 do i = lbx,ubx
                    jj = ey+j         !ey+1,ey+GN
                    var(i,jj,k,l) = bfr_r_N(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait(12)
        if(neigh(B)>=0) then
        !$acc parallel default(present) async(5)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = sz-GN+k-1    !sz-GN,sz-1
                    var(i,j,kk,l) = bfr_r_B(i,j,kk,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel 
        endif

        if(neigh(F)>=0) then
        !$acc parallel default(present) async(6)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = ez+k         !ez+1,ez+GN
                    var(i,j,kk,l) = bfr_r_F(i,j,kk,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel 
        endif

        if(neigh(W)>=0) then
          !$acc wait(1)
        endif
        if(neigh(E)>=0) then
          !$acc wait(2)
        endif
        if(neigh(S)>=0) then
          !$acc wait(3)
        endif
        if(neigh(N)>=0) then
          !$acc wait(4)
        endif
        if(neigh(B)>=0) then
          !$acc wait(5)
        endif
        if(neigh(F)>=0) then
          !$acc wait(6)
        endif


        return
end subroutine mpi_communicate4D_real8_bfr




subroutine mpi_communicate4D_real8_bfr_cuda_aware(comm,neigh,dims,&
                                       bfr_s_E, bfr_s_W, bfr_r_E, bfr_r_W, &
                                       bfr_s_N, bfr_s_S, bfr_r_N, bfr_r_S, &
                                       bfr_s_B, bfr_s_F, bfr_r_B, bfr_r_F, &
                                       var)

        implicit none
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: var     
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_E
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_W
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_E
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_W
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_N
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_S
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_N
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_S
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_B
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_F
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_B
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_F
        integer              , dimension(6)       , intent(in)    :: neigh
        integer                                   , intent(in)    :: dims
        integer                                   , intent(in)    :: comm

        ! local declaration
        integer, parameter    :: dtype = MPI_RP
        integer, dimension(4) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer, dimension(4) :: req_xx
        integer, dimension(4) :: req_yy
        integer, dimension(4) :: req_zz
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: i,j,k,l,err = 0, ii, jj, kk

        !
        ! === buffer dimensions
        !
        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)*bfr_dims_x(4)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)*bfr_dims_y(4)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)*bfr_dims_z(4)


        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN    
                    ii = ex-GN+i  !ex+1-GN,ex
                    bfr_s_E(ii,j,k,l) = var(ii,j,k,l)
                    ii = sx-1+i   !sx,sx+GN-1
                    bfr_s_W(ii,j,k,l) = var(ii,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN    
                 do i = lbx,ubx
                    jj = ey-GN+j !ey+1-GN,ey
                    bfr_s_N(i,jj,k,l) = var(i,jj,k,l)
                    jj = sy-1+j  !sy,sy+GN-1
                    bfr_s_S(i,jj,k,l) = var(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN 
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = ez-GN+k  !ez+1-GN,ez
                    bfr_s_F(i,j,kk,l) = var(i,j,kk,l)
                    kk = sz-1+k   !sz,sz+GN-1
                    bfr_s_B(i,j,kk,l) = var(i,j,kk,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel


        !$acc host_data use_device(bfr_s_E,bfr_r_W) &
        !$acc           use_device(bfr_s_W,bfr_r_E) &
        !$acc           use_device(bfr_s_S,bfr_r_N) &
        !$acc           use_device(bfr_s_N,bfr_r_S) &
        !$acc           use_device(bfr_s_B,bfr_r_F) &
        !$acc           use_device(bfr_s_F,bfr_r_B)

        !
        ! === east - west communications
        !
        call MPI_IRECV(bfr_r_W,bfr_size_x,dtype,neigh(W),tag+1,comm,req_xx(2),err)
        call MPI_IRECV(bfr_r_E,bfr_size_x,dtype,neigh(E),tag+2,comm,req_xx(4),err)
        
        !$acc wait(1)
        call MPI_ISEND(bfr_s_E,bfr_size_x,dtype,neigh(E),tag+1,comm,req_xx(1),err)
        call MPI_ISEND(bfr_s_W,bfr_size_x,dtype,neigh(W),tag+2,comm,req_xx(3),err)

        ! wait till communications end (corner needed)
        call MPI_waitall(size(req_xx), req_xx, MPI_STATUSES_IGNORE, err)

        !
        ! ==== nord-south communications
        !
        call MPI_IRECV(bfr_r_S,bfr_size_y,dtype,neigh(S),tag+3,comm,req_yy(2),err)
        call MPI_IRECV(bfr_r_N,bfr_size_y,dtype,neigh(N),tag+4,comm,req_yy(4),err)

        !$acc wait(2)
        call MPI_ISEND(bfr_s_N,bfr_size_y,dtype,neigh(N),tag+3,comm,req_yy(1),err)
        call MPI_ISEND(bfr_s_S,bfr_size_y,dtype,neigh(S),tag+4,comm,req_yy(3),err)

        call MPI_waitall(size(req_yy), req_yy, MPI_STATUSES_IGNORE, err)

        !
        ! ==== backward-forwvar communications
        !
        call MPI_IRECV(bfr_r_B,bfr_size_z,dtype,neigh(B),tag+5,comm,req_zz(2),err)
        call MPI_IRECV(bfr_r_F,bfr_size_z,dtype,neigh(F),tag+6,comm,req_zz(4),err)

        !$acc wait(3)
        call MPI_ISEND(bfr_s_F,bfr_size_z,dtype,neigh(F),tag+5,comm,req_zz(1),err)
        call MPI_ISEND(bfr_s_B,bfr_size_z,dtype,neigh(B),tag+6,comm,req_zz(3),err)

        call MPI_waitall(size(req_zz), req_zz, MPI_STATUSES_IGNORE, err)

        !$acc end host_data

        if(neigh(W)>=0) then
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN 
                    ii = sx-GN+i-1  !sx-GN,sx-1
                    var(ii,j,k,l) = bfr_r_W(ii,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(neigh(E)>=0) then
        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN 
                    ii = ex+i       !ex+1,ex+GN
                    var(ii,j,k,l) = bfr_r_E(ii,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(neigh(S)>=0) then
        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN 
                 do i = lbx,ubx
                    jj = sy-GN+j-1    !sy-GN,sy-1
                    var(i,jj,k,l) = bfr_r_S(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif
        
        if(neigh(N)>=0) then
        !$acc parallel default(present) async(4)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN 
                 do i = lbx,ubx
                    jj = ey+j         !ey+1,ey+GN
                    var(i,jj,k,l) = bfr_r_N(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(neigh(B)>=0) then
        !$acc parallel default(present) async(5)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = sz-GN+k-1    !sz-GN,sz-1
                    var(i,j,kk,l) = bfr_r_B(i,j,kk,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel 
        endif

        if(neigh(F)>=0) then
        !$acc parallel default(present) async(6)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = ez+k         !ez+1,ez+GN
                    var(i,j,kk,l) = bfr_r_F(i,j,kk,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel 
        endif

        if(neigh(W)>=0) then
          !$acc wait(1)
        endif
        if(neigh(E)>=0) then
          !$acc wait(2)
        endif
        if(neigh(S)>=0) then
          !$acc wait(3)
        endif
        if(neigh(N)>=0) then
          !$acc wait(4)
        endif
        if(neigh(B)>=0) then
          !$acc wait(5)
        endif
        if(neigh(F)>=0) then
          !$acc wait(6)
        endif


        return
end subroutine mpi_communicate4D_real8_bfr_cuda_aware


































subroutine mpi_communicate3D_real8_bfr(comm,neigh,dims,&
                                       bfr_s_E, bfr_s_W, bfr_r_E, bfr_r_W, &
                                       bfr_s_N, bfr_s_S, bfr_r_N, bfr_r_S, &
                                       bfr_s_B, bfr_s_F, bfr_r_B, bfr_r_F, &
                                       var)

        implicit none
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: var     
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_E
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_W
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_E
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_W
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_N
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_S
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_N
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_S
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_B
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_F
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_B
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_F
        integer              , dimension(6)     , intent(in)    :: neigh
        integer                                 , intent(in)    :: dims
        integer                                 , intent(in)    :: comm

        ! local declaration
        integer, parameter    :: dtype = MPI_RP
        integer, dimension(3) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer, dimension(4) :: req_xx
        integer, dimension(4) :: req_yy
        integer, dimension(4) :: req_zz
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: i,j,k,err = 0, ii, jj, kk

        !
        ! === buffer dimensions
        !
        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)


        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN    
                 ii = ex-GN+i  !ex+1-GN,ex
                 bfr_s_E(ii,j,k) = var(ii,j,k)
                 ii = sx-1+i   !sx,sx+GN-1
                 bfr_s_W(ii,j,k) = var(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN    
              do i = lbx,ubx
                 jj = ey-GN+j !ey+1-GN,ey
                 bfr_s_N(i,jj,k) = var(i,jj,k)
                 jj = sy-1+j  !sy,sy+GN-1
                 bfr_s_S(i,jj,k) = var(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN 
           do    j = lby,uby
              do i = lbx,ubx
                 kk = ez-GN+k  !ez+1-GN,ez
                 bfr_s_F(i,j,kk) = var(i,j,kk)
                 kk = sz-1+k   !sz,sz+GN-1
                 bfr_s_B(i,j,kk) = var(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc wait(1)
        !$acc update host(bfr_s_E) async(4)
        !$acc update host(bfr_s_W) async(5)

        !$acc wait(2)
        !$acc update host(bfr_s_N) async(6)
        !$acc update host(bfr_s_S) async(7)

        !$acc wait(3)
        !$acc update host(bfr_s_F) async(8)
        !$acc update host(bfr_s_B) async(9)

        !
        ! === east - west communications
        !
        call MPI_IRECV(bfr_r_W,bfr_size_x,dtype,neigh(W),tag+1,comm,req_xx(2),err)
        call MPI_IRECV(bfr_r_E,bfr_size_x,dtype,neigh(E),tag+2,comm,req_xx(4),err)

        !$acc wait(4)
        call MPI_ISEND(bfr_s_E,bfr_size_x,dtype,neigh(E),tag+1,comm,req_xx(1),err)
        !$acc wait(5)
        call MPI_ISEND(bfr_s_W,bfr_size_x,dtype,neigh(W),tag+2,comm,req_xx(3),err)

        ! wait till communications end (corner needed)
        call MPI_waitall(size(req_xx), req_xx, MPI_STATUSES_IGNORE, err)
        !$acc update device(bfr_r_W,bfr_r_E) async(10)

        !
        ! ==== nord-south communications
        !
        call MPI_IRECV(bfr_r_S,bfr_size_y,dtype,neigh(S),tag+3,comm,req_yy(2),err)
        call MPI_IRECV(bfr_r_N,bfr_size_y,dtype,neigh(N),tag+4,comm,req_yy(4),err)

        !$acc wait(6)
        call MPI_ISEND(bfr_s_N,bfr_size_y,dtype,neigh(N),tag+3,comm,req_yy(1),err)
        !$acc wait(7)
        call MPI_ISEND(bfr_s_S,bfr_size_y,dtype,neigh(S),tag+4,comm,req_yy(3),err)

        call MPI_waitall(size(req_yy), req_yy, MPI_STATUSES_IGNORE, err)
        !$acc update device(bfr_r_N,bfr_r_S) async(11)

        !
        ! ==== backward-forwvar communications
        !
        call MPI_IRECV(bfr_r_B,bfr_size_z,dtype,neigh(B),tag+5,comm,req_zz(2),err)
        call MPI_IRECV(bfr_r_F,bfr_size_z,dtype,neigh(F),tag+6,comm,req_zz(4),err)

        !$acc wait(8)
        call MPI_ISEND(bfr_s_F,bfr_size_z,dtype,neigh(F),tag+5,comm,req_zz(1),err)
        !$acc wait(9)
        call MPI_ISEND(bfr_s_B,bfr_size_z,dtype,neigh(B),tag+6,comm,req_zz(3),err)

        call MPI_waitall(size(req_zz), req_zz, MPI_STATUSES_IGNORE, err)
        !$acc update device(bfr_r_F,bfr_r_B) async(12)

        !$acc wait(10)
        if(neigh(W)>=0) then
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN 
                 ii = sx-GN+i-1  !sx-GN,sx-1
                 var(ii,j,k) = bfr_r_W(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(neigh(E)>=0) then
        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN 
                 ii = ex+i       !ex+1,ex+GN
                 var(ii,j,k) = bfr_r_E(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait(11)
        if(neigh(S)>=0) then
        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN 
              do i = lbx,ubx
                 jj = sy-GN+j-1    !sy-GN,sy-1
                 var(i,jj,k) = bfr_r_S(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(neigh(N)>=0) then
        !$acc parallel default(present) async(4)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN 
              do i = lbx,ubx
                 jj = ey+j         !ey+1,ey+GN
                 var(i,jj,k) = bfr_r_N(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait(12)
        if(neigh(B)>=0) then
        !$acc parallel default(present) async(5)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN
           do    j = lby,uby
              do i = lbx,ubx
                 kk = sz-GN+k-1    !sz-GN,sz-1
                 var(i,j,kk) = bfr_r_B(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel 
        endif

        if(neigh(F)>=0) then
        !$acc parallel default(present) async(6)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN
           do    j = lby,uby
              do i = lbx,ubx
                 kk = ez+k         !ez+1,ez+GN
                 var(i,j,kk) = bfr_r_F(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel 
        endif

        if(neigh(W)>=0) then
          !$acc wait(1)
        endif
        if(neigh(E)>=0) then
          !$acc wait(2)
        endif
        if(neigh(S)>=0) then
          !$acc wait(3)
        endif
        if(neigh(N)>=0) then
          !$acc wait(4)
        endif
        if(neigh(B)>=0) then
          !$acc wait(5)
        endif
        if(neigh(F)>=0) then
          !$acc wait(6)
        endif
        
        return
end subroutine mpi_communicate3D_real8_bfr




subroutine mpi_communicate3D_real8_bfr_cuda_aware(comm,neigh,dims,&
                                       bfr_s_E, bfr_s_W, bfr_r_E, bfr_r_W, &
                                       bfr_s_N, bfr_s_S, bfr_r_N, bfr_r_S, &
                                       bfr_s_B, bfr_s_F, bfr_r_B, bfr_r_F, &
                                       var)

        implicit none
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: var     
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_E
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_W
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_E
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_W
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_N
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_S
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_N
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_S
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_B
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_F
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_B
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_F
        integer              , dimension(6)     , intent(in)    :: neigh
        integer                                 , intent(in)    :: dims
        integer                                 , intent(in)    :: comm

        ! local declaration
        integer, parameter    :: dtype = MPI_RP
        integer, dimension(3) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer, dimension(4) :: req_xx
        integer, dimension(4) :: req_yy
        integer, dimension(4) :: req_zz
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: i,j,k,err = 0, ii, jj, kk

        !
        ! === buffer dimensions
        !
        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)


        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN    
                 ii = ex-GN+i  !ex+1-GN,ex
                 bfr_s_E(ii,j,k) = var(ii,j,k)
                 ii = sx-1+i   !sx,sx+GN-1
                 bfr_s_W(ii,j,k) = var(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN    
              do i = lbx,ubx
                 jj = ey-GN+j !ey+1-GN,ey
                 bfr_s_N(i,jj,k) = var(i,jj,k)
                 jj = sy-1+j  !sy,sy+GN-1
                 bfr_s_S(i,jj,k) = var(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN 
           do    j = lby,uby
              do i = lbx,ubx
                 kk = ez-GN+k  !ez+1-GN,ez
                 bfr_s_F(i,j,kk) = var(i,j,kk)
                 kk = sz-1+k   !sz,sz+GN-1
                 bfr_s_B(i,j,kk) = var(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc host_data use_device(bfr_s_E,bfr_r_W) &
        !$acc           use_device(bfr_s_W,bfr_r_E) &
        !$acc           use_device(bfr_s_S,bfr_r_N) &
        !$acc           use_device(bfr_s_N,bfr_r_S) &
        !$acc           use_device(bfr_s_B,bfr_r_F) &
        !$acc           use_device(bfr_s_F,bfr_r_B)

        !
        ! === east - west communications
        !
        call MPI_IRECV(bfr_r_W,bfr_size_x,dtype,neigh(W),tag+1,comm,req_xx(2),err)
        call MPI_IRECV(bfr_r_E,bfr_size_x,dtype,neigh(E),tag+2,comm,req_xx(4),err)

        !$acc wait(1)
        call MPI_ISEND(bfr_s_E,bfr_size_x,dtype,neigh(E),tag+1,comm,req_xx(1),err)
        call MPI_ISEND(bfr_s_W,bfr_size_x,dtype,neigh(W),tag+2,comm,req_xx(3),err)

        ! wait till communications end (corner needed)
        call MPI_waitall(size(req_xx), req_xx, MPI_STATUSES_IGNORE, err)

        !
        ! ==== nord-south communications
        !
        call MPI_IRECV(bfr_r_S,bfr_size_y,dtype,neigh(S),tag+3,comm,req_yy(2),err)
        call MPI_IRECV(bfr_r_N,bfr_size_y,dtype,neigh(N),tag+4,comm,req_yy(4),err)

        !$acc wait(2)
        call MPI_ISEND(bfr_s_N,bfr_size_y,dtype,neigh(N),tag+3,comm,req_yy(1),err)
        call MPI_ISEND(bfr_s_S,bfr_size_y,dtype,neigh(S),tag+4,comm,req_yy(3),err)

        call MPI_waitall(size(req_yy), req_yy, MPI_STATUSES_IGNORE, err)

        !
        ! ==== backward-forwvar communications
        !
        call MPI_IRECV(bfr_r_B,bfr_size_z,dtype,neigh(B),tag+5,comm,req_zz(2),err)
        call MPI_IRECV(bfr_r_F,bfr_size_z,dtype,neigh(F),tag+6,comm,req_zz(4),err)

        !$acc wait(3)
        call MPI_ISEND(bfr_s_F,bfr_size_z,dtype,neigh(F),tag+5,comm,req_zz(1),err)
        call MPI_ISEND(bfr_s_B,bfr_size_z,dtype,neigh(B),tag+6,comm,req_zz(3),err)

        call MPI_waitall(size(req_zz), req_zz, MPI_STATUSES_IGNORE, err)

        !$acc end host_data

        if(neigh(W)>=0) then
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN 
                 ii = sx-GN+i-1  !sx-GN,sx-1
                 var(ii,j,k) = bfr_r_W(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(neigh(E)>=0) then
        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN 
                 ii = ex+i       !ex+1,ex+GN
                 var(ii,j,k) = bfr_r_E(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(neigh(S)>=0) then
        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN 
              do i = lbx,ubx
                 jj = sy-GN+j-1    !sy-GN,sy-1
                 var(i,jj,k) = bfr_r_S(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(neigh(N)>=0) then
        !$acc parallel default(present) async(4)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN 
              do i = lbx,ubx
                 jj = ey+j         !ey+1,ey+GN
                 var(i,jj,k) = bfr_r_N(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(neigh(B)>=0) then
        !$acc parallel default(present) async(5)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN
           do    j = lby,uby
              do i = lbx,ubx
                 kk = sz-GN+k-1    !sz-GN,sz-1
                 var(i,j,kk) = bfr_r_B(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel 
        endif

        if(neigh(F)>=0) then
        !$acc parallel default(present) async(6)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN
           do    j = lby,uby
              do i = lbx,ubx
                 kk = ez+k         !ez+1,ez+GN
                 var(i,j,kk) = bfr_r_F(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel 
        endif

        if(neigh(W)>=0) then
          !$acc wait(1)
        endif
        if(neigh(E)>=0) then
          !$acc wait(2)
        endif
        if(neigh(S)>=0) then
          !$acc wait(3)
        endif
        if(neigh(N)>=0) then
          !$acc wait(4)
        endif
        if(neigh(B)>=0) then
          !$acc wait(5)
        endif
        if(neigh(F)>=0) then
          !$acc wait(6)
        endif

        return
end subroutine mpi_communicate3D_real8_bfr_cuda_aware























subroutine mpi_communicate4D_real8_sendrecv(comm,neigh,dims,&
                                       bfr_s_E, bfr_s_W, bfr_r_E, bfr_r_W, &
                                       bfr_s_N, bfr_s_S, bfr_r_N, bfr_r_S, &
                                       bfr_s_B, bfr_s_F, bfr_r_B, bfr_r_F, &
                                       var)


        implicit none
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: var     
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_E
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_W
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_E
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_W
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_N
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_S
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_N
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_S
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_B
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_F
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_B
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_F
        integer              , dimension(6)       , intent(in)    :: neigh
        integer                                   , intent(in)    :: dims
        integer                                   , intent(in)    :: comm

        ! local declaration
        integer, parameter    :: dtype = MPI_RP
        integer, dimension(4) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, dimension(MPI_STATUS_SIZE) :: status
        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0, i,j,k,l, ii, jj, kk

        !
        ! === buffer dimensions
        !
        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)*bfr_dims_x(4)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)*bfr_dims_y(4)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)*bfr_dims_z(4)

        !
        ! === east - west communications
        !
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN    
                    ii = ex-GN+i  !ex+1-GN,ex
                    bfr_s_E(ii,j,k,l) = var(ii,j,k,l)
                    ii = sx-1+i   !sx,sx+GN-1
                    bfr_s_W(ii,j,k,l) = var(ii,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN    
                 do i = lbx,ubx
                    jj = ey-GN+j !ey+1-GN,ey
                    bfr_s_N(i,jj,k,l) = var(i,jj,k,l)
                    jj = sy-1+j  !sy,sy+GN-1
                    bfr_s_S(i,jj,k,l) = var(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN 
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = ez-GN+k  !ez+1-GN,ez
                    bfr_s_F(i,j,kk,l) = var(i,j,kk,l)
                    kk = sz-1+k   !sz,sz+GN-1
                    bfr_s_B(i,j,kk,l) = var(i,j,kk,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc wait(1)
        !$acc update host(bfr_s_E) async(4)
        !$acc update host(bfr_s_W) async(5)

        !$acc wait(2)
        !$acc update host(bfr_s_S) async(6)
        !$acc update host(bfr_s_N) async(7)

        !$acc wait(3)
        !$acc update host(bfr_s_B) async(8)
        !$acc update host(bfr_s_F) async(9)

        !$acc wait(4)
        call MPI_SENDRECV(bfr_s_E,bfr_size_x,dtype,neigh(E),tag+1, &
                          bfr_r_W,bfr_size_x,dtype,neigh(W),tag+1, &
                          comm,status,err)
        !$acc update device(bfr_r_W) async(10)

        !$acc wait(5)
        call MPI_SENDRECV(bfr_s_W,bfr_size_x,dtype,neigh(W),tag+2, &
                          bfr_r_E,bfr_size_x,dtype,neigh(E),tag+2, &
                          comm,status,err)
        !$acc update device(bfr_r_E) async(11)
        
        !$acc wait(6)
        call MPI_SENDRECV(bfr_s_S,bfr_size_y,dtype,neigh(S),tag+4, &
                          bfr_r_N,bfr_size_y,dtype,neigh(N),tag+4, &
                          comm,status,err)
        !$acc update device(bfr_r_N) async(12)

        !$acc wait(7)
        call MPI_SENDRECV(bfr_s_N,bfr_size_y,dtype,neigh(N),tag+3, &
                          bfr_r_S,bfr_size_y,dtype,neigh(S),tag+3, & 
                          comm,status,err)
        !$acc update device(bfr_r_S) async(13)

        !$acc wait(8)
        call MPI_SENDRECV(bfr_s_B,bfr_size_z,dtype,neigh(B),tag+6, &
                          bfr_r_F,bfr_size_z,dtype,neigh(F),tag+6, &
                          comm,status,err)
        !$acc update device(bfr_r_F) async(14)

        !$acc wait(9)
        call MPI_SENDRECV(bfr_s_F,bfr_size_z,dtype,neigh(F),tag+5, &
                          bfr_r_B,bfr_size_z,dtype,neigh(B),tag+5, &
                          comm,status,err)
        !$acc update device(bfr_r_B) async(15)

        !$acc wait(10)
        if(neigh(W) >=0) then
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN 
                    ii = sx-GN+i-1  !sx-GN,sx-1
                    var(ii,j,k,l) = bfr_r_W(ii,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait(11)
        if(neigh(E)>=0) then
        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN 
                    ii = ex+i       !ex+1,ex+GN
                    var(ii,j,k,l) = bfr_r_E(ii,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif


        !$acc wait(12)
        if(neigh(N)>=0) then
        !$acc parallel default(present) async(4)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN 
                 do i = lbx,ubx
                    jj = ey+j         !ey+1,ey+GN
                    var(i,jj,k,l) = bfr_r_N(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait(13)
        if(neigh(S)>=0) then
        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN 
                 do i = lbx,ubx
                    jj = sy-GN+j-1    !sy-GN,sy-1
                    var(i,jj,k,l) = bfr_r_S(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait(14)
        if(neigh(F)>=0) then
        !$acc parallel default(present) async(6)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = ez+k         !ez+1,ez+GN
                    var(i,j,kk,l) = bfr_r_F(i,j,kk,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel 
        endif

        !$acc wait(15)
        if(neigh(B)>=0) then
        !$acc parallel default(present) async(5)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = sz-GN+k-1    !sz-GN,sz-1
                    var(i,j,kk,l) = bfr_r_B(i,j,kk,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel 
        endif

        !$acc wait

        return
end subroutine mpi_communicate4D_real8_sendrecv



subroutine mpi_communicate4D_real8_sendrecv_cuda_aware(comm,neigh,dims,&
                                       bfr_s_E, bfr_s_W, bfr_r_E, bfr_r_W, &
                                       bfr_s_N, bfr_s_S, bfr_r_N, bfr_r_S, &
                                       bfr_s_B, bfr_s_F, bfr_r_B, bfr_r_F, &
                                       var)


        implicit none
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: var     
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_E
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_W
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_E
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_W
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_N
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_S
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_N
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_S
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_B
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_s_F
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_B
        real(rp), allocatable, dimension(:,:,:,:) , intent(inout) :: bfr_r_F
        integer              , dimension(6)       , intent(in)    :: neigh
        integer                                   , intent(in)    :: dims
        integer                                   , intent(in)    :: comm

        ! local declaration
        integer, parameter    :: dtype = MPI_RP
        integer, dimension(4) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, dimension(MPI_STATUS_SIZE) :: status
        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0, i,j,k,l, ii, jj, kk

        !
        ! === buffer dimensions
        !
        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)*bfr_dims_x(4)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)*bfr_dims_y(4)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)*bfr_dims_z(4)

        !
        ! === east - west communications
        !
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN    
                    ii = ex-GN+i  !ex+1-GN,ex
                    bfr_s_E(ii,j,k,l) = var(ii,j,k,l)
                    ii = sx-1+i   !sx,sx+GN-1
                    bfr_s_W(ii,j,k,l) = var(ii,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN    
                 do i = lbx,ubx
                    jj = ey-GN+j !ey+1-GN,ey
                    bfr_s_N(i,jj,k,l) = var(i,jj,k,l)
                    jj = sy-1+j  !sy,sy+GN-1
                    bfr_s_S(i,jj,k,l) = var(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN 
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = ez-GN+k  !ez+1-GN,ez
                    bfr_s_F(i,j,kk,l) = var(i,j,kk,l)
                    kk = sz-1+k   !sz,sz+GN-1
                    bfr_s_B(i,j,kk,l) = var(i,j,kk,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
       
        !$acc host_data use_device(bfr_s_E,bfr_r_W) &
        !$acc           use_device(bfr_s_W,bfr_r_E) &
        !$acc           use_device(bfr_s_S,bfr_r_N) &
        !$acc           use_device(bfr_s_N,bfr_r_S) &
        !$acc           use_device(bfr_s_B,bfr_r_F) &
        !$acc           use_device(bfr_s_F,bfr_r_B)

        !$acc wait(1)
        call MPI_SENDRECV(bfr_s_E,bfr_size_x,dtype,neigh(E),tag+1, &
                          bfr_r_W,bfr_size_x,dtype,neigh(W),tag+1, &
                          comm,status,err)

        call MPI_SENDRECV(bfr_s_W,bfr_size_x,dtype,neigh(W),tag+2, &
                          bfr_r_E,bfr_size_x,dtype,neigh(E),tag+2, &
                          comm,status,err)
       
        !$acc wait(2)
        call MPI_SENDRECV(bfr_s_S,bfr_size_y,dtype,neigh(S),tag+4, &
                          bfr_r_N,bfr_size_y,dtype,neigh(N),tag+4, &
                          comm,status,err)

        call MPI_SENDRECV(bfr_s_N,bfr_size_y,dtype,neigh(N),tag+3, &
                          bfr_r_S,bfr_size_y,dtype,neigh(S),tag+3, & 
                          comm,status,err)

        !$acc wait(3)
        call MPI_SENDRECV(bfr_s_B,bfr_size_z,dtype,neigh(B),tag+6, &
                          bfr_r_F,bfr_size_z,dtype,neigh(F),tag+6, &
                          comm,status,err)

        call MPI_SENDRECV(bfr_s_F,bfr_size_z,dtype,neigh(F),tag+5, &
                          bfr_r_B,bfr_size_z,dtype,neigh(B),tag+5, &
                          comm,status,err)
        !$acc end host_data

        if(neigh(W)>=0) then
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN 
                    ii = sx-GN+i-1  !sx-GN,sx-1
                    var(ii,j,k,l) = bfr_r_W(ii,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(neigh(E)>=0) then
        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN 
                    ii = ex+i       !ex+1,ex+GN
                    var(ii,j,k,l) = bfr_r_E(ii,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(neigh(S)>=0) then
        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN 
                 do i = lbx,ubx
                    jj = sy-GN+j-1    !sy-GN,sy-1
                    var(i,jj,k,l) = bfr_r_S(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(neigh(N)>=0) then
        !$acc parallel default(present) async(4)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN 
                 do i = lbx,ubx
                    jj = ey+j         !ey+1,ey+GN
                    var(i,jj,k,l) = bfr_r_N(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(neigh(B)>=0) then
        !$acc parallel default(present) async(5)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = sz-GN+k-1    !sz-GN,sz-1
                    var(i,j,kk,l) = bfr_r_B(i,j,kk,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel 
        endif

        if(neigh(F)>=0) then
        !$acc parallel default(present) async(6)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = ez+k         !ez+1,ez+GN
                    var(i,j,kk,l) = bfr_r_F(i,j,kk,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel 
        endif

        if(neigh(W)>=0) then
          !$acc wait(1)
        endif
        if(neigh(E)>=0) then
          !$acc wait(2)
        endif
        if(neigh(S)>=0) then
          !$acc wait(3)
        endif
        if(neigh(N)>=0) then
          !$acc wait(4)
        endif
        if(neigh(B)>=0) then
          !$acc wait(5)
        endif
        if(neigh(F)>=0) then
          !$acc wait(6)
        endif

        return
end subroutine mpi_communicate4D_real8_sendrecv_cuda_aware



















subroutine mpi_communicate3D_real8_sendrecv(comm,neigh,dims,&
                                       bfr_s_E, bfr_s_W, bfr_r_E, bfr_r_W, &
                                       bfr_s_N, bfr_s_S, bfr_r_N, bfr_r_S, &
                                       bfr_s_B, bfr_s_F, bfr_r_B, bfr_r_F, &
                                       var)


        implicit none
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: var     
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_E
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_W
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_E
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_W
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_N
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_S
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_N
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_S
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_B
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_F
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_B
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_F
        integer              , dimension(6)     , intent(in)    :: neigh
        integer                                 , intent(in)    :: dims
        integer                                 , intent(in)    :: comm

        ! local declaration
        integer, parameter    :: dtype = MPI_RP
        integer, dimension(3) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, dimension(MPI_STATUS_SIZE) :: status
        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0, i,j,k, ii, jj, kk

        !
        ! === buffer dimensions
        !
        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)

        !
        ! === east - west communications
        !
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN    
                 ii = ex-GN+i  !ex+1-GN,ex
                 bfr_s_E(ii,j,k) = var(ii,j,k)
                 ii = sx-1+i   !sx,sx+GN-1
                 bfr_s_W(ii,j,k) = var(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN    
              do i = lbx,ubx
                 jj = ey-GN+j !ey+1-GN,ey
                 bfr_s_N(i,jj,k) = var(i,jj,k)
                 jj = sy-1+j  !sy,sy+GN-1
                 bfr_s_S(i,jj,k) = var(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN 
           do    j = lby,uby
              do i = lbx,ubx
                 kk = ez-GN+k  !ez+1-GN,ez
                 bfr_s_F(i,j,kk) = var(i,j,kk)
                 kk = sz-1+k   !sz,sz+GN-1
                 bfr_s_B(i,j,kk) = var(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc wait(1)
        !$acc update host(bfr_s_E) async(4)
        !$acc update host(bfr_s_W) async(5)

        !$acc wait(2)
        !$acc update host(bfr_s_S) async(6)
        !$acc update host(bfr_s_N) async(7)

        !$acc wait(3)
        !$acc update host(bfr_s_B) async(8)
        !$acc update host(bfr_s_F) async(9)

        !$acc wait(4)
        call MPI_SENDRECV(bfr_s_E,bfr_size_x,dtype,neigh(E),tag+1, &
                          bfr_r_W,bfr_size_x,dtype,neigh(W),tag+1, &
                          comm,status,err)
        !$acc update device(bfr_r_W) async(10)

        !$acc wait(5)
        call MPI_SENDRECV(bfr_s_W,bfr_size_x,dtype,neigh(W),tag+2, &
                          bfr_r_E,bfr_size_x,dtype,neigh(E),tag+2, &
                          comm,status,err)
        !$acc update device(bfr_r_E) async(11)
        
        !$acc wait(6)
        call MPI_SENDRECV(bfr_s_S,bfr_size_y,dtype,neigh(S),tag+4, &
                          bfr_r_N,bfr_size_y,dtype,neigh(N),tag+4, &
                          comm,status,err)
        !$acc update device(bfr_r_N) async(12)

        !$acc wait(7)
        call MPI_SENDRECV(bfr_s_N,bfr_size_y,dtype,neigh(N),tag+3, &
                          bfr_r_S,bfr_size_y,dtype,neigh(S),tag+3, & 
                          comm,status,err)
        !$acc update device(bfr_r_S) async(13)

        !$acc wait(8)
        call MPI_SENDRECV(bfr_s_B,bfr_size_z,dtype,neigh(B),tag+6, &
                          bfr_r_F,bfr_size_z,dtype,neigh(F),tag+6, &
                          comm,status,err)
        !$acc update device(bfr_r_F) async(14)

        !$acc wait(9)
        call MPI_SENDRECV(bfr_s_F,bfr_size_z,dtype,neigh(F),tag+5, &
                          bfr_r_B,bfr_size_z,dtype,neigh(B),tag+5, &
                          comm,status,err)
        !$acc update device(bfr_r_B) async(15)

        !$acc wait(10)
        if(neigh(W)>=0) then
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN 
                 ii = sx-GN+i-1  !sx-GN,sx-1
                 var(ii,j,k) = bfr_r_W(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait(11)
        if(neigh(E)>=0) then
        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN 
                 ii = ex+i       !ex+1,ex+GN
                 var(ii,j,k) = bfr_r_E(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait(12)
        if(neigh(N)>=0) then
        !$acc parallel default(present) async(4)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN 
              do i = lbx,ubx
                 jj = ey+j         !ey+1,ey+GN
                 var(i,jj,k) = bfr_r_N(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait(13)
        if(neigh(S)>=0) then
        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN 
              do i = lbx,ubx
                 jj = sy-GN+j-1    !sy-GN,sy-1
                 var(i,jj,k) = bfr_r_S(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait(14)
        if(neigh(F)>=0) then
        !$acc parallel default(present) async(6)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN
           do    j = lby,uby
              do i = lbx,ubx
                 kk = ez+k         !ez+1,ez+GN
                 var(i,j,kk) = bfr_r_F(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel 
        endif

        !$acc wait(15)
        if(neigh(B)>=0) then
        !$acc parallel default(present) async(5)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN
           do    j = lby,uby
              do i = lbx,ubx
                 kk = sz-GN+k-1    !sz-GN,sz-1
                 var(i,j,kk) = bfr_r_B(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel 
        endif

        !$acc wait

        return
end subroutine mpi_communicate3D_real8_sendrecv





subroutine mpi_communicate3D_real8_sendrecv_cuda_aware(comm,neigh,dims,&
                                       bfr_s_E, bfr_s_W, bfr_r_E, bfr_r_W, &
                                       bfr_s_N, bfr_s_S, bfr_r_N, bfr_r_S, &
                                       bfr_s_B, bfr_s_F, bfr_r_B, bfr_r_F, &
                                       var)


        implicit none
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: var     
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_E
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_W
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_E
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_W
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_N
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_S
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_N
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_S
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_B
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_F
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_B
        real(rp), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_F
        integer              , dimension(6)     , intent(in)    :: neigh
        integer                                 , intent(in)    :: dims
        integer                                 , intent(in)    :: comm

        ! local declaration
        integer, parameter    :: dtype = MPI_RP
        integer, dimension(3) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, dimension(MPI_STATUS_SIZE) :: status
        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0, i,j,k, ii, jj, kk

        !
        ! === buffer dimensions
        !
        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)

        !
        ! === east - west communications
        !
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN    
                 ii = ex-GN+i  !ex+1-GN,ex
                 bfr_s_E(ii,j,k) = var(ii,j,k)
                 ii = sx-1+i   !sx,sx+GN-1
                 bfr_s_W(ii,j,k) = var(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN    
              do i = lbx,ubx
                 jj = ey-GN+j !ey+1-GN,ey
                 bfr_s_N(i,jj,k) = var(i,jj,k)
                 jj = sy-1+j  !sy,sy+GN-1
                 bfr_s_S(i,jj,k) = var(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN 
           do    j = lby,uby
              do i = lbx,ubx
                 kk = ez-GN+k  !ez+1-GN,ez
                 bfr_s_F(i,j,kk) = var(i,j,kk)
                 kk = sz-1+k   !sz,sz+GN-1
                 bfr_s_B(i,j,kk) = var(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel
        
        !$acc host_data use_device(bfr_s_E,bfr_r_W) &
        !$acc           use_device(bfr_s_W,bfr_r_E) &
        !$acc           use_device(bfr_s_S,bfr_r_N) &
        !$acc           use_device(bfr_s_N,bfr_r_S) &
        !$acc           use_device(bfr_s_B,bfr_r_F) &
        !$acc           use_device(bfr_s_F,bfr_r_B)

        !$acc wait(1)
        call MPI_SENDRECV(bfr_s_E,bfr_size_x,dtype,neigh(E),tag+1, &
                          bfr_r_W,bfr_size_x,dtype,neigh(W),tag+1, &
                          comm,status,err)

        call MPI_SENDRECV(bfr_s_W,bfr_size_x,dtype,neigh(W),tag+2, &
                          bfr_r_E,bfr_size_x,dtype,neigh(E),tag+2, &
                          comm,status,err)
        
        !$acc wait(2)
        call MPI_SENDRECV(bfr_s_S,bfr_size_y,dtype,neigh(S),tag+4, &
                          bfr_r_N,bfr_size_y,dtype,neigh(N),tag+4, &
                          comm,status,err)

        call MPI_SENDRECV(bfr_s_N,bfr_size_y,dtype,neigh(N),tag+3, &
                          bfr_r_S,bfr_size_y,dtype,neigh(S),tag+3, & 
                          comm,status,err)

        !$acc wait(3)
        call MPI_SENDRECV(bfr_s_B,bfr_size_z,dtype,neigh(B),tag+6, &
                          bfr_r_F,bfr_size_z,dtype,neigh(F),tag+6, &
                          comm,status,err)

        call MPI_SENDRECV(bfr_s_F,bfr_size_z,dtype,neigh(F),tag+5, &
                          bfr_r_B,bfr_size_z,dtype,neigh(B),tag+5, &
                          comm,status,err)
        !$acc end host_data

        if(neigh(W)>=0) then
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN 
                 ii = sx-GN+i-1  !sx-GN,sx-1
                 var(ii,j,k) = bfr_r_W(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(neigh(E)>=0) then
        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN 
                 ii = ex+i       !ex+1,ex+GN
                 var(ii,j,k) = bfr_r_E(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(neigh(S)>=0) then
        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN 
              do i = lbx,ubx
                 jj = sy-GN+j-1    !sy-GN,sy-1
                 var(i,jj,k) = bfr_r_S(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(neigh(N)>=0) then
        !$acc parallel default(present) async(4)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN 
              do i = lbx,ubx
                 jj = ey+j         !ey+1,ey+GN
                 var(i,jj,k) = bfr_r_N(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(neigh(B)>=0) then
        !$acc parallel default(present) async(5)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN
           do    j = lby,uby
              do i = lbx,ubx
                 kk = sz-GN+k-1    !sz-GN,sz-1
                 var(i,j,kk) = bfr_r_B(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel 
        endif

        if(neigh(F)>=0) then
        !$acc parallel default(present) async(6)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN
           do    j = lby,uby
              do i = lbx,ubx
                 kk = ez+k         !ez+1,ez+GN
                 var(i,j,kk) = bfr_r_F(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel 
        endif

        !$acc wait

        return
end subroutine mpi_communicate3D_real8_sendrecv_cuda_aware

















!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! -------------------------INTEGERS----------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











subroutine mpi_communicateInt3D_isend_irecv_ddt(comm,tp_send,tp_recv,neigh,dims,var)
        implicit none
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: var
        integer                , dimension(6)    , intent(in)    :: tp_send
        integer                , dimension(6)    , intent(in)    :: tp_recv
        integer                , dimension(6)    , intent(in)    :: neigh
        integer                                  , intent(in)    :: dims
        integer                                  , intent(in)    :: comm

        ! local declaration
        integer, dimension(4) :: req_xx
        integer, dimension(4) :: req_yy
        integer, dimension(4) :: req_zz

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0

        !$acc update host(var)

        !
        ! === E/W communications
        !
        call MPI_ISEND(var,1,tp_send(E),neigh(E),tag+1,comm,req_xx(1),err)
        call MPI_IRECV(var,1,tp_recv(W),neigh(W),tag+1,comm,req_xx(2),err)

        call MPI_ISEND(var,1,tp_send(W),neigh(W),tag+2,comm,req_xx(3),err)
        call MPI_IRECV(var,1,tp_recv(E),neigh(E),tag+2,comm,req_xx(4),err)

        ! wait till communications end along x cause corners
        call MPI_waitall(size(req_xx), req_xx, MPI_STATUSES_IGNORE, err)

        ! === N/S communications
        !
        call MPI_ISEND(var,1,tp_send(N),neigh(N),tag+3,comm,req_yy(1),err)
        call MPI_IRECV(var,1,tp_recv(S),neigh(S),tag+3,comm,req_yy(2),err)

        call MPI_ISEND(var,1,tp_send(S),neigh(S),tag+4,comm,req_yy(3),err)
        call MPI_IRECV(var,1,tp_recv(N),neigh(N),tag+4,comm,req_yy(4),err)

        call MPI_waitall(size(req_yy), req_yy, MPI_STATUSES_IGNORE, err)
        !
        ! === B/F communications
        !
        call MPI_ISEND(var,1,tp_send(F),neigh(F),tag+5,comm,req_zz(1),err)
        call MPI_IRECV(var,1,tp_recv(B),neigh(B),tag+5,comm,req_zz(2),err)

        call MPI_ISEND(var,1,tp_send(B),neigh(B),tag+6,comm,req_zz(3),err)
        call MPI_IRECV(var,1,tp_recv(F),neigh(F),tag+6,comm,req_zz(4),err)
         
        ! wait till communications end
        call MPI_waitall(size(req_zz), req_zz, MPI_STATUSES_IGNORE, err)

        !$acc update device(var)

        return
end subroutine mpi_communicateInt3D_isend_irecv_ddt





subroutine mpi_communicateInt3D_isend_irecv_ddt_cuda_aware(comm,tp_send,tp_recv,neigh,dims,var)
        implicit none
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: var
        integer                , dimension(6)    , intent(in)    :: tp_send
        integer                , dimension(6)    , intent(in)    :: tp_recv
        integer                , dimension(6)    , intent(in)    :: neigh
        integer                                  , intent(in)    :: dims
        integer                                  , intent(in)    :: comm

        ! local declaration
        integer, dimension(4) :: req_xx
        integer, dimension(4) :: req_yy
        integer, dimension(4) :: req_zz

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0

        !$acc host_data use_device(var)

        !
        ! === E/W communications
        !
        call MPI_ISEND(var,1,tp_send(E),neigh(E),tag+1,comm,req_xx(1),err)
        call MPI_IRECV(var,1,tp_recv(W),neigh(W),tag+1,comm,req_xx(2),err)

        call MPI_ISEND(var,1,tp_send(W),neigh(W),tag+2,comm,req_xx(3),err)
        call MPI_IRECV(var,1,tp_recv(E),neigh(E),tag+2,comm,req_xx(4),err)

        ! wait till communications end along x cause corners
        call MPI_waitall(size(req_xx), req_xx, MPI_STATUSES_IGNORE, err)

        ! === N/S communications
        !
        call MPI_ISEND(var,1,tp_send(N),neigh(N),tag+3,comm,req_yy(1),err)
        call MPI_IRECV(var,1,tp_recv(S),neigh(S),tag+3,comm,req_yy(2),err)

        call MPI_ISEND(var,1,tp_send(S),neigh(S),tag+4,comm,req_yy(3),err)
        call MPI_IRECV(var,1,tp_recv(N),neigh(N),tag+4,comm,req_yy(4),err)

        call MPI_waitall(size(req_yy), req_yy, MPI_STATUSES_IGNORE, err)
        !
        ! === B/F communications
        !
        call MPI_ISEND(var,1,tp_send(F),neigh(F),tag+5,comm,req_zz(1),err)
        call MPI_IRECV(var,1,tp_recv(B),neigh(B),tag+5,comm,req_zz(2),err)

        call MPI_ISEND(var,1,tp_send(B),neigh(B),tag+6,comm,req_zz(3),err)
        call MPI_IRECV(var,1,tp_recv(F),neigh(F),tag+6,comm,req_zz(4),err)
         
        ! wait till communications end
        call MPI_waitall(size(req_zz), req_zz, MPI_STATUSES_IGNORE, err)

        !$acc end host_data

        return
end subroutine mpi_communicateInt3D_isend_irecv_ddt_cuda_aware































subroutine mpi_communicateInt3D_sendrecv(comm,neigh,dims,&
                                       bfr_s_E, bfr_s_W, bfr_r_E, bfr_r_W, &
                                       bfr_s_N, bfr_s_S, bfr_r_N, bfr_r_S, &
                                       bfr_s_B, bfr_s_F, bfr_r_B, bfr_r_F, &
                                       var)


        implicit none
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: var     
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_E
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_W
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_E
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_W
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_N
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_S
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_N
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_S
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_B
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_F
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_B
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_F
        integer                , dimension(6)     , intent(in)    :: neigh
        integer                                   , intent(in)    :: dims
        integer                                   , intent(in)    :: comm

        ! local declaration
        integer, parameter    :: dtype = MPI_INTEGER1
        integer, dimension(3) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, dimension(MPI_STATUS_SIZE) :: status
        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0, i,j,k, ii, jj, kk

        !
        ! === buffer dimensions
        !
        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)

        !
        ! === east - west communications
        !
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN    
                 ii = ex-GN+i  !ex+1-GN,ex
                 bfr_s_E(ii,j,k) = var(ii,j,k)
                 ii = sx-1+i   !sx,sx+GN-1
                 bfr_s_W(ii,j,k) = var(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN    
              do i = lbx,ubx
                 jj = ey-GN+j !ey+1-GN,ey
                 bfr_s_N(i,jj,k) = var(i,jj,k)
                 jj = sy-1+j  !sy,sy+GN-1
                 bfr_s_S(i,jj,k) = var(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN 
           do    j = lby,uby
              do i = lbx,ubx
                 kk = ez-GN+k  !ez+1-GN,ez
                 bfr_s_F(i,j,kk) = var(i,j,kk)
                 kk = sz-1+k   !sz,sz+GN-1
                 bfr_s_B(i,j,kk) = var(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc wait(1)
        !$acc update host(bfr_s_E) async(4)
        !$acc update host(bfr_s_W) async(5)

        !$acc wait(2)
        !$acc update host(bfr_s_S) async(6)
        !$acc update host(bfr_s_N) async(7)

        !$acc wait(3)
        !$acc update host(bfr_s_B) async(8)
        !$acc update host(bfr_s_F) async(9)

        !$acc wait(4)
        call MPI_SENDRECV(bfr_s_E,bfr_size_x,dtype,neigh(E),tag+1, &
                          bfr_r_W,bfr_size_x,dtype,neigh(W),tag+1, &
                          comm,status,err)
        !$acc update device(bfr_r_W) async(10)

        !$acc wait(5)
        call MPI_SENDRECV(bfr_s_W,bfr_size_x,dtype,neigh(W),tag+2, &
                          bfr_r_E,bfr_size_x,dtype,neigh(E),tag+2, &
                          comm,status,err)
        !$acc update device(bfr_r_E) async(11)
        
        !$acc wait(6)
        call MPI_SENDRECV(bfr_s_S,bfr_size_y,dtype,neigh(S),tag+4, &
                          bfr_r_N,bfr_size_y,dtype,neigh(N),tag+4, &
                          comm,status,err)
        !$acc update device(bfr_r_N) async(12)

        !$acc wait(7)
        call MPI_SENDRECV(bfr_s_N,bfr_size_y,dtype,neigh(N),tag+3, &
                          bfr_r_S,bfr_size_y,dtype,neigh(S),tag+3, & 
                          comm,status,err)
        !$acc update device(bfr_r_S) async(13)

        !$acc wait(8)
        call MPI_SENDRECV(bfr_s_B,bfr_size_z,dtype,neigh(B),tag+6, &
                          bfr_r_F,bfr_size_z,dtype,neigh(F),tag+6, &
                          comm,status,err)
        !$acc update device(bfr_r_F) async(14)

        !$acc wait(9)
        call MPI_SENDRECV(bfr_s_F,bfr_size_z,dtype,neigh(F),tag+5, &
                          bfr_r_B,bfr_size_z,dtype,neigh(B),tag+5, &
                          comm,status,err)
        !$acc update device(bfr_r_B) async(15)


        !$acc wait(10)
        if(neigh(W)>=0) then
          !$acc parallel default(present) 
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN 
                   ii = sx-GN+i-1  !sx-GN,sx-1
                   var(ii,j,k) = bfr_r_W(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait(11)
        if(neigh(E)>=0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN 
                   ii = ex+i       !ex+1,ex+GN
                   var(ii,j,k) = bfr_r_E(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait(12)
        if(neigh(N)>=0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN 
                do i = lbx,ubx
                   jj = ey+j         !ey+1,ey+GN
                   var(i,jj,k) = bfr_r_N(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait(13)
        if(neigh(S)>=0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN 
                do i = lbx,ubx
                   jj = sy-GN+j-1    !sy-GN,sy-1
                   var(i,jj,k) = bfr_r_S(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait(14)
        if(neigh(F)>=0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = ez+k         !ez+1,ez+GN
                   var(i,j,kk) = bfr_r_F(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel 
        endif

        !$acc wait(15)
        if(neigh(B)>=0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = sz-GN+k-1    !sz-GN,sz-1
                   var(i,j,kk) = bfr_r_B(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel 
        endif


        return
end subroutine mpi_communicateInt3D_sendrecv





subroutine mpi_communicateInt3D_sendrecv_cuda_aware(comm,neigh,dims,&
                                       bfr_s_E, bfr_s_W, bfr_r_E, bfr_r_W, &
                                       bfr_s_N, bfr_s_S, bfr_r_N, bfr_r_S, &
                                       bfr_s_B, bfr_s_F, bfr_r_B, bfr_r_F, &
                                       var)


        implicit none
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: var     
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_E
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_W
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_E
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_W
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_N
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_S
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_N
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_S
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_B
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_F
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_B
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_F
        integer                , dimension(6)     , intent(in)    :: neigh
        integer                                   , intent(in)    :: dims
        integer                                   , intent(in)    :: comm

        ! local declaration
        integer, parameter    :: dtype = MPI_INTEGER1
        integer, dimension(3) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, dimension(MPI_STATUS_SIZE) :: status
        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0, i,j,k, ii, jj, kk

        !
        ! === buffer dimensions
        !
        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)

        !
        ! === east - west communications
        !
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN    
                 ii = ex-GN+i  !ex+1-GN,ex
                 bfr_s_E(ii,j,k) = var(ii,j,k)
                 ii = sx-1+i   !sx,sx+GN-1
                 bfr_s_W(ii,j,k) = var(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN    
              do i = lbx,ubx
                 jj = ey-GN+j !ey+1-GN,ey
                 bfr_s_N(i,jj,k) = var(i,jj,k)
                 jj = sy-1+j  !sy,sy+GN-1
                 bfr_s_S(i,jj,k) = var(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN 
           do    j = lby,uby
              do i = lbx,ubx
                 kk = ez-GN+k  !ez+1-GN,ez
                 bfr_s_F(i,j,kk) = var(i,j,kk)
                 kk = sz-1+k   !sz,sz+GN-1
                 bfr_s_B(i,j,kk) = var(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc host_data use_device(bfr_s_E,bfr_r_W) &
        !$acc           use_device(bfr_s_W,bfr_r_E) &
        !$acc           use_device(bfr_s_S,bfr_r_N) &
        !$acc           use_device(bfr_s_N,bfr_r_S) &
        !$acc           use_device(bfr_s_B,bfr_r_F) &
        !$acc           use_device(bfr_s_F,bfr_r_B)

        !$acc wait(1)
        call MPI_SENDRECV(bfr_s_E,bfr_size_x,dtype,neigh(E),tag+1, &
                          bfr_r_W,bfr_size_x,dtype,neigh(W),tag+1, &
                          comm,status,err)

        call MPI_SENDRECV(bfr_s_W,bfr_size_x,dtype,neigh(W),tag+2, &
                          bfr_r_E,bfr_size_x,dtype,neigh(E),tag+2, &
                          comm,status,err)
        
        !$acc wait(2)
        call MPI_SENDRECV(bfr_s_S,bfr_size_y,dtype,neigh(S),tag+4, &
                          bfr_r_N,bfr_size_y,dtype,neigh(N),tag+4, &
                          comm,status,err)

        call MPI_SENDRECV(bfr_s_N,bfr_size_y,dtype,neigh(N),tag+3, &
                          bfr_r_S,bfr_size_y,dtype,neigh(S),tag+3, & 
                          comm,status,err)

        !$acc wait(3)
        call MPI_SENDRECV(bfr_s_B,bfr_size_z,dtype,neigh(B),tag+6, &
                          bfr_r_F,bfr_size_z,dtype,neigh(F),tag+6, &
                          comm,status,err)

        call MPI_SENDRECV(bfr_s_F,bfr_size_z,dtype,neigh(F),tag+5, &
                          bfr_r_B,bfr_size_z,dtype,neigh(B),tag+5, &
                          comm,status,err)

        !$acc end host_data


        if(neigh(W)>=0) then
          !$acc parallel default(present) async(1)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN 
                   ii = sx-GN+i-1  !sx-GN,sx-1
                   var(ii,j,k) = bfr_r_W(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(E)>=0) then
          !$acc parallel default(present) async(2)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN 
                   ii = ex+i       !ex+1,ex+GN
                   var(ii,j,k) = bfr_r_E(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(N)>=0) then
          !$acc parallel default(present) async(3)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN 
                do i = lbx,ubx
                   jj = ey+j         !ey+1,ey+GN
                   var(i,jj,k) = bfr_r_N(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(S)>=0) then
          !$acc parallel default(present) async(4)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN 
                do i = lbx,ubx
                   jj = sy-GN+j-1    !sy-GN,sy-1
                   var(i,jj,k) = bfr_r_S(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(F)>=0) then
          !$acc parallel default(present) async(5)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = ez+k         !ez+1,ez+GN
                   var(i,j,kk) = bfr_r_F(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel 
        endif

        if(neigh(B)>=0) then
          !$acc parallel default(present) async(6)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = sz-GN+k-1    !sz-GN,sz-1
                   var(i,j,kk) = bfr_r_B(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel 
        endif

        if(neigh(W)>=0) then
          !$acc wait(1)
        endif
        if(neigh(E)>=0) then
          !$acc wait(2)
        endif
        if(neigh(N)>=0) then
          !$acc wait(3)
        endif
        if(neigh(S)>=0) then
          !$acc wait(4)
        endif
        if(neigh(F)>=0) then
          !$acc wait(5)
        endif
        if(neigh(B)>=0) then
          !$acc wait(6)
        endif

        return
end subroutine mpi_communicateInt3D_sendrecv_cuda_aware


















subroutine mpi_communicateInt3D_isend_irecv_bfr(comm,neigh,dims,&
                                       bfr_s_E, bfr_s_W, bfr_r_E, bfr_r_W, &
                                       bfr_s_N, bfr_s_S, bfr_r_N, bfr_r_S, &
                                       bfr_s_B, bfr_s_F, bfr_r_B, bfr_r_F, &
                                       var)

        implicit none
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: var     
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_E
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_W
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_E
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_W
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_N
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_S
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_N
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_S
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_B
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_F
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_B
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_F
        integer                , dimension(6)     , intent(in)    :: neigh
        integer                                   , intent(in)    :: dims
        integer                                   , intent(in)    :: comm

        ! local declaration
        integer, parameter    :: dtype = MPI_INTEGER1
        integer, dimension(3) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer, dimension(4) :: req_xx
        integer, dimension(4) :: req_yy
        integer, dimension(4) :: req_zz
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: i,j,k,err = 0, ii, jj, kk

        !
        ! === buffer dimensions
        !
        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)


        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN    
                 ii = ex-GN+i  !ex+1-GN,ex
                 bfr_s_E(ii,j,k) = var(ii,j,k)
                 ii = sx-1+i   !sx,sx+GN-1
                 bfr_s_W(ii,j,k) = var(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN    
              do i = lbx,ubx
                 jj = ey-GN+j !ey+1-GN,ey
                 bfr_s_N(i,jj,k) = var(i,jj,k)
                 jj = sy-1+j  !sy,sy+GN-1
                 bfr_s_S(i,jj,k) = var(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN 
           do    j = lby,uby
              do i = lbx,ubx
                 kk = ez-GN+k  !ez+1-GN,ez
                 bfr_s_F(i,j,kk) = var(i,j,kk)
                 kk = sz-1+k   !sz,sz+GN-1
                 bfr_s_B(i,j,kk) = var(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc wait(1)
        !$acc update host(bfr_s_E) async(4)
        !$acc update host(bfr_s_W) async(5)

        !$acc wait(2)
        !$acc update host(bfr_s_N) async(6)
        !$acc update host(bfr_s_S) async(7)

        !$acc wait(3)
        !$acc update host(bfr_s_F) async(8)
        !$acc update host(bfr_s_B) async(9)

        !
        ! === east - west communications
        !
        call MPI_IRECV(bfr_r_W,bfr_size_x,dtype,neigh(W),tag+1,comm,req_xx(2),err)
        call MPI_IRECV(bfr_r_E,bfr_size_x,dtype,neigh(E),tag+2,comm,req_xx(4),err)

        !$acc wait(4)
        call MPI_ISEND(bfr_s_E,bfr_size_x,dtype,neigh(E),tag+1,comm,req_xx(1),err)
        !$acc wait(5)
        call MPI_ISEND(bfr_s_W,bfr_size_x,dtype,neigh(W),tag+2,comm,req_xx(3),err)

        ! wait till communications end (corner needed)
        call MPI_waitall(size(req_xx), req_xx, MPI_STATUSES_IGNORE, err)
        !$acc update device(bfr_r_W,bfr_r_E) async(10)

        !
        ! ==== nord-south communications
        !
        call MPI_IRECV(bfr_r_S,bfr_size_y,dtype,neigh(S),tag+3,comm,req_yy(2),err)
        call MPI_IRECV(bfr_r_N,bfr_size_y,dtype,neigh(N),tag+4,comm,req_yy(4),err)

        !$acc wait(6)
        call MPI_ISEND(bfr_s_N,bfr_size_y,dtype,neigh(N),tag+3,comm,req_yy(1),err)
        !$acc wait(7)
        call MPI_ISEND(bfr_s_S,bfr_size_y,dtype,neigh(S),tag+4,comm,req_yy(3),err)

        call MPI_waitall(size(req_yy), req_yy, MPI_STATUSES_IGNORE, err)
        !$acc update device(bfr_r_N,bfr_r_S) async(11)

        !
        ! ==== backward-forwvar communications
        !
        call MPI_IRECV(bfr_r_B,bfr_size_z,dtype,neigh(B),tag+5,comm,req_zz(2),err)
        call MPI_IRECV(bfr_r_F,bfr_size_z,dtype,neigh(F),tag+6,comm,req_zz(4),err)

        !$acc wait(8)
        call MPI_ISEND(bfr_s_F,bfr_size_z,dtype,neigh(F),tag+5,comm,req_zz(1),err)
        !$acc wait(9)
        call MPI_ISEND(bfr_s_B,bfr_size_z,dtype,neigh(B),tag+6,comm,req_zz(3),err)

        call MPI_waitall(size(req_zz), req_zz, MPI_STATUSES_IGNORE, err)
        !$acc update device(bfr_r_F,bfr_r_B) async(12)

        !$acc wait(10)
        if(neigh(W)>=0) then
          !$acc parallel default(present) 
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN 
                   ii = sx-GN+i-1  !sx-GN,sx-1
                   var(ii,j,k) = bfr_r_W(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(E)>=0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN 
                   ii = ex+i       !ex+1,ex+GN
                   var(ii,j,k) = bfr_r_E(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait(11)
        if(neigh(N)>=0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN 
                do i = lbx,ubx
                   jj = ey+j         !ey+1,ey+GN
                   var(i,jj,k) = bfr_r_N(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(S)>=0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN 
                do i = lbx,ubx
                   jj = sy-GN+j-1    !sy-GN,sy-1
                   var(i,jj,k) = bfr_r_S(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait(12)
        if(neigh(F)>=0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = ez+k         !ez+1,ez+GN
                   var(i,j,kk) = bfr_r_F(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel 
        endif

        if(neigh(B)>=0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = sz-GN+k-1    !sz-GN,sz-1
                   var(i,j,kk) = bfr_r_B(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel 
        endif


        return
end subroutine mpi_communicateInt3D_isend_irecv_bfr








subroutine mpi_communicateInt3D_isend_irecv_bfr_cuda_aware(comm,neigh,dims,&
                                       bfr_s_E, bfr_s_W, bfr_r_E, bfr_r_W, &
                                       bfr_s_N, bfr_s_S, bfr_r_N, bfr_r_S, &
                                       bfr_s_B, bfr_s_F, bfr_r_B, bfr_r_F, &
                                       var)

        implicit none
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: var     
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_E
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_W
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_E
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_W
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_N
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_S
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_N
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_S
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_B
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_s_F
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_B
        integer(1), allocatable, dimension(:,:,:) , intent(inout) :: bfr_r_F
        integer                , dimension(6)     , intent(in)    :: neigh
        integer                                   , intent(in)    :: dims
        integer                                   , intent(in)    :: comm

        ! local declaration
        integer, parameter    :: dtype = MPI_INTEGER1
        integer, dimension(3) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer, dimension(4) :: req_xx
        integer, dimension(4) :: req_yy
        integer, dimension(4) :: req_zz
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: i,j,k,err = 0, ii, jj, kk

        !
        ! === buffer dimensions
        !
        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)


        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN    
                 ii = ex-GN+i  !ex+1-GN,ex
                 bfr_s_E(ii,j,k) = var(ii,j,k)
                 ii = sx-1+i   !sx,sx+GN-1
                 bfr_s_W(ii,j,k) = var(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN    
              do i = lbx,ubx
                 jj = ey-GN+j !ey+1-GN,ey
                 bfr_s_N(i,jj,k) = var(i,jj,k)
                 jj = sy-1+j  !sy,sy+GN-1
                 bfr_s_S(i,jj,k) = var(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN 
           do    j = lby,uby
              do i = lbx,ubx
                 kk = ez-GN+k  !ez+1-GN,ez
                 bfr_s_F(i,j,kk) = var(i,j,kk)
                 kk = sz-1+k   !sz,sz+GN-1
                 bfr_s_B(i,j,kk) = var(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc host_data use_device(bfr_s_E,bfr_r_W) &
        !$acc           use_device(bfr_s_W,bfr_r_E) &
        !$acc           use_device(bfr_s_S,bfr_r_N) &
        !$acc           use_device(bfr_s_N,bfr_r_S) &
        !$acc           use_device(bfr_s_B,bfr_r_F) &
        !$acc           use_device(bfr_s_F,bfr_r_B)

        !
        ! === east - west communications
        !
        !$acc wait(1)
        call MPI_IRECV(bfr_r_W,bfr_size_x,dtype,neigh(W),tag+1,comm,req_xx(2),err)
        call MPI_IRECV(bfr_r_E,bfr_size_x,dtype,neigh(E),tag+2,comm,req_xx(4),err)

        call MPI_ISEND(bfr_s_E,bfr_size_x,dtype,neigh(E),tag+1,comm,req_xx(1),err)
        call MPI_ISEND(bfr_s_W,bfr_size_x,dtype,neigh(W),tag+2,comm,req_xx(3),err)

        ! wait till communications end (corner needed)
        call MPI_waitall(size(req_xx), req_xx, MPI_STATUSES_IGNORE, err)

        !
        ! ==== nord-south communications
        !
        !$acc wait(2)
        call MPI_IRECV(bfr_r_S,bfr_size_y,dtype,neigh(S),tag+3,comm,req_yy(2),err)
        call MPI_IRECV(bfr_r_N,bfr_size_y,dtype,neigh(N),tag+4,comm,req_yy(4),err)

        call MPI_ISEND(bfr_s_N,bfr_size_y,dtype,neigh(N),tag+3,comm,req_yy(1),err)
        call MPI_ISEND(bfr_s_S,bfr_size_y,dtype,neigh(S),tag+4,comm,req_yy(3),err)

        call MPI_waitall(size(req_yy), req_yy, MPI_STATUSES_IGNORE, err)

        !
        ! ==== backward-forwvar communications
        !
        !$acc wait(3)
        call MPI_IRECV(bfr_r_B,bfr_size_z,dtype,neigh(B),tag+5,comm,req_zz(2),err)
        call MPI_IRECV(bfr_r_F,bfr_size_z,dtype,neigh(F),tag+6,comm,req_zz(4),err)

        call MPI_ISEND(bfr_s_F,bfr_size_z,dtype,neigh(F),tag+5,comm,req_zz(1),err)
        call MPI_ISEND(bfr_s_B,bfr_size_z,dtype,neigh(B),tag+6,comm,req_zz(3),err)

        call MPI_waitall(size(req_zz), req_zz, MPI_STATUSES_IGNORE, err)

        !$acc end host_data

        if(neigh(W)>=0) then
          !$acc parallel default(present) async(1)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN 
                   ii = sx-GN+i-1  !sx-GN,sx-1
                   var(ii,j,k) = bfr_r_W(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(E)>=0) then
          !$acc parallel default(present) async(2)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN 
                   ii = ex+i       !ex+1,ex+GN
                   var(ii,j,k) = bfr_r_E(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(N)>=0) then
          !$acc parallel default(present) async(3)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN 
                do i = lbx,ubx
                   jj = ey+j         !ey+1,ey+GN
                   var(i,jj,k) = bfr_r_N(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(S)>=0) then
          !$acc parallel default(present) async(4)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN 
                do i = lbx,ubx
                   jj = sy-GN+j-1    !sy-GN,sy-1
                   var(i,jj,k) = bfr_r_S(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(F)>=0) then
          !$acc parallel default(present) async(5)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = ez+k         !ez+1,ez+GN
                   var(i,j,kk) = bfr_r_F(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel 
        endif

        if(neigh(B)>=0) then
          !$acc parallel default(present) async(6)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = sz-GN+k-1    !sz-GN,sz-1
                   var(i,j,kk) = bfr_r_B(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel 
        endif


        if(neigh(W)>=0) then
          !$acc wait(1)
        endif
        if(neigh(E)>=0) then
          !$acc wait(2)
        endif
        if(neigh(N)>=0) then
          !$acc wait(3)
        endif
        if(neigh(S)>=0) then
          !$acc wait(4)
        endif
        if(neigh(F)>=0) then
          !$acc wait(5)
        endif
        if(neigh(B)>=0) then
          !$acc wait(6)
        endif


        return
end subroutine mpi_communicateInt3D_isend_irecv_bfr_cuda_aware























































end module mpi_comm_module
