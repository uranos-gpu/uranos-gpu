module mpi_comm_module
use parameters_module, only: rp, MPI_RP, mpi_opt_level, cuda_aware
use profiling_module,  only: StartProfRange, EndProfRange
use mpi_module
use mpi

implicit none
private
public mpi_share, mpi_share_int1, mpi_share_planes

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
        case(1)

          if(cuda_aware) then
          call mpi_communicate4D_real8_sendrecv_cuda_aware(mpi_comm_cart,&
                my_neighbour,&
                phi_bfr_send_E, phi_bfr_send_W, phi_bfr_recv_E, phi_bfr_recv_W, &
                phi_bfr_send_N, phi_bfr_send_S, phi_bfr_recv_N, phi_bfr_recv_S, &
                phi_bfr_send_B, phi_bfr_send_F, phi_bfr_recv_B, phi_bfr_recv_F, &
                phi)

          else
          call mpi_communicate4D_real8_sendrecv(mpi_comm_cart,my_neighbour,&
                phi_bfr_send_E, phi_bfr_send_W, phi_bfr_recv_E, phi_bfr_recv_W, &
                phi_bfr_send_N, phi_bfr_send_S, phi_bfr_recv_N, phi_bfr_recv_S, &
                phi_bfr_send_B, phi_bfr_send_F, phi_bfr_recv_B, phi_bfr_recv_F, &
                phi)
          endif

        case(2)

          if(cuda_aware) then
          call mpi_communicate4D_real8_bfr_cuda_aware(mpi_comm_cart,&
                my_neighbour, &
                phi_bfr_send_E, phi_bfr_send_W, phi_bfr_recv_E, phi_bfr_recv_W, &
                phi_bfr_send_N, phi_bfr_send_S, phi_bfr_recv_N, phi_bfr_recv_S, &
                phi_bfr_send_B, phi_bfr_send_F, phi_bfr_recv_B, phi_bfr_recv_F, &
                phi)
          else
          call mpi_communicate4D_real8_bfr(mpi_comm_cart,my_neighbour, &
                phi_bfr_send_E, phi_bfr_send_W, phi_bfr_recv_E, phi_bfr_recv_W, &
                phi_bfr_send_N, phi_bfr_send_S, phi_bfr_recv_N, phi_bfr_recv_S, &
                phi_bfr_send_B, phi_bfr_send_F, phi_bfr_recv_B, phi_bfr_recv_F, &
                phi)
          endif

        case(3)

          if(cuda_aware) then
          call mpi_communicate4D_real8_cuda_aware(mpi_comm_cart,&
                     type_send_cons,type_recv_cons,my_neighbour,dims,phi)
          else
          call mpi_communicate4D_real8(mpi_comm_cart,&
                     type_send_cons,type_recv_cons,my_neighbour,dims,phi)
          endif

        case(4)
         call mpi_communicate4D_real8_bfr_one_sided(mpi_comm_cart,&
               my_neighbour, &
               phi_bfr_send_E, phi_bfr_send_W, phi_bfr_recv_E, phi_bfr_recv_W, &
               phi_bfr_send_N, phi_bfr_send_S, phi_bfr_recv_N, phi_bfr_recv_S, &
               phi_bfr_send_B, phi_bfr_send_F, phi_bfr_recv_B, phi_bfr_recv_F, &
               phi)
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
        case(1)

          if(cuda_aware) then
          call mpi_communicate3D_real8_sendrecv_cuda_aware(mpi_comm_cart,&
                my_neighbour,&
                bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                var)

          else
          call mpi_communicate3D_real8_sendrecv(mpi_comm_cart,my_neighbour,&
                bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                var)
          endif

        case(2)

          if(cuda_aware) then
          call mpi_communicate3D_real8_bfr_cuda_aware(mpi_comm_cart, &
                my_neighbour, &
                bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                var)
          else
          call mpi_communicate3D_real8_bfr(mpi_comm_cart,my_neighbour, &
                bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                var)
          endif

        case(3)
          if(cuda_aware) then
          call mpi_communicate3D_real8_cuda_aware(mpi_comm_cart,&
                     type_send_prim,type_recv_prim,my_neighbour,dims,var)
          else
          call mpi_communicate3D_real8(mpi_comm_cart,&
                     type_send_prim,type_recv_prim,my_neighbour,dims,var)
          endif

        case(4)
         call mpi_communicate3D_real8_bfr_one_sided(mpi_comm_cart, &
               my_neighbour, &
               bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
               bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
               bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
               var)
        case default

          print*, 'mpi_opt_level equal to ', mpi_opt_level, ' is not implemented'
          stop

        endselect

        return
end subroutine mpi_share3D

subroutine mpi_shareInt3D(mpi_comm_cart,type_send_prim,type_recv_prim,my_neighbour, &
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
        integer                                  , intent(in)    :: mpi_comm_cart

        selectcase(mpi_opt_level)
        case(1)
        if(cuda_aware) then
        call mpi_communicateInt3D_sendrecv_cuda_aware(mpi_comm_cart,&
                        my_neighbour,&
                        bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                        bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                        bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                        var)
        else
        call mpi_communicateInt3D_sendrecv(mpi_comm_cart,my_neighbour,&
                        bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                        bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                        bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                        var)
        endif

        case(2)
        if(cuda_aware) then
        call mpi_communicateInt3D_isend_irecv_bfr_cuda_aware(mpi_comm_cart,&
                        my_neighbour ,&
                        bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                        bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                        bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                        var)
        else
        call mpi_communicateInt3D_isend_irecv_bfr(mpi_comm_cart,my_neighbour,&
                        bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                        bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                        bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                        var)
        endif

        case(3)
         if(cuda_aware) then
         call mpi_communicateInt3D_isend_irecv_ddt_cuda_aware(mpi_comm_cart,&
                    type_send_prim,type_recv_prim,my_neighbour,var)
         else
         call mpi_communicateInt3D_isend_irecv_ddt(mpi_comm_cart,&
                    type_send_prim,type_recv_prim,my_neighbour,var)
         endif
        case(4)
         call mpi_communicateInt3D_bfr_one_sided(mpi_comm_cart,&
                        my_neighbour ,&
                        bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W, &
                        bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S, &
                        bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F, &
                        var)
        case default

          print*, 'mpi_opt_level equal to ', mpi_opt_level, ' is not implemented'
          stop

        endselect

        return
end subroutine mpi_shareInt3D

subroutine mpi_share_planes(comm,neigh,dims, &
                pl_sl_x, pl_sh_x, pl_lo_x, pl_hi_x, &
                pl_sl_y, pl_sh_y, pl_lo_y, pl_hi_y, &
                pl_sl_z, pl_sh_z, pl_lo_z, pl_hi_z)
        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_sl_x, pl_sh_x
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_lo_x, pl_hi_x
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_sl_y, pl_sh_y
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_lo_y, pl_hi_y
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_sl_z, pl_sh_z
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_lo_z, pl_hi_z
        integer , dimension(6)                 , intent(in)    :: neigh
        integer                                , intent(in)    :: dims
        integer                                , intent(in)    :: comm

        selectcase(mpi_opt_level)
        case(1)

          if(cuda_aware) then
          call mpi_planes_sendrecv_cuda_aware(comm,neigh,dims, &
                pl_sl_x, pl_sh_x, pl_lo_x, pl_hi_x, &
                pl_sl_y, pl_sh_y, pl_lo_y, pl_hi_y, &
                pl_sl_z, pl_sh_z, pl_lo_z, pl_hi_z)
          else
          call mpi_planes_sendrecv(comm,neigh,dims, &
                pl_sl_x, pl_sh_x, pl_lo_x, pl_hi_x, &
                pl_sl_y, pl_sh_y, pl_lo_y, pl_hi_y, &
                pl_sl_z, pl_sh_z, pl_lo_z, pl_hi_z)
          endif

        case(2,3,4)

          if(cuda_aware) then
          call mpi_planes_isend_irecv_cuda_aware(comm,neigh,dims, &
                pl_sl_x, pl_sh_x, pl_lo_x, pl_hi_x, &
                pl_sl_y, pl_sh_y, pl_lo_y, pl_hi_y, &
                pl_sl_z, pl_sh_z, pl_lo_z, pl_hi_z)
          else
          call mpi_planes_isend_irecv(comm,neigh,dims, &
                pl_sl_x, pl_sh_x, pl_lo_x, pl_hi_x, &
                pl_sl_y, pl_sh_y, pl_lo_y, pl_hi_y, &
                pl_sl_z, pl_sh_z, pl_lo_z, pl_hi_z)
          endif

        case default

          print*, 'mpi_opt_level equal to ', mpi_opt_level, ' is not implemented'
          stop

        endselect

        return
end subroutine mpi_share_planes

subroutine mpi_planes_sendrecv(comm,neigh,dims, &
                pl_sl_x, pl_sh_x, pl_lo_x, pl_hi_x, &
                pl_sl_y, pl_sh_y, pl_lo_y, pl_hi_y, &
                pl_sl_z, pl_sh_z, pl_lo_z, pl_hi_z)
        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_sl_x, pl_sh_x
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_lo_x, pl_hi_x
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_sl_y, pl_sh_y
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_lo_y, pl_hi_y
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_sl_z, pl_sh_z
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_lo_z, pl_hi_z
        integer , dimension(6)                 , intent(in)    :: neigh
        integer                                , intent(in)    :: dims
        integer                                , intent(in)    :: comm

        integer, dimension(MPI_STATUS_SIZE) :: istat
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: np_x, np_y, np_z, err

        np_x = size(pl_sl_x)
        np_y = size(pl_sl_y)
        np_z = size(pl_sl_z)

        !$acc update host(pl_sl_x,pl_sh_x,pl_lo_x,pl_hi_x, &
        !$acc             pl_sl_y,pl_sh_y,pl_lo_y,pl_hi_y, &
        !$acc             pl_sl_z,pl_sh_z,pl_lo_z,pl_hi_z)

        call MPI_SENDRECV(pl_sh_x, np_x, MPI_RP, neigh(E), 41, &
                          pl_lo_x, np_x, MPI_RP, neigh(W), 41, comm, istat, err)
        call MPI_SENDRECV(pl_sl_x, np_x, MPI_RP, neigh(W), 42, &
                          pl_hi_x, np_x, MPI_RP, neigh(E), 42, comm, istat, err)
        call MPI_SENDRECV(pl_sh_y, np_y, MPI_RP, neigh(N), 43, &
                          pl_lo_y, np_y, MPI_RP, neigh(S), 43, comm, istat, err)
        call MPI_SENDRECV(pl_sl_y, np_y, MPI_RP, neigh(S), 44, &
                          pl_hi_y, np_y, MPI_RP, neigh(N), 44, comm, istat, err)
        if(dims == 3) then
          call MPI_SENDRECV(pl_sh_z, np_z, MPI_RP, neigh(F), 45, &
                            pl_lo_z, np_z, MPI_RP, neigh(B), 45, comm, istat, err)
          call MPI_SENDRECV(pl_sl_z, np_z, MPI_RP, neigh(B), 46, &
                            pl_hi_z, np_z, MPI_RP, neigh(F), 46, comm, istat, err)
        endif

        !$acc update device(pl_lo_x,pl_hi_x,pl_lo_y,pl_hi_y,pl_lo_z,pl_hi_z)

        return
end subroutine mpi_planes_sendrecv

subroutine mpi_planes_sendrecv_cuda_aware(comm,neigh,dims, &
                pl_sl_x, pl_sh_x, pl_lo_x, pl_hi_x, &
                pl_sl_y, pl_sh_y, pl_lo_y, pl_hi_y, &
                pl_sl_z, pl_sh_z, pl_lo_z, pl_hi_z)
        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_sl_x, pl_sh_x
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_lo_x, pl_hi_x
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_sl_y, pl_sh_y
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_lo_y, pl_hi_y
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_sl_z, pl_sh_z
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_lo_z, pl_hi_z
        integer , dimension(6)                 , intent(in)    :: neigh
        integer                                , intent(in)    :: dims
        integer                                , intent(in)    :: comm

        integer, dimension(MPI_STATUS_SIZE) :: istat
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: np_x, np_y, np_z, err

        np_x = size(pl_sl_x)
        np_y = size(pl_sl_y)
        np_z = size(pl_sl_z)

        !$acc host_data use_device(pl_sl_x,pl_sh_x,pl_lo_x,pl_hi_x, &
        !$acc                      pl_sl_y,pl_sh_y,pl_lo_y,pl_hi_y, &
        !$acc                      pl_sl_z,pl_sh_z,pl_lo_z,pl_hi_z)

        call MPI_SENDRECV(pl_sh_x, np_x, MPI_RP, neigh(E), 41, &
                          pl_lo_x, np_x, MPI_RP, neigh(W), 41, comm, istat, err)
        call MPI_SENDRECV(pl_sl_x, np_x, MPI_RP, neigh(W), 42, &
                          pl_hi_x, np_x, MPI_RP, neigh(E), 42, comm, istat, err)
        call MPI_SENDRECV(pl_sh_y, np_y, MPI_RP, neigh(N), 43, &
                          pl_lo_y, np_y, MPI_RP, neigh(S), 43, comm, istat, err)
        call MPI_SENDRECV(pl_sl_y, np_y, MPI_RP, neigh(S), 44, &
                          pl_hi_y, np_y, MPI_RP, neigh(N), 44, comm, istat, err)
        if(dims == 3) then
          call MPI_SENDRECV(pl_sh_z, np_z, MPI_RP, neigh(F), 45, &
                            pl_lo_z, np_z, MPI_RP, neigh(B), 45, comm, istat, err)
          call MPI_SENDRECV(pl_sl_z, np_z, MPI_RP, neigh(B), 46, &
                            pl_hi_z, np_z, MPI_RP, neigh(F), 46, comm, istat, err)
        endif

        !$acc end host_data

        return
end subroutine mpi_planes_sendrecv_cuda_aware

subroutine mpi_planes_isend_irecv(comm,neigh,dims, &
                pl_sl_x, pl_sh_x, pl_lo_x, pl_hi_x, &
                pl_sl_y, pl_sh_y, pl_lo_y, pl_hi_y, &
                pl_sl_z, pl_sh_z, pl_lo_z, pl_hi_z)
        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_sl_x, pl_sh_x
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_lo_x, pl_hi_x
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_sl_y, pl_sh_y
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_lo_y, pl_hi_y
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_sl_z, pl_sh_z
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_lo_z, pl_hi_z
        integer , dimension(6)                 , intent(in)    :: neigh
        integer                                , intent(in)    :: dims
        integer                                , intent(in)    :: comm

        integer, dimension(12) :: req
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: np_x, np_y, np_z, nrq, err

        np_x = size(pl_sl_x)
        np_y = size(pl_sl_y)
        np_z = size(pl_sl_z)

        !$acc update host(pl_sl_x,pl_sh_x,pl_lo_x,pl_hi_x, &
        !$acc             pl_sl_y,pl_sh_y,pl_lo_y,pl_hi_y, &
        !$acc             pl_sl_z,pl_sh_z,pl_lo_z,pl_hi_z)

        call MPI_IRECV(pl_lo_x, np_x, MPI_RP, neigh(W), 41, comm, req(1), err)
        call MPI_IRECV(pl_hi_x, np_x, MPI_RP, neigh(E), 42, comm, req(2), err)
        call MPI_IRECV(pl_lo_y, np_y, MPI_RP, neigh(S), 43, comm, req(3), err)
        call MPI_IRECV(pl_hi_y, np_y, MPI_RP, neigh(N), 44, comm, req(4), err)
        call MPI_ISEND(pl_sh_x, np_x, MPI_RP, neigh(E), 41, comm, req(5), err)
        call MPI_ISEND(pl_sl_x, np_x, MPI_RP, neigh(W), 42, comm, req(6), err)
        call MPI_ISEND(pl_sh_y, np_y, MPI_RP, neigh(N), 43, comm, req(7), err)
        call MPI_ISEND(pl_sl_y, np_y, MPI_RP, neigh(S), 44, comm, req(8), err)
        nrq = 8
        if(dims == 3) then
          call MPI_IRECV(pl_lo_z, np_z, MPI_RP, neigh(B), 45, comm, req( 9), err)
          call MPI_IRECV(pl_hi_z, np_z, MPI_RP, neigh(F), 46, comm, req(10), err)
          call MPI_ISEND(pl_sh_z, np_z, MPI_RP, neigh(F), 45, comm, req(11), err)
          call MPI_ISEND(pl_sl_z, np_z, MPI_RP, neigh(B), 46, comm, req(12), err)
          nrq = 12
        endif

        call MPI_waitall(nrq, req(1:nrq), MPI_STATUSES_IGNORE, err)

        !$acc update device(pl_lo_x,pl_hi_x,pl_lo_y,pl_hi_y,pl_lo_z,pl_hi_z)

        return
end subroutine mpi_planes_isend_irecv

subroutine mpi_planes_isend_irecv_cuda_aware(comm,neigh,dims, &
                pl_sl_x, pl_sh_x, pl_lo_x, pl_hi_x, &
                pl_sl_y, pl_sh_y, pl_lo_y, pl_hi_y, &
                pl_sl_z, pl_sh_z, pl_lo_z, pl_hi_z)
        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_sl_x, pl_sh_x
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_lo_x, pl_hi_x
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_sl_y, pl_sh_y
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_lo_y, pl_hi_y
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_sl_z, pl_sh_z
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: pl_lo_z, pl_hi_z
        integer , dimension(6)                 , intent(in)    :: neigh
        integer                                , intent(in)    :: dims
        integer                                , intent(in)    :: comm

        integer, dimension(12) :: req
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: np_x, np_y, np_z, nrq, err

        np_x = size(pl_sl_x)
        np_y = size(pl_sl_y)
        np_z = size(pl_sl_z)

        !$acc host_data use_device(pl_sl_x,pl_sh_x,pl_lo_x,pl_hi_x, &
        !$acc                      pl_sl_y,pl_sh_y,pl_lo_y,pl_hi_y, &
        !$acc                      pl_sl_z,pl_sh_z,pl_lo_z,pl_hi_z)

        call MPI_IRECV(pl_lo_x, np_x, MPI_RP, neigh(W), 41, comm, req(1), err)
        call MPI_IRECV(pl_hi_x, np_x, MPI_RP, neigh(E), 42, comm, req(2), err)
        call MPI_IRECV(pl_lo_y, np_y, MPI_RP, neigh(S), 43, comm, req(3), err)
        call MPI_IRECV(pl_hi_y, np_y, MPI_RP, neigh(N), 44, comm, req(4), err)
        call MPI_ISEND(pl_sh_x, np_x, MPI_RP, neigh(E), 41, comm, req(5), err)
        call MPI_ISEND(pl_sl_x, np_x, MPI_RP, neigh(W), 42, comm, req(6), err)
        call MPI_ISEND(pl_sh_y, np_y, MPI_RP, neigh(N), 43, comm, req(7), err)
        call MPI_ISEND(pl_sl_y, np_y, MPI_RP, neigh(S), 44, comm, req(8), err)
        nrq = 8
        if(dims == 3) then
          call MPI_IRECV(pl_lo_z, np_z, MPI_RP, neigh(B), 45, comm, req( 9), err)
          call MPI_IRECV(pl_hi_z, np_z, MPI_RP, neigh(F), 46, comm, req(10), err)
          call MPI_ISEND(pl_sh_z, np_z, MPI_RP, neigh(F), 45, comm, req(11), err)
          call MPI_ISEND(pl_sl_z, np_z, MPI_RP, neigh(B), 46, comm, req(12), err)
          nrq = 12
        endif

        call MPI_waitall(nrq, req(1:nrq), MPI_STATUSES_IGNORE, err)

        !$acc end host_data

        return
end subroutine mpi_planes_isend_irecv_cuda_aware

subroutine mpi_communicate3D_real8(comm,tp_send,tp_recv,neigh,dims,var)

        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: var
        integer              , dimension(6)    , intent(in)    :: tp_send
        integer              , dimension(6)    , intent(in)    :: tp_recv
        integer              , dimension(6)    , intent(in)    :: neigh
        integer                                , intent(in)    :: dims
        integer                                , intent(in)    :: comm

        integer, dimension(4)       :: req_x
        integer, dimension(4)       :: req_y
        integer, dimension(4)       :: req_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0

        !$acc update host(var)
        call MPI_IRECV(var,1,tp_recv(W),neigh(W),tag+1,comm,req_x(1),err)
        call MPI_IRECV(var,1,tp_recv(E),neigh(E),tag+2,comm,req_x(2),err)

        call MPI_ISEND(var,1,tp_send(E),neigh(E),tag+1,comm,req_x(3),err)
        call MPI_ISEND(var,1,tp_send(W),neigh(W),tag+2,comm,req_x(4),err)

        call MPI_waitall(size(req_x), req_x, MPI_STATUSES_IGNORE, err)

        call MPI_IRECV(var,1,tp_recv(S),neigh(S),tag+3,comm,req_y(1),err)
        call MPI_IRECV(var,1,tp_recv(N),neigh(N),tag+4,comm,req_y(2),err)

        call MPI_ISEND(var,1,tp_send(N),neigh(N),tag+3,comm,req_y(3),err)
        call MPI_ISEND(var,1,tp_send(S),neigh(S),tag+4,comm,req_y(4),err)

        call MPI_waitall(size(req_y), req_y, MPI_STATUSES_IGNORE, err)

        if(dims == 3) then
        call MPI_IRECV(var,1,tp_recv(B),neigh(B),tag+5,comm,req_z(1),err)
        call MPI_IRECV(var,1,tp_recv(F),neigh(F),tag+6,comm,req_z(2),err)

        call MPI_ISEND(var,1,tp_send(F),neigh(F),tag+5,comm,req_z(3),err)
        call MPI_ISEND(var,1,tp_send(B),neigh(B),tag+6,comm,req_z(4),err)

        call MPI_waitall(size(req_z), req_z, MPI_STATUSES_IGNORE, err)
        endif

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

        integer, dimension(4)       :: req_x
        integer, dimension(4)       :: req_y
        integer, dimension(4)       :: req_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0

        !$acc host_data use_device(var)
        call MPI_IRECV(var,1,tp_recv(W),neigh(W),tag+1,comm,req_x(1),err)
        call MPI_IRECV(var,1,tp_recv(E),neigh(E),tag+2,comm,req_x(2),err)

        call MPI_ISEND(var,1,tp_send(E),neigh(E),tag+1,comm,req_x(3),err)
        call MPI_ISEND(var,1,tp_send(W),neigh(W),tag+2,comm,req_x(4),err)

        call MPI_waitall(size(req_x), req_x, MPI_STATUSES_IGNORE, err)

        call MPI_IRECV(var,1,tp_recv(S),neigh(S),tag+3,comm,req_y(1),err)
        call MPI_IRECV(var,1,tp_recv(N),neigh(N),tag+4,comm,req_y(2),err)

        call MPI_ISEND(var,1,tp_send(N),neigh(N),tag+3,comm,req_y(3),err)
        call MPI_ISEND(var,1,tp_send(S),neigh(S),tag+4,comm,req_y(4),err)

        call MPI_waitall(size(req_y), req_y, MPI_STATUSES_IGNORE, err)

        if(dims == 3) then
        call MPI_IRECV(var,1,tp_recv(B),neigh(B),tag+5,comm,req_z(1),err)
        call MPI_IRECV(var,1,tp_recv(F),neigh(F),tag+6,comm,req_z(2),err)

        call MPI_ISEND(var,1,tp_send(F),neigh(F),tag+5,comm,req_z(3),err)
        call MPI_ISEND(var,1,tp_send(B),neigh(B),tag+6,comm,req_z(4),err)

        call MPI_waitall(size(req_z), req_z, MPI_STATUSES_IGNORE, err)
        endif

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

        integer, dimension(4)       :: req_x
        integer, dimension(4)       :: req_y
        integer, dimension(4)       :: req_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0

        !$acc update host(var)
        call MPI_IRECV(var,1,tp_recv(W),neigh(W),tag+1,comm,req_x(1),err)
        call MPI_IRECV(var,1,tp_recv(E),neigh(E),tag+2,comm,req_x(2),err)

        call MPI_ISEND(var,1,tp_send(E),neigh(E),tag+1,comm,req_x(3),err)
        call MPI_ISEND(var,1,tp_send(W),neigh(W),tag+2,comm,req_x(4),err)

        call MPI_waitall(size(req_x), req_x, MPI_STATUSES_IGNORE, err)

        call MPI_IRECV(var,1,tp_recv(S),neigh(S),tag+3,comm,req_y(1),err)
        call MPI_IRECV(var,1,tp_recv(N),neigh(N),tag+4,comm,req_y(2),err)

        call MPI_ISEND(var,1,tp_send(N),neigh(N),tag+3,comm,req_y(3),err)
        call MPI_ISEND(var,1,tp_send(S),neigh(S),tag+4,comm,req_y(4),err)

        call MPI_waitall(size(req_y), req_y, MPI_STATUSES_IGNORE, err)

        if(dims == 3) then
        call MPI_IRECV(var,1,tp_recv(B),neigh(B),tag+5,comm,req_z(1),err)
        call MPI_IRECV(var,1,tp_recv(F),neigh(F),tag+6,comm,req_z(2),err)

        call MPI_ISEND(var,1,tp_send(F),neigh(F),tag+5,comm,req_z(3),err)
        call MPI_ISEND(var,1,tp_send(B),neigh(B),tag+6,comm,req_z(4),err)

        call MPI_waitall(size(req_z), req_z, MPI_STATUSES_IGNORE, err)
        endif

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

        integer, dimension(4)       :: req_x
        integer, dimension(4)       :: req_y
        integer, dimension(4)       :: req_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0

        !$acc host_data use_device(var)
        call MPI_IRECV(var,1,tp_recv(W),neigh(W),tag+1,comm,req_x(1),err)
        call MPI_IRECV(var,1,tp_recv(E),neigh(E),tag+2,comm,req_x(2),err)

        call MPI_ISEND(var,1,tp_send(E),neigh(E),tag+1,comm,req_x(3),err)
        call MPI_ISEND(var,1,tp_send(W),neigh(W),tag+2,comm,req_x(4),err)

        call MPI_waitall(size(req_x), req_x, MPI_STATUSES_IGNORE, err)

        call MPI_IRECV(var,1,tp_recv(S),neigh(S),tag+3,comm,req_y(1),err)
        call MPI_IRECV(var,1,tp_recv(N),neigh(N),tag+4,comm,req_y(2),err)

        call MPI_ISEND(var,1,tp_send(N),neigh(N),tag+3,comm,req_y(3),err)
        call MPI_ISEND(var,1,tp_send(S),neigh(S),tag+4,comm,req_y(4),err)

        call MPI_waitall(size(req_y), req_y, MPI_STATUSES_IGNORE, err)

        if(dims == 3) then
        call MPI_IRECV(var,1,tp_recv(B),neigh(B),tag+5,comm,req_z(1),err)
        call MPI_IRECV(var,1,tp_recv(F),neigh(F),tag+6,comm,req_z(2),err)

        call MPI_ISEND(var,1,tp_send(F),neigh(F),tag+5,comm,req_z(3),err)
        call MPI_ISEND(var,1,tp_send(B),neigh(B),tag+6,comm,req_z(4),err)

        call MPI_waitall(size(req_z), req_z, MPI_STATUSES_IGNORE, err)
        endif

        !$acc end host_data

        return
end subroutine mpi_communicate4D_real8_cuda_aware

subroutine mpi_communicate4D_real8_bfr(comm,neigh,&
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
        integer                                   , intent(in)    :: comm

        integer, parameter    :: dtype = MPI_RP
        integer, dimension(4) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer, dimension(4) :: req_x
        integer, dimension(4) :: req_y
        integer, dimension(4) :: req_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: i,j,k,l,err = 0, ii, jj, kk, q

        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)*bfr_dims_x(4)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)*bfr_dims_y(4)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)*bfr_dims_z(4)

        if(neigh(E) >= 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(4)
          do          l = 1,5
             do       k = lbz,ubz
                do    j = lby,uby
                   do i = 1,GN
                      ii = ex-GN+i
                      bfr_s_E(ii,j,k,l) = var(ii,j,k,l)
                   enddo
                enddo
             enddo
          enddo
          !$acc end parallel
          !$acc update host(bfr_s_E)
        endif

        if(neigh(W) >= 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(4)
          do          l = 1,5
             do       k = lbz,ubz
                do    j = lby,uby
                   do i = 1,GN
                      ii = sx-1+i
                      bfr_s_W(ii,j,k,l) = var(ii,j,k,l)
                   enddo
                enddo
             enddo
          enddo
          !$acc end parallel
          !$acc update host(bfr_s_W)
        endif

        q = 0
        if(neigh(E) >= 0) then
        call MPI_IRECV(bfr_r_E,bfr_size_x,dtype,neigh(E),tag+1,comm,req_x(q+1),err)
        call MPI_ISEND(bfr_s_E,bfr_size_x,dtype,neigh(E),tag+2,comm,req_x(q+2),err)
        q = q+2
        endif

        if(neigh(W) >= 0) then
        call MPI_IRECV(bfr_r_W,bfr_size_x,dtype,neigh(W),tag+2,comm,req_x(q+1),err)
        call MPI_ISEND(bfr_s_W,bfr_size_x,dtype,neigh(W),tag+1,comm,req_x(q+2),err)
        q = q+2
        endif

        if(q > 0) then
           call MPI_waitall(q, req_x, MPI_STATUSES_IGNORE, err)
           !$acc update device(bfr_r_W,bfr_r_E)
        endif
        if(neigh(W)>=0) then
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN
                    ii = sx-GN+i-1
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
                    ii = ex+i
                    var(ii,j,k,l) = bfr_r_E(ii,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        if(neigh(N) >= 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(4)
          do          l = 1,5
             do       k = lbz,ubz
                do    j = 1,GN
                   do i = lbx,ubx
                      jj = ey-GN+j
                      bfr_s_N(i,jj,k,l) = var(i,jj,k,l)
                   enddo
                enddo
             enddo
          enddo
          !$acc end parallel
          !$acc update host(bfr_s_N)
        endif

        if(neigh(S) >= 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(4)
          do          l = 1,5
             do       k = lbz,ubz
                do    j = 1,GN
                   do i = lbx,ubx
                      jj = sy-1+j
                      bfr_s_S(i,jj,k,l) = var(i,jj,k,l)
                   enddo
                enddo
             enddo
          enddo
          !$acc end parallel
          !$acc update host(bfr_s_S)
        endif

        q = 0
        if(neigh(N) >= 0) then
        call MPI_IRECV(bfr_r_N,bfr_size_y,dtype,neigh(N),tag+1,comm,req_y(q+1),err)
        call MPI_ISEND(bfr_s_N,bfr_size_y,dtype,neigh(N),tag+2,comm,req_y(q+2),err)
        q = q+2
        endif

        if(neigh(S) >= 0) then
        call MPI_IRECV(bfr_r_S,bfr_size_y,dtype,neigh(S),tag+2,comm,req_y(q+1),err)
        call MPI_ISEND(bfr_s_S,bfr_size_y,dtype,neigh(S),tag+1,comm,req_y(q+2),err)
        q = q+2
        endif

        if(q > 0) then
          call MPI_waitall(q, req_y, MPI_STATUSES_IGNORE, err)
          !$acc update device(bfr_r_N,bfr_r_S)
        endif
        if(neigh(S)>=0) then
        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN
                 do i = lbx,ubx
                    jj = sy-GN+j-1
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
                    jj = ey+j
                    var(i,jj,k,l) = bfr_r_N(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        if(neigh(F) >= 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(4)
          do          l = 1,5
             do       k = 1,GN
                do    j = lby,uby
                   do i = lbx,ubx
                      kk = ez-GN+k
                      bfr_s_F(i,j,kk,l) = var(i,j,kk,l)
                   enddo
                enddo
             enddo
          enddo
          !$acc end parallel
          !$acc update host(bfr_s_F)
        endif

        if(neigh(B) >= 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(4)
          do          l = 1,5
             do       k = 1,GN
                do    j = lby,uby
                   do i = lbx,ubx
                      kk = sz-1+k
                      bfr_s_B(i,j,kk,l) = var(i,j,kk,l)
                   enddo
                enddo
             enddo
          enddo
          !$acc end parallel
          !$acc update host(bfr_s_B)
        endif

        q = 0
        if(neigh(F) >= 0) then
        call MPI_IRECV(bfr_r_F,bfr_size_z,dtype,neigh(F),tag+1,comm,req_z(q+1),err)
        call MPI_ISEND(bfr_s_F,bfr_size_z,dtype,neigh(F),tag+2,comm,req_z(q+2),err)
        q = q+2
        endif

        if(neigh(B) >= 0) then
        call MPI_IRECV(bfr_r_B,bfr_size_z,dtype,neigh(B),tag+2,comm,req_z(q+1),err)
        call MPI_ISEND(bfr_s_B,bfr_size_z,dtype,neigh(B),tag+1,comm,req_z(q+2),err)
        q = q+2
        endif

        if(q > 0) then
          call MPI_waitall(q, req_z, MPI_STATUSES_IGNORE, err)
          !$acc update device(bfr_r_F,bfr_r_B)
        endif
        if(neigh(B)>=0) then
        !$acc parallel default(present) async(5)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = sz-GN+k-1
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
                    kk = ez+k
                    var(i,j,kk,l) = bfr_r_F(i,j,kk,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        return
end subroutine mpi_communicate4D_real8_bfr

subroutine mpi_communicate4D_real8_bfr_cuda_aware(comm,neigh,&
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
        integer                                   , intent(in)    :: comm

        integer, parameter    :: dtype = MPI_RP
        integer, dimension(4) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer, dimension(4) :: req_x
        integer, dimension(4) :: req_y
        integer, dimension(4) :: req_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: i,j,k,l,err = 0, ii, jj, kk, q
        integer, dimension(6) :: nb

        nb = neigh
        where(nb == rank) nb = MPI_PROC_NULL

        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)*bfr_dims_x(4)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)*bfr_dims_y(4)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)*bfr_dims_z(4)

        call StartProfRange('halo_x_pack')
        if(neigh(E) == rank) then
        !$acc parallel default(present)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN
                    var(ex+i     ,j,k,l) = var(sx-1+i ,j,k,l)
                    var(sx-GN+i-1,j,k,l) = var(ex-GN+i,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(nb(E) >= 0) then
          !$acc parallel default(present) async(1)
          !$acc loop gang, vector collapse(4)
          do          l = 1,5
             do       k = lbz,ubz
                do    j = lby,uby
                   do i = 1,GN
                      ii = ex-GN+i
                      bfr_s_E(ii,j,k,l) = var(ii,j,k,l)
                   enddo
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(nb(W) >= 0) then
          !$acc parallel default(present) async(2)
          !$acc loop gang, vector collapse(4)
          do          l = 1,5
             do       k = lbz,ubz
                do    j = lby,uby
                   do i = 1,GN
                      ii = sx-1+i
                      bfr_s_W(ii,j,k,l) = var(ii,j,k,l)
                   enddo
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait

        call EndProfRange
        call StartProfRange('halo_x_mpi')
        !$acc host_data use_device(bfr_s_E,bfr_r_W) &
        !$acc           use_device(bfr_s_W,bfr_r_E)

        q = 0
        if(nb(E) >= 0) then
        call MPI_IRECV(bfr_r_E,bfr_size_x,dtype,nb(E),tag+1,comm,req_x(q+1),err)
        call MPI_ISEND(bfr_s_E,bfr_size_x,dtype,nb(E),tag+2,comm,req_x(q+2),err)
        q = q+2
        endif

        if(nb(W) >= 0) then
        call MPI_IRECV(bfr_r_W,bfr_size_x,dtype,nb(W),tag+2,comm,req_x(q+1),err)
        call MPI_ISEND(bfr_s_W,bfr_size_x,dtype,nb(W),tag+1,comm,req_x(q+2),err)
        q = q+2
        endif

        if(q > 0) call MPI_waitall(q, req_x, MPI_STATUSES_IGNORE, err)

        !$acc end host_data
        call EndProfRange
        call StartProfRange('halo_x_unpack')

        if(nb(W)>=0) then
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN
                    ii = sx-GN+i-1
                    var(ii,j,k,l) = bfr_r_W(ii,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(nb(E)>=0) then
        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN
                    ii = ex+i
                    var(ii,j,k,l) = bfr_r_E(ii,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        call EndProfRange
        call StartProfRange('halo_y_pack')
        if(neigh(N) == rank) then
        !$acc parallel default(present)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN
                 do i = lbx,ubx
                    var(i,ey+j     ,k,l) = var(i,sy-1+j ,k,l)
                    var(i,sy-GN+j-1,k,l) = var(i,ey-GN+j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(nb(N) >= 0) then
          !$acc parallel default(present) async(3)
          !$acc loop gang, vector collapse(4)
          do          l = 1,5
             do       k = lbz,ubz
                do    j = 1,GN
                   do i = lbx,ubx
                      jj = ey-GN+j
                      bfr_s_N(i,jj,k,l) = var(i,jj,k,l)
                   enddo
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(nb(S) >= 0) then
          !$acc parallel default(present) async(4)
          !$acc loop gang, vector collapse(4)
          do          l = 1,5
             do       k = lbz,ubz
                do    j = 1,GN
                   do i = lbx,ubx
                      jj = sy-1+j
                      bfr_s_S(i,jj,k,l) = var(i,jj,k,l)
                   enddo
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait

        call EndProfRange
        call StartProfRange('halo_y_mpi')
        !$acc host_data use_device(bfr_s_S,bfr_r_N) &
        !$acc           use_device(bfr_s_N,bfr_r_S)

        q = 0
        if(nb(S) >= 0) then
        call MPI_IRECV(bfr_r_S,bfr_size_y,dtype,nb(S),tag+1,comm,req_y(q+1),err)
        call MPI_ISEND(bfr_s_S,bfr_size_y,dtype,nb(S),tag+2,comm,req_y(q+2),err)
        q = q+2
        endif

        if(nb(N) >= 0) then
        call MPI_IRECV(bfr_r_N,bfr_size_y,dtype,nb(N),tag+2,comm,req_y(q+1),err)
        call MPI_ISEND(bfr_s_N,bfr_size_y,dtype,nb(N),tag+1,comm,req_y(q+2),err)
        q = q+2
        endif

        if(q > 0) call MPI_waitall(q, req_y, MPI_STATUSES_IGNORE, err)

        !$acc end host_data
        call EndProfRange
        call StartProfRange('halo_y_unpack')

        if(nb(S)>=0) then
        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN
                 do i = lbx,ubx
                    jj = sy-GN+j-1
                    var(i,jj,k,l) = bfr_r_S(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(nb(N)>=0) then
        !$acc parallel default(present) async(4)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN
                 do i = lbx,ubx
                    jj = ey+j
                    var(i,jj,k,l) = bfr_r_N(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        call EndProfRange
        call StartProfRange('halo_z_pack')
        if(neigh(F) == rank) then
        !$acc parallel default(present)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN
              do    j = lby,uby
                 do i = lbx,ubx
                    var(i,j,ez+k     ,l) = var(i,j,sz-1+k ,l)
                    var(i,j,sz-GN+k-1,l) = var(i,j,ez-GN+k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(nb(F) >= 0) then
          !$acc parallel default(present) async(5)
          !$acc loop gang, vector collapse(4)
          do          l = 1,5
             do       k = 1,GN
                do    j = lby,uby
                   do i = lbx,ubx
                      kk = ez-GN+k
                      bfr_s_F(i,j,kk,l) = var(i,j,kk,l)
                   enddo
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(nb(B) >= 0) then
          !$acc parallel default(present) async(6)
          !$acc loop gang, vector collapse(4)
          do          l = 1,5
             do       k = 1,GN
                do    j = lby,uby
                   do i = lbx,ubx
                      kk = sz-1+k
                      bfr_s_B(i,j,kk,l) = var(i,j,kk,l)
                   enddo
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait

        call EndProfRange
        call StartProfRange('halo_z_mpi')
        !$acc host_data use_device(bfr_s_B,bfr_r_F) &
        !$acc           use_device(bfr_s_F,bfr_r_B)

        q = 0
        if(nb(B) >= 0) then
        call MPI_IRECV(bfr_r_B,bfr_size_z,dtype,nb(B),tag+1,comm,req_z(q+1),err)
        call MPI_ISEND(bfr_s_B,bfr_size_z,dtype,nb(B),tag+2,comm,req_z(q+2),err)
        q = q+2
        endif

        if(nb(F) >= 0) then
        call MPI_IRECV(bfr_r_F,bfr_size_z,dtype,nb(F),tag+2,comm,req_z(q+1),err)
        call MPI_ISEND(bfr_s_F,bfr_size_z,dtype,nb(F),tag+1,comm,req_z(q+2),err)
        q = q+2
        endif

        if(q > 0) call MPI_waitall(q, req_z, MPI_STATUSES_IGNORE, err)

        !$acc end host_data
        call EndProfRange
        call StartProfRange('halo_z_unpack')

        if(nb(B)>=0) then
        !$acc parallel default(present) async(5)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = sz-GN+k-1
                    var(i,j,kk,l) = bfr_r_B(i,j,kk,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(nb(F)>=0) then
        !$acc parallel default(present) async(6)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = ez+k
                    var(i,j,kk,l) = bfr_r_F(i,j,kk,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        call EndProfRange
        return
end subroutine mpi_communicate4D_real8_bfr_cuda_aware

subroutine mpi_communicate4D_real8_bfr_one_sided(comm,neigh,&
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
        integer                                   , intent(in)    :: comm

        integer, parameter    :: dtype = MPI_RP
        integer, dimension(4) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer, dimension(4) :: req_x
        integer, dimension(4) :: req_y
        integer, dimension(4) :: req_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: i,j,k,l,err, ii, jj, kk, q

        integer(KIND=MPI_ADDRESS_KIND) :: target_disp

        target_disp = 0
        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)*bfr_dims_x(4)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)*bfr_dims_y(4)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)*bfr_dims_z(4)

        if(neigh(E) >= 0) then
          !$acc parallel default(present) async(1)
          !$acc loop gang, vector collapse(4)
          do          l = 1,5
             do       k = lbz,ubz
                do    j = lby,uby
                   do i = 1,GN
                      ii = ex-GN+i
                      bfr_s_E(ii,j,k,l) = var(ii,j,k,l)
                   enddo
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(W) >= 0) then
          !$acc parallel default(present) async(2)
          !$acc loop gang, vector collapse(4)
          do          l = 1,5
             do       k = lbz,ubz
                do    j = lby,uby
                   do i = 1,GN
                      ii = sx-1+i
                      bfr_s_W(ii,j,k,l) = var(ii,j,k,l)
                   enddo
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait
        !$acc host_data use_device(bfr_s_E) &
        !$acc           use_device(bfr_s_W)

        if (neigh(W) >= 0) call MPI_WIN_POST(mpi_group_W, 0, win_phi_bfr_W, err)
        if (neigh(E) >= 0) call MPI_WIN_POST(mpi_group_E, 0, win_phi_bfr_E, err)

        if (neigh(W) >= 0) then
         call MPI_WIN_START(mpi_group_W, 0, win_phi_bfr_E, err)
         call MPI_PUT(bfr_s_W,bfr_size_x,dtype,neigh(W),target_disp,bfr_size_x,dtype,win_phi_bfr_E,err)
         !$acc wait
         call MPI_WIN_COMPLETE(win_phi_bfr_E, err)
        endif
        if (neigh(E) >= 0) then
         call MPI_WIN_START(mpi_group_E, 0, win_phi_bfr_W, err)
         call MPI_PUT(bfr_s_E,bfr_size_x,dtype,neigh(E),target_disp,bfr_size_x,dtype,win_phi_bfr_W,err)
         !$acc wait
         call MPI_WIN_COMPLETE(win_phi_bfr_W, err)
        endif

        if (neigh(W) >= 0) call MPI_WIN_WAIT(win_phi_bfr_W, err)
        if (neigh(E) >= 0) call MPI_WIN_WAIT(win_phi_bfr_E, err)
        !$acc end host_data

        if(neigh(W)>=0) then
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN
                    ii = sx-GN+i-1
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
                    ii = ex+i
                    var(ii,j,k,l) = bfr_r_E(ii,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        if(neigh(N) >= 0) then
          !$acc parallel default(present) async(3)
          !$acc loop gang, vector collapse(4)
          do          l = 1,5
             do       k = lbz,ubz
                do    j = 1,GN
                   do i = lbx,ubx
                      jj = ey-GN+j
                      bfr_s_N(i,jj,k,l) = var(i,jj,k,l)
                   enddo
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(S) >= 0) then
          !$acc parallel default(present) async(4)
          !$acc loop gang, vector collapse(4)
          do          l = 1,5
             do       k = lbz,ubz
                do    j = 1,GN
                   do i = lbx,ubx
                      jj = sy-1+j
                      bfr_s_S(i,jj,k,l) = var(i,jj,k,l)
                   enddo
                enddo
             enddo
          enddo
          !$acc end parallel
        endif
        !$acc wait
        !$acc host_data use_device(bfr_s_S) &
        !$acc           use_device(bfr_s_N)

        if (neigh(S) >= 0) call MPI_WIN_POST(mpi_group_S, 0, win_phi_bfr_S, err)
        if (neigh(N) >= 0) call MPI_WIN_POST(mpi_group_N, 0, win_phi_bfr_N, err)

        if (neigh(S) >= 0) then
         call MPI_WIN_START(mpi_group_S, 0, win_phi_bfr_N, err)
         call MPI_PUT(bfr_s_S,bfr_size_y,dtype,neigh(S),target_disp,bfr_size_y,dtype,win_phi_bfr_N,err)
         !$acc wait
         call MPI_WIN_COMPLETE(win_phi_bfr_N, err)
        endif
        if (neigh(N) >= 0) then
         call MPI_WIN_START(mpi_group_N, 0, win_phi_bfr_S, err)
         call MPI_PUT(bfr_s_N,bfr_size_y,dtype,neigh(N),target_disp,bfr_size_y,dtype,win_phi_bfr_S,err)
         !$acc wait
         call MPI_WIN_COMPLETE(win_phi_bfr_S, err)
        endif

        if (neigh(S) >= 0) call MPI_WIN_WAIT(win_phi_bfr_S, err)
        if (neigh(N) >= 0) call MPI_WIN_WAIT(win_phi_bfr_N, err)
        !$acc end host_data

        if(neigh(S)>=0) then
        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN
                 do i = lbx,ubx
                    jj = sy-GN+j-1
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
                    jj = ey+j
                    var(i,jj,k,l) = bfr_r_N(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        if(neigh(F) >= 0) then
          !$acc parallel default(present) async(5)
          !$acc loop gang, vector collapse(4)
          do          l = 1,5
             do       k = 1,GN
                do    j = lby,uby
                   do i = lbx,ubx
                      kk = ez-GN+k
                      bfr_s_F(i,j,kk,l) = var(i,j,kk,l)
                   enddo
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(B) >= 0) then
          !$acc parallel default(present) async(6)
          !$acc loop gang, vector collapse(4)
          do          l = 1,5
             do       k = 1,GN
                do    j = lby,uby
                   do i = lbx,ubx
                      kk = sz-1+k
                      bfr_s_B(i,j,kk,l) = var(i,j,kk,l)
                   enddo
                enddo
             enddo
          enddo
          !$acc end parallel
        endif
        !$acc wait
        !$acc host_data use_device(bfr_s_B) &
        !$acc           use_device(bfr_s_F)

        if (neigh(B) >= 0) call MPI_WIN_POST(mpi_group_B, 0, win_phi_bfr_B, err)
        if (neigh(F) >= 0) call MPI_WIN_POST(mpi_group_F, 0, win_phi_bfr_F, err)

        if (neigh(B) >= 0) then
         call MPI_WIN_START(mpi_group_B, 0, win_phi_bfr_F, err)
         call MPI_PUT(bfr_s_B,bfr_size_z,dtype,neigh(B),target_disp,bfr_size_z,dtype,win_phi_bfr_F,err)
         !$acc wait
         call MPI_WIN_COMPLETE(win_phi_bfr_F, err)
        endif
        if (neigh(F) >= 0) then
         call MPI_WIN_START(mpi_group_F, 0, win_phi_bfr_B, err)
         call MPI_PUT(bfr_s_F,bfr_size_z,dtype,neigh(F),target_disp,bfr_size_z,dtype,win_phi_bfr_B,err)
         !$acc wait
         call MPI_WIN_COMPLETE(win_phi_bfr_B, err)
        endif

        if (neigh(B) >= 0) call MPI_WIN_WAIT(win_phi_bfr_B, err)
        if (neigh(F) >= 0) call MPI_WIN_WAIT(win_phi_bfr_F, err)
        !$acc end host_data

        if(neigh(B)>=0) then
        !$acc parallel default(present) async(5)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = sz-GN+k-1
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
                    kk = ez+k
                    var(i,j,kk,l) = bfr_r_F(i,j,kk,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        return
end subroutine mpi_communicate4D_real8_bfr_one_sided

subroutine mpi_communicate3D_real8_bfr(comm,neigh,&
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
        integer                                 , intent(in)    :: comm

        integer, parameter    :: dtype = MPI_RP
        integer, dimension(3) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer, dimension(4) :: req_x
        integer, dimension(4) :: req_y
        integer, dimension(4) :: req_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: i,j,k,err = 0, ii, jj, kk, q

        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)

        if(neigh(E) >= 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN
                   ii = ex-GN+i
                   bfr_s_E(ii,j,k) = var(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
          !$acc update host(bfr_s_E)
        endif

        if(neigh(W) >= 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN
                   ii = sx-1+i
                   bfr_s_W(ii,j,k) = var(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
          !$acc update host(bfr_s_W)
        endif

        q = 0
        if(neigh(E) >= 0) then
        call MPI_IRECV(bfr_r_E,bfr_size_x,dtype,neigh(E),tag+1,comm,req_x(q+1),err)
        call MPI_ISEND(bfr_s_E,bfr_size_x,dtype,neigh(E),tag+2,comm,req_x(q+2),err)
        q = q+2
        endif

        if(neigh(W) >= 0) then
        call MPI_IRECV(bfr_r_W,bfr_size_x,dtype,neigh(W),tag+2,comm,req_x(q+1),err)
        call MPI_ISEND(bfr_s_W,bfr_size_x,dtype,neigh(W),tag+1,comm,req_x(q+2),err)
        q = q+2
        endif

        if(q > 0) then
          call MPI_waitall(q, req_x, MPI_STATUSES_IGNORE, err)
          !$acc update device(bfr_r_W,bfr_r_E)
        endif

        if(neigh(W)>=0) then
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN
                 ii = sx-GN+i-1
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
                 ii = ex+i
                 var(ii,j,k) = bfr_r_E(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        if(neigh(N) >= 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN
                do i = lbx,ubx
                   jj = ey-GN+j
                   bfr_s_N(i,jj,k) = var(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
          !$acc update host(bfr_s_N)
        endif

        if(neigh(S) >= 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN
                do i = lbx,ubx
                   jj = sy-1+j
                   bfr_s_S(i,jj,k) = var(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
          !$acc update host(bfr_s_S)
        endif

        q = 0
        if(neigh(N) >= 0) then
        call MPI_IRECV(bfr_r_N,bfr_size_y,dtype,neigh(N),tag+1,comm,req_y(q+1),err)
        call MPI_ISEND(bfr_s_N,bfr_size_y,dtype,neigh(N),tag+2,comm,req_y(q+2),err)
        q = q+2
        endif

        if(neigh(S) >= 0) then
        call MPI_IRECV(bfr_r_S,bfr_size_y,dtype,neigh(S),tag+2,comm,req_y(q+1),err)
        call MPI_ISEND(bfr_s_S,bfr_size_y,dtype,neigh(S),tag+1,comm,req_y(q+2),err)
        q = q+2
        endif

        if(q > 0) then
          call MPI_waitall(q, req_y, MPI_STATUSES_IGNORE, err)
          !$acc update device(bfr_r_N,bfr_r_S)
        endif

        if(neigh(S)>=0) then
        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN
              do i = lbx,ubx
                 jj = sy-GN+j-1
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
                 jj = ey+j
                 var(i,jj,k) = bfr_r_N(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        if(neigh(F) >= 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = ez-GN+k
                   bfr_s_F(i,j,kk) = var(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel
          !$acc update host(bfr_s_F)
        endif

        if(neigh(B) >= 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = sz-1+k
                   bfr_s_B(i,j,kk) = var(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel
          !$acc update host(bfr_s_B)
        endif

        q = 0
        if(neigh(F) >= 0) then
        call MPI_IRECV(bfr_r_F,bfr_size_z,dtype,neigh(F),tag+1,comm,req_z(q+1),err)
        call MPI_ISEND(bfr_s_F,bfr_size_z,dtype,neigh(F),tag+2,comm,req_z(q+2),err)
        q = q+2
        endif

        if(neigh(B) >= 0) then
        call MPI_IRECV(bfr_r_B,bfr_size_z,dtype,neigh(B),tag+2,comm,req_z(q+1),err)
        call MPI_ISEND(bfr_s_B,bfr_size_z,dtype,neigh(B),tag+1,comm,req_z(q+2),err)
        q = q+2
        endif

        if(q > 0) then
          call MPI_waitall(q, req_z, MPI_STATUSES_IGNORE, err)
          !$acc update device(bfr_r_F,bfr_r_B)
        endif

        if(neigh(B)>=0) then
        !$acc parallel default(present) async(5)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN
           do    j = lby,uby
              do i = lbx,ubx
                 kk = sz-GN+k-1
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
                 kk = ez+k
                 var(i,j,kk) = bfr_r_F(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        return
end subroutine mpi_communicate3D_real8_bfr

subroutine mpi_communicate3D_real8_bfr_cuda_aware(comm,neigh,&
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
        integer                                 , intent(in)    :: comm

        integer, parameter    :: dtype = MPI_RP
        integer, dimension(3) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer, dimension(4) :: req_x
        integer, dimension(4) :: req_y
        integer, dimension(4) :: req_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: i,j,k,err = 0, ii, jj, kk, q
        integer, dimension(6) :: nb

        nb = neigh
        where(nb == rank) nb = MPI_PROC_NULL

        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)

        if(neigh(E) == rank) then
        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN
                    var(ex+i     ,j,k) = var(sx-1+i ,j,k)
                    var(sx-GN+i-1,j,k) = var(ex-GN+i,j,k)
                 enddo
              enddo
           enddo
        !$acc end parallel
        endif

        if(nb(E) >= 0) then
          !$acc parallel default(present) async(1)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN
                   ii = ex-GN+i
                   bfr_s_E(ii,j,k) = var(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(nb(W) >= 0) then
          !$acc parallel default(present) async(2)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN
                   ii = sx-1+i
                   bfr_s_W(ii,j,k) = var(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait
        !$acc host_data use_device(bfr_s_E,bfr_r_W) &
        !$acc           use_device(bfr_s_W,bfr_r_E)

        q = 0
        if(nb(W) >= 0) then
        call MPI_IRECV(bfr_r_W,bfr_size_x,dtype,nb(W),tag+1,comm,req_x(q+1),err)
        call MPI_ISEND(bfr_s_W,bfr_size_x,dtype,nb(W),tag+2,comm,req_x(q+2),err)
        q = q+2
        endif

        if(nb(E) >= 0) then
        call MPI_IRECV(bfr_r_E,bfr_size_x,dtype,nb(E),tag+2,comm,req_x(q+1),err)
        call MPI_ISEND(bfr_s_E,bfr_size_x,dtype,nb(E),tag+1,comm,req_x(q+2),err)
        q = q+2
        endif

        if(q > 0) call MPI_waitall(q, req_x, MPI_STATUSES_IGNORE, err)
        !$acc end host_data

        if(nb(W)>=0) then
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN
                 ii = sx-GN+i-1
                 var(ii,j,k) = bfr_r_W(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(nb(E)>=0) then
        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN
                 ii = ex+i
                 var(ii,j,k) = bfr_r_E(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        if(neigh(N) == rank) then
        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
           do       k = lbz,ubz
              do    j = 1,GN
                 do i = lbx,ubx
                    var(i,ey+j     ,k) = var(i,sy-1+j ,k)
                    var(i,sy-GN+j-1,k) = var(i,ey-GN+j,k)
                 enddo
              enddo
           enddo
        !$acc end parallel
        endif

        if(nb(N) >= 0) then
          !$acc parallel default(present) async(3)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN
                do i = lbx,ubx
                   jj = ey-GN+j
                   bfr_s_N(i,jj,k) = var(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(nb(S) >= 0) then
          !$acc parallel default(present) async(4)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN
                do i = lbx,ubx
                   jj = sy-1+j
                   bfr_s_S(i,jj,k) = var(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif
        !$acc wait
        !$acc host_data use_device(bfr_s_S,bfr_r_N) &
        !$acc           use_device(bfr_s_N,bfr_r_S)

        q = 0
        if(nb(S) >= 0) then
        call MPI_IRECV(bfr_r_S,bfr_size_y,dtype,nb(S),tag+1,comm,req_y(q+1),err)
        call MPI_ISEND(bfr_s_S,bfr_size_y,dtype,nb(S),tag+2,comm,req_y(q+2),err)
        q = q+2
        endif

        if(nb(N) >= 0) then
        call MPI_IRECV(bfr_r_N,bfr_size_y,dtype,nb(N),tag+2,comm,req_y(q+1),err)
        call MPI_ISEND(bfr_s_N,bfr_size_y,dtype,nb(N),tag+1,comm,req_y(q+2),err)
        q = q+2
        endif

        if(q > 0) call MPI_waitall(q, req_y, MPI_STATUSES_IGNORE, err)
        !$acc end host_data

        if(nb(S)>=0) then
        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN
              do i = lbx,ubx
                 jj = sy-GN+j-1
                 var(i,jj,k) = bfr_r_S(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(nb(N)>=0) then
        !$acc parallel default(present) async(4)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN
              do i = lbx,ubx
                 jj = ey+j
                 var(i,jj,k) = bfr_r_N(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        if(neigh(F) == rank) then
        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
           do       k = 1,GN
              do    j = lby,uby
                 do i = lbx,ubx
                    var(i,j,ez+k     ) = var(i,j,sz-1+k )
                    var(i,j,sz-GN+k-1) = var(i,j,ez-GN+k)
                 enddo
              enddo
           enddo
        !$acc end parallel
        endif

        if(nb(F) >= 0) then
          !$acc parallel default(present) async(5)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = ez-GN+k
                   bfr_s_F(i,j,kk) = var(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(nb(B) >= 0) then
          !$acc parallel default(present) async(6)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = sz-1+k
                   bfr_s_B(i,j,kk) = var(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif
        !$acc wait
        !$acc host_data use_device(bfr_s_B,bfr_r_F) &
        !$acc           use_device(bfr_s_F,bfr_r_B)

        q = 0
        if(nb(B) >= 0) then
        call MPI_IRECV(bfr_r_B,bfr_size_z,dtype,nb(B),tag+1,comm,req_z(q+1),err)
        call MPI_ISEND(bfr_s_B,bfr_size_z,dtype,nb(B),tag+2,comm,req_z(q+2),err)
        q = q+2
        endif

        if(nb(F) >= 0) then
        call MPI_IRECV(bfr_r_F,bfr_size_z,dtype,nb(F),tag+2,comm,req_z(q+1),err)
        call MPI_ISEND(bfr_s_F,bfr_size_z,dtype,nb(F),tag+1,comm,req_z(q+2),err)
        q = q+2
        endif

        if(q > 0) call MPI_waitall(q, req_z, MPI_STATUSES_IGNORE, err)
        !$acc end host_data

        if(nb(B)>=0) then
        !$acc parallel default(present) async(5)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN
           do    j = lby,uby
              do i = lbx,ubx
                 kk = sz-GN+k-1
                 var(i,j,kk) = bfr_r_B(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        if(nb(F)>=0) then
        !$acc parallel default(present) async(6)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN
           do    j = lby,uby
              do i = lbx,ubx
                 kk = ez+k
                 var(i,j,kk) = bfr_r_F(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        return
end subroutine mpi_communicate3D_real8_bfr_cuda_aware

subroutine mpi_communicate3D_real8_bfr_one_sided(comm,neigh,&
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
        integer                                 , intent(in)    :: comm

        integer, parameter    :: dtype = MPI_RP
        integer, dimension(3) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer, dimension(4) :: req_x
        integer, dimension(4) :: req_y
        integer, dimension(4) :: req_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: i,j,k,err = 0, ii, jj, kk, q

        integer(KIND=MPI_ADDRESS_KIND), parameter :: target_disp = 0

        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)

        if(neigh(E) >= 0) then
          !$acc parallel default(present) async(1)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN
                   ii = ex-GN+i
                   bfr_s_E(ii,j,k) = var(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(W) >= 0) then
          !$acc parallel default(present) async(2)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN
                   ii = sx-1+i
                   bfr_s_W(ii,j,k) = var(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait
        !$acc host_data use_device(bfr_s_E,bfr_r_W) &
        !$acc           use_device(bfr_s_W,bfr_r_E)

        if (neigh(W) >= 0) call MPI_WIN_POST(mpi_group_W, 0, win_bfr_W, err)
        if (neigh(E) >= 0) call MPI_WIN_POST(mpi_group_E, 0, win_bfr_E, err)

        if (neigh(W) >= 0) then
         call MPI_WIN_START(mpi_group_W, 0, win_bfr_E, err)
         call MPI_PUT(bfr_s_W,bfr_size_x,dtype,neigh(W),target_disp,bfr_size_x,dtype,win_bfr_E,err)
         !$acc wait
         call MPI_WIN_COMPLETE(win_bfr_E, err)
        endif
        if (neigh(E) >= 0) then
         call MPI_WIN_START(mpi_group_E, 0, win_bfr_W, err)
         call MPI_PUT(bfr_s_E,bfr_size_x,dtype,neigh(E),target_disp,bfr_size_x,dtype,win_bfr_W,err)
         !$acc wait
         call MPI_WIN_COMPLETE(win_bfr_W, err)
        endif

        if (neigh(W) >= 0) call MPI_WIN_WAIT(win_bfr_W, err)
        if (neigh(E) >= 0) call MPI_WIN_WAIT(win_bfr_E, err)
        !$acc end host_data

        if(neigh(W)>=0) then
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN
                 ii = sx-GN+i-1
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
                 ii = ex+i
                 var(ii,j,k) = bfr_r_E(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        if(neigh(N) >= 0) then
          !$acc parallel default(present) async(3)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN
                do i = lbx,ubx
                   jj = ey-GN+j
                   bfr_s_N(i,jj,k) = var(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(S) >= 0) then
          !$acc parallel default(present) async(4)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN
                do i = lbx,ubx
                   jj = sy-1+j
                   bfr_s_S(i,jj,k) = var(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif
        !$acc wait
        !$acc host_data use_device(bfr_s_S,bfr_r_N) &
        !$acc           use_device(bfr_s_N,bfr_r_S)

        if (neigh(S) >= 0) call MPI_WIN_POST(mpi_group_S, 0, win_bfr_S, err)
        if (neigh(N) >= 0) call MPI_WIN_POST(mpi_group_N, 0, win_bfr_N, err)

        if (neigh(S) >= 0) then
         call MPI_WIN_START(mpi_group_S, 0, win_bfr_N, err)
         call MPI_PUT(bfr_s_S,bfr_size_y,dtype,neigh(S),target_disp,bfr_size_y,dtype,win_bfr_N,err)
         !$acc wait
         call MPI_WIN_COMPLETE(win_bfr_N, err)
        endif
        if (neigh(N) >= 0) then
         call MPI_WIN_START(mpi_group_N, 0, win_bfr_S, err)
         call MPI_PUT(bfr_s_N,bfr_size_y,dtype,neigh(N),target_disp,bfr_size_y,dtype,win_bfr_S,err)
         !$acc wait
         call MPI_WIN_COMPLETE(win_bfr_S, err)
        endif

        if (neigh(S) >= 0) call MPI_WIN_WAIT(win_bfr_S, err)
        if (neigh(N) >= 0) call MPI_WIN_WAIT(win_bfr_N, err)
        !$acc end host_data

        if(neigh(S)>=0) then
        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN
              do i = lbx,ubx
                 jj = sy-GN+j-1
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
                 jj = ey+j
                 var(i,jj,k) = bfr_r_N(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        if(neigh(F) >= 0) then
          !$acc parallel default(present) async(5)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = ez-GN+k
                   bfr_s_F(i,j,kk) = var(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(B) >= 0) then
          !$acc parallel default(present) async(6)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = sz-1+k
                   bfr_s_B(i,j,kk) = var(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif
        !$acc wait
        !$acc host_data use_device(bfr_s_B,bfr_r_F) &
        !$acc           use_device(bfr_s_F,bfr_r_B)

        if (neigh(B) >= 0) call MPI_WIN_POST(mpi_group_B, 0, win_bfr_B, err)
        if (neigh(F) >= 0) call MPI_WIN_POST(mpi_group_F, 0, win_bfr_F, err)

        if (neigh(B) >= 0) then
         call MPI_WIN_START(mpi_group_B, 0, win_bfr_F, err)
         call MPI_PUT(bfr_s_B,bfr_size_z,dtype,neigh(B),target_disp,bfr_size_z,dtype,win_bfr_F,err)
         !$acc wait
         call MPI_WIN_COMPLETE(win_bfr_F, err)
        endif
        if (neigh(F) >= 0) then
         call MPI_WIN_START(mpi_group_F, 0, win_bfr_B, err)
         call MPI_PUT(bfr_s_F,bfr_size_z,dtype,neigh(F),target_disp,bfr_size_z,dtype,win_bfr_B,err)
         !$acc wait
         call MPI_WIN_COMPLETE(win_bfr_B, err)
        endif

        if (neigh(B) >= 0) call MPI_WIN_WAIT(win_bfr_B, err)
        if (neigh(F) >= 0) call MPI_WIN_WAIT(win_bfr_F, err)
        !$acc end host_data

        if(neigh(B)>=0) then
        !$acc parallel default(present) async(5)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN
           do    j = lby,uby
              do i = lbx,ubx
                 kk = sz-GN+k-1
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
                 kk = ez+k
                 var(i,j,kk) = bfr_r_F(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        return
end subroutine mpi_communicate3D_real8_bfr_one_sided

subroutine mpi_communicate4D_real8_sendrecv(comm,neigh,&
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
        integer                                   , intent(in)    :: comm

        integer, parameter    :: dtype = MPI_RP
        integer, dimension(4) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, dimension(MPI_STATUS_SIZE) :: status
        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0, i,j,k,l, ii, jj, kk

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
                    ii = ex-GN+i
                    bfr_s_E(ii,j,k,l) = var(ii,j,k,l)
                    ii = sx-1+i
                    bfr_s_W(ii,j,k,l) = var(ii,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc wait(1)
        !$acc update host(bfr_s_E) async(4)
        !$acc update host(bfr_s_W) async(5)

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

        !$acc wait(10)
        if(neigh(W) >=0) then
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN
                    ii = sx-GN+i-1
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
                    ii = ex+i
                    var(ii,j,k,l) = bfr_r_E(ii,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN
                 do i = lbx,ubx
                    jj = ey-GN+j
                    bfr_s_N(i,jj,k,l) = var(i,jj,k,l)
                    jj = sy-1+j
                    bfr_s_S(i,jj,k,l) = var(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc wait(2)
        !$acc update host(bfr_s_S) async(6)
        !$acc update host(bfr_s_N) async(7)

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

        !$acc wait(12)
        if(neigh(N)>=0) then
        !$acc parallel default(present) async(4)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN
                 do i = lbx,ubx
                    jj = ey+j
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
                    jj = sy-GN+j-1
                    var(i,jj,k,l) = bfr_r_S(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = ez-GN+k
                    bfr_s_F(i,j,kk,l) = var(i,j,kk,l)
                    kk = sz-1+k
                    bfr_s_B(i,j,kk,l) = var(i,j,kk,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc wait(3)
        !$acc update host(bfr_s_B) async(8)
        !$acc update host(bfr_s_F) async(9)

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

        !$acc wait(14)
        if(neigh(F)>=0) then
        !$acc parallel default(present) async(6)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = ez+k
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
                    kk = sz-GN+k-1
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

subroutine mpi_communicate4D_real8_sendrecv_cuda_aware(comm,neigh,&
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
        integer                                   , intent(in)    :: comm

        integer, parameter    :: dtype = MPI_RP
        integer, dimension(4) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, dimension(MPI_STATUS_SIZE) :: status
        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0, i,j,k,l, ii, jj, kk

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
                    ii = ex-GN+i
                    bfr_s_E(ii,j,k,l) = var(ii,j,k,l)
                    ii = sx-1+i
                    bfr_s_W(ii,j,k,l) = var(ii,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc wait
        !$acc host_data use_device(bfr_s_E,bfr_r_W) &
        !$acc           use_device(bfr_s_W,bfr_r_E)

        !$acc wait(1)
        call MPI_SENDRECV(bfr_s_E,bfr_size_x,dtype,neigh(E),tag+1, &
                          bfr_r_W,bfr_size_x,dtype,neigh(W),tag+1, &
                          comm,status,err)

        call MPI_SENDRECV(bfr_s_W,bfr_size_x,dtype,neigh(W),tag+2, &
                          bfr_r_E,bfr_size_x,dtype,neigh(E),tag+2, &
                          comm,status,err)
        !$acc end host_data

        if(neigh(W)>=0) then
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = 1,GN
                    ii = sx-GN+i-1
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
                    ii = ex+i
                    var(ii,j,k,l) = bfr_r_E(ii,j,k,l)
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

        !$acc wait

        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN
                 do i = lbx,ubx
                    jj = ey-GN+j
                    bfr_s_N(i,jj,k,l) = var(i,jj,k,l)
                    jj = sy-1+j
                    bfr_s_S(i,jj,k,l) = var(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        !$acc wait
        !$acc host_data use_device(bfr_s_S,bfr_r_N) &
        !$acc           use_device(bfr_s_N,bfr_r_S)

        !$acc wait(2)
        call MPI_SENDRECV(bfr_s_S,bfr_size_y,dtype,neigh(S),tag+4, &
                          bfr_r_N,bfr_size_y,dtype,neigh(N),tag+4, &
                          comm,status,err)

        call MPI_SENDRECV(bfr_s_N,bfr_size_y,dtype,neigh(N),tag+3, &
                          bfr_r_S,bfr_size_y,dtype,neigh(S),tag+3, &
                          comm,status,err)
        !$acc end host_data

        if(neigh(S)>=0) then
        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = lbz,ubz
              do    j = 1,GN
                 do i = lbx,ubx
                    jj = sy-GN+j-1
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
                    jj = ey+j
                    var(i,jj,k,l) = bfr_r_N(i,jj,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif
        if(neigh(S)>=0) then
          !$acc wait(3)
        endif
        if(neigh(N)>=0) then
          !$acc wait(4)
        endif

        !$acc wait

        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = ez-GN+k
                    bfr_s_F(i,j,kk,l) = var(i,j,kk,l)
                    kk = sz-1+k
                    bfr_s_B(i,j,kk,l) = var(i,j,kk,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        !$acc wait
        !$acc host_data use_device(bfr_s_B,bfr_r_F) &
        !$acc           use_device(bfr_s_F,bfr_r_B)

        !$acc wait(3)
        call MPI_SENDRECV(bfr_s_B,bfr_size_z,dtype,neigh(B),tag+6, &
                          bfr_r_F,bfr_size_z,dtype,neigh(F),tag+6, &
                          comm,status,err)

        call MPI_SENDRECV(bfr_s_F,bfr_size_z,dtype,neigh(F),tag+5, &
                          bfr_r_B,bfr_size_z,dtype,neigh(B),tag+5, &
                          comm,status,err)
        !$acc end host_data

        if(neigh(B)>=0) then
        !$acc parallel default(present) async(5)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = 1,GN
              do    j = lby,uby
                 do i = lbx,ubx
                    kk = sz-GN+k-1
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
                    kk = ez+k
                    var(i,j,kk,l) = bfr_r_F(i,j,kk,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel
        endif
        if(neigh(B)>=0) then
          !$acc wait(5)
        endif
        if(neigh(F)>=0) then
          !$acc wait(6)
        endif

        !$acc wait

        return
end subroutine mpi_communicate4D_real8_sendrecv_cuda_aware

subroutine mpi_communicate3D_real8_sendrecv(comm,neigh,&
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
        integer                                 , intent(in)    :: comm

        integer, parameter    :: dtype = MPI_RP
        integer, dimension(3) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, dimension(MPI_STATUS_SIZE) :: status
        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0, i,j,k, ii, jj, kk

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
                 ii = ex-GN+i
                 bfr_s_E(ii,j,k) = var(ii,j,k)
                 ii = sx-1+i
                 bfr_s_W(ii,j,k) = var(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc wait(1)
        !$acc update host(bfr_s_E) async(4)
        !$acc update host(bfr_s_W) async(5)

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

        !$acc wait(10)
        if(neigh(W)>=0) then
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN
                 ii = sx-GN+i-1
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
                 ii = ex+i
                 var(ii,j,k) = bfr_r_E(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN
              do i = lbx,ubx
                 jj = ey-GN+j
                 bfr_s_N(i,jj,k) = var(i,jj,k)
                 jj = sy-1+j
                 bfr_s_S(i,jj,k) = var(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc wait(2)
        !$acc update host(bfr_s_S) async(6)
        !$acc update host(bfr_s_N) async(7)

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

        !$acc wait(12)
        if(neigh(N)>=0) then
        !$acc parallel default(present) async(4)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN
              do i = lbx,ubx
                 jj = ey+j
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
                 jj = sy-GN+j-1
                 var(i,jj,k) = bfr_r_S(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN
           do    j = lby,uby
              do i = lbx,ubx
                 kk = ez-GN+k
                 bfr_s_F(i,j,kk) = var(i,j,kk)
                 kk = sz-1+k
                 bfr_s_B(i,j,kk) = var(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc wait(3)
        !$acc update host(bfr_s_B) async(8)
        !$acc update host(bfr_s_F) async(9)

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

        !$acc wait(14)
        if(neigh(F)>=0) then
        !$acc parallel default(present) async(6)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN
           do    j = lby,uby
              do i = lbx,ubx
                 kk = ez+k
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
                 kk = sz-GN+k-1
                 var(i,j,kk) = bfr_r_B(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        return
end subroutine mpi_communicate3D_real8_sendrecv

subroutine mpi_communicate3D_real8_sendrecv_cuda_aware(comm,neigh,&
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
        integer                                 , intent(in)    :: comm

        integer, parameter    :: dtype = MPI_RP
        integer, dimension(3) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, dimension(MPI_STATUS_SIZE) :: status
        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0, i,j,k, ii, jj, kk

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
                 ii = ex-GN+i
                 bfr_s_E(ii,j,k) = var(ii,j,k)
                 ii = sx-1+i
                 bfr_s_W(ii,j,k) = var(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc wait
        !$acc host_data use_device(bfr_s_E,bfr_r_W) &
        !$acc           use_device(bfr_s_W,bfr_r_E)

        !$acc wait(1)
        call MPI_SENDRECV(bfr_s_E,bfr_size_x,dtype,neigh(E),tag+1, &
                          bfr_r_W,bfr_size_x,dtype,neigh(W),tag+1, &
                          comm,status,err)

        call MPI_SENDRECV(bfr_s_W,bfr_size_x,dtype,neigh(W),tag+2, &
                          bfr_r_E,bfr_size_x,dtype,neigh(E),tag+2, &
                          comm,status,err)
        !$acc end host_data

        if(neigh(W)>=0) then
        !$acc parallel default(present) async(1)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = 1,GN
                 ii = sx-GN+i-1
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
                 ii = ex+i
                 var(ii,j,k) = bfr_r_E(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN
              do i = lbx,ubx
                 jj = ey-GN+j
                 bfr_s_N(i,jj,k) = var(i,jj,k)
                 jj = sy-1+j
                 bfr_s_S(i,jj,k) = var(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        !$acc wait
        !$acc host_data use_device(bfr_s_S,bfr_r_N) &
        !$acc           use_device(bfr_s_N,bfr_r_S)

        !$acc wait(2)
        call MPI_SENDRECV(bfr_s_S,bfr_size_y,dtype,neigh(S),tag+4, &
                          bfr_r_N,bfr_size_y,dtype,neigh(N),tag+4, &
                          comm,status,err)

        call MPI_SENDRECV(bfr_s_N,bfr_size_y,dtype,neigh(N),tag+3, &
                          bfr_r_S,bfr_size_y,dtype,neigh(S),tag+3, &
                          comm,status,err)
        !$acc end host_data

        if(neigh(S)>=0) then
        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN
              do i = lbx,ubx
                 jj = sy-GN+j-1
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
                 jj = ey+j
                 var(i,jj,k) = bfr_r_N(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN
           do    j = lby,uby
              do i = lbx,ubx
                 kk = ez-GN+k
                 bfr_s_F(i,j,kk) = var(i,j,kk)
                 kk = sz-1+k
                 bfr_s_B(i,j,kk) = var(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel
        !$acc wait
        !$acc host_data use_device(bfr_s_B,bfr_r_F) &
        !$acc           use_device(bfr_s_F,bfr_r_B)

        !$acc wait(3)
        call MPI_SENDRECV(bfr_s_B,bfr_size_z,dtype,neigh(B),tag+6, &
                          bfr_r_F,bfr_size_z,dtype,neigh(F),tag+6, &
                          comm,status,err)

        call MPI_SENDRECV(bfr_s_F,bfr_size_z,dtype,neigh(F),tag+5, &
                          bfr_r_B,bfr_size_z,dtype,neigh(B),tag+5, &
                          comm,status,err)
        !$acc end host_data

        if(neigh(B)>=0) then
        !$acc parallel default(present) async(5)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN
           do    j = lby,uby
              do i = lbx,ubx
                 kk = sz-GN+k-1
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
                 kk = ez+k
                 var(i,j,kk) = bfr_r_F(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel
        endif

        !$acc wait

        return
end subroutine mpi_communicate3D_real8_sendrecv_cuda_aware

subroutine mpi_communicateInt3D_isend_irecv_ddt(comm,tp_send,tp_recv,neigh,var)
        implicit none
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: var
        integer                , dimension(6)    , intent(in)    :: tp_send
        integer                , dimension(6)    , intent(in)    :: tp_recv
        integer                , dimension(6)    , intent(in)    :: neigh
        integer                                  , intent(in)    :: comm

        integer, dimension(4) :: req_x
        integer, dimension(4) :: req_y
        integer, dimension(4) :: req_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0

        !$acc update host(var)

        call MPI_ISEND(var,1,tp_send(E),neigh(E),tag+1,comm,req_x(1),err)
        call MPI_IRECV(var,1,tp_recv(W),neigh(W),tag+1,comm,req_x(2),err)

        call MPI_ISEND(var,1,tp_send(W),neigh(W),tag+2,comm,req_x(3),err)
        call MPI_IRECV(var,1,tp_recv(E),neigh(E),tag+2,comm,req_x(4),err)

        call MPI_waitall(size(req_x), req_x, MPI_STATUSES_IGNORE, err)

        call MPI_ISEND(var,1,tp_send(N),neigh(N),tag+3,comm,req_y(1),err)
        call MPI_IRECV(var,1,tp_recv(S),neigh(S),tag+3,comm,req_y(2),err)

        call MPI_ISEND(var,1,tp_send(S),neigh(S),tag+4,comm,req_y(3),err)
        call MPI_IRECV(var,1,tp_recv(N),neigh(N),tag+4,comm,req_y(4),err)

        call MPI_waitall(size(req_y), req_y, MPI_STATUSES_IGNORE, err)
        call MPI_ISEND(var,1,tp_send(F),neigh(F),tag+5,comm,req_z(1),err)
        call MPI_IRECV(var,1,tp_recv(B),neigh(B),tag+5,comm,req_z(2),err)

        call MPI_ISEND(var,1,tp_send(B),neigh(B),tag+6,comm,req_z(3),err)
        call MPI_IRECV(var,1,tp_recv(F),neigh(F),tag+6,comm,req_z(4),err)

        call MPI_waitall(size(req_z), req_z, MPI_STATUSES_IGNORE, err)

        !$acc update device(var)

        return
end subroutine mpi_communicateInt3D_isend_irecv_ddt

subroutine mpi_communicateInt3D_isend_irecv_ddt_cuda_aware(comm,tp_send,tp_recv,neigh,var)
        implicit none
        integer(1), allocatable, dimension(:,:,:), intent(inout) :: var
        integer                , dimension(6)    , intent(in)    :: tp_send
        integer                , dimension(6)    , intent(in)    :: tp_recv
        integer                , dimension(6)    , intent(in)    :: neigh
        integer                                  , intent(in)    :: comm

        integer, dimension(4) :: req_x
        integer, dimension(4) :: req_y
        integer, dimension(4) :: req_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0

        !$acc host_data use_device(var)

        call MPI_ISEND(var,1,tp_send(E),neigh(E),tag+1,comm,req_x(1),err)
        call MPI_IRECV(var,1,tp_recv(W),neigh(W),tag+1,comm,req_x(2),err)

        call MPI_ISEND(var,1,tp_send(W),neigh(W),tag+2,comm,req_x(3),err)
        call MPI_IRECV(var,1,tp_recv(E),neigh(E),tag+2,comm,req_x(4),err)

        call MPI_waitall(size(req_x), req_x, MPI_STATUSES_IGNORE, err)

        call MPI_ISEND(var,1,tp_send(N),neigh(N),tag+3,comm,req_y(1),err)
        call MPI_IRECV(var,1,tp_recv(S),neigh(S),tag+3,comm,req_y(2),err)

        call MPI_ISEND(var,1,tp_send(S),neigh(S),tag+4,comm,req_y(3),err)
        call MPI_IRECV(var,1,tp_recv(N),neigh(N),tag+4,comm,req_y(4),err)

        call MPI_waitall(size(req_y), req_y, MPI_STATUSES_IGNORE, err)
        call MPI_ISEND(var,1,tp_send(F),neigh(F),tag+5,comm,req_z(1),err)
        call MPI_IRECV(var,1,tp_recv(B),neigh(B),tag+5,comm,req_z(2),err)

        call MPI_ISEND(var,1,tp_send(B),neigh(B),tag+6,comm,req_z(3),err)
        call MPI_IRECV(var,1,tp_recv(F),neigh(F),tag+6,comm,req_z(4),err)

        call MPI_waitall(size(req_z), req_z, MPI_STATUSES_IGNORE, err)

        !$acc end host_data

        return
end subroutine mpi_communicateInt3D_isend_irecv_ddt_cuda_aware

subroutine mpi_communicateInt3D_sendrecv(comm,neigh,&
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
        integer                                   , intent(in)    :: comm

        integer, parameter    :: dtype = MPI_INTEGER1
        integer, dimension(3) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, dimension(MPI_STATUS_SIZE) :: status
        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0, i,j,k, ii, jj, kk

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
                 ii = ex-GN+i
                 bfr_s_E(ii,j,k) = var(ii,j,k)
                 ii = sx-1+i
                 bfr_s_W(ii,j,k) = var(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc wait(1)
        !$acc update host(bfr_s_E) async(4)
        !$acc update host(bfr_s_W) async(5)

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

        !$acc wait(10)
        if(neigh(W)>=0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN
                   ii = sx-GN+i-1
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
                   ii = ex+i
                   var(ii,j,k) = bfr_r_E(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait

        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN
              do i = lbx,ubx
                 jj = ey-GN+j
                 bfr_s_N(i,jj,k) = var(i,jj,k)
                 jj = sy-1+j
                 bfr_s_S(i,jj,k) = var(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc wait(2)
        !$acc update host(bfr_s_S) async(6)
        !$acc update host(bfr_s_N) async(7)

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

        !$acc wait(12)
        if(neigh(N)>=0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN
                do i = lbx,ubx
                   jj = ey+j
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
                   jj = sy-GN+j-1
                   var(i,jj,k) = bfr_r_S(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait

        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN
           do    j = lby,uby
              do i = lbx,ubx
                 kk = ez-GN+k
                 bfr_s_F(i,j,kk) = var(i,j,kk)
                 kk = sz-1+k
                 bfr_s_B(i,j,kk) = var(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc wait(3)
        !$acc update host(bfr_s_B) async(8)
        !$acc update host(bfr_s_F) async(9)

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

        !$acc wait(14)
        if(neigh(F)>=0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = ez+k
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
                   kk = sz-GN+k-1
                   var(i,j,kk) = bfr_r_B(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait

        return
end subroutine mpi_communicateInt3D_sendrecv

subroutine mpi_communicateInt3D_sendrecv_cuda_aware(comm,neigh,&
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
        integer                                   , intent(in)    :: comm

        integer, parameter    :: dtype = MPI_INTEGER1
        integer, dimension(3) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, dimension(MPI_STATUS_SIZE) :: status
        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: err = 0, i,j,k, ii, jj, kk

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
                 ii = ex-GN+i
                 bfr_s_E(ii,j,k) = var(ii,j,k)
                 ii = sx-1+i
                 bfr_s_W(ii,j,k) = var(ii,j,k)
              enddo
           enddo
        enddo
        !$acc end parallel

        !$acc wait
        !$acc host_data use_device(bfr_s_E,bfr_r_W) &
        !$acc           use_device(bfr_s_W,bfr_r_E)

        !$acc wait(1)
        call MPI_SENDRECV(bfr_s_E,bfr_size_x,dtype,neigh(E),tag+1, &
                          bfr_r_W,bfr_size_x,dtype,neigh(W),tag+1, &
                          comm,status,err)

        call MPI_SENDRECV(bfr_s_W,bfr_size_x,dtype,neigh(W),tag+2, &
                          bfr_r_E,bfr_size_x,dtype,neigh(E),tag+2, &
                          comm,status,err)
        !$acc end host_data

        if(neigh(W)>=0) then
          !$acc parallel default(present) async(1)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN
                   ii = sx-GN+i-1
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
                   ii = ex+i
                   var(ii,j,k) = bfr_r_E(ii,j,k)
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

        !$acc wait

        !$acc parallel default(present) async(2)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = 1,GN
              do i = lbx,ubx
                 jj = ey-GN+j
                 bfr_s_N(i,jj,k) = var(i,jj,k)
                 jj = sy-1+j
                 bfr_s_S(i,jj,k) = var(i,jj,k)
              enddo
           enddo
        enddo
        !$acc end parallel
        !$acc wait
        !$acc host_data use_device(bfr_s_S,bfr_r_N) &
        !$acc           use_device(bfr_s_N,bfr_r_S)

        !$acc wait(2)
        call MPI_SENDRECV(bfr_s_S,bfr_size_y,dtype,neigh(S),tag+4, &
                          bfr_r_N,bfr_size_y,dtype,neigh(N),tag+4, &
                          comm,status,err)

        call MPI_SENDRECV(bfr_s_N,bfr_size_y,dtype,neigh(N),tag+3, &
                          bfr_r_S,bfr_size_y,dtype,neigh(S),tag+3, &
                          comm,status,err)
        !$acc end host_data

        if(neigh(N)>=0) then
          !$acc parallel default(present) async(3)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN
                do i = lbx,ubx
                   jj = ey+j
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
                   jj = sy-GN+j-1
                   var(i,jj,k) = bfr_r_S(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif
        if(neigh(N)>=0) then
          !$acc wait(3)
        endif
        if(neigh(S)>=0) then
          !$acc wait(4)
        endif

        !$acc wait

        !$acc parallel default(present) async(3)
        !$acc loop gang, vector collapse(3)
        do       k = 1,GN
           do    j = lby,uby
              do i = lbx,ubx
                 kk = ez-GN+k
                 bfr_s_F(i,j,kk) = var(i,j,kk)
                 kk = sz-1+k
                 bfr_s_B(i,j,kk) = var(i,j,kk)
              enddo
           enddo
        enddo
        !$acc end parallel
        !$acc wait
        !$acc host_data use_device(bfr_s_B,bfr_r_F) &
        !$acc           use_device(bfr_s_F,bfr_r_B)

        !$acc wait(3)
        call MPI_SENDRECV(bfr_s_B,bfr_size_z,dtype,neigh(B),tag+6, &
                          bfr_r_F,bfr_size_z,dtype,neigh(F),tag+6, &
                          comm,status,err)

        call MPI_SENDRECV(bfr_s_F,bfr_size_z,dtype,neigh(F),tag+5, &
                          bfr_r_B,bfr_size_z,dtype,neigh(B),tag+5, &
                          comm,status,err)
        !$acc end host_data

        if(neigh(F)>=0) then
          !$acc parallel default(present) async(5)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = ez+k
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
                   kk = sz-GN+k-1
                   var(i,j,kk) = bfr_r_B(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif
        if(neigh(F)>=0) then
          !$acc wait(5)
        endif
        if(neigh(B)>=0) then
          !$acc wait(6)
        endif

        !$acc wait

        return
end subroutine mpi_communicateInt3D_sendrecv_cuda_aware

subroutine mpi_communicateInt3D_isend_irecv_bfr(comm,neigh,&
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
        integer                                   , intent(in)    :: comm

        integer, parameter    :: dtype = MPI_INTEGER1
        integer, dimension(3) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer, dimension(4) :: req_x
        integer, dimension(4) :: req_y
        integer, dimension(4) :: req_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: i,j,k,err = 0, ii, jj, kk, q

        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)

        if(neigh(E) >= 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN
                   ii = ex-GN+i
                   bfr_s_E(ii,j,k) = var(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
          !$acc update host(bfr_s_E)
        endif

        if(neigh(W) >= 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN
                   ii = sx-1+i
                   bfr_s_W(ii,j,k) = var(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
          !$acc update host(bfr_s_W)
        endif

        q = 0
        if(neigh(E) >= 0) then
        call MPI_IRECV(bfr_r_E,bfr_size_x,dtype,neigh(E),tag+1,comm,req_x(q+1),err)
        call MPI_ISEND(bfr_s_E,bfr_size_x,dtype,neigh(E),tag+2,comm,req_x(q+2),err)
        q = q+2
        endif

        if(neigh(W) >= 0) then
        call MPI_IRECV(bfr_r_W,bfr_size_x,dtype,neigh(W),tag+2,comm,req_x(q+1),err)
        call MPI_ISEND(bfr_s_W,bfr_size_x,dtype,neigh(W),tag+1,comm,req_x(q+2),err)
        q = q+2
        endif

        if(q > 0) then
          call MPI_waitall(q, req_x, MPI_STATUSES_IGNORE, err)
          !$acc update device(bfr_r_W,bfr_r_E)
        endif

        if(neigh(W)>=0) then
          !$acc parallel default(present) async(1)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN
                   ii = sx-GN+i-1
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
                   ii = ex+i
                   var(ii,j,k) = bfr_r_E(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait

        if(neigh(N) >= 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN
                do i = lbx,ubx
                   jj = ey-GN+j
                   bfr_s_N(i,jj,k) = var(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
          !$acc update host(bfr_s_N)
        endif

        if(neigh(S) >= 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN
                do i = lbx,ubx
                   jj = sy-1+j
                   bfr_s_S(i,jj,k) = var(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
          !$acc update host(bfr_s_S)
        endif

        q = 0
        if(neigh(N) >= 0) then
        call MPI_IRECV(bfr_r_N,bfr_size_y,dtype,neigh(N),tag+1,comm,req_y(q+1),err)
        call MPI_ISEND(bfr_s_N,bfr_size_y,dtype,neigh(N),tag+2,comm,req_y(q+2),err)
        q = q+2
        endif

        if(neigh(S) >= 0) then
        call MPI_IRECV(bfr_r_S,bfr_size_y,dtype,neigh(S),tag+2,comm,req_y(q+1),err)
        call MPI_ISEND(bfr_s_S,bfr_size_y,dtype,neigh(S),tag+1,comm,req_y(q+2),err)
        q = q+2
        endif

        if(q > 0) then
          call MPI_waitall(q, req_y, MPI_STATUSES_IGNORE, err)
          !$acc update device(bfr_r_N,bfr_r_S)
        endif

        if(neigh(N)>=0) then
          !$acc parallel default(present) async(3)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN
                do i = lbx,ubx
                   jj = ey+j
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
                   jj = sy-GN+j-1
                   var(i,jj,k) = bfr_r_S(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait

        if(neigh(F) >= 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = ez-GN+k
                   bfr_s_F(i,j,kk) = var(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel
          !$acc update host(bfr_s_F)
        endif

        if(neigh(B) >= 0) then
          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = sz-1+k
                   bfr_s_B(i,j,kk) = var(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel
          !$acc update host(bfr_s_B)
        endif

        q = 0
        if(neigh(F) >= 0) then
        call MPI_IRECV(bfr_r_F,bfr_size_z,dtype,neigh(F),tag+1,comm,req_z(q+1),err)
        call MPI_ISEND(bfr_s_F,bfr_size_z,dtype,neigh(F),tag+2,comm,req_z(q+2),err)
        q = q+2
        endif

        if(neigh(B) >= 0) then
        call MPI_IRECV(bfr_r_B,bfr_size_z,dtype,neigh(B),tag+2,comm,req_z(q+1),err)
        call MPI_ISEND(bfr_s_B,bfr_size_z,dtype,neigh(B),tag+1,comm,req_z(q+2),err)
        q = q+2
        endif

        if(q > 0) then
          call MPI_waitall(q, req_z, MPI_STATUSES_IGNORE, err)
          !$acc update device(bfr_r_F,bfr_r_B)
        endif

        if(neigh(F)>=0) then
          !$acc parallel default(present) async(5)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = ez+k
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
                   kk = sz-GN+k-1
                   var(i,j,kk) = bfr_r_B(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait

        return
end subroutine mpi_communicateInt3D_isend_irecv_bfr

subroutine mpi_communicateInt3D_isend_irecv_bfr_cuda_aware(comm,neigh     ,&
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
        integer                                   , intent(in)    :: comm

        integer, parameter    :: dtype = MPI_INTEGER1
        integer, dimension(3) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer, dimension(4) :: req_x
        integer, dimension(4) :: req_y
        integer, dimension(4) :: req_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: i,j,k,err = 0, ii, jj, kk, q

        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)

        if(neigh(E) >= 0) then
          !$acc parallel default(present) async(1)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN
                   ii = ex-GN+i
                   bfr_s_E(ii,j,k) = var(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(W) >= 0) then
          !$acc parallel default(present) async(2)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN
                   ii = sx-1+i
                   bfr_s_W(ii,j,k) = var(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait
        !$acc host_data use_device(bfr_s_E,bfr_r_W) &
        !$acc           use_device(bfr_s_W,bfr_r_E)

        q = 0
        if(neigh(W) >= 0) then
        call MPI_IRECV(bfr_r_W,bfr_size_x,dtype,neigh(W),tag+1,comm,req_x(q+1),err)
        call MPI_ISEND(bfr_s_W,bfr_size_x,dtype,neigh(W),tag+2,comm,req_x(q+2),err)
        q = q+2
        endif

        if(neigh(E) >= 0) then
        call MPI_IRECV(bfr_r_E,bfr_size_x,dtype,neigh(E),tag+2,comm,req_x(q+1),err)
        call MPI_ISEND(bfr_s_E,bfr_size_x,dtype,neigh(E),tag+1,comm,req_x(q+2),err)
        q = q+2
        endif

        if(q > 0) call MPI_waitall(q, req_x, MPI_STATUSES_IGNORE, err)
        !$acc end host_data

        if(neigh(W)>=0) then
          !$acc parallel default(present) async(1)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN
                   ii = sx-GN+i-1
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
                   ii = ex+i
                   var(ii,j,k) = bfr_r_E(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait

        if(neigh(N) >= 0) then
          !$acc parallel default(present) async(3)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN
                do i = lbx,ubx
                   jj = ey-GN+j
                   bfr_s_N(i,jj,k) = var(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(S) >= 0) then
          !$acc parallel default(present) async(4)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN
                do i = lbx,ubx
                   jj = sy-1+j
                   bfr_s_S(i,jj,k) = var(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif
        !$acc wait
        !$acc host_data use_device(bfr_s_S,bfr_r_N) &
        !$acc           use_device(bfr_s_N,bfr_r_S)

        q = 0
        if(neigh(S) >= 0) then
        call MPI_IRECV(bfr_r_S,bfr_size_y,dtype,neigh(S),tag+1,comm,req_y(q+1),err)
        call MPI_ISEND(bfr_s_S,bfr_size_y,dtype,neigh(S),tag+2,comm,req_y(q+2),err)
        q = q+2
        endif

        if(neigh(N) >= 0) then
        call MPI_IRECV(bfr_r_N,bfr_size_y,dtype,neigh(N),tag+2,comm,req_y(q+1),err)
        call MPI_ISEND(bfr_s_N,bfr_size_y,dtype,neigh(N),tag+1,comm,req_y(q+2),err)
        q = q+2
        endif

        if(q > 0) call MPI_waitall(q, req_y, MPI_STATUSES_IGNORE, err)
        !$acc end host_data

        if(neigh(N)>=0) then
          !$acc parallel default(present) async(3)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN
                do i = lbx,ubx
                   jj = ey+j
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
                   jj = sy-GN+j-1
                   var(i,jj,k) = bfr_r_S(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait

        if(neigh(F) >= 0) then
          !$acc parallel default(present) async(5)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = ez-GN+k
                   bfr_s_F(i,j,kk) = var(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(B) >= 0) then
          !$acc parallel default(present) async(6)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = sz-1+k
                   bfr_s_B(i,j,kk) = var(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif
        !$acc wait
        !$acc host_data use_device(bfr_s_B,bfr_r_F) &
        !$acc           use_device(bfr_s_F,bfr_r_B)

        q = 0
        if(neigh(B) >= 0) then
        call MPI_IRECV(bfr_r_B,bfr_size_z,dtype,neigh(B),tag+1,comm,req_z(q+1),err)
        call MPI_ISEND(bfr_s_B,bfr_size_z,dtype,neigh(B),tag+2,comm,req_z(q+2),err)
        q = q+2
        endif

        if(neigh(F) >= 0) then
        call MPI_IRECV(bfr_r_F,bfr_size_z,dtype,neigh(F),tag+2,comm,req_z(q+1),err)
        call MPI_ISEND(bfr_s_F,bfr_size_z,dtype,neigh(F),tag+1,comm,req_z(q+2),err)
        q = q+2
        endif

        if(q > 0) call MPI_waitall(q, req_z, MPI_STATUSES_IGNORE, err)
        !$acc end host_data

        if(neigh(F)>=0) then
          !$acc parallel default(present) async(5)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = ez+k
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
                   kk = sz-GN+k-1
                   var(i,j,kk) = bfr_r_B(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait

        return
end subroutine mpi_communicateInt3D_isend_irecv_bfr_cuda_aware

subroutine mpi_communicateInt3D_bfr_one_sided(comm,neigh     ,&
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
        integer                                   , intent(in)    :: comm

        integer, parameter    :: dtype = MPI_INTEGER1
        integer, dimension(3) :: bfr_dims_x, bfr_dims_y, bfr_dims_z
        integer, dimension(4) :: req_x
        integer, dimension(4) :: req_y
        integer, dimension(4) :: req_z
        integer               :: bfr_size_x, bfr_size_y, bfr_size_z

        integer, parameter :: tag = 0
        integer, parameter :: W = 1, E = 2, S = 3, N = 4, B = 5, F = 6
        integer            :: i,j,k,err = 0, ii, jj, kk, q

        integer(KIND=MPI_ADDRESS_KIND), parameter :: target_disp = 0

        bfr_dims_x = shape(bfr_s_E)
        bfr_dims_y = shape(bfr_s_N)
        bfr_dims_z = shape(bfr_s_F)

        bfr_size_x = bfr_dims_x(1)*bfr_dims_x(2)*bfr_dims_x(3)
        bfr_size_y = bfr_dims_y(1)*bfr_dims_y(2)*bfr_dims_y(3)
        bfr_size_z = bfr_dims_z(1)*bfr_dims_z(2)*bfr_dims_z(3)

        if(neigh(E) >= 0) then
          !$acc parallel default(present) async(1)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN
                   ii = ex-GN+i
                   bfr_s_E(ii,j,k) = var(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(W) >= 0) then
          !$acc parallel default(present) async(2)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN
                   ii = sx-1+i
                   bfr_s_W(ii,j,k) = var(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait
        !$acc host_data use_device(bfr_s_E,bfr_r_W) &
        !$acc           use_device(bfr_s_W,bfr_r_E)

        if (neigh(W) >= 0) call MPI_WIN_POST(mpi_group_W, 0, win_i13D_bfr_W, err)
        if (neigh(E) >= 0) call MPI_WIN_POST(mpi_group_E, 0, win_i13D_bfr_E, err)

        if (neigh(W) >= 0) then
         call MPI_WIN_START(mpi_group_W, 0, win_i13D_bfr_E, err)
         call MPI_PUT(bfr_s_W,bfr_size_x,dtype,neigh(W),target_disp,bfr_size_x,dtype,win_i13D_bfr_E,err)
         !$acc wait
         call MPI_WIN_COMPLETE(win_i13D_bfr_E, err)
        endif
        if (neigh(E) >= 0) then
         call MPI_WIN_START(mpi_group_E, 0, win_i13D_bfr_W, err)
         call MPI_PUT(bfr_s_E,bfr_size_x,dtype,neigh(E),target_disp,bfr_size_x,dtype,win_i13D_bfr_W,err)
         !$acc wait
         call MPI_WIN_COMPLETE(win_i13D_bfr_W, err)
        endif

        if (neigh(W) >= 0) call MPI_WIN_WAIT(win_i13D_bfr_W, err)
        if (neigh(E) >= 0) call MPI_WIN_WAIT(win_i13D_bfr_E, err)
        !$acc end host_data

        if(neigh(W)>=0) then
          !$acc parallel default(present) async(1)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = lby,uby
                do i = 1,GN
                   ii = sx-GN+i-1
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
                   ii = ex+i
                   var(ii,j,k) = bfr_r_E(ii,j,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait

        if(neigh(N) >= 0) then
          !$acc parallel default(present) async(3)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN
                do i = lbx,ubx
                   jj = ey-GN+j
                   bfr_s_N(i,jj,k) = var(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(S) >= 0) then
          !$acc parallel default(present) async(4)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN
                do i = lbx,ubx
                   jj = sy-1+j
                   bfr_s_S(i,jj,k) = var(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif
        !$acc wait
        !$acc host_data use_device(bfr_s_S,bfr_r_N) &
        !$acc           use_device(bfr_s_N,bfr_r_S)

        if (neigh(S) >= 0) call MPI_WIN_POST(mpi_group_S, 0, win_i13D_bfr_S, err)
        if (neigh(N) >= 0) call MPI_WIN_POST(mpi_group_N, 0, win_i13D_bfr_N, err)

        if (neigh(S) >= 0) then
         call MPI_WIN_START(mpi_group_S, 0, win_i13D_bfr_N, err)
         call MPI_PUT(bfr_s_S,bfr_size_y,dtype,neigh(S),target_disp,bfr_size_y,dtype,win_i13D_bfr_N,err)
         !$acc wait
         call MPI_WIN_COMPLETE(win_i13D_bfr_N, err)
        endif
        if (neigh(N) >= 0) then
         call MPI_WIN_START(mpi_group_N, 0, win_i13D_bfr_S, err)
         call MPI_PUT(bfr_s_N,bfr_size_y,dtype,neigh(N),target_disp,bfr_size_y,dtype,win_i13D_bfr_S,err)
         !$acc wait
         call MPI_WIN_COMPLETE(win_i13D_bfr_S, err)
        endif

        if (neigh(S) >= 0) call MPI_WIN_WAIT(win_i13D_bfr_S, err)
        if (neigh(N) >= 0) call MPI_WIN_WAIT(win_i13D_bfr_N, err)
        !$acc end host_data

        if(neigh(N)>=0) then
          !$acc parallel default(present) async(3)
          !$acc loop gang, vector collapse(3)
          do       k = lbz,ubz
             do    j = 1,GN
                do i = lbx,ubx
                   jj = ey+j
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
                   jj = sy-GN+j-1
                   var(i,jj,k) = bfr_r_S(i,jj,k)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait

        if(neigh(F) >= 0) then
          !$acc parallel default(present) async(5)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = ez-GN+k
                   bfr_s_F(i,j,kk) = var(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        if(neigh(B) >= 0) then
          !$acc parallel default(present) async(6)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = sz-1+k
                   bfr_s_B(i,j,kk) = var(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif
        !$acc wait
        !$acc host_data use_device(bfr_s_B,bfr_r_F) &
        !$acc           use_device(bfr_s_F,bfr_r_B)

        if (neigh(B) >= 0) call MPI_WIN_POST(mpi_group_B, 0, win_i13D_bfr_B, err)
        if (neigh(F) >= 0) call MPI_WIN_POST(mpi_group_F, 0, win_i13D_bfr_F, err)

        if (neigh(B) >= 0) then
         call MPI_WIN_START(mpi_group_B, 0, win_i13D_bfr_F, err)
         call MPI_PUT(bfr_s_B,bfr_size_z,dtype,neigh(B),target_disp,bfr_size_z,dtype,win_i13D_bfr_F,err)
         !$acc wait
         call MPI_WIN_COMPLETE(win_i13D_bfr_F, err)
        endif
        if (neigh(F) >= 0) then
         call MPI_WIN_START(mpi_group_F, 0, win_i13D_bfr_B, err)
         call MPI_PUT(bfr_s_F,bfr_size_z,dtype,neigh(F),target_disp,bfr_size_z,dtype,win_i13D_bfr_B,err)
         !$acc wait
         call MPI_WIN_COMPLETE(win_i13D_bfr_B, err)
        endif

        if (neigh(B) >= 0) call MPI_WIN_WAIT(win_i13D_bfr_B, err)
        if (neigh(F) >= 0) call MPI_WIN_WAIT(win_i13D_bfr_F, err)
        !$acc end host_data

        if(neigh(F)>=0) then
          !$acc parallel default(present) async(5)
          !$acc loop gang, vector collapse(3)
          do       k = 1,GN
             do    j = lby,uby
                do i = lbx,ubx
                   kk = ez+k
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
                   kk = sz-GN+k-1
                   var(i,j,kk) = bfr_r_B(i,j,kk)
                enddo
             enddo
          enddo
          !$acc end parallel
        endif

        !$acc wait

        return
end subroutine mpi_communicateInt3D_bfr_one_sided

end module mpi_comm_module
