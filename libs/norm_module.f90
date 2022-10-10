! -------------------------------------------------------------------------        
!       This module contains useful routines in order to compute norm of 
!       field of various dimensions both serial and mpi.
! -------------------------------------------------------------------------        
module norm_module
use parameters_module, only: rp
use mpi

implicit none
private
public p_norm, inf_norm, mean

interface p_norm
  module procedure p_norm1D, p_norm2D, p_norm3D
end interface p_norm

interface inf_norm
  module procedure inf_norm1D, inf_norm2D, inf_norm3D
end interface inf_norm

interface mean
  module procedure mean_1D, mean_3D
endinterface mean

contains
function p_norm1D(A,p,lb,ub,mpi) result(p_norm)
! ----------------------------------------------------------------
!       Computation of the p norm of a 3D field.
! ----------------------------------------------------------------
        implicit none
        real(rp), allocatable, dimension(:), intent(in) :: A
        real(rp)                           , intent(in) :: p
        integer                            , intent(in) :: lb, ub
        logical , optional                 , intent(in) :: mpi
        real(rp)                                        :: p_norm
        
        ! local declarations
        real(rp) :: lcl_sum, gbl_sum
        integer  :: lcl_pts, gbl_pts
        integer  :: err = 0
        
        lcl_pts = (ub-lb+1)
        
        lcl_sum = sum(abs(A(lb:ub))**p)

        ! MPI
        if(present(mpi)) then
          call MPI_allreduce(lcl_sum, gbl_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, err)
          call MPI_allreduce(lcl_pts, gbl_pts, 1, MPI_INTEGER         , MPI_SUM, mpi_comm_world, err)
        else
          gbl_sum = lcl_sum
          gbl_pts = lcl_pts
        endif

        p_norm = (gbl_sum/gbl_pts)**(1.0_rp/p)

        return
end function p_norm1D



function inf_norm1D(A,lb,ub,mpi) result(inf_norm)
! ----------------------------------------------------------------
!       Computation of the infinite norm of a 1D field.
! ----------------------------------------------------------------
        implicit none
        real(rp), allocatable, dimension(:), intent(in) :: A
        integer , dimension(1)             , intent(in) :: lb, ub
        logical , optional                 , intent(in) :: mpi 
        real(rp)                                        :: inf_norm

        ! local declarations
        real(rp) :: lcl_inf_norm
        integer  :: err = 0
        
        ! compute local inf_norm
        lcl_inf_norm = maxval(abs(A(lb(1):ub(1))))

        ! reduce MPI
        if(present(mpi)) then
          call MPI_allreduce(lcl_inf_norm, inf_norm, 1, MPI_DOUBLE_PRECISION, MPI_MAX, mpi_comm_world, err)
        else
          inf_norm = lcl_inf_norm
        endif
        
        return
end function inf_norm1D



function p_norm2D(A,p,lb,ub,mpi) result(p_norm)
! ----------------------------------------------------------------
!       Computation of the p norm of a 2D field.
! ----------------------------------------------------------------
        implicit none
        real(rp), allocatable, dimension(:,:), intent(in) :: A
        real(rp)                             , intent(in) :: p
        integer              , dimension(2)  , intent(in) :: lb, ub 
        logical, optional                    , intent(in) :: mpi
        real(rp)                                          :: p_norm

        ! local declarations
        real(rp) :: lcl_sum, gbl_sum
        integer  :: lcl_pts, gbl_pts
        integer  :: err = 0

        lcl_pts = (ub(1)-lb(1)+1)*(ub(2)-lb(2)+1)
        
        lcl_sum = sum(abs(A(lb(1):ub(1), lb(2):ub(2)))**p)

        ! MPI
        if(present(mpi)) then
          call MPI_allreduce(lcl_sum, gbl_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, err)
          call MPI_allreduce(lcl_pts, gbl_pts, 1, MPI_INTEGER         , MPI_SUM, mpi_comm_world, err)
        else
          gbl_sum = lcl_sum
          gbl_pts = lcl_pts
        endif

        p_norm = (gbl_sum/gbl_pts)**(1.0_rp/p)

        return
end function p_norm2D



function inf_norm2D(A,lb,ub,mpi) result(inf_norm)
! ----------------------------------------------------------------
!       Computation of the infinite norm of a 2D field.
! ----------------------------------------------------------------
        implicit none
        real(rp), allocatable, dimension(:,:), intent(in) :: A
        integer              , dimension(2)  , intent(in) :: lb, ub
        logical, optional                    , intent(in) :: mpi
        real(rp)                                          :: inf_norm

        ! local declarations
        real(rp) :: lcl_inf_norm
        integer  :: err = 0
        
        lcl_inf_norm = maxval(abs(A(lb(1):ub(1), lb(2):ub(2))))
        
        ! MPI
        if(present(mpi)) then
          call MPI_allreduce(lcl_inf_norm, inf_norm, 1, MPI_DOUBLE_PRECISION, MPI_MAX, mpi_comm_world, err)
        else
          inf_norm = lcl_inf_norm
        endif

        return
end function inf_norm2D



function p_norm3D(A,p,lb,ub,mpi) result(p_norm)
! ----------------------------------------------------------------
!       Computation of the p norm of a 3D field.
! ----------------------------------------------------------------
        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(in) :: A
        real(rp)                               , intent(in) :: p
        integer , dimension(3)                 , intent(in) :: lb, ub
        logical , optional                     , intent(in) :: mpi
        real(rp)                                            :: p_norm
        
        ! local declarations
        real(rp) :: lcl_sum, gbl_sum
        integer  :: lcl_pts, gbl_pts
        integer  :: err = 0

        lcl_pts = (ub(1)-lb(1)+1)*(ub(2)-lb(2)+1)*(ub(3)-lb(3)+1)

        lcl_sum = sum(abs(A(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)))**p)

        ! MPI
        if(present(mpi)) then
          call MPI_allreduce(lcl_sum, gbl_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, err)
          call MPI_allreduce(lcl_pts, gbl_pts, 1, MPI_INTEGER         , MPI_SUM, mpi_comm_world, err)
        else
          gbl_sum = lcl_sum
          gbl_pts = lcl_pts
        endif

        p_norm = (gbl_sum/gbl_pts)**(1.0_rp/p)

        return
end function p_norm3D



function inf_norm3D(A,lb,ub,mpi) result(inf_norm)
! ----------------------------------------------------------------
!       Computation of the infinite norm of a 3D field.
! ----------------------------------------------------------------
        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(in) :: A
        integer , dimension(3)                 , intent(in) :: lb, ub
        logical , optional                     , intent(in) :: mpi 
        real(rp)                                            :: inf_norm

        ! local declarations
        real(rp) :: lcl_inf_norm
        integer  :: err = 0
        
        lcl_inf_norm = maxval(abs(A(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3))))

        ! MPI
        if(present(mpi)) then
          call MPI_allreduce(lcl_inf_norm, inf_norm, 1, MPI_DOUBLE_PRECISION, MPI_MAX, mpi_comm_world, err)
        else
          inf_norm = lcl_inf_norm
        endif

        return
end function inf_norm3D


function mean_1D(A,lb,ub,mpi) result(a_mean)
! ----------------------------------------------------------------
!       Computation of the mean of a 1D field A
! ----------------------------------------------------------------
        implicit none
        real(rp), allocatable, dimension(:), intent(in) :: A
        integer                            , intent(in) :: lb, ub
        logical, optional                  , intent(in) :: mpi
        real(rp)                                        :: a_mean

        ! local declarations
        real(rp) :: lcl_sum, gbl_sum
        integer  :: lcl_npt, gbl_npt
        integer  :: err = 0

        lcl_sum = sum(a(lb:ub))

        lcl_npt = ub - lb + 1

        if(present(mpi)) then

          call MPI_allreduce(lcl_sum, gbl_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, err)
          call MPI_allreduce(lcl_npt, gbl_npt, 1, MPI_INTEGER         , MPI_SUM, mpi_comm_world, err)

        else

          gbl_sum = lcl_sum
          gbl_npt = lcl_npt

        endif

        a_mean = gbl_sum/real(gbl_npt,rp)

        return
end function mean_1D



function mean_3D(A,lb,ub,mpi) result(a_mean)
! ----------------------------------------------------------------
!       Computation of the mean of a 1D field A
! ----------------------------------------------------------------
        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(in) :: A
        integer              , dimension(3)    , intent(in) :: lb, ub
        logical, optional                      , intent(in) :: mpi
        real(rp)                                            :: a_mean

        ! local declarations
        real(rp) :: lcl_sum, gbl_sum
        integer  :: lcl_npt, gbl_npt
        integer  :: err = 0

        lcl_sum = sum(a(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)))

        lcl_npt = (ub(1) - lb(1) + 1)*(ub(2) - lb(2) + 1)*(ub(3) - lb(3) + 1)

        if(present(mpi)) then

          call MPI_allreduce(lcl_sum, gbl_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, err)
          call MPI_allreduce(lcl_npt, gbl_npt, 1, MPI_INTEGER         , MPI_SUM, mpi_comm_world, err)

        else

          gbl_sum = lcl_sum
          gbl_npt = lcl_npt

        endif

        a_mean = gbl_sum/real(gbl_npt,rp)

        return
end function mean_3D











end module norm_module
