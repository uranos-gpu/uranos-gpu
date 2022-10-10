module integration_module
! -------------------------------------------------------------
!
!       This module produces a standalone interface in order to 
!       integrate fields of general shape.
!       
!       USAGE:
!       call integrate field(arg1, arg2, ...argN) 
!       
!       does everthing you need both in serial and in parallel!
!
! -------------------------------------------------------------
  use parameters_module, only: rp
  use mpi
  !$ use omp_lib
  
  implicit none
  private
  public integrate_field, mean_field, rms_field
  
  interface integrate_field
    module procedure int_1Dfield, int_2Dfield, int_3Dfield, int_4Dfield
  end interface

  interface mean_field
    module procedure mean_3Dfield, mean_4Dfield
  end interface

  interface rms_field
    module procedure rms_3Dfield
  end interface

contains
subroutine int_1Dfield(field,sx,ex,dx,int,mpi) 
! -----------------------------------------------------------
!       This function compute the integral of a 1D field
! -----------------------------------------------------------
        implicit none
        real(rp), dimension(:)    , allocatable, intent(in)  :: field       !< field we want to integrate
        real(rp), dimension(:)    , allocatable, intent(in)  :: dx          !< space coordinates
        integer                                , intent(in)  :: sx,ex       !< integer init and end positions
        logical , optional                     , intent(in)  :: mpi         !< handle for mpi calculation
        real(rp)                               , intent(out) :: Int         !< integral value

        ! local declaration
        real(rp) :: my_int
        integer  :: i, err = 0

        my_int = 0.0_rp

        !$omp parallel do default(private) &
        !$omp shared(sx,ex,field,dx) &
        !$omp reduction(+:my_int)
        do i = sx,ex

           my_int = my_int + field(i) * dx(i)

        enddo
        !$omp end parallel do


        if(present(mpi)) then
           call MPI_allreduce(my_int, int, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, err)
        else
           int = my_int
        endif


        return
end subroutine int_1Dfield



subroutine int_2Dfield(field,sx,ex,dx,sy,ey,dy,int,mpi) 
! -----------------------------------------------------------
!       This function compute the integral of a 2D field
! -----------------------------------------------------------
        implicit none
        real(rp), dimension(:,:)  , allocatable, intent(in)  :: field       !< field we want to integrate
        real(rp), dimension(:)    , allocatable, intent(in)  :: dx, dy      !< space coordinates
        integer                                , intent(in)  :: sx,ex,sy,ey !< integer init and end positions
        logical , optional                     , intent(in)  :: mpi         !< handle for mpi calculation
        real(rp)                               , intent(out) :: Int         !< integral value

        ! local declaration
        real(rp), dimension(2) :: step
        real(rp)               :: dv, my_int
        integer                :: i,j, err = 0

        my_int = 0.0_rp

        !$omp parallel do default(private) &
        !$omp shared(sx,sy,ex,ey,field,dx,dy) &
        !$omp reduction(+:my_int)
        do j = sy,ey
           do i = sx,ex

              step(1) = dx(i)
              step(2) = dy(j)

              dv = step(1) * step(2)

              my_int = my_int + field(i,j) * dv

           enddo
        enddo
        !$omp end parallel do


        if(present(mpi)) then
           call MPI_allreduce(my_int, int, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, err)
        else
           int = my_int
        endif


        return
end subroutine int_2Dfield



subroutine int_3Dfield(field,sx,ex,dx,sy,ey,dy,sz,ez,dz,int,mpi) 
! -----------------------------------------------------------
!       This function compute the integral of a 3D field
! -----------------------------------------------------------
        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(in)  :: field       !< field we want to integrate
        real(rp), dimension(:)    , allocatable, intent(in)  :: dx, dy, dz  !< space coordinates
        integer                                , intent(in)  :: sx,ex       !< integer init and end positions
        integer                                , intent(in)  :: sy,ey       !< integer init and end positions
        integer                                , intent(in)  :: sz,ez       !< integer init and end positions
        logical , optional                     , intent(in)  :: mpi         !< handle for mpi calculation
        real(rp)                               , intent(out) :: Int         !< integral value

        ! local declaration
        real(rp), dimension(3) :: step
        real(rp)               :: dv, my_int
        integer                :: i,j,k, err = 0

        my_int = 0.0_rp

        !$omp parallel do default(private) &
        !$omp shared(sx,sy,sz,ex,ey,ez,field,dx,dy,dz) &
        !$omp reduction(+:my_int)
        do k = sz,ez
           do j = sy,ey
              do i = sx,ex

                 step(1) = dx(i)
                 step(2) = dy(j)
                 step(3) = dz(k)

                 dv = step(1) * step(2) * step(3)

                 my_int = my_int + field(i,j,k) * dv

              enddo
           enddo
        enddo
        !$omp end parallel do


        if(present(mpi)) then
           call MPI_allreduce(my_int, int, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, err)
        else
           int = my_int
        endif


        return
end subroutine int_3Dfield


subroutine int_4Dfield(field,id,sx,ex,dx,sy,ey,dy,sz,ez,dz,int,mpi) 
! -----------------------------------------------------------
!       This function compute the integral of a 4D field
! -----------------------------------------------------------
        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(in)  :: field       !< field we want to integrate
        integer                                  , intent(in)  :: id          !< field index
        integer                                  , intent(in)  :: sx,ex       !< integer init and end positions
        integer                                  , intent(in)  :: sy,ey       !< integer init and end positions
        integer                                  , intent(in)  :: sz,ez       !< integer init and end positions
        real(rp), dimension(:)      , allocatable, intent(in)  :: dx, dy, dz  !< space coordinates
        logical , optional                       , intent(in)  :: mpi         !< handle for mpi calculation
        real(rp)                                 , intent(out) :: Int         !< integral value

        ! local declaration
        real(rp), dimension(3) :: step
        real(rp)               :: dv, my_int
        integer                :: i,j,k, err = 0

        my_int = 0.0_rp

        !$omp parallel do default(private) &
        !$omp shared(sx,sy,sz,ex,ey,ez,field,id,dx,dy,dz) &
        !$omp reduction(+:my_int)
        do k = sz,ez
           do j = sy,ey
              do i = sx,ex

                 step(1) = dx(i)
                 step(2) = dy(j)
                 step(3) = dz(k)

                 dv = step(1) * step(2) * step(3)

                 my_int = my_int + field(i,j,k,id) * dv

              enddo
           enddo
        enddo
        !$omp end parallel do


        if(present(mpi)) then
           call MPI_allreduce(my_int, int, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, err)
        else
           int = my_int
        endif


        return
end subroutine int_4Dfield






subroutine mean_3Dfield(field,sx,ex,dx,sy,ey,dy,sz,ez,dz,mean,mpi) 
! -----------------------------------------------------------
!       This function compute the integral of a 3D field
! -----------------------------------------------------------
        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(in)  :: field       !< field we want to integrate
        real(rp), dimension(:)    , allocatable, intent(in)  :: dx, dy, dz  !< space coordinates
        integer                                , intent(in)  :: sx,ex       !< integer init and end positions
        integer                                , intent(in)  :: sy,ey       !< integer init and end positions
        integer                                , intent(in)  :: sz,ez       !< integer init and end positions
        logical , optional                     , intent(in)  :: mpi         !< handle for mpi calculation
        real(rp)                               , intent(out) :: mean        !< integral mean 

        ! local declaration
        real(rp), dimension(3) :: step
        real(rp)               :: dv, my_int, my_vol, int, vol
        integer                :: i,j,k, err = 0

        my_int = 0.0_rp
        my_vol = 0.0_rp

        !$omp parallel do default(private) &
        !$omp shared(sx,sy,sz,ex,ey,ez,field,dx,dy,dz) &
        !$omp reduction(+:my_int) reduction(+:my_vol)
        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex

                 step(1) = dx(i)
                 step(2) = dy(j)
                 step(3) = dz(k)

                 dv = step(1) * step(2) * step(3)

                 my_int = my_int + field(i,j,k) * dv
                 my_vol = my_vol + dv

              enddo
           enddo
        enddo
        !$omp end parallel do


        if(present(mpi)) then
           call MPI_allreduce(my_int, int, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, err)
           call MPI_allreduce(my_vol, vol, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, err)

        else
           int = my_int
           vol = my_vol
        endif

        mean = int/vol

        return
end subroutine mean_3Dfield




subroutine mean_4Dfield(field,sx,ex,dx,sy,ey,dy,sz,ez,dz,sv,ev,mean,mpi) 
! -----------------------------------------------------------
!       This function compute the integral of a 4D field
! -----------------------------------------------------------
        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(in)  :: field       !< field we want to integrate
        real(rp), dimension(:)      , allocatable, intent(in)  :: dx, dy, dz  !< space coordinates
        integer                                  , intent(in)  :: sx,ex       !< integer init and end positions of x dimension
        integer                                  , intent(in)  :: sy,ey       !< integer init and end positions of y dimension
        integer                                  , intent(in)  :: sz,ez       !< integer init and end positions of z dimension
        integer                                  , intent(in)  :: sv,ev       !> integer init and end positions of codimension
        logical , optional                       , intent(in)  :: mpi         !< handle for mpi calculation
        real(rp), dimension(sv:ev)               , intent(out) :: mean        !< integral mean 

        ! local declaration
        real(rp), dimension(3)     :: step
        real(rp), dimension(sv:ev) :: lcl_int, gbl_int
        integer , dimension(2)     :: err = 0
        real(rp)                   :: lcl_vol, gbl_vol, dv
        integer                    :: id, i,j,k


        do id = sv,ev
           lcl_int(id) = 0.0_rp
           lcl_vol     = 0.0_rp
           !$omp parallel do default(private) &
           !$omp shared(sx,sy,sz,ex,ey,ez,field,dx,dy,dz,id) &
           !$omp reduction(+:lcl_int) reduction(+:lcl_vol)
           do k       = sz,ez
              do j    = sy,ey
                 do i = sx,ex

                    step(1) = dx(i)
                    step(2) = dy(j)
                    step(3) = dz(k)

                    dv = step(1) * step(2) * step(3)

                    lcl_int(id) = lcl_int(id) + field(i,j,k,id) * dv
                    lcl_vol     = lcl_vol + dv

                 enddo
              enddo
           enddo
           !$omp end parallel do
        enddo


        if(present(mpi)) then
           call MPI_allreduce(lcl_int, gbl_int, size(lcl_int), MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, err(1))
           call MPI_allreduce(lcl_vol, gbl_vol, 1            , MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, err(2))

           if(sum(err) .ne. 0) stop ' ERROR in mean_4Dfield!'
        else
           gbl_int = lcl_int
           gbl_vol = lcl_vol
        endif

        mean = gbl_int/gbl_vol

        return
end subroutine mean_4Dfield







subroutine rms_3Dfield(field,sx,ex,dx,sy,ey,dy,sz,ez,dz,rms,mpi) 
! -----------------------------------------------------------
!       This function compute the rms of a 3D field
! -----------------------------------------------------------
        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(in)  :: field       !< field we want to integrate
        real(rp), dimension(:)    , allocatable, intent(in)  :: dx, dy, dz  !< space coordinates
        integer                                , intent(in)  :: sx,ex       !< integer init and end positions
        integer                                , intent(in)  :: sy,ey       !< integer init and end positions
        integer                                , intent(in)  :: sz,ez       !< integer init and end positions
        logical , optional                     , intent(in)  :: mpi         !< handle for mpi calculation
        real(rp)                               , intent(out) :: rms         !< integral rms 

        ! local declaration
        real(rp), dimension(3) :: step
        real(rp)               :: lcl_rms, lcl_vol
        real(rp)               :: gbl_rms, gbl_vol
        real(rp)               :: dv, mean
        integer                :: i,j,k, err = 0
        
        ! 
        ! === compute the integral mean of the field
        !
        call mean_3Dfield(field,sx,ex,dx,sy,ey,dy,sz,ez,dz,mean,mpi)
        !
        ! === compute the field root mean square
        !
        lcl_rms = 0.0_rp
        lcl_vol = 0.0_rp

        !$omp parallel do default(private) &
        !$omp shared(sx,sy,sz,ex,ey,ez,field,dx,dy,dz,mean) &
        !$omp reduction(+:lcl_rms) reduction(+:lcl_vol)
        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex

                 step(1) = dx(i)
                 step(2) = dy(j)
                 step(3) = dz(k)

                 dv = step(1) * step(2) * step(3)

                 lcl_rms = lcl_rms + (field(i,j,k) - mean)**2 * dv
                 lcl_vol = lcl_vol + dv

              enddo
           enddo
        enddo
        !$omp end parallel do
        !
        ! === MPI staff
        !
        if(present(mpi)) then
           call MPI_allreduce(lcl_rms, gbl_rms, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, err)
           call MPI_allreduce(lcl_vol, gbl_vol, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, err)

        else
           gbl_rms = lcl_rms
           gbl_vol = lcl_vol

        endif

        rms = sqrt(gbl_rms/gbl_vol)

        return
end subroutine rms_3Dfield








end module integration_module
