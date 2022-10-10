module statistics_module
implicit none

private
public mean0D, mean1D, mean2D, rey_averaged, gbl_mean!, ySymmetry, yAntMetry

interface gbl_mean
  module procedure gbl_mean3Dto1D, gbl_mean3Dto2D, gbl_mean4Dto1D, gbl_mean4Dto2D
end interface gbl_mean

interface mean2D
  module procedure MeanField4Dto2D, MeanField3Dto2D
end interface mean2D

interface mean1D
  module procedure MeanField4Dto1D, MeanField3Dto1D
end interface mean1D

interface mean0D
  module procedure MeanField4Dto0D, MeanField3Dto0D
end interface mean0D

interface rey_averaged
  module procedure ReynoldsAverage0D, ReynoldsAverage1D, ReynoldsAverage2D, ReynoldsAverage3D, ReynoldsAverage4D
end interface rey_averaged

!interface ySymmetry
!  module procedure ySymmetry1D, ySymmetry2D
!end interface ySymmetry
!
!interface yAntMetry
!  module procedure yAntMetry1D, yAntMetry2D
!end interface yAntMetry


contains
subroutine gbl_mean3Dto1D(v,dir,n1,n2,time,trestart,dt,sx,ex,sy,ey,sz,ez,comm,vmean,mpif)      

        use parameters_module, only: rp

        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(in)    :: v
        real(rp), allocatable, dimension(:)    , intent(inout) :: vmean
        real(rp)                               , intent(in)    :: time, tRestart,dt
        integer                                , intent(in)    :: n1, n2, sx, ex, sy, ey, sz, ez
        integer                                , intent(in)    :: comm
        character(*)                           , intent(in)    :: dir
        logical                                , intent(in)    :: mpif

        !local declaration
        real(rp), allocatable, dimension(:) :: tmp
        integer , dimension(2)              :: b
        integer                             :: err = 0
       
        b = 0
        if(dir == 'xy') b = (/sz,ez/)
        if(dir == 'xz') b = (/sy,ey/)
        if(dir == 'yz') b = (/sx,ex/)

        allocate(tmp(b(1):b(2)), stat = err)
        if(err .ne. 0) stop ' Allocatation error in gbl_mean'

        call MeanField3Dto1D(v,dir,n1,n2,sx,ex,sy,ey,sz,ez,comm,tmp,mpif)

        call rey_averaged(tmp,time,tRestart,dt,vmean)

        deallocate(tmp)

        return
end subroutine gbl_mean3Dto1D


subroutine gbl_mean3Dto2D(v,dir,nSlice,time,trestart,dt,sx,ex,sy,ey,sz,ez,comm,vmean,mpif)      

        use parameters_module, only: rp

        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(in)    :: v
        real(rp), allocatable, dimension(:,:)  , intent(inout) :: vmean
        real(rp)                               , intent(in)    :: time, tRestart,dt
        integer                                , intent(in)    :: nSlice, sx, ex, sy, ey, sz, ez
        integer                                , intent(in)    :: comm
        character(*)                           , intent(in)    :: dir
        logical                                , intent(in)    :: mpif

        !local declaration
        real(rp), allocatable, dimension(:,:) :: tmp
        integer , dimension(2)                :: lb, ub
        integer                               :: err = 0
       
        lb = 0; ub = 0
        if    (dir == 'x') then
          lb = (/sy,sz/); ub = (/ey,ez/)
        elseif(dir == 'y') then
          lb = (/sx,sz/); ub = (/ex,ez/)
        elseif(dir == 'z') then 
          lb = (/sx,sy/); ub = (/ex,ey/)
        endif
        
        allocate(tmp(lb(1):ub(1), lb(2):ub(2)), stat=err)
        if(err .ne. 0) stop ' Allocatation error in gbl_mean'

        call MeanField3Dto2D(v,dir,nSlice,sx,ex,sy,ey,sz,ez,comm,tmp,mpif)

        call rey_averaged(tmp,time,tRestart,dt,vmean)

        deallocate(tmp)

        return
end subroutine gbl_mean3Dto2D


subroutine gbl_mean4Dto1D(v,dir,n1,n2,time,trestart,dt,sx,ex,sy,ey,sz,ez,comm,vmean,mpif)      

        use parameters_module, only: rp
        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(in)    :: v
        real(rp), allocatable, dimension(:,:)    , intent(inout) :: vmean
        real(rp)                                 , intent(in)    :: time, tRestart,dt
        integer                                  , intent(in)    :: n1, n2, sx, ex, sy, ey, sz, ez
        integer                                  , intent(in)    :: comm
        character(*)                             , intent(in)    :: dir
        logical                                  , intent(in)    :: mpif

        !local declaration
        real(rp), allocatable, dimension(:,:) :: tmp
        integer , dimension(2)                :: b
        integer                               :: err = 0
       
        b = 0
        if(dir == 'xy') b = (/sz,ez/)
        if(dir == 'xz') b = (/sy,ey/)
        if(dir == 'yz') b = (/sx,ex/)

        allocate(tmp(b(1):b(2), lbound(v,4):ubound(v,4)), stat = err)
        if(err .ne. 0) stop ' Allocatation error in gbl_mean'
        
        call MeanField4Dto1D(v,dir,n1,n2,sx,ex,sy,ey,sz,ez,comm,tmp,mpif)

        call rey_averaged(tmp,time,tRestart,dt,vmean)

        deallocate(tmp)

        return
end subroutine gbl_mean4Dto1D



subroutine gbl_mean4Dto2D(v,dir,nSlice,time,trestart,dt,sx,ex,sy,ey,sz,ez,comm,vmean,mpif)      

        use parameters_module, only: rp 
        
        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(in)    :: v
        real(rp), allocatable, dimension(:,:,:)  , intent(inout) :: vmean
        real(rp)                                 , intent(in)    :: time, tRestart,dt
        integer                                  , intent(in)    :: nSlice, sx, ex, sy, ey, sz, ez
        integer                                  , intent(in)    :: comm
        character(*)                             , intent(in)    :: dir
        logical                                  , intent(in)    :: mpif

        !local declaration
        real(rp), allocatable, dimension(:,:,:) :: tmp
        integer , dimension(2)                  :: lb, ub
        integer                                 :: err = 0
       
        lb = 0; ub = 0
        if    (dir == 'x') then
          lb = (/sy,sz/); ub = (/ey,ez/)
        elseif(dir == 'y') then
          lb = (/sx,sz/); ub = (/ex,ez/)
        elseif(dir == 'z') then 
          lb = (/sx,sy/); ub = (/ex,ey/)
        endif
        
        allocate(tmp(lb(1):ub(1), lb(2):ub(2),lbound(v,4):ubound(v,4)), stat=err)
        if(err .ne. 0) stop ' Allocatation error in gbl_mean'

        call MeanField4Dto2D(v,dir,nSlice,sx,ex,sy,ey,sz,ez,comm,tmp,mpif)

        call rey_averaged(tmp,time,tRestart,dt,vmean)

        deallocate(tmp)

        return
end subroutine gbl_mean4Dto2D



subroutine MeanField4Dto2D(v,dir,nSlice,sx,ex,sy,ey,sz,ez,comm,vmean,mpif)

        use parameters_module, only: rp, MPI_RP, msm
        use mpi
        
        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(in)    :: v
        real(rp), allocatable, dimension(:,:,:)  , intent(inout) :: vmean
        integer                                  , intent(in)    :: nSlice
        integer                                  , intent(in)    :: sx, ex, sy, ey, sz, ez
        integer                                  , intent(in)    :: comm
        character(*)                             , intent(in)    :: dir
        logical                                  , intent(in)    :: mpif

        ! local declaration

        real(rp) :: lcl_sum, irSlice
        integer  :: i,j,k,l, lb,ub,err = 0

        lb = lbound(v,4)
        ub = ubound(v,4)
        
        selectcase(dir)

          case('x') ! x mean field
                
            do l       = lb,ub
               do k    = sz,ez
                  do j = sy,ey
                     lcl_sum = sum(v(sx:ex,j,k,l))
                     if(mpif) then
                       call MPI_ALLREDUCE(lcl_sum, vmean(j,k,l), 1, MPI_RP, msm, comm,err)
                     else
                       vmean(j,k,l) = lcl_sum
                     endif
                  enddo
               enddo
            enddo
          
          case('y') ! y mean field

            do l       = lb,ub
               do k    = sz,ez
                  do i = sx,ex
                     lcl_sum = sum(v(i,sy:ey,k,l))
                     if(mpif) then
                       call MPI_ALLREDUCE(lcl_sum, vmean(i,k,l), 1, MPI_RP, msm, comm,err)
                     else
                       vmean(i,k,l) = lcl_sum
                     endif
                  enddo
               enddo
            enddo
          
          case('z') ! z mean field

            do l       = lb,ub
               do j    = sy,ey
                  do i = sx,ex
                     lcl_sum = sum(v(i,j,sz:ez,l))
                     if(mpif) then
                       call MPI_ALLREDUCE(lcl_sum, vmean(i,j,l), 1, MPI_RP, msm, comm,err)
                     else
                       vmean(i,j,l) = lcl_sum
                     endif
                  enddo
               enddo
            enddo
          
        endselect

        irSlice = 1.0_rp/real(nSlice,rp)
        vmean = vmean*irSlice

        return
end subroutine MeanField4Dto2D



subroutine MeanField4Dto1D(v,dir,nSlice1,nSlice2,sx,ex,sy,ey,sz,ez,comm,vmean,mpif)

        use parameters_module, only: rp, MPI_RP, msm
        use mpi

        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(in)    :: v
        real(rp), allocatable, dimension(:,:)    , intent(inout) :: vmean
        integer                                  , intent(in)    :: nSlice1, nSlice2
        integer                                  , intent(in)    :: sx, ex, sy, ey, sz, ez
        integer                                  , intent(in)    :: comm
        character(*)                             , intent(in)    :: dir
        logical                                  , intent(in)    :: mpif

        ! local declarations
        real(rp)           :: lcl_sum, irSlice
        integer            :: i,j,k,l,lb,ub,err=0

        lb = lbound(v,4)
        ub = ubound(v,4)
        
        selectcase(dir)

          case('xy') ! xy mean field

            do l    = lb,ub
               do k = sz,ez
                  lcl_sum = sum(v(sx:ex,sy:ey,k,l))
                  if(mpif) then
                    call MPI_ALLREDUCE(lcl_sum, vmean(k,l), 1, MPI_RP, msm, comm,err)
                  else
                    vmean(k,l) = lcl_sum
                  endif
               enddo
            enddo

          case('xz') ! xz mean field

            do l    = lb,ub
               do j = sy,ey
                  lcl_sum = sum(v(sx:ex,j,sz:ez,l))
                  if(mpif) then
                    call MPI_ALLREDUCE(lcl_sum, vmean(j,l), 1, MPI_RP, msm, comm,err)
                  else
                    vmean(j,l) = lcl_sum
                  endif
               enddo
            enddo

          case('yz') ! yz mean field

            do l    = lb,ub
               do i = sx,ex
                  lcl_sum = sum(v(i,sy:ey,sz:ez,l))
                  if(mpif) then
                    call MPI_ALLREDUCE(lcl_sum, vmean(i,l), 1, MPI_RP, msm, comm,err)
                  else
                    vmean(i,l) = lcl_sum
                  endif
               enddo
            enddo

        end select

        irSlice = 1.0_rp/(real(nSlice1,rp)*real(nSlice2,rp))
        vmean = vmean*irSlice

        return
end subroutine MeanField4Dto1D




subroutine MeanField4Dto0D(v,nSlice1,nSlice2,nSlice3,sx,ex,sy,ey,sz,ez,comm,vmean,mpif)

        use parameters_module, only: rp, MPI_RP, msm
        use mpi

        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(in)    :: v
        real(rp), allocatable, dimension(:)      , intent(inout) :: vmean
        integer                                  , intent(in)    :: nSlice1, nSlice2, nSlice3
        integer                                  , intent(in)    :: sx, ex, sy, ey, sz, ez
        integer                                  , intent(in)    :: comm
        logical                                  , intent(in)    :: mpif

        ! local declarations
        real(rp)           :: lcl_sum, irSlice
        integer            :: l,lb,ub,err=0

        lb = lbound(v,4)
        ub = ubound(v,4)
        
        do l    = lb,ub
           lcl_sum = sum(v(sx:ex,sy:ey,sz:ez,l))
           if(mpif) then
             call MPI_ALLREDUCE(lcl_sum, vmean(l), 1, MPI_RP, msm, comm,err)
           else
             vmean(l) = lcl_sum
           endif
        enddo

        irSlice = 1.0_rp/(real(nSlice1,rp)*real(nSlice2,rp)*real(nSlice3,rp))
        vmean = vmean*irSlice

        return
end subroutine MeanField4Dto0D
! --------------------------------------------------------------------------------
! --------------------------------------------------------------------------------
! --------------------------------------------------------------------------------


subroutine MeanField3Dto2D(v,dir,nSlice,sx,ex,sy,ey,sz,ez,comm,vmean,mpif)

        use parameters_module, only: rp, MPI_RP, msm
        use mpi
        
        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(in)    :: v
        real(rp), allocatable, dimension(:,:)  , intent(inout) :: vmean
        integer                                , intent(in)    :: nSlice
        integer                                , intent(in)    :: sx, ex, sy, ey, sz, ez
        integer                                , intent(in)    :: comm
        character(*)                           , intent(in)    :: dir
        logical                                , intent(in)    :: mpif

        ! local declaration

        real(rp) :: lcl_sum, irSlice
        integer  :: i,j,k, err = 0

        selectcase(dir)

          case('x') ! x mean field
                
            do k    = sz,ez
               do j = sy,ey
                  lcl_sum = sum(v(sx:ex,j,k))
                  if(mpif) then
                    call MPI_ALLREDUCE(lcl_sum, vmean(j,k), 1, MPI_RP, msm, comm,err)
                  else
                    vmean(j,k) = lcl_sum
                  endif
               enddo
            enddo
          
          case('y') ! y mean field

            do k    = sz,ez
               do i = sx,ex
                  lcl_sum = sum(v(i,sy:ey,k))
                  if(mpif) then
                    call MPI_ALLREDUCE(lcl_sum, vmean(i,k), 1, MPI_RP, msm, comm,err)
                  else
                    vmean(i,k) = lcl_sum
                  endif
               enddo
            enddo
          
          case('z') ! z mean field

            do j    = sy,ey
               do i = sx,ex
                  lcl_sum = sum(v(i,j,sz:ez))
                  if(mpif) then
                    call MPI_ALLREDUCE(lcl_sum, vmean(i,j), 1, MPI_RP, msm, comm,err)
                  else
                    vmean(i,j) = lcl_sum
                  endif
               enddo
            enddo
          
        endselect

        irSlice = 1.0_rp/real(nSlice,rp)
        vmean = vmean*irSlice

        return
end subroutine MeanField3Dto2D



subroutine MeanField3Dto1D(v,dir,nSlice1,nSlice2,sx,ex,sy,ey,sz,ez,comm,vmean,mpif)

        use parameters_module, only: rp, MPI_RP, msm
        use mpi

        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(in)    :: v
        real(rp), allocatable, dimension(:)    , intent(inout) :: vmean
        integer                                , intent(in)    :: nSlice1, nSlice2
        integer                                , intent(in)    :: sx, ex, sy, ey, sz, ez
        integer                                , intent(in)    :: comm
        character(*)                           , intent(in)    :: dir
        logical                                , intent(in)    :: mpif

        ! local declarations
        real(rp)           :: lcl_sum, irSlice
        integer            :: i,j,k,err=0

        selectcase(dir)

          case('xy') ! xy mean field

            do k = sz,ez
               lcl_sum = sum(v(sx:ex,sy:ey,k))
               if(mpif) then
                 call MPI_ALLREDUCE(lcl_sum, vmean(k), 1, MPI_RP, msm, comm,err)
               else
                 vmean(k) = lcl_sum
               endif
            enddo

          case('xz') ! xz mean field

            do j = sy,ey
               lcl_sum = sum(v(sx:ex,j,sz:ez))
               if(mpif) then
                 call MPI_ALLREDUCE(lcl_sum, vmean(j), 1, MPI_RP, msm, comm,err)
               else
                 vmean(j) = lcl_sum
               endif
            enddo

          case('yz') ! yz mean field

            do i = sx,ex
               lcl_sum = sum(v(i,sy:ey,sz:ez))
               if(mpif) then
                 call MPI_ALLREDUCE(lcl_sum, vmean(i), 1, MPI_RP, msm, comm,err)
               else
                 vmean(i) = lcl_sum
               endif
            enddo

        end select

        irSlice = 1.0_rp/(real(nSlice1,rp)*real(nSlice2,rp))
        vmean = vmean*irSlice

        return
end subroutine MeanField3Dto1D


subroutine MeanField3Dto0D(v,nSlice1,nSlice2,nSlice3,sx,ex,sy,ey,sz,ez,comm,vmean,mpif)

        use parameters_module, only: rp, MPI_RP, msm
        use mpi

        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(in)    :: v
        real(rp), allocatable                  , intent(inout) :: vmean
        integer                                , intent(in)    :: nSlice1, nSlice2, nSlice3
        integer                                , intent(in)    :: sx, ex, sy, ey, sz, ez
        integer                                , intent(in)    :: comm
        logical                                , intent(in)    :: mpif

        ! local declarations
        real(rp)           :: lcl_sum, irSlice
        integer            :: err=0

        lcl_sum = sum(v(sx:ex,sy:ey,sz:ez))
        if(mpif) then
          call MPI_ALLREDUCE(lcl_sum, vmean, 1, MPI_RP, msm, comm,err)
        else
          vmean = lcl_sum
        endif

        irSlice = 1.0_rp/(real(nSlice1,rp)*real(nSlice2,rp)*real(nSlice3,rp))
        vmean = vmean*irSlice

        return
end subroutine MeanField3Dto0D



subroutine ReynoldsAverage0D(v,time,time_restart,dt,ReyAv)

        use parameters_module, only: rp

        implicit none
        real(rp), intent(in)    :: v
        real(rp), intent(inout) :: ReyAv
        real(rp), intent(in)    :: time
        real(rp), intent(in)    :: time_restart
        real(rp), intent(in)    :: dt

        if(abs(time - time_restart) < 1.0E-14_rp) then
          ReyAv = v
        else
          ReyAv = ((time-time_restart - dt) * ReyAv + v * dt) / (time-time_restart)
        endif

        return
end subroutine ReynoldsAverage0D




subroutine ReynoldsAverage1D(v,time,time_restart,dt,ReyAv)

        use parameters_module, only: rp

        implicit none
        real(rp), allocatable, dimension(:), intent(in)    :: v
        real(rp), allocatable, dimension(:), intent(inout) :: ReyAv
        real(rp)                           , intent(in)    :: time
        real(rp)                           , intent(in)    :: time_restart
        real(rp)                           , intent(in)    :: dt

        integer :: s1,e1
        s1 = lbound(v,1)
        e1 = ubound(v,1)
        
        if(abs(time - time_restart) < 1.0E-14_rp) then
          ReyAv(s1:e1) = v(s1:e1)
        else
          ReyAv(s1:e1) = ((time-time_restart - dt) * ReyAv(s1:e1) + v(s1:e1) * dt) / (time-time_restart)
        endif

        return
end subroutine ReynoldsAverage1D

subroutine ReynoldsAverage2D(v,time,time_restart,dt,ReyAv)

        use parameters_module, only: rp

        implicit none
        real(rp), allocatable, dimension(:,:), intent(in)    :: v
        real(rp), allocatable, dimension(:,:), intent(inout) :: ReyAv
        real(rp)                             , intent(in)    :: time
        real(rp)                             , intent(in)    :: time_restart
        real(rp)                             , intent(in)    :: dt

        integer :: s1,e1,s2,e2
        s1 = lbound(v,1); s2 = lbound(v,2)
        e1 = ubound(v,1); e2 = ubound(v,2)
        
        if(abs(time - time_restart) < 1.0E-14_rp) then
          ReyAv(s1:e1,s2:e2) = v(s1:e1,s2:e2)
        else
          ReyAv(s1:e1,s2:e2) = ((time - time_restart - dt) &
                  * ReyAv(s1:e1,s2:e2) + v(s1:e1,s2:e2) * dt) / (time - time_restart)
        endif

        return
end subroutine ReynoldsAverage2D

subroutine ReynoldsAverage3D(v,time,time_restart,dt,ReyAv)

        use parameters_module, only: rp

        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(in)    :: v
        real(rp), allocatable, dimension(:,:,:), intent(inout) :: ReyAv
        real(rp)                               , intent(in)    :: time
        real(rp)                               , intent(in)    :: time_restart
        real(rp)                               , intent(in)    :: dt

        integer :: s1,e1,s2,e2,s3,e3
        s1 = lbound(v,1); s2 = lbound(v,2); s3 = lbound(v,3)
        e1 = ubound(v,1); e2 = ubound(v,2); e3 = ubound(v,3)
        
        if(abs(time - time_restart) < 1.0E-14_rp) then
          ReyAv(s1:e1,s2:e2,s3:e3) = v(s1:e1,s2:e2,s3:e3)
        else
          ReyAv(s1:e1,s2:e2,s3:e3) = ((time - time_restart - dt) &
                  * ReyAv(s1:e1,s2:e2,s3:e3) + v(s1:e1,s2:e2,s3:e3) * dt) / (time - time_restart)
        endif

        return
end subroutine ReynoldsAverage3D

subroutine ReynoldsAverage4D(v,time,time_restart,dt,ReyAv)

        use parameters_module, only: rp

        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(in)    :: v
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: ReyAv
        real(rp)                                 , intent(in)    :: time
        real(rp)                                 , intent(in)    :: time_restart
        real(rp)                                 , intent(in)    :: dt

        integer :: s1,e1,s2,e2,s3,e3,s4,e4
        s1 = lbound(v,1); s2 = lbound(v,2); s3 = lbound(v,3); s4 = lbound(v,4)
        e1 = ubound(v,1); e2 = ubound(v,2); e3 = ubound(v,3); e4 = ubound(v,4)
        
        if(abs(time - time_restart) < 1.0E-14_rp) then
          ReyAv(s1:e1,s2:e2,s3:e3,s4:e4) = v(s1:e1,s2:e2,s3:e3,s4:e4)
        else
          ReyAv(s1:e1,s2:e2,s3:e3,s4:e4) = ((time - time_restart - dt) &
                  * ReyAv(s1:e1,s2:e2,s3:e3,s4:e4) + v(s1:e1,s2:e2,s3:e3,s4:e4) * dt) / (time - time_restart)
        endif

        return
end subroutine ReynoldsAverage4D































end module statistics_module
