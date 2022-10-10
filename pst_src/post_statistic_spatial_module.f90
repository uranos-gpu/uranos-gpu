module post_statistic_spatial_module
implicit none

private
public global_spatial_average, global_spatial_standard_deviation

contains
subroutine global_spatial_average(v,periodic_dir,symmtric_dir,v_spatial_average)

        use parameters_module, only: rp
        use mpi_module       , only: lbx,ubx, lby,uby, lbz,ubz

        implicit none
        character(1), dimension(3)             , intent(in)    :: periodic_dir
        character(1), dimension(3)             , intent(in)    :: symmtric_dir
        real(rp), dimension(:,:,:), allocatable, intent(in)    :: v
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: v_spatial_average

        ! local declarations
        real(rp), dimension(:,:,:), allocatable :: tmp1, tmp2
        integer                                 :: err = 0

        ! Exploing periodic directions in order to compute the average
        ! ---------------------------
        
        ! only x is periodic
        if    (periodic_dir(1)//periodic_dir(2)//periodic_dir(3) == 'TFF') then

          call s_averaged_x(v,v_spatial_average)
        
        ! only y is periodic
        elseif(periodic_dir(1)//periodic_dir(2)//periodic_dir(3) == 'FTF') then

          call s_averaged_y(v,v_spatial_average)

        ! only z is periodic
        elseif(periodic_dir(1)//periodic_dir(2)//periodic_dir(3) == 'FFT') then

          call s_averaged_z(v,v_spatial_average)
        
        ! x and y are periodic
        elseif(periodic_dir(1)//periodic_dir(2)//periodic_dir(3) == 'TTF') then

          allocate(tmp1(lbx:ubx, lby:uby, lbz:ubz), stat = err)
          if(err .ne. 0) stop ' Allocation error in global averaged.'; tmp1 = 0.0_rp

          call s_averaged_x(v   ,tmp1)
          call s_averaged_y(tmp1,v_spatial_average)

          deallocate(tmp1)

        ! x and z are periodic
        elseif(periodic_dir(1)//periodic_dir(2)//periodic_dir(3) == 'TFT') then

          allocate(tmp1(lbx:ubx, lby:uby, lbz:ubz), stat = err)
          if(err .ne. 0) stop ' Allocation error in global averaged.'; tmp1 = 0.0_rp

          call s_averaged_x(v   ,tmp1)
          call s_averaged_z(tmp1,v_spatial_average)

          deallocate(tmp1)

        ! y and z are periodic
        elseif(periodic_dir(1)//periodic_dir(2)//periodic_dir(3) == 'FTT') then

          allocate(tmp1(lbx:ubx, lby:uby, lbz:ubz), stat = err)
          if(err .ne. 0) stop ' Allocation error in global averaged.'; tmp1 = 0.0_rp

          call s_averaged_y(v   ,tmp1)
          call s_averaged_z(tmp1,v_spatial_average)

          deallocate(tmp1)

        ! x, y and z are periodic
        elseif(periodic_dir(1)//periodic_dir(2)//periodic_dir(3) == 'TTT') then

          allocate(tmp1(lbx:ubx, lby:uby, lbz:ubz), &
                   tmp2(lbx:ubx, lby:uby, lbz:ubz), stat = err)
          if(err .ne. 0) stop ' Allocation error in global averaged.'

          tmp1 = 0.0_rp
          tmp2 = 0.0_rp

          call s_averaged_x(v   ,tmp1)
          call s_averaged_y(tmp1,tmp2)
          call s_averaged_y(tmp2,v_spatial_average)

          deallocate(tmp1, tmp2)

        ! DEFAULT: there are no periodic directions
        else
        
          !$omp workshare
          v_spatial_average = v
          !$omp end workshare

        endif
        
        ! exploiting symmetrical directions
        ! --------------------------------------------------------------------
        if    (symmtric_dir(1)//symmtric_dir(2)//symmtric_dir(3) == 'TFF') then
          !call symmetry_x(v_spatial_average)

        elseif(symmtric_dir(1)//symmtric_dir(2)//symmtric_dir(3) == 'FTF') then
!          call symmetry_y(v_spatial_average)

        elseif(symmtric_dir(1)//symmtric_dir(2)//symmtric_dir(3) == 'FFT') then
          !call symmetry_z(v_spatial_average)

        endif

        return
end subroutine global_spatial_average


subroutine s_averaged_x(v,s_mean_v)
! -------------------------------------------------------------------------
!       This subroutine computes the spatial mean of a field v along x
! -------------------------------------------------------------------------
        use parameters_module, only: rp
        use mpi_module       , only: lbx,ubx,lby,uby,lbz,ubz,sx,ex
        use mesh_module      , only: xstep

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(in)    :: v          !< field we want to average
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: s_mean_v   !< averaged field

        ! local declaration
        real(rp), dimension(:), allocatable :: tmp             !< temporay array
        real(rp)                            :: i_sum_dx        !< grid spaces
        integer                             :: j,k,il,err = 0
        
        allocate(tmp((uby-lby+1) * (ubz-lbz+1)), stat = err)
        if(err.ne.0) stop ' Allocation error in s_averaged_x'
        tmp = 0.0_rp

        ! storing the mean value in the tmp array
        il     = 1
        do k = lbz,ubz
           do j = lby,uby
                
              tmp(il) = sum(v(sx:ex,j,k) * xstep(sx:ex))

              il = il + 1

           enddo
        enddo
        
        ! riscaling the average
        i_sum_dx = 1.0_rp/sum(xstep(sx:ex))
        il       = 1
        do k = lbz,ubz
           do j = lby,uby

              s_mean_v(sx:ex,j,k) = i_sum_dx * tmp(il)

              ! extrapolate the mean in the ghost region
              s_mean_v(lbx:sx-1,j,k) = s_mean_v(sx,j,k)
              s_mean_v(ex+1:ubx,j,k) = s_mean_v(ex,j,k)
        
              ! update il
              il = il +1

           enddo
        enddo
        
        deallocate(tmp)

        return
end subroutine s_averaged_x


subroutine s_averaged_y(v,s_mean_v)
! -------------------------------------------------------------------------
!       This subroutine computes the spatial mean of a field v along y
! -------------------------------------------------------------------------
        use parameters_module, only: rp
        use mpi_module       , only: lbx,ubx,lby,uby,lbz,ubz,sy,ey
        use mesh_module      , only: ystep

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(in)    :: v          !< field we want to average
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: s_mean_v   !< averaged field

        ! local declaration
        real(rp), dimension(:), allocatable :: tmp             !< temporay array
        real(rp)                            :: i_sum_dy        !< grid spaces
        integer                             :: i,k,jl, err = 0
        
        allocate(tmp((ubx-lbx+1) * (ubz-lbz+1)), stat = err)
        if(err.ne.0) stop ' Allocation error in s_averaged_y'
        tmp = 0.0_rp
       
        ! storing the mean value in the tmp array
        jl     = 1
        do k = lbz,ubz
           do i = lbx,ubx
                
              tmp(jl) = sum(v(i,sy:ey,k)*ystep(sy:ey))

              jl = jl + 1

           enddo
        enddo
        
        ! riscaling the average
        i_sum_dy = 1.0_rp/sum(ystep(sy:ey))
        jl       = 1
        do k = lbz,ubz
           do i = lbx,ubx

              s_mean_v(i,sy:ey,k) = i_sum_dy * tmp(jl)

              ! extrapolate the mean in the ghost region
              s_mean_v(i,lby:sy-1,k) = s_mean_v(i,sy,k)
              s_mean_v(i,ey+1:uby,k) = s_mean_v(i,ey,k)

              ! update jl
              jl = jl +1

           enddo
        enddo
        deallocate(tmp)

        return
end subroutine s_averaged_y


subroutine s_averaged_z(v,s_mean_v)
! -------------------------------------------------------------------------
!       This subroutine computes the spatial mean of a field v along z
! -------------------------------------------------------------------------
        use parameters_module, only: rp
        use mpi_module       , only: lbx,ubx,lby,uby,lbz,ubz,sz,ez
        use mesh_module      , only: zstep

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(in)    :: v          !< field we want to average
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: s_mean_v   !< averaged field

        ! local declaration
        real(rp), dimension(:), allocatable :: tmp             !< temporay array
        real(rp)                            :: i_sum_dz        !< grid spaces
        integer                             :: i,j,kl, err = 0                  
        
        allocate(tmp((ubx-lbx+1) * (uby-lby+1)), stat = err)
        if(err.ne.0) stop ' Allocation error in s_averaged_z'
        tmp = 0.0_rp
       
        ! storing the mean value in the tmp array
        kl     = 1
        do j = lby,uby
           do i = lbx,ubx
                
              tmp(kl) = sum(v(i,j,sz:ez)*zstep(sz:ez))

              kl = kl + 1

           enddo
        enddo
        
        ! riscaling the average
        i_sum_dz = 1.0_rp/sum(zstep(sz:ez))
        kl       = 1
        do j = lby,uby
           do i = lbx,ubx

              s_mean_v(i,j,sz:ez) = i_sum_dz * tmp(kl)

              ! extrapolate the mean in the ghost region
              s_mean_v(i,j,lbz:sz-1) = s_mean_v(i,j,sz)
              s_mean_v(i,j,ez+1:ubz) = s_mean_v(i,j,ez)

              ! update kl
              kl = kl +1

           enddo
        enddo
        deallocate(tmp)

        return
end subroutine s_averaged_z


subroutine global_spatial_standard_deviation(v,v_mean,periodic_dir,v_rms)

        use parameters_module, only: rp

        implicit none
        real(rp)    , dimension(:,:,:), allocatable, intent(in)    :: v
        real(rp)    , dimension(:,:,:), allocatable, intent(in)    :: v_mean
        real(rp)    , dimension(:,:,:), allocatable, intent(inout) :: v_rms
        character(1), dimension(3)                 , intent(in)    :: periodic_dir


        ! only x is periodic
        if    (periodic_dir(1)//periodic_dir(2)//periodic_dir(3) == 'TFF') then
          call standard_deviation_x(v,v_mean,v_rms)

        ! only y is periodic
        elseif(periodic_dir(1)//periodic_dir(2)//periodic_dir(3) == 'FTF') then
          call standard_deviation_y(v,v_mean,v_rms)

        ! only z is periodic
        elseif(periodic_dir(1)//periodic_dir(2)//periodic_dir(3) == 'FFT') then
          call standard_deviation_z(v,v_mean,v_rms)

        ! xy are periodic
        elseif(periodic_dir(1)//periodic_dir(2)//periodic_dir(3) == 'TTF') then
          continue

        ! yz are periodic
        elseif(periodic_dir(1)//periodic_dir(2)//periodic_dir(3) == 'FTT') then
          continue

        ! xz are periodic
        elseif(periodic_dir(1)//periodic_dir(2)//periodic_dir(3) == 'TFT') then
          continue

        ! xyz are periodic
        elseif(periodic_dir(1)//periodic_dir(2)//periodic_dir(3) == 'TTT') then
          continue
        
        ! no periodic directions
        else
          continue

        endif




        return
end subroutine global_spatial_standard_deviation




subroutine standard_deviation_x(v,v_mean,v_rms)

        use parameters_module, only: rp
        use mpi_module       , only: lby, uby, lbz, ubz, sx, ex
        use mesh_module      , only: xstep

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(in)    :: v
        real(rp), dimension(:,:,:), allocatable, intent(in)    :: v_mean
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: v_rms

        ! local declaration
        real(rp), dimension(:), allocatable :: tmp          !< temporay array
        real(rp)                            :: i_sum        !< grid spaces
        integer                             :: j,k,il,err = 0
        
        allocate(tmp((uby-lby+1) * (ubz-lbz+1)), stat = err)
        if(err /= 0) stop ' Allocation error in standard_deviation_x'

        il = 1
        do k    = lbz,ubz
           do j = lby,uby
              tmp(il) = sum( (v(sx:ex,j,k) - v_mean(sx:ex,j,k))**2 * xstep(sx:ex))
              il = il + 1
           enddo
        enddo
        
        il = 1
        i_sum = 1.0_rp/sum(xstep(sx:ex))
        do k    = lbz,ubz
           do j = lby,uby
              v_rms(sx:ex,j,k) = sqrt(i_sum*tmp(il))
              il = il + 1
           enddo
        enddo

        deallocate(tmp)
        
        return
end subroutine standard_deviation_x



subroutine standard_deviation_y(v,v_mean,v_rms)
        
        use parameters_module, only: rp
        use mpi_module       , only: lbx, ubx, lbz, ubz, sy, ey
        use mesh_module      , only: ystep
                
        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(in)    :: v
        real(rp), dimension(:,:,:), allocatable, intent(in)    :: v_mean
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: v_rms

        ! local declarations
        real(rp), dimension(:), allocatable :: tmp
        real(rp)                            :: i_sum
        integer                             :: i,k, il, err = 0
        
        allocate(tmp((ubx-lbx+1) * (ubz-lbz+1)), stat = err)
        if(err.ne.0) stop ' Allocation error in standard_deviation_y'

        il = 1
        do k    = lbz,ubz
           do i = lbx,ubx
              tmp(il) = sum( (v(i,sy:ey,k) - v_mean(i,sy:ey,k))**2 * ystep(sy:ey))
              il = il + 1
           enddo
        enddo

        il = 1
        i_sum = 1.0_rp/sum(ystep(sy:ey))
        do k    = lbz,ubz
           do i = lbx,ubx
              v_rms(i,sy:ey,k) = sqrt(tmp(il)*i_sum)
              il = il + 1
           enddo
        enddo
        
        deallocate(tmp)

        return
end subroutine standard_deviation_y







subroutine standard_deviation_z(v,v_mean,v_rms)
        
        use parameters_module, only: rp
        use mpi_module       , only: lbx, ubx, lby, uby, sz, ez
        use mesh_module      , only: zstep
                
        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(in)    :: v
        real(rp), dimension(:,:,:), allocatable, intent(in)    :: v_mean
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: v_rms

        ! local declarations
        real(rp), dimension(:), allocatable :: tmp
        real(rp)                            :: i_sum
        integer                             :: i,j, il, err = 0
        
        allocate(tmp((ubx-lbx+1) * (uby-lby+1)), stat = err)
        if(err.ne.0) stop ' Allocation error in standard_deviation_z'

        il = 1
        do j    = lby,uby
           do i = lbx,ubx
              tmp(il) = sum( (v(i,j,sz:ez) - v_mean(i,j,sz:ez))**2 * zstep(sz:ez))
              il = il + 1
           enddo
        enddo

        il = 1
        i_sum = 1.0_rp/sum(zstep(sz:ez))
        do j    = lby,uby
           do i = lbx,ubx
              v_rms(i,j,sz:ez) = sqrt(tmp(il)*i_sum)
              il = il + 1
           enddo
        enddo
        
        deallocate(tmp)

        return
end subroutine standard_deviation_z































end module post_statistic_spatial_module
