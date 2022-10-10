module post_statistic_module
use storage_module
use post_storage_module

implicit none
private
public compute_statistic
contains

subroutine compute_statistic
        implicit none
        integer :: err = 0
        
        if(Rey_average) then

          allocate(wrk(lbx:ubx, lby:uby, lbz:ubz), stat = err)
          if(err.ne.0) stop 'Allocation error in compute_statistic'

          wrk = phi(:,:,:,1)
          call global_averaged(wrk,Re_av%rho)

          wrk = phi(:,:,:,2)
          call global_averaged(wrk,Re_av%rhu)

          wrk = phi(:,:,:,3)
          call global_averaged(wrk,Re_av%rhv)

          wrk = phi(:,:,:,4)
          call global_averaged(wrk,Re_av%rhw)

          wrk = phi(:,:,:,2)*phi(:,:,:,2)/phi(:,:,:,1)
          call global_averaged(wrk,Re_av%ruu)

          wrk = phi(:,:,:,3)*phi(:,:,:,3)/phi(:,:,:,1)
          call global_averaged(wrk,Re_av%rvv)

          wrk = phi(:,:,:,4)*phi(:,:,:,4)/phi(:,:,:,1)
          call global_averaged(wrk,Re_av%rww)

          wrk = phi(:,:,:,2)*phi(:,:,:,3)/phi(:,:,:,1)
          call global_averaged(wrk,Re_av%ruv)

          ! Mach favre
          wrk = sqrt((U**2 + V**2 + W**2)/(gamma0*T))
          call global_averaged(wrk,Re_av%Mac)
          deallocate(wrk)


          call global_averaged(P,Re_av%prs)
          call global_averaged(T,Re_av%Tmp)
          if(viscous) call global_averaged(VIS,Re_av%VIS)

        endif


        return
end subroutine compute_statistic


subroutine global_averaged(v,v_av)
! -------------------------------------------------------------------
!       This subroutine performs the spatial averages in respect of 
!       the periodic directions and the time average if it is needed
! -------------------------------------------------------------------
        use parameters_module            , only: rp
        use mpi_module                   , only: lbx,ubx, lby,uby, lbz,ubz
        use post_storage_module          , only: f, file_ini, file_dt, periodic_dir, symmtric_dir, time_from_restart
        use post_statistic_time_module   , only: Rey_averaged
        use post_statistic_spatial_module, only: global_spatial_average
        
        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(in)    :: v
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: v_av

        ! local declarations
        real(rp), dimension(:,:,:), allocatable :: tmp_av
        integer                                 :: err = 0

        allocate(tmp_av(lbx:ubx, lby:uby, lbz:ubz), stat = err)
        if(err.ne.0) stop 'Allocation error in global averaged.'
        
        ! ========== Compute global space average ============= !
        call global_spatial_average(v,periodic_dir,symmtric_dir,tmp_av)

        ! ============== Compute time average ================= !
        call Rey_averaged(f,file_ini,file_dt,time_from_restart,tmp_av,v_av)

        deallocate(tmp_av)
        
        return
end subroutine global_averaged



!subroutine global_rootmsqr(v,v_mean,v_rms,type_time_average)
!
!        use parameters_module            , only: rp
!        use post_storage_module          , only: time_from_restart
!        use mpi_module                   , only: lbx,ubx, lby,uby, lbz,ubz
!        use post_statistic_spatial_module, only: global_spatial_standard_deviation
!        use post_statistic_time_module   , only: Rey_averaged, Fav_averaged
!
!        implicit none
!        character(*)                           , intent(in)    :: type_time_average
!        real(rp), dimension(:,:,:), allocatable, intent(in)    :: v
!        real(rp), dimension(:,:,:), allocatable, intent(in)    :: v_mean
!        real(rp), dimension(:,:,:), allocatable, intent(inout) :: v_rms
!
!        ! local declaration
!        real(rp), dimension(:,:,:), allocatable :: tmp1
!        integer                                 :: err = 0
!
!        allocate(tmp1(lbx:ubx, lby:uby, lbz:ubz), stat = err)
!        if(err /= 0) stop ' Allocation error in space_and_time_rms.'
!        
!        ! ============= Compute standard deviation of v ================ !
!        call global_spatial_standard_deviation(v,v_mean,periodic_dir,tmp1)
!
!        ! ============= mean in time the global space rms ============== !
!        if    (type_time_average == 'Rey_average') then
!          call Rey_averaged(f,file_ini,file_dt,time_from_restart,tmp1,v_rms)
!
!        elseif(type_time_average == 'Fav_average') then
!          call Fav_averaged(f,file_ini,file_dt,time_from_restart,tmp1,v_rms)
!
!        endif
!
!        deallocate(tmp1)
!        return
!end subroutine global_rootmsqr



end module post_statistic_module

