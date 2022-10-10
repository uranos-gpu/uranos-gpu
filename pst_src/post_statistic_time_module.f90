module post_statistic_time_module

implicit none
private
public Rey_averaged, Fav_averaged, root_mean_square

interface Rey_averaged
  module procedure Rey_averaged1D, Rey_averaged2D, Rey_averaged3D
end interface Rey_averaged

interface Fav_averaged
  module procedure Fav_averaged3D
end interface Fav_averaged

interface root_mean_square
  module procedure v_root_mean_square_0D
end interface root_mean_square

contains
subroutine Rey_averaged1D(file_id,file_ini,dt,delta_time,v,Re_av_v)
! ----------------------------------------------------------------------
!       This subroutine compute the Reynolds avagered of a scalar variable v
! ----------------------------------------------------------------------
        use parameters_module, only: rp

        implicit none
        integer , intent(in)    :: file_id
        integer , intent(in)    :: file_ini
        real(rp), intent(in)    :: dt
        real(rp), intent(in)    :: delta_time
        real(rp), intent(in)    :: v
        real(rp), intent(inout) :: Re_av_v

        if(file_id == file_ini) then 
          Re_av_v = v

        else
          Re_av_v = ((delta_time-dt) * Re_av_v + v * dt) / delta_time

        endif
        return
end subroutine Rey_averaged1D


subroutine Rey_averaged2D(file_id,file_ini,dt,delta_time,v,Re_av_v)
! ----------------------------------------------------------------------
!       This subroutine compute the Reynolds avagered of a 3D variable v
! ----------------------------------------------------------------------
        use parameters_module, only: rp

        implicit none
        integer                              , intent(in)    :: file_id
        integer                              , intent(in)    :: file_ini
        real(rp)                             , intent(in)    :: dt
        real(rp)                             , intent(in)    :: delta_time
        real(rp), dimension(:,:), allocatable, intent(in)    :: v
        real(rp), dimension(:,:), allocatable, intent(inout) :: Re_av_v

        if(file_id == file_ini) then 
                
                !$omp workshare
                Re_av_v(:,:) = v(:,:)
                !$omp end workshare

        else

                !$omp workshare
                Re_av_v(:,:) = ((delta_time-dt) * Re_av_v(:,:) + v(:,:) * dt) / delta_time
                !$omp end workshare

        endif

        return
end subroutine Rey_averaged2D


subroutine Rey_averaged3D(file_id,file_ini,dt,delta_time,v,Re_av_v)
! ----------------------------------------------------------------------
!       This subroutine compute the Reynolds avagered of a 3D variable v
! ----------------------------------------------------------------------
        use parameters_module, only: rp

        implicit none
        integer                                , intent(in)    :: file_id
        integer                                , intent(in)    :: file_ini
        real(rp)                               , intent(in)    :: dt
        real(rp)                               , intent(in)    :: delta_time
        real(rp), dimension(:,:,:), allocatable, intent(in)    :: v
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: Re_av_v

        if(file_id == file_ini) then 
                
                !$omp workshare
                Re_av_v(:,:,:) = v(:,:,:)
                !$omp end workshare

        else

                !$omp workshare
                Re_av_v(:,:,:) = ((delta_time-dt) * Re_av_v(:,:,:) + v(:,:,:) * dt) / delta_time
                !$omp end workshare

        endif

        return
end subroutine Rey_averaged3D

subroutine Fav_averaged3D(file_id,file_ini,dt,delta_time,v,Fa_av_v)
! ------------------------------------------------------
!       This subroutine compute the Favre average of a
!       variable v
! ------------------------------------------------------
        use parameters_module, only: rp
        use storage_module   , only: phi

        implicit none
        integer                                , intent(in)    :: file_id
        integer                                , intent(in)    :: file_ini
        real(rp)                               , intent(in)    :: dt
        real(rp)                               , intent(in)    :: delta_time
        real(rp), dimension(:,:,:), allocatable, intent(in)    :: v
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: Fa_av_v

        if(file_id == file_ini) then 
                        
                !$omp workshare
                Fa_av_v = phi(:,:,:,1)*v(:,:,:)
                !$omp end workshare

        else

                !$omp workshare
                Fa_av_v(:,:,:) = ((delta_time-dt) * Fa_av_v(:,:,:) + phi(:,:,:,1) * v(:,:,:) * dt) / delta_time
                !$omp end workshare

        endif
        
        return
end subroutine Fav_averaged3D


subroutine v_root_mean_square_0D(file_id,file_ini,dt,delta_time,v,v_mean,v_rms)
! ------------------------------------------------------------
!       Comnputation of the RMS of a scalar variable v
! ------------------------------------------------------------

        use parameters_module, only: rp

        implicit none
        integer , intent(in)    :: file_id
        integer , intent(in)    :: file_ini
        real(rp), intent(in)    :: dt
        real(rp), intent(in)    :: delta_time
        real(rp), intent(in)    :: v
        real(rp), intent(in)    :: v_mean
        real(rp), intent(inout) :: v_rms

        !local
        real(rp) :: v_rms_old, v_rms_new

        if(file_id == file_ini) then
          v_rms = 0.0_rp

        else

          v_rms_old = (delta_time-dt) * v_rms**2
          v_rms_new = (v-v_mean)**2 * dt
        
          v_rms = sqrt((v_rms_old + v_rms_new)/delta_time)

        endif

        return
end subroutine v_root_mean_square_0D



end module post_statistic_time_module
