module post_stat_computation_module
use parameters_module
use post_stat_storage_module
implicit none

private
public ComputeAdditionalStats



contains
subroutine ComputeAdditionalStats
        implicit none

        selectcase(pcharAll)
        case('FFT')
          call ComputeAdditionalStats2D(vmean2D)

        case('TFT')

        endselect

        return
end subroutine ComputeAdditionalStats


subroutine ComputeAdditionalStats2D(vmean2D)
        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(in) :: vmean2D

        integer  :: i,j
        real(rp) :: r_, ir, ru, rv, rw
        real(rp) :: r_x, r_y, r_z
        real(rp) :: u_x, u_y, u_z
        real(rp) :: v_x, v_y, v_z
        real(rp) :: w_x, w_y, w_z
        real(rp) :: omx, omy, omz

        do j    = sy,ey
           do i = sx,ex
        
              r_ = vmean2D(i,j,1)
              ir = 1.0_rp/r_
              ru = vmean2D(i,j,2)
              rv = vmean2D(i,j,3)
              rw = vmean2D(i,j,4)

              ! Favre Velocity
              FavreVelU2D(i,j) = ru*ir
              FavreVelV2D(i,j) = rv*ir
              FavreVelW2D(i,j) = rw*ir

              ! Favre Mach number
              FavreMach2D(i,j) = vmean2D(i,j,19)*ir

              ! Total Viscosity
              MeanMuTot2D(i,j) = vmean2D(i,j,15) + vmean2D(i,j,16)

              ! Reynolds stress
              FavreRUU_2D(i,j) = (vmean2D(i,j,6) - ru*ru*ir)*ir
              FavreRVV_2D(i,j) = (vmean2D(i,j,7) - rv*rv*ir)*ir
              FavreRWW_2D(i,j) = (vmean2D(i,j,8) - rw*rw*ir)*ir
              FavreRUV_2D(i,j) = (vmean2D(i,j,9) - ru*rv*ir)*ir

              ! Velocity RMS
              FavreUrms2D(i,j)  = sqrt(FavreRUU_2D(i,j))
              FavreVrms2D(i,j)  = sqrt(FavreRVV_2D(i,j))
              FavreWrms2D(i,j)  = sqrt(FavreRWW_2D(i,j))

              FavreEkT_2D(i,j) = 0.5_rp*r_*(FavreUrms2D(i,j)**2 + &
                                            FavreVrms2D(i,j)**2 + &
                                            FavreWrms2D(i,j)**2)

              ! rms weno sensor
              RMS_Sens_2D(i,j) = sqrt(vmean2D(i,j,27) - vmean2D(i,j,26)**2)
              RMS_Wflg_2D(i,j) = sqrt(vmean2D(i,j,29) - vmean2D(i,j,28)**2)

              ! Vorticity components
              r_x = vmean2D(i,j,20)
              r_y = vmean2D(i,j,21)
              r_z = vmean2D(i,j,22)

              u_x = (vmean2D(i,j,23)*r_ - r_x*ru)*ir**2
              u_y = (vmean2D(i,j,24)*r_ - r_y*ru)*ir**2
              u_z = (vmean2D(i,j,25)*r_ - r_z*ru)*ir**2

              v_x = (vmean2D(i,j,26)*r_ - r_x*rv)*ir**2
              v_y = (vmean2D(i,j,27)*r_ - r_y*rv)*ir**2
              v_z = (vmean2D(i,j,28)*r_ - r_z*rv)*ir**2

              w_x = (vmean2D(i,j,29)*r_ - r_x*rw)*ir**2
              w_y = (vmean2D(i,j,30)*r_ - r_y*rw)*ir**2
              w_z = (vmean2D(i,j,31)*r_ - r_z*rw)*ir**2

              omx = w_y - v_z
              omy = u_z - w_x
              omz = v_x - u_y

              FavreOMX_2D(i,j) = omx
              FavreOMY_2D(i,j) = omy
              FavreOMZ_2D(i,j) = omz

              ! velocity gradient components
              Favre_Ux_2D(i,j) = u_x
              Favre_Vx_2D(i,j) = v_x
              Favre_Uy_2D(i,j) = u_y
              Favre_Vy_2D(i,j) = v_y

              Favre_absOM(i,j) = sqrt(omx**2 + omy**2 + omz**2)

           enddo
        enddo

        return
end subroutine ComputeAdditionalStats2D











end module post_stat_computation_module
