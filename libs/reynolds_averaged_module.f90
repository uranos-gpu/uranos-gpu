module reynolds_averaged_module
implicit none

private
public reynolds_averaged

interface reynolds_averaged
  module procedure reynolds_averaged_0D
end interface

contains
subroutine reynolds_averaged_0D(v,time,time_restart,dt,Rey_av_v)
        implicit none

        integer , parameter :: rp = 8
        real(rp), parameter :: toll = 1.0E-14_rp

        real(rp), intent(in)    :: v
        real(rp), intent(in)    :: time
        real(rp), intent(in)    :: time_restart
        real(rp), intent(in)    :: dt
        real(rp), intent(inout) :: Rey_av_v
        

        if(abs(time - time_restart) < toll) then
          Rey_av_v = v

        else
          Rey_av_v = ((time-dt) * Rey_av_v + v * dt) / time

        endif

        return
end subroutine reynolds_averaged_0D




end module reynolds_averaged_module
