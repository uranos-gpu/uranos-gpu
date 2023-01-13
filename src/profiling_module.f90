module profiling_module
use nvtx
use parameters_module, only: profiler
implicit none


contains

subroutine StartProfRange(char)
        implicit none
        character(*), intent(in) :: char

        select case(profiler)
        case('NVTX')

          call nvtxStartRange(char)

        case('rocTX')

          ! ...

        case default

          print*, 'Profiler type ', trim(profiler), ' is not implemented'
          stop

        endselect


        return
end subroutine StartProfRange


subroutine EndProfRange
        implicit none

        select case(profiler)
        case('NVTX')
          call nvtxEndRange

        case('rocTX')

          ! ...

        endselect

        return
end subroutine EndProfRange


end module profiling_module
