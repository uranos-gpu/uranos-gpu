module profiling_module
use nvtx
use parameters_module, only: profiler
implicit none


contains

subroutine StartProfRange(char)
        use nvtx
        implicit none
        character(*), intent(in) :: char

        select case(profiler)
        case('NVTX')
          #ifdef NVTX
          call nvtxStartRange('char')
          #endif

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
          #ifdef NVTX
          call nvtxEndRange
          #endif

        case('rocTX')

          ! ...

        endselect

        return
end subroutine EndProfRange


end module profiling_module
