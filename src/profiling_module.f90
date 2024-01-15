module profiling_module
use nvtx
use roctx
implicit none


contains

subroutine StartProfRange(char)
        use nvtx
        use roctx
        implicit none
        character(*), intent(in) :: char

#ifdef NVTX
        call nvtxStartRange(char)
#endif

#ifdef ROCTX
        call roctxStartRange(char)
#endif

        return
end subroutine StartProfRange


subroutine EndProfRange
        implicit none

#ifdef NVTX
        call nvtxEndRange
#endif

#ifdef ROCTX
        call roctxEndRange
#endif


        return
end subroutine EndProfRange


end module profiling_module
