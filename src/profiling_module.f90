module profiling_module
use nvtx
implicit none


contains

subroutine StartProfRange(char)
        use nvtx
        implicit none
        character(*), intent(in) :: char

        #ifdef NVTX
        call nvtxStartRange(char)
        #endif

        #ifdef rocTX

        ! ...

        #endif

        return
end subroutine StartProfRange


subroutine EndProfRange
        implicit none

        #ifdef NVTX
        call nvtxEndRange
        #endif

        #ifdef rocTX

        ! ...

        #endif


        return
end subroutine EndProfRange


end module profiling_module
