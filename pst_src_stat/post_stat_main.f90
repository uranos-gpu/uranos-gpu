program post_stat_main
!================================================================================================
! Starting development: 25/05/2021, Venice, Italy.
!
!██████╗  ██████╗ ███████╗████████╗   ██╗   ██╗██████╗  █████╗ ███╗   ██╗ ██████╗ ███████╗    ███████╗████████╗ █████╗ ████████╗███████╗
!██╔══██╗██╔═══██╗██╔════╝╚══██╔══╝   ██║   ██║██╔══██╗██╔══██╗████╗  ██║██╔═══██╗██╔════╝    ██╔════╝╚══██╔══╝██╔══██╗╚══██╔══╝██╔════╝
!██████╔╝██║   ██║███████╗   ██║█████╗██║   ██║██████╔╝███████║██╔██╗ ██║██║   ██║███████╗    ███████╗   ██║   ███████║   ██║   ███████╗
!██╔═══╝ ██║   ██║╚════██║   ██║╚════╝██║   ██║██╔══██╗██╔══██║██║╚██╗██║██║   ██║╚════██║    ╚════██║   ██║   ██╔══██║   ██║   ╚════██║
!██║     ╚██████╔╝███████║   ██║      ╚██████╔╝██║  ██║██║  ██║██║ ╚████║╚██████╔╝███████║    ███████║   ██║   ██║  ██║   ██║   ███████║
!╚═╝      ╚═════╝ ╚══════╝   ╚═╝       ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝ ╚═════╝ ╚══════╝    ╚══════╝   ╚═╝   ╚═╝  ╚═╝   ╚═╝   ╚══════╝
!
! Post-Uranos Stat is a post-treatment code expecially developed in order to post treat 
! the solutions obtained with URANOS and get statistics.
! The code is entirely written in Fortran90 by Francesco De Vanna. 
!
! To run the code ask to: fra.devanna@gmail.com
!================================================================================================

use parameters_module
use post_stat_input_module
use GetRetau_module
use post_stat_input_module
use post_stat_storage_module
use post_stat_output_module
use post_stat_computation_module
use m_npy

        implicit none
        character(100) :: path

        ! read input file and main directory
        call ReadInputFile
        call init_Reynolds
        call ReadInputDirectoryContent('BINARY_STAT',MainDir)
        if(wmles) call ReadInputDirectoryContent('BINARY_WMLES_STAT',WmlesDir)
        call GetRestartFile(MainDir,file_ini,file_end)
        call GetPeriodicDirs
        call InitStatsFields
        call init_grid

        if(NPY) then
            path = 'DATA/'//trim(data_dir)//'/NPY/' 
            call save_npy(trim(path)//'x.npy',x(sx:ex))
            call save_npy(trim(path)//'y.npy',y(sy:ey))
            call save_npy(trim(path)//'z.npy',z(sz:ez))
        endif

        ! init variable
        do f = file_ini, file_end, skip_file
          
           ifile = trim(MainDir%name)//'/'//trim(MainDir%content(f))
           call ReadBinaryStatFile(ifile)
        
           if(wmles) then
             ifile_wmles = trim(WmlesDir%name)//'/'//trim(WmlesDir%content(f))
             call ReadBinaryWmlesStatFile(ifile_wmles)
           endif

           call ComputeAdditionalStats

           call WriteScreenOutput

           call WriteVTKStats

           if(NPY) call WriteNPYStats

           call WriteTecplotStats

           call WriteTXTStats

        enddo

        call DestroyStatsFields

end program post_stat_main
