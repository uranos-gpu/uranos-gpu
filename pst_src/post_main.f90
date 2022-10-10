!================================================================================================
! Starting development: 12/10/2018, Venice, Italy.
!
!  ██████╗  ██████╗ ███████╗████████╗   ██╗   ██╗██████╗  █████╗ ███╗   ██╗ ██████╗ ███████╗
!  ██╔══██╗██╔═══██╗██╔════╝╚══██╔══╝   ██║   ██║██╔══██╗██╔══██╗████╗  ██║██╔═══██╗██╔════╝
!  ██████╔╝██║   ██║███████╗   ██║█████╗██║   ██║██████╔╝███████║██╔██╗ ██║██║   ██║███████╗
!  ██╔═══╝ ██║   ██║╚════██║   ██║╚════╝██║   ██║██╔══██╗██╔══██║██║╚██╗██║██║   ██║╚════██║
!  ██║     ╚██████╔╝███████║   ██║      ╚██████╔╝██║  ██║██║  ██║██║ ╚████║╚██████╔╝███████║
!  ╚═╝      ╚═════╝ ╚══════╝   ╚═╝       ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝ ╚═════╝ ╚══════╝
!
! Post-Uranos is a post-treatment code expecially developed in order to post treat 
! the solutions obtained with URANOS
! The code is entirely written in Fortran90 by Francesco De Vanna. 
!
! To run the code ask to: fra.devanna@gmail.com
!================================================================================================
program main_post
!$ use omp_lib
use parameters_module
use mesh_module
use inflow_module
use post_input_module
use post_storage_module
use post_computation_module
use post_output_module
use post_bc_module
use post_solutions_module
use post_statistic_module
use sgs_module
use GetRetau_module
use m_npy

        implicit none
        character(100) :: path
       
        ! read input file and main dir
        call read_input_file
        call read_input_dir
        call make_saving_directories
        call init_Reynolds

        ! look for symmetrical direction
        call look_for_symmetry

        ! init variables
        call init_cons_fields
        call init_inst_fields
        call init_exct_fields
        call init_stat_fields
        call init_FD_coefficients
        call init_grid
        if(inflow) call init_inflow_profile

        if(NPY) then 
            path = 'DATA/'//trim(data_dir)//'/NPY/' 
            call save_npy(trim(path)//'x.npy',x(sx:ex))
            call save_npy(trim(path)//'y.npy',y(sy:ey))
            call save_npy(trim(path)//'z.npy',z(sz:ez))
        endif
        
#ifdef OPENMP
        s_cpu_time = omp_get_wtime()
#else
        call cpu_time(s_cpu_time)
#endif
        do f = file_ini, size(main_dir%content), skip_file
        
                ifile = trim(main_dir%name)//'/'//trim(main_dir%content(f))

                call read_binary_file(ifile)
                call apply_pst_bc_conditions

                call compute_fields
                call compute_exact_fields
                if(LES) call compute_subgrid_model
                call compute_statistic
                call write_post_treatments
        
        enddo
#ifdef OPENMP
        e_cpu_time = omp_get_wtime()
#else
        call cpu_time(e_cpu_time)
#endif
        print*, '-----------------------------------------'
        write(*,'(A,f18.6,A)') " cpu time = ", e_cpu_time-s_cpu_time, ' s'
        print*, '-----------------------------------------'

        call destroy_variables

end program main_post
