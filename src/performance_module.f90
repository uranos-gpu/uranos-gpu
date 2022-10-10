module performance_module
use mpi
use parameters_module, only: rp, data_dir, e_cpu_time, s_cpu_time

implicit none
real(rp) :: s_cfl_time,  t_cfl_time
real(rp) :: s_out_time,  t_out_time
real(rp) :: s_irk_time,  t_irk_time
real(rp) :: s_rhs_time,  t_rhs_time
real(rp) :: s_crk_time,  t_crk_time
real(rp) :: s_bcs_time,  t_bcs_time
real(rp) :: s_mpi_time,  t_mpi_time
real(rp) :: s_upd_time,  t_upd_time
real(rp) :: s_upg_time,  t_upg_time
real(rp) :: s_wcs_time,  t_wcs_time
real(rp) :: s_les_time,  t_les_time
real(rp) :: s_sts_time,  t_sts_time

real(rp) :: s_adv_time,  t_adv_time
real(rp) :: s_vis_time,  t_vis_time
real(rp) :: s_flx_time,  t_flx_time
real(rp) :: s_fly_time,  t_fly_time
real(rp) :: s_flz_time,  t_flz_time
real(rp) :: s_LFr_time,  t_LFr_time
real(rp) :: s_eix_time,  t_eix_time
real(rp) :: s_eiy_time,  t_eiy_time
real(rp) :: s_eiz_time,  t_eiz_time

real(rp) :: s_wrc_time,  t_wrc_time
real(rp) :: s_roe_time,  t_roe_time
real(rp) :: s_SUr_time,  t_SUr_time
real(rp) :: s_tOp_time,  t_tOp_time
real(rp) :: s_prA_time,  t_prA_time

real(rp) :: s_adx_time,  t_adx_time
real(rp) :: s_ady_time,  t_ady_time
real(rp) :: s_adz_time,  t_adz_time

real(rp) :: s_ibm_time,  t_ibm_time
real(rp) :: s_ibc_time,  t_ibc_time
real(rp) :: s_Mib_time,  t_Mib_time
real(rp) :: s_dlt_time,  t_dlt_time
real(rp) :: s_dld_time,  t_dld_time
real(rp) :: s_bpt_time,  t_bpt_time
real(rp) :: s_itp_time,  t_itp_time
real(rp) :: s_flr_time,  t_flr_time
real(rp) :: s_nei_time,  t_nei_time
real(rp) :: s_ucv_time,  t_ucv_time

integer(8) :: t_cfl_calls
integer(8) :: t_out_calls
integer(8) :: t_irk_calls
integer(8) :: t_rhs_calls
integer(8) :: t_crk_calls
integer(8) :: t_bcs_calls
integer(8) :: t_mpi_calls
integer(8) :: t_upd_calls
integer(8) :: t_upg_calls
integer(8) :: t_wcs_calls
integer(8) :: t_adv_calls
integer(8) :: t_vis_calls
integer(8) :: t_flx_calls
integer(8) :: t_fly_calls
integer(8) :: t_flz_calls
integer(8) :: t_LFr_calls
integer(8) :: t_wrc_calls
integer(8) :: t_roe_calls
integer(8) :: t_SUr_calls
integer(8) :: t_adx_calls
integer(8) :: t_ady_calls
integer(8) :: t_adz_calls
integer(8) :: t_tOp_calls
integer(8) :: t_prA_calls
integer(8) :: t_eix_calls
integer(8) :: t_eiy_calls
integer(8) :: t_eiz_calls
integer(8) :: t_les_calls
integer(8) :: t_sts_calls

integer(8) :: t_ibm_calls
integer(8) :: t_ibc_calls
integer(8) :: t_Mib_calls
integer(8) :: t_dlt_calls
integer(8) :: t_dld_calls
integer(8) :: t_bpt_calls
integer(8) :: t_itp_calls
integer(8) :: t_flr_calls
integer(8) :: t_nei_calls
integer(8) :: t_ucv_calls

contains

subroutine mpi_stime(stime)
        implicit none
        real(rp), intent(out) :: stime

        stime = MPI_WTIME()

        return
end subroutine mpi_stime


subroutine mpi_etime(stime,tcalls,ttime)
        implicit none
        real(rp)  , intent(in )   :: stime
        integer(8), intent(inout) :: tcalls
        real(rp)  , intent(inout) :: ttime
        real(rp)                  :: dtime

        dtime = MPI_WTIME() - stime

        ttime = ttime + dtime

        tcalls = tcalls + 1
        
        return
end subroutine mpi_etime


subroutine resume_time(rank,nprocs)
        implicit none
        integer, intent(in) :: rank
        integer, intent(in) :: nprocs


        integer             :: err = 0, funit = 99
        character(50)       :: filename
        real(rp)            :: ref_time
        
        ! ITR DATA
        real(rp)     , dimension(13) :: my_itr_time_array
        integer(8)   , dimension(13) :: my_itr_call_array
        character(60), dimension(13) :: itr_subname_array

        ! flx DATA
        real(rp)     , dimension(4) :: my_flx_time_array
        integer(8)   , dimension(4) :: my_flx_call_array
        character(60), dimension(4) :: flx_subname_array

        ! RHS DATA
        real(rp)     , dimension(12) :: my_rhs_time_array
        integer(8)   , dimension(12) :: my_rhs_call_array
        character(60), dimension(12) :: rhs_subname_array

        ! IBM DATA
        real(rp)     , dimension(9) :: my_ibm_time_array
        integer(8)   , dimension(9) :: my_ibm_call_array
        character(60), dimension(9) :: ibm_subname_array

        filename = "DATA/"//trim(data_dir)//'/performance.txt'

        if(rank == 0) then
          open(unit = funit, file = filename, iostat = err)
          if(err /= 0) stop 'Error in writing performance file'
        endif

        ref_time = (e_cpu_time - s_cpu_time)

        ! ============================================
        ! ITERATION TIME 
        ! ============================================
        my_itr_time_array( 1) = t_cfl_time
        my_itr_time_array( 2) = t_out_time
        my_itr_time_array( 3) = t_irk_time
        my_itr_time_array( 4) = t_rhs_time
        my_itr_time_array( 5) = t_crk_time
        my_itr_time_array( 6) = t_bcs_time
        my_itr_time_array( 7) = t_mpi_time
        my_itr_time_array( 8) = t_upd_time
        my_itr_time_array( 9) = t_ibm_time
        my_itr_time_array(10) = t_upg_time
        my_itr_time_array(11) = t_wcs_time
        my_itr_time_array(12) = t_sts_time
        my_itr_time_array(13) = t_les_time

        my_itr_call_array( 1) = t_cfl_calls
        my_itr_call_array( 2) = t_out_calls
        my_itr_call_array( 3) = t_irk_calls
        my_itr_call_array( 4) = t_rhs_calls
        my_itr_call_array( 5) = t_crk_calls
        my_itr_call_array( 6) = t_bcs_calls
        my_itr_call_array( 7) = t_mpi_calls
        my_itr_call_array( 8) = t_upd_calls
        my_itr_call_array( 9) = t_ibm_calls
        my_itr_call_array(10) = t_upg_calls
        my_itr_call_array(11) = t_wcs_calls
        my_itr_call_array(12) = t_sts_calls
        my_itr_call_array(13) = t_les_calls

        itr_subname_array( 1) = '!< compute_dt'
        itr_subname_array( 2) = '!< write_all'
        itr_subname_array( 3) = '!< init_runge_kutta'
        itr_subname_array( 4) = '!< rhs_navier_stokes'
        itr_subname_array( 5) = '!< runge_kutta'
        itr_subname_array( 6) = '!< set_bc_conditions'
        itr_subname_array( 7) = '!< mpi_bc_communications'
        itr_subname_array( 8) = '!< update_all'
        itr_subname_array( 9) = '!< ibm_correction'
        itr_subname_array(10) = '!< update_all_ghost'
        itr_subname_array(11) = '!< mpi_wait_conservatives'
        itr_subname_array(12) = '!< write_statistics'
        itr_subname_array(13) = '!< compute_subgrid_model'

        call write_performances_reordered(funit,rank,nprocs,'ITR',&
                                          ref_time               ,&
                                          my_itr_time_array      ,&
                                          my_itr_call_array      ,&
                                          itr_subname_array)

        ! ============================================
        ! FLUX COMPUTATION TIME 
        ! ============================================
        my_flx_time_array(1) = t_adx_time
        my_flx_time_array(2) = t_ady_time
        my_flx_time_array(3) = t_adz_time
        my_flx_time_array(4) = t_vis_time

        my_flx_call_array(1) = t_adx_calls
        my_flx_call_array(2) = t_ady_calls
        my_flx_call_array(3) = t_adz_calls
        my_flx_call_array(4) = t_vis_calls

        flx_subname_array(1) = '!< advection x'
        flx_subname_array(2) = '!< advection y'
        flx_subname_array(3) = '!< advection z'
        flx_subname_array(4) = '!< viscous_therms'

        call write_performances_reordered(funit,rank,nprocs,'FLX',&
                                          ref_time               ,&
                                          my_flx_time_array      ,&
                                          my_flx_call_array      ,&
                                          flx_subname_array)


        ! ============================================
        ! RHS COMPUTATION TIME 
        ! ============================================
        my_rhs_time_array( 1) = t_LFr_time
        my_rhs_time_array( 2) = t_flx_time
        my_rhs_time_array( 3) = t_fly_time
        my_rhs_time_array( 4) = t_flz_time
        my_rhs_time_array( 5) = t_wrc_time
        my_rhs_time_array( 6) = t_roe_time
        my_rhs_time_array( 7) = t_SUr_time
        my_rhs_time_array( 8) = t_top_time
        my_rhs_time_array( 9) = t_eix_time
        my_rhs_time_array(10) = t_eiy_time
        my_rhs_time_array(11) = t_eiz_time
        my_rhs_time_array(12) = t_prA_time

        my_rhs_call_array( 1) = t_LFr_calls
        my_rhs_call_array( 2) = t_flx_calls
        my_rhs_call_array( 3) = t_fly_calls
        my_rhs_call_array( 4) = t_flz_calls
        my_rhs_call_array( 5) = t_wrc_calls
        my_rhs_call_array( 6) = t_roe_calls
        my_rhs_call_array( 7) = t_Sur_calls
        my_rhs_call_array( 8) = t_top_calls
        my_rhs_call_array( 9) = t_eix_calls
        my_rhs_call_array(10) = t_eiy_calls
        my_rhs_call_array(11) = t_eiz_calls
        my_rhs_call_array(12) = t_prA_calls

        rhs_subname_array( 1) = '!< Lax_friedrichs'
        rhs_subname_array( 2) = '!< flux_x'
        rhs_subname_array( 3) = '!< flux_y'
        rhs_subname_array( 4) = '!< flux_z'
        rhs_subname_array( 5) = '!< weno_reconstruction'
        rhs_subname_array( 6) = '!< roe_mean_state'
        rhs_subname_array( 7) = '!< semi_upwind_reconstruction'
        rhs_subname_array( 8) = '!< energy_preserving_tilde_op'
        rhs_subname_array( 9) = '!< eigenmatrix_x'
        rhs_subname_array(10) = '!< eigenmatrix_y'
        rhs_subname_array(11) = '!< eigenmatrix_z'
        rhs_subname_array(12) = '!< compute_pressure_array'

        call write_performances_reordered(funit,rank,nprocs,'RHS',&
                                          ref_time               ,&
                                          my_rhs_time_array      ,&
                                          my_rhs_call_array      ,&
                                          rhs_subname_array)

        
        if(rank == 0) then
          close(funit)
          call execute_command_line('cat '//trim(filename))
        endif

        return
end subroutine resume_time


subroutine write_performances_reordered(funit,rank,nprocs,header,adim,time_array,call_array,name_array)
        implicit none

        real(rp)     , dimension(:), intent(in) :: time_array
        integer(8)   , dimension(:), intent(in) :: call_array
        character(60), dimension(:), intent(in) :: name_array
        real(rp)                   , intent(in) :: adim
        integer                    , intent(in) :: funit
        integer                    , intent(in) :: rank
        integer                    , intent(in) :: nprocs
        character(3)               , intent(in) :: header
        
        real(rp)  , dimension(1:size(time_array)) :: gbl_time_array
        integer(8), dimension(1:size(time_array)) :: gbl_call_array
        logical   , dimension(1:size(time_array)) :: gbl_mask_array

        integer, parameter  :: real8 = MPI_DOUBLE_PRECISION
        integer, parameter  :: int_8 = MPI_INTEGER8
        integer, parameter  :: sum   = MPI_SUM
        integer, parameter  :: comm  = mpi_comm_world

        integer             :: err = 0, index, i
        real(rp)            :: sub_time, adi_time, cumulatv, gbl_adim, sum_times
        integer(8)          :: sub_call
        character(60)       :: sub_name

        call MPI_allreduce(time_array, gbl_time_array, size(time_array), real8, sum, comm, err)
        call MPI_allreduce(call_array, gbl_call_array, size(call_array), int_8, sum, comm, err)
        call MPI_allreduce(adim      , gbl_adim     , 1                , real8, sum, comm, err)

        gbl_time_array = gbl_time_array/real(nprocs,rp)

        if(rank == 0) then
        
          write(funit,*) ' '//header//' PERFORMANCES'
          write(funit,*) ' ----------------------- '
          write(funit,*) '     TIME [s]', '     #CALLS ', '   CUM/ITR%', '  TIME/ITR%', '   SUB NAME'

          gbl_mask_array(:) = .true.
          cumulatv          = 0.0_rp
          do i = 1,size(gbl_time_array,1)

             index = maxloc(gbl_time_array,1,gbl_mask_array)
               
             ! subroutine time
             sub_time = gbl_time_array(index)

             ! subroutine calls
             sub_call = int(real(gbl_call_array(index))/real(nprocs,rp))

             ! adimensional time
             adi_time = sub_time * 100 * nprocs/ gbl_adim

             ! cumulative percentage
             cumulatv = cumulatv + adi_time

             ! subrotutine name
             sub_name = '       '//trim(name_array(index))

             write(funit,11) sub_time, sub_call, cumulatv, adi_time, sub_name

             gbl_mask_array(index) = .false.

          enddo

          write(funit,'(A)') '      -------------'

          sum_times = 0.0_rp
          do i = 1, size(gbl_time_array,1)
             sum_times = sum_times + gbl_time_array(i)
          enddo

          write(funit,'(e18.6)') sum_times
          write(funit,*)

        endif

        11 format(e18.6,I10,2f8.2,A)
        return
end subroutine write_performances_reordered


end module performance_module
