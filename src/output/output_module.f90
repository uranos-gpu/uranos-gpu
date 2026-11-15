module output_module
use parameters_module
use mpi_module
use storage_module
use profiling_module    , only: StartProfRange, EndProfRange
use output_restart      , only: mpi_write, mpi_write_stat1D, mpi_write_stat2D
use output_screen       , only: screen_output, compute_average_iteration_time, screen_elapse_time
use output_residuals    , only: dump_residuals
use output_probes       , only: init_probes, sample_probes, flush_probes, write_wall_map, end_probes
use output_stats_channel, only: stat_turbulent_channel
use output_stats_tgv    , only: stat_taylor_green_vortex
implicit none

private
public output_init, output_start, output_step, output_finalize

contains
subroutine output_init
        implicit none

        call init_probes

        return
end subroutine output_init

subroutine output_start
        implicit none

        str_it_time = MPI_WTIME()

        if(it /= 0 .or. .not. printres) return

        if(due(itPrb))        call sample_probes
        if(due(itWallMap))    call write_wall_map
        call write_checkpoint

        return
end subroutine output_start

subroutine output_step
        implicit none

        if(.not. printres) return

        if(mod(it-1,rsOut) == 0) call dump_residuals(it-1)

        if(due(itout))             call compute_average_iteration_time
        if(due(itPrb))             call sample_probes
        if(due(itWallMap))         call write_wall_map
        if(due(itout))             call write_checkpoint
        if(due(StOut))             call write_statistics

        return
end subroutine output_step

subroutine output_finalize
        implicit none

        if(printres .and. istop == 1 .and. mod(it,itout) /= 0) call write_checkpoint

        call flush_probes
        call end_probes

        e_cpu_time = MPI_WTIME()
        call screen_elapse_time

        return
end subroutine output_finalize

subroutine write_checkpoint
        use inflow_module, only: inflow_turbulent_dump_state, iflow_vf_old

        implicit none

        call StartProfRange('write_checkpoint')

        call flush_probes

        call screen_output

        !$acc update host(phi)
        call mpi_write

        if(turb_inflow) then
          !$acc update host(iflow_vf_old)
          call inflow_turbulent_dump_state
        endif

        call EndProfRange

        return
end subroutine write_checkpoint

subroutine write_statistics
        use onlineStats    , only: compute_1DStats, compute_2DStats

        implicit none
        integer :: funit, ios

        call StartProfRange('write_statistics')

        itStat = itStat + 1

        select case(ic)
        case('turbulent_channel')
          call stat_turbulent_channel
        case('TGV')
          call stat_taylor_green_vortex
        end select

        select case(stat_kind)
        case('TFT')
          call compute_1DStats(phi,T,GAMMA_T,VIS,nx,nz,nVAve1D,itStat, &
                               sx,ex,sy,ey,sz,ez,comm2Dxz,vmean1D)
          !$acc update host(vmean1D)
          call mpi_write_stat1D(sy,ey,ny,nvAve1D,'BINARY_STAT',vmean1D, &
                                comm1Dy, (cart_coord(1)==0 .and. cart_coord(3)==0))

        case('FFT','FTT')
          call compute_2DStats(phi,U,V,W,T,GAMMA_T,VIS,LMD,isoSGS,nz,nVAve2D,itStat, &
                               sx,ex,sy,ey,sz,ez,comm1Dz,vmean2D)
          !$acc update host(vmean2D)
          call mpi_write_stat2D(sx,ex,sy,ey,nx,ny,vmean2D, &
                                comm2Dxy, (cart_coord(3)==0))
        end select

        if(rank == root .and. (stat_kind == 'TFT' .or. stat_kind == 'FFT' .or. stat_kind == 'FTT')) then
          open(newunit=funit, file='DATA/'//trim(data_dir)//'/BINARY_STAT/stats_index.txt', &
               position='append', action='write', iostat=ios)
          if(ios == 0) then
            write(funit,'(2i12)') itStat, it
            close(funit)
          endif
        endif

        call EndProfRange

        return
end subroutine write_statistics

logical function due(n)
        implicit none
        integer, intent(in) :: n

        due = n > 0
        if(due) due = mod(it,n) == 0

        return
end function due

end module output_module
