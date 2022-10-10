module post_stat_output_module
use parameters_module
use post_stat_storage_module
use FileModule
use vtk_utils_module
use Allocate_module
use m_npy

implicit none
contains
subroutine WriteScreenOutput
        implicit none
        
        integer :: i,j
        real(rp) :: r_, ir, p_, T_, M_

        real(rp) :: rmin, pmin, Tmin, Mmin
        real(rp) :: rmax, pmax, Tmax, Mmax

        rmin = huge(0.0_rp)
        pmin = huge(0.0_rp)
        Tmin = huge(0.0_rp)
        Mmin = huge(0.0_rp)

        rmax = 0.0_rp
        pmax = 0.0_rp
        Tmax = 0.0_rp
        Mmax = 0.0_rp

        selectcase(pcharAll)
        case('FFT')
        do    j = sy,ey
           do i = sx,ex

              r_ = vmean2D(i,j,1)
              ir = 1.0_rp/r_
              p_ = vmean2D(i,j,11)
              T_ = vmean2D(i,j,13)
              M_ = vmean2D(i,j,19)*ir

              rmax = max(r_, rmax)
              pmax = max(p_, pmax)
              Tmax = max(T_, Tmax)
              Mmax = max(M_, Mmax)

              rmin = min(r_, rmin)
              pmin = min(p_, pmin)
              Tmin = min(T_, Tmin)
              Mmin = min(M_, Mmin)

           enddo
        enddo

        case('TFT')



        endselect

        if(mod(itstat,5) == 0) then
          write(*,*) '-----------------------------------------------------------',&
                     '-----------------------------------------------------------',&
                     '----------------- post Uranos Stat'
          write(*,*) '   itStat  file %      Mmin           Mmax           rmin',&
        '           rmax           pmin           pmax           Tmin           Tmax'

        endif

        write(*,10) itStat, real(f,rp)/real(file_end,rp)*100, &
                   Mmin, Mmax, &
                   rmin, rmax, &
                   pmin, pmax, &
                   Tmin, Tmax

        10 format(I10,f8.3,10(e15.3))

        return
end subroutine WriteScreenOutput



subroutine WriteVTKStats
        implicit none

        integer, dimension(2) :: lo, up
        type(FileType)        :: vtkFile
        integer               :: n

        real(rp), allocatable, dimension(:,:) :: tmp_2D

        vtkFile%name = trim(output_file_name)
        vtkFile%dir  = trim(data_dir)//'/VTK_STATS'
        call OpenVTKFile(vtkFile,itStat)

        selectcase(pcharAll)

          case('FFT')

          lo = (/1,1/)
          up = (/nx,ny/)

          call AllocateReal(tmp_2D,lbx,ubx,lby,uby)
        
          call WriteVTKHeader(vtkFile,x,y,lo,up)
       
          call scalarVTK(vtkFile, lo, up, 'mean_density', vmean2D, 1)

          call scalarVTK(vtkFile, lo, up, 'mean_pressure', vmean2D, 11)

          tmp_2D = FavreVelU2D/u_inf
          call scalarVTK(vtkFile, lo, up, 'favre_x-velocity', tmp_2D)
          !
          tmp_2D = FavreVelV2D/u_inf
          call scalarVTK(vtkFile, lo, up, 'favre_y-velocity', tmp_2D)
          !
          tmp_2D = FavreVelW2D/u_inf
          call scalarVTK(vtkFile, lo, up, 'favre_z-velocity', tmp_2D)
          !
          tmp_2D = sqrt(FavreVelU2D**2 + FavreVelV2D**2 + FavreVelW2D**2)/u_inf
          call scalarVTK(vtkFile, lo, up, 'favre_velocity_module', tmp_2D)

          call scalarVTK(vtkFile, lo, up, 'favre_mach', FavreMach2D)
          call scalarVTK(vtkFile, lo, up, 'mean_total_viscosity', MeanMuTot2D)
          call scalarVTK(vtkFile, lo, up, 'mean_laminar_viscosity', vmean2D,15)
          call scalarVTK(vtkFile, lo, up, 'mean_turbulent_viscosity', vmean2D,16)

          call scalarVTK(vtkFile, lo, up, 'mean_rhou"u"', FavreRUU_2D)
          call scalarVTK(vtkFile, lo, up, 'mean_rhov"v"', FavreRVV_2D)
          call scalarVTK(vtkFile, lo, up, 'mean_rhow"w"', FavreRWW_2D)
          call scalarVTK(vtkFile, lo, up, 'mean_rhou"v"', FavreRUV_2D)
        
          tmp_2D = FavreUrms2D/u_inf
          call scalarVTK(vtkFile, lo, up, 'u_rms', tmp_2D)
          tmp_2D = FavreVrms2D/u_inf
          call scalarVTK(vtkFile, lo, up, 'v_rms', tmp_2D)
          tmp_2D = FavreWrms2D/u_inf
          call scalarVTK(vtkFile, lo, up, 'w_rms', tmp_2D)
          tmp_2D = FavreRUV_2D/u_inf**2
          call scalarVTK(vtkFile, lo, up, 'uvrms', tmp_2D)

          tmp_2D = FavreEkT_2D/(0.5_rp*u_inf**2)
          call scalarVTK(vtkFile, lo, up, 'mean_turbulent_kin_energy', tmp_2D)

          call scalarVTK(vtkFile, lo, up, 'mean_sensor', vmean2D,26)
          call scalarVTK(vtkFile, lo, up, 'mean_weno_flag', vmean2D,28)

          call scalarVTK(vtkFile, lo, up, 'rms_sensor', RMS_Sens_2D)
          call scalarVTK(vtkFile, lo, up, 'rms_weno_flag', RMS_Wflg_2D)

          call CloseFile(vtkFile)
        
          call DeallocateReal(tmp_2D)

        endselect


        return
end subroutine WriteVTKStats



subroutine WriteNPYStats
        implicit none
        character(100) :: path

        path = 'DATA/'//trim(data_dir)//'/NPY/' 
        call save_npy(trim(path)//'vmean2D'//trim(str(itStat))//'.npy',vmean2D(sx:ex,sy:ey,:))

        return
end subroutine WriteNPYStats










subroutine WriteTecplotStats
        implicit none
        type(FileType) :: tecPlotFile
        integer :: i,j

        tecPlotFile%name = 'wavplot'
        tecPlotFile%dir  = trim(data_dir)//'/TECPLOT_STATS'
        call OpenNewFile(tecplotFile,itStat)

        write(tecPlotFile%unit,*) 'zone i=',ex,', j=',ey
        do j = sy,ey
           do i = sx,ex
              write(tecPlotFile%unit,*) x(i),y(j), &
                vmean2D(i,j,1)          , &
                FavreVelU2D(i,j)/u_inf  , &
                FavreVelV2D(i,j)/u_inf  , &
                FavreVelW2D(i,j)/u_inf  , &
                !
                sqrt(FavreRUU_2D(i,j))/u_inf, &
                sqrt(FavreRVV_2D(i,j))/u_inf, &
                sqrt(FavreRWW_2D(i,j))/u_inf
           enddo
        enddo
        call CloseFile(tecPlotFile)


       return
end subroutine WriteTecplotStats




subroutine WriteTXTStats

        use post_stat_output_bl_module

        implicit none

        select case(ic)
          case('turbulent_BL','swbli')
          call write_bl_streamwise_stats
          call write_bl_wall_normal_stats
          if(wmles) call write_bl_wall_normal_stats_wmles

        endselect

        return
end subroutine WriteTXTStats















end module post_stat_output_module
