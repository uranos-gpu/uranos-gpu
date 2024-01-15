module post_output_vtk_module
use storage_module
use post_storage_module
use post_computation_module
use vtk_utils_module
use FileModule
implicit none
private

public write_vtk, write_xz_plane_vtk
contains

subroutine write_vtk(dims,dir)
! -----------------------------------------------------------------------
!
!       This subroutine write a file in binary.vtk format
!
! -----------------------------------------------------------------------

        implicit none
        integer     , intent(in) :: dims
        character(*), intent(in) :: dir

        ! indices
        integer , dimension(3) :: id   !< indices
        integer , dimension(3) :: lo   !< lower indices
        integer , dimension(3) :: up   !< upper indices
        
        ! variabili per i file
        integer :: ifich

        type(FileType) :: vtkFile

        ! inizializzo gli indici
        id(1) = nx+2
        id(2) = ny+2
        id(3) = nz+2

        lo(:) = (/sx-1, sy-1, sz-1/)
        up(:) = (/ex+1, ey+1, ez+1/)
        
        if(dims == 2) then
          id(3) = 1
          lo(3) = 1
          up(3) = 1
        endif

        vtkFile%name = trim(output_file_name)
        vtkFile%dir  = trim(data_dir)//'/'//trim(dir)
        call OpenVTKFile(vtkFile,it)
        ifich = vtkFile%unit

        call WriteVTKHeader(vtkFile,x,y,z,id,lo,up)
                
        call print_field_vtk(lo,up,vtkFile)

        call CloseFile(vtkFile)
        return
end subroutine write_vtk


subroutine write_xz_plane_vtk(plane)
        implicit none

        integer, intent(in) :: plane

        ! indices
        integer , dimension(3) :: id   !< indices
        integer , dimension(3) :: lo   !< lower indices
        integer , dimension(3) :: up   !< upper indices

        character(len=1), parameter :: newline = achar(10)
        character(len=100)          :: s_buffer
        character(len=30)           :: buf2
        character(len=256)          :: GRID_header
        type(FileType) :: vtkFile
        integer        :: funit, nx, nz

        ! inizializzo gli indici
        id(1) = nx+2
        id(2) = plane
        id(3) = nz+2

        lo(:) = (/sx-1, plane, sz-1/)
        up(:) = (/ex+1, plane, ez+1/)

        vtkFile%name = trim(output_file_name)
        vtkFile%dir  = trim(data_dir)//'/VTK_PLANE_'//trim(str(plane))
        call OpenVTKFile(vtkFile,it)

        funit = vtkFile%unit
        nx = up(1) - lo(1) + 1
        nz = up(3) - lo(3) + 1

        write(unit = funit, iostat = err) '# vtk DataFile Version 3.0' // newline
        write(unit = funit, iostat = err) 'test file' // newline
        write(unit = funit, iostat = err) 'BINARY' // newline
        write(unit = funit, iostat = err) newline

        ! === grid header
        write(s_buffer, fmt = '(A)'     , iostat = err) 'DATASET RECTILINEAR_GRID'
        write(unit = funit, iostat = err) trim(s_buffer) // newline
        write(s_buffer, FMT = '(A,3I12)', iostat = err) 'DIMENSIONS', nx, 1, nz
        write(unit = funit, iostat = err) trim(s_buffer) // newline

        ! === x coordinate
        write(buf2,'(i8)') nx
        GRID_header="X_COORDINATES "//trim(adjustl(buf2))//" double"//newline
        write(unit = funit) trim(GRID_header),real(x(lo(1):up(1)),rp),newline

        ! === y coordinate
        write(buf2,'(i8)') 1
        GRID_header="Y_COORDINATES "//trim(adjustl(buf2))//" double"//newline
        write(unit = funit) trim(GRID_header),0.0_rp,newline

        ! === y coordinate
        write(buf2,'(i8)') nz
        GRID_header="Z_COORDINATES "//trim(adjustl(buf2))//" double"//newline
        write(unit = funit) trim(GRID_header),real(z(lo(3):up(3)),rp),newline

        ! ==== writing the data
        write(unit = funit, iostat = err) newline
        write(s_buffer, FMT = '(A,I12)', iostat = err) 'POINT_DATA', nx*nz
        write(unit = funit, iostat = err) trim(s_buffer) // newline
                
        call print_field_vtk(lo,up,vtkFile)

        call CloseFile(vtkFile)

        return
end subroutine write_xz_plane_vtk



subroutine print_field_vtk(lo,up,vtkFile)
        implicit none
        integer, dimension(3), intent(in) :: lo, up
        type(FileType)       , intent(in) :: vtkFile


        if(density)             call scalarVTK(vtkFile, lo,up,'Density'    , phi,1)
        if(velocity) then
           allocate(wrk(lbx:ubx, lby:uby, lbz:ubz))
           wrk = U/u_inf
           call scalarVTK(vtkFile, lo,up,'X-Velocity'   , wrk)
           wrk = V/u_inf
           call scalarVTK(vtkFile, lo,up,'Y-Velocity'   , wrk)
           wrk = W/u_inf
           call scalarVTK(vtkFile, lo,up,'Z-Velocity'   , wrk)
           wrk = sqrt(U**2 + V**2 + W**2)/u_inf
           call scalarVTK(vtkFile, lo,up,'Velocity-magnitude'   , wrk)
           deallocate(wrk)
        endif
        if(pressure)            call scalarVTK(vtkFile, lo,up,'Pressure'   , P)
        if(temperature)         call scalarVTK(vtkFile, lo,up,'Temperature', T)
        
        if(mach_)               call vectorVTK(vtkFile, lo,up,'Mach'             , Mach_x, Mach_y, Mach_z)
        if(speed_div)           call scalarVTK(vtkFile, lo,up,'Divergency'       , div)
        if(hybrid_weno)         call scalarVTK(vtkFile, lo,up,'Shock_sensor'     , SSENSOR)
        if(sdensity)            call scalarVTK(vtkFile, lo,up,'Schlieren_density', s_density)

        if(vorticity)           call scalarVTK(vtkFile, lo,up,'Vorticity_x'        , vor_x)
        if(vorticity)           call scalarVTK(vtkFile, lo,up,'Vorticity_y'        , vor_y)
        if(vorticity)           call scalarVTK(vtkFile, lo,up,'Vorticity_z'        , vor_z)
        if(vorticity_magnitude) call scalarVTK(vtkFile, lo,up,'Vorticity_magnitude', abs_vor)

        if(hybrid_weno)         call scalarVTK(vtkFile, lo,up,'Hybrid_weno_flag', weno_flag)

        if(qcriterion)          call scalarVTK(vtkFile, lo,up,'Q_Criterion', q_criterion)

        if(MixedVelocity)       call ScalarVTK(vtkFile, lo, up,'UV',UV)
        if(MixedVelocity)       call ScalarVTK(vtkFile, lo, up,'UW',UW)


        if(charactheristic) then

           call scalarVTK(vtkFile, lo,up,'L1', L1)
           call scalarVTK(vtkFile, lo,up,'L2', L2)
           call scalarVTK(vtkFile, lo,up,'L3', L3)
           call scalarVTK(vtkFile, lo,up,'L4', L4)
           call scalarVTK(vtkFile, lo,up,'L5', L5)

        endif

        ! statistics
        if(Rey_average) then

          allocate(wrk(lbx:ubx, lby:uby, lbz:ubz))

          call scalarVTK(vtkFile, lo,up,'r_favre', Re_av%rho)

          wrk = Re_av%rhu/Re_av%rho
          call scalarVTK(vtkFile, lo,up,'u_favre', wrk)
          wrk = Re_av%rhv/Re_av%rho
          call scalarVTK(vtkFile, lo,up,'v_favre', wrk)
          wrk = Re_av%rhw/Re_av%rho
          call scalarVTK(vtkFile, lo,up,'w_favre', wrk)

          call scalarVTK(vtkFile, lo,up,'p_favre', Re_av%prs)

          call compute_urms(wrk)
          call scalarVTK(vtkFile, lo,up,'urms', wrk)

          call compute_vrms(wrk)
          call scalarVTK(vtkFile, lo,up,'vrms', wrk)

          call compute_wrms(wrk)
          call scalarVTK(vtkFile, lo,up,'wrms', wrk)

          call compute_ekt(wrk)
          call scalarVTK(vtkFile, lo,up,'ekt', wrk)

          call scalarVTK(vtkFile, lo,up,'Mach_favre', Re_av%Mac)

          deallocate(wrk)
        endif

        return
end subroutine print_field_vtk




      

end module post_output_vtk_module
