module output_module
use parameters_module
use mpi_module
use storage_module
use norm_module
implicit none

private
public write_restart, write_restart_stat, write_statistics, &
       screen_elapse_time, mpi_plot, screen_output, &
       compute_average_iteration_time

interface mpi_plot
  module procedure mpi_plot_scl, mpi_plot_vec 
end interface 


contains 
subroutine write_restart

        implicit none

        !$acc update host(phi,U,V,W,P,T,LMD,VIS)

        call screen_output
        call mpi_write

        return
end subroutine write_restart

subroutine write_restart_stat

        use onlineStats

        implicit none
        character(1), dimension(3) :: pchar
        character(3)               :: pcharAll

        !$acc update host(phi,U,V,W,P,T,LMD,VIS)

        itStat = itStat + 1

        pchar = 'F'
        if(periods(1)) pchar(1) = 'T'
        if(periods(2)) pchar(2) = 'T'
        if(periods(3)) pchar(3) = 'T'
        pcharAll = pchar(1)//pchar(2)//pchar(3)
        
        ! write binary statistics
        selectcase(pcharAll)

        case('TFT') ! xz are periodic
          call compute_1DStats(phi,VIS,nx,nz,nVAve1D,itStat, &
                            sx,ex,sy,ey,sz,ez,comm2Dxz,vmean1D,mpi_flag)

          call mpi_write_stat1D(sy,ey,ny,nvAve1D,'BINARY_STAT',vmean1D)

        case('FFT') ! z  is periodid
          call compute_2DStats(phi,VIS,nz,nVAve2D,itStat, &
                            sx,ex,sy,ey,sz,ez,comm1Dz,vmean2D,mpi_flag)

          call mpi_write_stat2D(sx,ex,sy,ey,nx,ny,vmean2D)
          if(wmles) then
            call mpi_write_stat1D(sx,ex,nx,nVWMLESData,&
                    'BINARY_WMLES_STAT',vmean1D_wmles)
          endif
        case('FTT') ! yz are periodic
          call compute_2DStats(phi,VIS,nz,nVAve2D,itStat, &
                            sx,ex,sy,ey,sz,ez,comm1Dz,vmean2D,mpi_flag)
          call mpi_write_stat2D(sx,ex,sy,ey,nx,ny,vmean2D)

        endselect

        return
end subroutine write_restart_stat



subroutine mpi_write
        implicit none
        
        ! local declarations
        character(dl)                       :: filename, iteration
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer                             :: fh
        integer(kind = MPI_OFFSET_KIND)     :: initial_displacement
        integer, parameter                  :: array_rank=3
        integer, dimension(array_rank)      :: shape_array, shape_sub_array, start_coord
        integer, dimension(array_rank)      :: shape_view_array, shape_sub_view_array, start_view_coord
        integer                             :: type_sub_array, type_sub_view_array

        ! init filename
        write(iteration, '(i7.7)') it
        filename = "DATA/"//trim(data_dir)//"/BINARY/"//trim(output_file_name)//trim(iteration)//'.bin'

        CALL MPI_FILE_OPEN(mpi_comm_cart, filename, &
                 MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                 MPI_INFO_NULL, fh, mpi_err)

        if(mpi_err .ne. 0) then
          if(rank == root) write(*,'(A,A)') 'MPI_IO error in writing file ', trim(filename)
          call secure_stop
        endif
        
        ! subarrays shape
        shape_array     = SHAPE(phi(:,:,:,1))
        shape_sub_array = SHAPE(phi(sx:ex, sy:ey, sz:ez,1))
        start_coord     = (/GN, GN, GN/)

        call mpi_type_create_subarray(array_rank, shape_array, shape_sub_array, start_coord, &
                                      MPI_ORDER_FORTRAN, MPI_RP              , &
                                      type_sub_array, mpi_err)
        call mpi_type_commit(type_sub_array, mpi_err)

        ! view array shape
        shape_view_array     = (/nx, ny, nz/)
        shape_sub_view_array = SHAPE(phi(sx:ex, sy:ey, sz:ez,1))
        start_view_coord     =  (/sx-1 , sy-1 , sz-1/)

        call mpi_type_create_subarray(array_rank, shape_view_array, shape_sub_view_array, start_view_coord, &
                                      MPI_ORDER_FORTRAN, MPI_RP                             , &
                                      type_sub_view_array, mpi_err)
        call mpi_type_commit(type_sub_view_array, mpi_err)

        ! write domain size
        CALL MPI_FILE_WRITE_ALL(fh, nx , 1, MPI_INTEGER, status, mpi_err)
        CALL MPI_FILE_WRITE_ALL(fh, ny , 1, MPI_INTEGER, status, mpi_err)
        CALL MPI_FILE_WRITE_ALL(fh, nz , 1, MPI_INTEGER, status, mpi_err)

        ! write time, dt and iteration
        CALL MPI_FILE_WRITE_ALL(fh, it   , 1, MPI_INTEGER         , status, mpi_err)
        CALL MPI_FILE_WRITE_ALL(fh, time , 1, MPI_RP, status, mpi_err)
        CALL MPI_FILE_WRITE_ALL(fh, dt   , 1, MPI_RP, status, mpi_err)

        ! start writing conservative variables
        initial_displacement = 4*ip + 2*rp
        CALL MPI_FILE_SET_VIEW(fh, initial_displacement, MPI_RP, &
             type_sub_view_array, "native", MPI_INFO_NULL, mpi_err)
        
        CALL MPI_FILE_WRITE_ALL(fh, phi(:,:,:,1), 1, type_sub_array, status, mpi_err)
        CALL MPI_FILE_WRITE_ALL(fh, phi(:,:,:,2), 1, type_sub_array, status, mpi_err)
        CALL MPI_FILE_WRITE_ALL(fh, phi(:,:,:,3), 1, type_sub_array, status, mpi_err)
        CALL MPI_FILE_WRITE_ALL(fh, phi(:,:,:,4), 1, type_sub_array, status, mpi_err)
        CALL MPI_FILE_WRITE_ALL(fh, phi(:,:,:,5), 1, type_sub_array, status, mpi_err)

        CALL MPI_FILE_CLOSE(fh, mpi_err)

        call check_mpi('mpi_write', mpi_err)

        return
end subroutine mpi_write




subroutine compute_average_iteration_time
        implicit none

        end_it_time = MPI_WTIME()
        end_it_time = (end_it_time - str_it_time)/real(itOut,rp)
        str_it_time = MPI_WTIME()

        return
end subroutine compute_average_iteration_time






subroutine screen_output
! -------------------------------------------------------------
!
!       This subroutine prints screen informations
!
! -------------------------------------------------------------
        implicit none
        real(rp), parameter :: DivCheckToll = 1.0E-14_rp
        real(rp)            :: gbl_r_max, gbl_T_max, gbl_p_max, gbl_M_max
        real(rp)            :: lcl_r_max, lcl_T_max, lcl_p_max, lcl_M_max
        real(rp)            :: gbl_r_min, gbl_T_min, gbl_p_min
        real(rp)            :: lcl_r_min, lcl_T_min, lcl_p_min
        real(rp)            :: r_, ir, u_, v_, w_, p_, T_
        integer             :: i,j,k, err

        err = 0

        lcl_M_max = 0.0_rp
        lcl_r_max = 0.0_rp
        lcl_T_max = 0.0_rp
        lcl_p_max = 0.0_rp

        lcl_r_min = huge(0.0_rp)
        lcl_T_min = huge(0.0_rp)
        lcl_p_min = huge(0.0_rp)

        selectcase(dims)
          !
          !  === 2D CASE
          !
          case(2)
          k = 1
          do j    = sy,ey
             do i = sx,ex

                r_ = phi (i,j,k,1)
                ir = 1.0_rp/r_
                u_ = phi(i,j,k,2)*ir
                v_ = phi(i,j,k,3)*ir

                p_ = cv_i * (phi(i,j,k,5) - 0.5_rp*r_*(u_*u_ + v_*v_))
                T_ = p_*ir
                
                ! max quantities
                lcl_r_max = max(r_,lcl_r_max)
                lcl_p_max = max(p_,lcl_p_max)
                lcl_T_max = max(T_,lcl_T_max) 
                lcl_M_max = max(((u_*u_ + v_*v_)/(gamma0*T_)), lcl_M_max)

                ! min quantities
                lcl_r_min = min(r_,lcl_r_min)
                lcl_p_min = min(p_,lcl_p_min)
                lcl_T_min = min(T_,lcl_T_min) 

             enddo
          enddo
          !
          !  === 3D CASE
          !
          case(3)
!          !$acc parallel default(present)
!          !$acc loop gang, vector collapse(3) &
!          !$acc reduction(max:lcl_r_max,lcl_p_max,lcl_T_max,lcl_M_max) &
!          !$acc reduction(min:lcl_r_min,lcl_p_min,lcl_T_min)
          do k       = sz,ez
             do j    = sy,ey
                do i = sx,ex

                   r_ = phi (i,j,k,1)
                   ir = 1.0_rp/r_
                   u_ = phi(i,j,k,2)*ir
                   v_ = phi(i,j,k,3)*ir
                   w_ = phi(i,j,k,4)*ir

                   p_ = cv_i * (phi(i,j,k,5) - 0.5_rp*r_*(u_*u_ + v_*v_ + w_*w_))
                   T_ = p_*ir
                
                   ! max quantities
                   lcl_r_max = max(r_,lcl_r_max)
                   lcl_p_max = max(p_,lcl_p_max)
                   lcl_T_max = max(T_,lcl_T_max) 
                   lcl_M_max = max(((u_*u_ + v_*v_ + w_*w_)/(gamma0*T_)), lcl_M_max)

                   ! min quantities
                   lcl_r_min = min(r_,lcl_r_min)
                   lcl_p_min = min(p_,lcl_p_min)
                   lcl_T_min = min(T_,lcl_T_min) 

                enddo
             enddo
          enddo
!          !$acc end parallel

        endselect

        lcl_M_max = sqrt(lcl_M_max)

        call MPI_allreduce(lcl_T_max, gbl_T_max, 1, MPI_RP, MPI_MAX, mpi_comm_cart, err)
        call MPI_allreduce(lcl_r_max, gbl_r_max, 1, MPI_RP, MPI_MAX, mpi_comm_cart, err)
        call MPI_allreduce(lcl_M_max, gbl_M_max, 1, MPI_RP, MPI_MAX, mpi_comm_cart, err)
        call MPI_allreduce(lcl_p_max, gbl_p_max, 1, MPI_RP, MPI_MAX, mpi_comm_cart, err)

        call MPI_allreduce(lcl_T_min, gbl_T_min, 1, MPI_RP, MPI_MIN, mpi_comm_cart, err)
        call MPI_allreduce(lcl_r_min, gbl_r_min, 1, MPI_RP, MPI_MIN, mpi_comm_cart, err)
        call MPI_allreduce(lcl_p_min, gbl_p_min, 1, MPI_RP, MPI_MIN, mpi_comm_cart, err)

        ! ---- write screen output ---- !
        if(rank == root) then
          if(mod(it,10*itout)==0) then
            write(*,'(4x,11(A,5x))') ' it', 'time %', ' dt ', '       M max',              &
                                                              '     r max'  ,'     r min', &
                                                              '     T max'  ,'     T min', &
                                                              '     p max'  ,'     p min', &
                                                              '     titer'
            write(*,*) ' -----------------------------------------------------------------------------------', &
                       '---------------------------------------------------------------------------'
          endif

          write(*,10) it , time/tmax*100 , dt, gbl_M_max,            &
                                               gbl_r_max, gbl_r_min, &
                                               gbl_T_max, gbl_T_min, &
                                               gbl_p_max, gbl_p_min, &
                                               end_it_time
        endif

        if(gbl_r_min < DivCheckToll .or. &
           gbl_T_min < DivCheckToll .or. &
           dt        < DivCheckToll) then
           if(rank == root) print*, ' Divergency has been detected!'
           istop = 1
        endif

        10 format(I10,f8.3,10(e15.3))

        return                
end subroutine screen_output

subroutine screen_elapse_time
        implicit none
        
        if(rank == root) then
          print*, ' -----------------------------------------'
          write(*,'(A,f18.6,A)') " cpu time = ", e_cpu_time-s_cpu_time, ' s'
          print*, ' -----------------------------------------'
        endif
        
        return
end subroutine screen_elapse_time





subroutine write_statistics

        use output_tch

        implicit none

        !$acc update host(phi,U,V,W,P,T,LMD,VIS)
        !$acc update host(WMLES_DATA)

        ! updata statistical iterations
        !itStat = itStat + 1

        selectcase(ic)

          case('turbulent_channel','poiseuille_x')
          call stat_turbulent_channel

          case('turbulent_BL','swbli')
          call stat_turbulent_boundary_layer

        endselect

        return
end subroutine write_statistics


subroutine mpi_write_stat1D(sid,eid,n,nstats,dirname,vmean)

        use FileModule

        implicit none
        integer                              , intent(in) :: sid, eid, n, nstats
        character(*)                         , intent(in) :: dirname
        real(rp), dimension(:,:), allocatable, intent(in) :: vmean

        integer, parameter :: array_rank=1
        integer, parameter :: intprec = MPI_INTEGER

        integer, dimension(array_rank)      :: shape_array
        integer, dimension(array_rank)      :: shape_sub_array
        integer, dimension(array_rank)      :: start_coord

        integer, dimension(array_rank)      :: shape_view_array
        integer, dimension(array_rank)      :: shape_sub_view_array
        integer, dimension(array_rank)      :: start_view_coord

        integer(kind = MPI_OFFSET_KIND)     :: initial_displacement
        integer, dimension(MPI_STATUS_SIZE) :: status

        character(dl) :: iteration, filename, saving_dir
        integer       :: fh, err = 0, l
        integer       :: type_sub_array, type_sub_view_array
        !
        ! === open the file
        !
        write(iteration, '(i7.7)') itStat
        saving_dir = "DATA/"//trim(data_dir)//"/"//trim(dirname)
        filename = trim(saving_dir)//"/stats"//trim(iteration)//'.bin'

        CALL MPI_FILE_OPEN(mpi_comm_cart, filename, &
                 MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                 MPI_INFO_NULL, fh, err)

        if(err.ne.0) then
          if(rank == root) then
            write(*,'(A,A)') 'MPI_IO error in writing file ', trim(filename)
          endif
          stop
        endif
        !
        ! === subarrays 
        !
        shape_array     = SHAPE(vmean(:,1))
        shape_sub_array = SHAPE(vmean(sid:eid,1))
        start_coord     = (/GN/)

        call mpi_type_create_subarray(array_rank, shape_array, shape_sub_array, &
                                      start_coord, MPI_ORDER_FORTRAN          , &
                                      MPI_RP, type_sub_array, err)
        call mpi_type_commit(type_sub_array, err)
        !
        ! === view array 
        !
        shape_view_array     = (/n/)
        shape_sub_view_array = SHAPE(vmean(sid:eid,1))
        start_view_coord     = (/sid-1/)

        call mpi_type_create_subarray(array_rank, shape_view_array           , &
                                      shape_sub_view_array, start_view_coord , &
                                      MPI_ORDER_FORTRAN, MPI_RP, &
                                      type_sub_view_array, err)
        call mpi_type_commit(type_sub_view_array, mpi_err)
        !
        ! === start writing to the file
        !
        CALL MPI_FILE_WRITE_ALL(fh, itStat , 1, intprec, status, err)

        initial_displacement = 1*ip
        CALL MPI_FILE_SET_VIEW(fh, initial_displacement, MPI_RP, &
             type_sub_view_array, "native", MPI_INFO_NULL, err)

        do l = 1,nstats
           CALL MPI_FILE_WRITE_ALL(fh,vmean(:,l),1,type_sub_array,status,err)
        enddo

        CALL MPI_FILE_CLOSE(fh, err)

        return
end subroutine mpi_write_stat1D





subroutine mpi_write_stat2D(sx,ex,sy,ey,nx,ny,vmean)

        implicit none
        integer                                , intent(in) :: sx,ex,sy,ey
        integer                                , intent(in) :: nx,ny
        real(rp), dimension(:,:,:), allocatable, intent(in) :: vmean

        integer, parameter :: array_rank=2
        integer, parameter :: intprec = MPI_INTEGER

        integer, dimension(array_rank)      :: shape_array
        integer, dimension(array_rank)      :: shape_sub_array
        integer, dimension(array_rank)      :: start_coord

        integer, dimension(array_rank)      :: shape_view_array
        integer, dimension(array_rank)      :: shape_sub_view_array
        integer, dimension(array_rank)      :: start_view_coord

        integer(kind = MPI_OFFSET_KIND)     :: initial_displacement
        integer, dimension(MPI_STATUS_SIZE) :: status

        character(dl) :: iteration, filename
        integer       :: fh, err = 0, l
        integer       :: type_sub_array, type_sub_view_array
        !
        ! === open the file
        !
        write(iteration, '(i7.7)') itStat
        filename = "DATA/"//trim(data_dir)//"/BINARY_STAT/stats"&
                //trim(iteration)//'.bin'

        CALL MPI_FILE_OPEN(mpi_comm_cart, filename, &
                 MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                 MPI_INFO_NULL, fh, err)

        if(err.ne.0) then
          if(rank == root) then
            write(*,'(A,A)') 'MPI_IO error in writing file ', trim(filename)
          endif
          stop
        endif
        !
        ! === subarrays 
        !
        shape_array     = SHAPE(vmean(:,:,1))
        shape_sub_array = SHAPE(vmean(sx:ex,sy:ey,1))
        start_coord     = (/GN,GN/)

        call mpi_type_create_subarray(array_rank, shape_array, shape_sub_array, &
                                      start_coord, MPI_ORDER_FORTRAN          , &
                                      MPI_RP, type_sub_array, err)
        call mpi_type_commit(type_sub_array, err)
        !
        ! === view array 
        !
        shape_view_array     = (/nx,ny/)
        shape_sub_view_array = SHAPE(vmean(sx:ex,sy:ey,1))
        start_view_coord     = (/sx-1,sy-1/)

        call mpi_type_create_subarray(array_rank, shape_view_array           , &
                                      shape_sub_view_array, start_view_coord , &
                                      MPI_ORDER_FORTRAN, MPI_RP, &
                                      type_sub_view_array, err)
        call mpi_type_commit(type_sub_view_array, err)
        !
        ! === start writing to the file
        !
        CALL MPI_FILE_WRITE_ALL(fh, itStat , 1, intprec, status, err)

        initial_displacement = 1*ip
        CALL MPI_FILE_SET_VIEW(fh, initial_displacement, MPI_RP, &
             type_sub_view_array, "native", MPI_INFO_NULL, err)

        do l = 1,nvAve2D
           CALL MPI_FILE_WRITE_ALL(fh,vmean(:,:,l),1,type_sub_array,status,err)
        enddo

        CALL MPI_FILE_CLOSE(fh, err)

        return
end subroutine mpi_write_stat2D










subroutine stat_turbulent_boundary_layer
        
        use statistics_module
        use FileModule
        use real_to_integer_module
        use interpolation_module
        use fluid_functions_module, only: laminar_viscosity, KarmanSchoenherr, ReThetaIncSpalding
        use onlineStats
        use wmles_module!, only: compute_tau_wall_wmles2D, compute_tau_wall_wmles2D_static

        implicit none

        ! local declaration
        type(FileType)                          :: statFile
        type(FileType)                          :: cf_File
        integer                                 :: i, s, is
        integer , parameter                     :: nSt = 20      !< number of stations
        integer , dimension(nSt)                :: StPid
        real(rp), dimension(nsT)                :: xTgT
       
        real(rp), allocatable, dimension(:,:) :: frictVEC
        real(rp), dimension(0:ny) :: yPlus, uPlus, uPlusVD!, uTemp

        real(rp), dimension(lby:uby) :: rFav, uFav, Tfav, mFav
        real(rp) :: errR, errU, errT, errM

        real(rp) :: rWall, uWall, tmpWl, mWall, u_yWl, tauWl
        real(rp) :: u_tau, ReyTu, idvu, Matau
        real(rp) :: r, du, uc, d99, dOmega
        real(rp) :: omega, omTh, omeg099, omeg100
        real(rp) :: r_, ir, uu, vv, ww, uv, mMol, mTrb, mRat
        
        real(rp) :: R11, R22, R33, R12
        real(rp) :: pRms, rRms
        real(rp) :: Machj, deltaMachRey, deltaMachFav, Mach099, Mach100
        integer  :: j99, jMa, jj99, m, jOme, l

        real(rp) :: deltaSt
        real(rp) :: dely, unum, uden, udel, ufav099, ufav100, TtotFavre, ek2
        real(rp) :: dstar, theta, dstarinc, thetainc, shapef, shapefinc, ReyTheta
        real(rp) :: r_E, u_E, T_E, m_E, u_, r_p, u_p, dy_j, tmpAD
        real(rp) :: aa, Fc, cfInc, ReyThetaInc
  
        integer , parameter :: nw = wmles_npts
        real(rp), dimension(1:nw) :: yc, u_wm, T_wm, u_wmVD
        real(rp) :: sensor_med, sensor_rms
        real(rp) :: wenofl_med, wenofl_rms
        real(rp) :: uep, ReThetaSpalding
        logical  :: wmles_stat

        type(FileType) :: wmFile


        if(cart_dims(2) > 1) &
        stop ' Turbulent boundary layer online statistics supports only 2Ddecomp!'


        if(wmles) call compute_xz_plane_stats(WMLES_DATA_LW,nz,nVWMLESData,itStat,&
                                              comm1Dz,vmean1D_wmles,mpi_flag)

        !
        ! compute friction reynolds number and friction coefficient
        !
        call AllocateReal(frictVEC,sx,ex,2,21)
        if(wmles) then
          call DeallocateReal(frictVEC)
          call AllocateReal(frictVec,sx,ex,2,21+nvWmlesData)
        endif

        do i = sx,ex
        
           ! compute wall properties
           r     = Prandtl**(1.0_rp/3.0_rp)
           tmpAD = 1.0_rp + r*0.5_rp*(gamma0-1.0_rp)*Mach**2
           tmpWl = Trat*tmpAD
           rWall = vmean2D(i,sy,11)/tmpWl
           mWall = laminar_viscosity(tmpWl,Tref,vis_flag)

           uWall = vmean2D(i,sy,2)/vmean2D(i,sy,1)
           u_yWl = (uWall)/y(sy)
           tauWl = mu_inf*mWall*u_yWl
           if    (bc(3) == 'dws_adiabatic') then
             call compute_tau_wall_wmles2D(i,vmean2D,tmpWl,u_wm,T_wm,tauWl)
           elseif(bc(3) == 'static_dws_adiabatic') then
             tauWl = vmean1D_wmles(i,3)
           elseif(bc(3) == 'istatic_dws_adiabatic') then
             tauWl = vmean1D_wmles(i,3)
           elseif(bc(3) == 'vorticity_dws_adiabatic') then
             call compute_tau_wall_wmles2D_vorticity(i,vmean2D,vmean2D_aux,tmpWl,u_wm,T_wm,tauWl)
           elseif(bc(3) == 'idws_adiabatic') then
             tauWl = vmean1D_wmles(i,3)
           elseif(bc(3) == 'idws_isothermal') then
             tauWl = vmean1D_wmles(i,3)
           endif
           u_Tau = sqrt(abs(tauWl)/rWall)

           udel = 0.99_rp*u_inf
           j99  = 1
           do j=1,ny-1
            uu = vmean2D(i,j,2)/vmean2D(i,j,1)
            if (uu>udel) then
             j99 = j-1
             exit
            endif
           enddo
           ufav099 = vmean2D(i,j99  ,2)/vmean2D(i,j99  ,1)
           ufav100 = vmean2D(i,j99+1,2)/vmean2D(i,j99+1,1)
           dely = y(j99+1)-y(j99)
           unum = udel-ufav099
           uden = ufav100-ufav099
           d99 = y(j99)+dely*(unum/uden)
        
           ! get vorticity BL thickness
           omTh = 0.05_rp
           jOme = 1
           do j=1,ny-1
            omega = vmean2D_aux(i,j,19) 
            if (omega<omTh) then
             jOme = j-1
             exit
            endif
           enddo
           omeg099 = vmean2D_aux(i,jOme  ,19)
           omeg100 = vmean2D_aux(i,jOme-1,19)
           dely = y(jOme)-y(jOme-1)
           unum = omTh-omeg099
           uden = omeg099 - omeg100
           dOmega = y(jOme)+abs(dely*(unum/uden))

           ! compute the local shear reynolds number
           ReyTu = rWall*u_tau*d99/(mu_inf*mWall)
           !
           ! === integral boundary layer thickness
           !
           rFav = vmean2D(i,:,1)
           uFav = vmean2D(i,:,2)/vmean2D(i,:,1)
           TFav = vmean2D(i,:,13)
           mFav = vmean2D(i,:,15) + vmean2D(i,:,16)

           m = 4
           jj99 = min(max(j99-(m-1)/2,1),ny+1-m)
           call polint(y(jj99),rFav(jj99),m,d99,r_E,errR)
           call polint(y(jj99),uFav(jj99),m,d99,u_E,errU)
           call polint(y(jj99),TFav(jj99),m,d99,T_E,errT)
           call polint(y(jj99),mFav(jj99),m,d99,m_E,errM)

           dstar    = 0.0_rp
           theta    = 0.0_rp
           dstarinc = 0.0_rp
           thetainc = 0.0_rp
           do j = 1,j99

              r_ =  vmean2D(i,j,1)/r_E
              u_ = (vmean2D(i,j,2)/r_)/u_E

              r_p =  vmean2D(i,j+1,1)/r_E
              u_p = (vmean2D(i,j+1,2)/r_p)/u_E

              dy_j = y(j+1) - y(j)

              dstar = dstar + 0.5_rp*dy_j*((1.0_rp-r_*u_) + (1.0_rp-r_p*u_p))
              theta = theta + 0.5_rp*dy_j*((r_*u_*(1.0_rp-u_) + r_p*u_p*(1.0_rp-u_p)))

              dstarinc = dstarinc + 0.5_rp*dy_j*((1.0_rp-u_) + (1.0_rp-u_p))
              thetainc = thetainc + 0.5_rp*dy_j*((u_*(1.0_rp-u_) + u_p*(1.0_rp-u_p)))

           enddo
           shapef    = dstar/theta
           shapefinc = dstarinc/thetainc
           ReyTheta  = r_E*theta*u_E/(mu_inf*m_E)
        
           !
           ! === VAN DRIEST II
           !
           aa = (tmpWl/T_E-1.0_rp)/sqrt(tmpWl/T_E*(tmpWl/T_E-1))
           Fc = (tmpWl/T_E-1.0_rp)/(asin(aa)**2)

           cfInc       = tauWl/q_inf*Fc
           ReyThetaInc = m_E/mWall*ReyTheta

           ! compute Mach Line Reynolds
           jMa = 1
           do j = 1,ny-1
              machj = vmean2D(i,j,18)
              if(machj > 1.0_rp) then
                jMa = j-1
                exit
              endif
           enddo
           Mach099 = vmean2D(i,jMa  ,18)
           Mach100 = vmean2D(i,jMa+1,18)
           dely = y(jMa+1)-y(jMa)
           unum = 1.0_rp-Mach099
           uden = Mach100-Mach099
           deltaMachRey = y(jMa)+dely*(unum/uden)

           ! compute Mach Line favre
           jMa = 1
           do j = 1,ny-1
              machj = vmean2D(i,j,19)/vmean2D(i,j,1)
              if(machj > 1.0_rp) then
                jMa = j-1
                exit
              endif
           enddo
           Mach099 = vmean2D(i,jMa  ,19)/vmean2D(i,jMa+1,1)
           Mach100 = vmean2D(i,jMa+1,19)/vmean2D(i,jMa+1,1)
           dely = y(jMa+1)-y(jMa)
           unum = 1.0_rp-Mach099
           uden = Mach100-Mach099
           deltaMachFav = y(jMa)+dely*(unum/uden)

           ! compute spalding and chi formula
           uep = sqrt(2.0_rp/(tauWl/q_inf))
           ReThetaSpalding = ReThetaIncSpalding(uep)

           FrictVEC(i,2) = tauWl/q_inf
           FrictVEC(i,3) = ReyTu
           FrictVEC(i,4) = shapef
           FrictVEC(i,5) = shapefinc
           FrictVEC(i,6) = d99
           FrictVEC(i,7) = dstar
           FrictVEC(i,8) = theta
           FrictVEC(i,9) = u_tau
           FrictVEC(i,10)= ReyTheta
           FrictVEC(i,11)= ReyThetaInc
           FrictVEC(i,12)= cfInc
           !FrictVEC(i,13)= 0.024_rp*ReyThetaInc**(-0.25_rp)
           !Karman-Schoenherr
           FrictVEC(i,13) = KarmanSchoenherr(ReyThetaInc)
           FrictVEC(i,14) = deltaMachRey   ! sonic line Reynolds
           FRICTVEC(i,15) = vmean2D(i,sy,11) ! pwall
           FrictVEC(i,16) = y(wmles_strI)/d99
           FrictVEC(i,17) = rWall
           FrictVEC(i,18) = vmean2D(i,sy,11)
           FrictVEC(i,19) = dOmega
           FrictVEC(i,20) = ReThetaSpalding
           FrictVEC(i,21) = KarmanSchoenherr(ReThetaSpalding)

           if(wmles) then
             do l = 1,nvWmlesData
                FrictVEC(i,21+l) = vmean1D_wmles(i,l)
             enddo
           endif

        enddo
        
        ! plot cf coefficient
        cf_File%name = 'friction'
        cf_File%dir  = trim(data_dir)//'/FRICTION'
        call mpi_plot_vec(x,frictVEC,nx,sx,ex,rank,&
             mpi_comm_world,nprocs,root,cf_File,it)

        ! determine the ranks which own the control stations
        deltaSt = Lx/real(nSt,rp)
        xTgT(1) = 10.0_rp
        do is = 1,nSt
           xTgT(is) = is*deltaSt
        enddo

        StPid(:) = MPI_PROC_NULL
        do is = 1, nSt
           if(xTgT(is) > x(sx) .and. xTgT(is) < x(ex) .and. & 
             abs(0.5_rp*(z(sz) + z(sz-1)) - zmin) < 1.0E-14_rp) then
             StPid(is) = rank
           endif
        enddo

        !uTemp(:) = 0.0_rp

        StLoop:do is = 1, nSt
           if(rank == StPid(is)) then

           ! open files
           statFile%name = 'stat'
           statFile%dir  = trim(data_dir)//'/STATS_'//trim(str(int(xTgT(is))))
           call OpenNewFile(statFile,it)

           ! find integer of the StLoc station
           s = nearest_integer_opt(x,sx,ex,xTgt(is))
        
           ! compute wall properties
           r     = Prandtl**(1.0_rp/3.0_rp)
           tmpWl = Trat*(1.0_rp + r*0.5_rp*(gamma0-1.0_rp)*Mach**2)
           rWall = vmean2D(s,sy,11)/tmpWl
           mWall = laminar_viscosity(tmpWl,Tref,vis_flag)

           uWall = vmean2D(s,sy,2)/vmean2D(s,sy,1)
           u_yWl = (uWall)/y(sy)
           tauWl = mu_inf*mWall*u_yWl

           wmles_stat=.false.
           if(bc(3) == 'dws_adiabatic') then
             call compute_tau_wall_wmles2D(s,vmean2D,tmpWl,u_wm,T_wm,tauWl)
             wmles_stat=.true.
           elseif(bc(3) == 'static_dws_adiabatic') then
             tauWl = vmean1D_wmles(i,3)
             wmles_stat=.true.
           elseif(bc(3) == 'istatic_dws_adiabatic') then
             tauWl = vmean1D_wmles(i,3)
             wmles_stat=.true.
           elseif(bc(3) == 'vorticity_dws_adiabatic') then
             call compute_tau_wall_wmles2D_vorticity(s,vmean2D,vmean2D_aux,tmpWl,u_wm,T_wm,tauWl)
             wmles_stat=.true.
           elseif(bc(3) == 'idws_adiabatic') then
             tauWl = vmean1D_wmles(i,3)
             wmles_stat=.true.
           elseif(bc(3) == 'idws_isothermal') then
             tauWl = vmean1D_wmles(i,3)
             wmles_stat=.true.
           endif
        
           if(wmles_stat) then
             wmFile%dir = trim(data_dir)//'/WMLES_STATS_'//trim(str(int(xTgT(is))))             
             wmFile%name= 'wmles_stat'
             call OpenNewFile(wmFile,it)
             u_Tau = sqrt(tauWl/rWall)
        
             j = 1
             du = u_wm(j)
             r_ = vmean2D(s,wmles_strI,11)/T_wm(j)
             uc = sqrt(r_/rWall)
             u_wmVD(j) = uc*du
             do j = 2,nw
                du = u_wm(j) - u_wm(j-1)
                r_ = vmean2D(s,wmles_strI,11)/T_wm(j)
                uc = sqrt(r_/rWall)
                u_wmVD(j) = u_wmVD(j-1) + uc*du
             enddo
        
             do j = 1,nw
                write(wmFile%unit,'(20e18.9)')      &
                  rWall*yc(j)*u_tau/(mu_inf*mWall), &
                  u_wm(j)/u_tau                   , &
                  u_wmVD(j)/u_tau                 , &
                  T_wm(j)/tmpWl                   , &
                  yc(j)                           , &                            
                  u_wm(j)/u_inf                  
             enddo
             call CloseFile(wmFile)
           endif

           u_Tau = sqrt(tauWl/rWall)
           MaTau = u_tau/sqrt(gamma0*tmpWl)
        
           ! compute the local shear reynolds number
           idvu = rWall*u_tau/(mu_inf*mWall)

           ! compute uPlus and yPlus
           uPlus = 0.0_rp
           yPlus = 0.0_rp
           do j = sy,ey
              yPlus(j) = y(j)*idvu
              uPlus(j) = (vmean2D(s,j,2)/vmean2D(s,j,1))/u_tau
           enddo
        
           ! trasform with VAN DRIEST
           uPlusVD = 0.0_rp
           do j = sy,ey
              du = uPlus(j) - uPlus(j-1)
              uc = sqrt(vmean2D(s,j,1)/rWall)
              uPlusVD(j) = uPlusVD(j-1) + uc*du
           enddo

           write(statFile%unit,'(A)')       '# *************************************************'
           write(statFile%unit,'(A,f18.6)') '#  TURBULENT BOUNDARY LAYER ', Mach 
           write(statFile%unit,'(A)')       '#'
           write(statFile%unit,'(A)')       '#  author: Francesco De Vanna '
           write(statFile%unit,'(A)')       '#  e-mail: fra.devanna@gmail.com'
           write(statFile%unit,'(A)')       '#'
           write(statFile%unit,'(A,f18.6)') '#  LOCATION%                 ', x(s)/Lx*100
           write(statFile%unit,'(A,f18.6)') '#  FRICTION COEFFICIENT      ', frictVEC(s,2)
           write(statFile%unit,'(A,f18.6)') '#  FRICTION REYNOLDS NUMBER  ', frictVEC(s,3)
           write(statFile%unit,'(A,f18.6)') '#  THETA REYNOLDS NUMBER     ', frictVEC(s,10)
           write(statFile%unit,'(A,f18.6)') '#  SHAPE FACTOR              ', frictVEC(s,4)
           write(statFile%unit,'(A)')       '# *************************************************'
           write(statFile%unit,'(A)')       '# Column 1  : y '
           write(statFile%unit,'(A)')       '# Column 2  : y+'
           write(statFile%unit,'(A)')       '# Column 3  : u+'
           write(statFile%unit,'(A)')       '# Column 4  : u_VD+'
           write(statFile%unit,'(A)')       '# Column 5  : urms+'
           write(statFile%unit,'(A)')       '# Column 6  : vrms+'
           write(statFile%unit,'(A)')       '# Column 7  : wrms+'
           write(statFile%unit,'(A)')       '# Column 8  : uv+'
           write(statFile%unit,'(A)')       '# Column 9  : sqrt(rho/rhoWall)'
           write(statFile%unit,'(A)')       '# Column 10 : temp/tWall'
           write(statFile%unit,'(A)')       '# Column 11 : p/tauW'
           write(statFile%unit,'(A)')       '# Column 12 : prsrms+'
           write(statFile%unit,'(A)')       '# Column 13 : rhorms+'
           write(statFile%unit,'(A)')       '# Column 14 : MuMol/MuMolWall'
           write(statFile%unit,'(A)')       '# Column 14 : MuTrb/MuMolWall'
           write(statFile%unit,'(A)')       '# Column 16 : MuTrb/MuMol'
           write(statFile%unit,'(A)')       '# Column 17 : Ttot/Twall'
           write(statFile%unit,'(A)')       '# *************************************************'
        
           do j = sy,ey
                
              r_ = vmean2D(s,j,1)
              ir = 1.0_rp/r_
                
              uu = vmean2D_aux(s,j,3)
              vv = vmean2D_aux(s,j,4)
              ww = vmean2D_aux(s,j,5)
              uv = vmean2D_aux(s,j,6)

              R11 = (r_/rWall)*uu/u_tau**2
              R22 = (r_/rWall)*vv/u_tau**2
              R33 = (r_/rWall)*ww/u_tau**2
              R12 = (r_/rWall)*uv/u_tau**2

              ! prsRms, rhoRms, TmpRms
              rRms = sqrt(vmean2D(s,j,10)    - vmean2D(s,j,1 )**2)
              pRms = sqrt(vmean2D(s,j,12)    - vmean2D(s,j,11)**2)

              pRms = pRms/tauWl
              rRms = rRms/(gamma0*rWall*Matau**2)

              ! molecular viscosity, turbulent viscosity and viscosity ratio
              mMol = vmean2D(s,j,15)
              mTrb = vmean2D(s,j,16)
              mRat = mTrb/mMol

              ek2 = 0.5_rp*(vmean2D(s,j,6) + vmean2D(s,j,7) + vmean2D(s,j,8))
              TtotFavre = ir*(vmean2D(s,j,17) + (gamma0-1.0_rp)/gamma0 * ek2)
        
              ! WENO sensor
              sensor_med = vmean2D(s,j,32)
              sensor_rms = sqrt(vmean2D(s,j,33) - sensor_med**2)

              wenofl_med = vmean2D(s,j,34)
              wenofl_rms = sqrt(vmean2D(s,j,35) - wenofl_med**2)

              write(statFile%unit,'(100e18.9)') &
                      y(j)                   , &  !1
                      yPlus(j)               , &  !2
                      uPlus(j)               , &  !3
                      uPlusVD(j)             , &  !4
                      sqrt(R11)              , &  !5
                      sqrt(R22)              , &  !6
                      sqrt(R33)              , &  !7
                      R12                    , &  !8
                      sqrt(r_/rWall)         , &  !9
                      vmean2D(s,j,13)/tmpWl  , &  !10
                      vmean2D(s,j,11)/tauWl  , &  !11
                      pRms                   , &  !12
                      rRms                   , &  !13
                      mMol/mWall             , &  !14
                      mTrb/mWall             , &  !15
                      mRat                   , &  !16
                      TtotFavre/tmpWl        , &  !17
                      sensor_med             , &  !18
                      sensor_rms             , &  !19
                      wenofl_med             , &  !20
                      wenofl_rms             , &  !21
                      uPlus(j)*u_tau/u_inf   , &  !22
                      R11*rWall*u_tau**2/u_inf**2  , &  !23
                      R22*rWall*u_tau**2/u_inf**2  , &  !23
                      R33*rWall*u_tau**2/u_inf**2  , &  !23
                      R12*rWall*u_tau**2/u_inf**2  
                      
           enddo

           call CloseFile(statFile)

           endif
        enddo StLoop
        
        deallocate(frictVEC)
        
        return
end subroutine stat_turbulent_boundary_layer






subroutine WriteWallPressure

        implicit none
        real(rp), dimension(:,:), allocatable :: wallPressure
        integer                               :: i,k, funit = 99
        character(dl)                         :: fdir, fnom, ff
        character(20)                         :: chx, chz
        logical                               :: dir_exist

        if(my_neighbour(3) == MPI_PROC_NULL) then
        
          !
          ! === move pressure to wallPressure array
          !
          call AllocateReal(wallPressure,sx,ex,sz,ez)
          do    k = sz,ez
             do i = sx,ex
                wallPressure(i,k) = P(i,1,k)
             enddo
          enddo
          !
          ! === set the directory
          !
          fdir = 'DATA/'//trim(data_dir)//'/WALL_PRESSURE/'
          inquire(file = trim(fdir), exist= dir_exist)
          if(.not. dir_exist) call execute_command_line('mkdir -p '//trim(fdir))
          !
          ! === set the file
          ! 
          write(chx,'(I3.3)') cart_coord(1)
          write(chz,'(I3.3)') cart_coord(3)

          fnom = 'wallpressure_'//trim(chx)//'_'//trim(chz)//'.bin'
          ff   = trim(fdir)//trim(fnom)
          open(unit = funit, file=ff, form='unformatted', position='append')
          !
          ! === write down
          !
          write(funit) it, time, ex-sx+1, ez-sz+1
          write(funit) wallPressure
          !
          ! === close
          !
          close(funit)
          call DeallocateReal(wallPressure)
        
        endif

        return
end subroutine WriteWallPressure




















end module output_module
