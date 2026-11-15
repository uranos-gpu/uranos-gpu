module inflow_module
use parameters_module
use storage_module
use profiling_module

implicit none
private

integer, parameter :: DF_STATE_VERSION = 1

real(rp), allocatable, dimension(:,:)  , public :: iflow_mean
real(rp), allocatable, dimension(:,:,:), public :: iflow_turb
real(rp), allocatable, dimension(:)    , public :: iflow_dTdu
real(rp), allocatable, dimension(:)    , public :: iflow_urms
real(rp), allocatable, dimension(:,:,:), public :: iflow_vf_old
real(rp), allocatable, dimension(:,:,:), public :: iflow_vf_new

public :: init_inflow_profile
public :: inflow_turbulent_update
public :: inflow_turbulent_dump_state
public :: inflow_turbulent_load_state

private :: inflow_TurbulentBoundaryLayer
private :: inflow_turbulent_init
private :: dump_inflow_mean_profile

contains

subroutine init_inflow_profile
        implicit none

        if(trim(bc(1)) /= 'nscbc_inflow_generalized') return

        call AllocateReal(iflow_mean, sy,ey, 1,5)
        call AllocateReal(iflow_turb, 1,3,    sy,ey, sz,ez)
        call AllocateReal(iflow_dTdu, sy,ey)
        call AllocateReal(iflow_urms, sy,ey)

        selectcase(trim(inflow_profile))

          case('TurbulentBoundaryLayerInflow')
             call inflow_TurbulentBoundaryLayer()

          case default
             if(rank == root) write(*,'(A)') &
                " init_inflow_profile: inflow_profile '"//trim(inflow_profile)//"' is not implemented"
             call secure_stop

        end select

        call dump_inflow_mean_profile()

        if(turb_inflow) call inflow_turbulent_init

        return
end subroutine init_inflow_profile

subroutine inflow_TurbulentBoundaryLayer()
        use fluid_functions_module, only: MuskerProfile, CompressibleCorrection, GRATemperatureSlope
        use mesh_module           , only: west_wall_normal_distance, nx_xi, ny_xi

        implicit none

        real(rp), dimension(0:ny)           :: UPlusInc, UPlusCmp, RhoY, TmpY, velY
        real(rp), dimension(0:ny)           :: wdist, vbl
        real(rp), dimension(:), allocatable :: y_tmp
        real(rp) :: tWall, rWall, ue, cf, theta, dddx, rui
        integer  :: j

        call AllocateReal(y_tmp, 0, ny)
        call west_wall_normal_distance(wdist)
        y_tmp(0:ny) = ymin + wdist(0:ny)

        call MuskerProfile         (ny, y_tmp, ReTau, UPlusInc)
        call CompressibleCorrection(ny, UPlusInc, UPlusCmp, RhoY, TmpY, tWall, rWall)

        ue   = maxval(UPlusCmp)
        velY = UPlusCmp/ue

        cf    = 2.0_rp*rWall/(RhoY(ny)*ue*ue)
        theta = 0.0_rp
        do j = 1,ny
           theta = theta + 0.5_rp*(wdist(j) - wdist(j-1))                    &
                 * ( RhoY(j  )*velY(j  )*(1.0_rp - velY(j  ))               &
                   + RhoY(j-1)*velY(j-1)*(1.0_rp - velY(j-1)) )
        enddo
        theta = theta/RhoY(ny)
        dddx  = 0.5_rp*cf/theta

        rui    = 0.0_rp
        vbl(0) = 0.0_rp
        do j = 1,ny
           rui    = rui + 0.5_rp*(wdist(j) - wdist(j-1))                     &
                  * (RhoY(j)*velY(j) + RhoY(j-1)*velY(j-1))
           vbl(j) = dddx*u_inf*(wdist(j)*RhoY(j)*velY(j) - rui)/RhoY(j)
        enddo

        do j = sy, ey
           iflow_mean(j,1) = RhoY(j)
           iflow_mean(j,2) = u_inf*velY(j)*nx_xi(sx,j) - vbl(j)*ny_xi(sx,j)
           iflow_mean(j,3) = u_inf*velY(j)*ny_xi(sx,j) + vbl(j)*nx_xi(sx,j)
           iflow_mean(j,4) = 0.0_rp
           iflow_mean(j,5) = TmpY(j)
           iflow_dTdu(j)   = GRATemperatureSlope(velY(j), TmpY(j))/u_inf
        enddo

        call DeallocateReal(y_tmp)

        return
end subroutine inflow_TurbulentBoundaryLayer

subroutine inflow_turbulent_init

        use df_module
        use fileModule
        use mesh_module, only: west_wall_normal_distance

        implicit none
        real(rp), dimension(:), allocatable :: y_tmp
        real(rp), dimension(0:ny)           :: wdist
        integer                             :: j
        integer                             :: err = 0
#ifdef DEBUG
        type(FileType)         :: reyStressFile
        real(rp), dimension(3) :: rms, rmean, rmeansq
        real(rp)               :: rey12, rey13, rey23
        integer                :: m
#endif

        allocate(y_tmp(0:ny), stat = err)
        if(err .ne. 0) stop ' Allocation error in inflow_turbulent_init'

        allocate(DF_ylen(3,0:ny), DF_zlen(3,0:ny), DF_uconv(0:ny), stat = err)
        if(err .ne. 0) stop ' Allocation error in inflow_turbulent_init'
        DF_uconv = 1.0_rp

        allocate(DF_By(3,-DF_N:DF_N, 0:ny), DF_Bz(3,-DF_N:DF_N, 0:ny), stat = err)
        if(err .ne. 0) stop ' Allocation error in inflow_turbulent_init'

        allocate(DF_Rnd2D(3, 1-DF_N:ny+DF_N, 1-DF_N:nz+DF_N), stat = err)
        if(err .ne. 0) stop ' Allocation error in inflow_turbulent_init'

        allocate(iflow_vf_old(3, sy:ey, sz:ez), stat = err)
        if(err .ne. 0) stop ' Allocation error in inflow_turbulent_init'

        allocate(iflow_vf_new(3, sy:ey, sz:ez), stat = err)
        if(err .ne. 0) stop ' Allocation error in inflow_turbulent_init'

        call west_wall_normal_distance(wdist)
        y_tmp(0:ny) = ymin + wdist(0:ny)

        selectcase(inflow_profile)

          case('TurbulentBoundaryLayerInflow')
          call DFInputData_TBL(ReTau,Mach,Prandtl,y_tmp,ny)

          call DFIntegralLenght_TBL(y_tmp,ny)

        case default
          print*, 'Turbulent inflow is not implemented for inflow "', trim(inflow_profile), '"'
          stop
        endselect

        call DFCoefficients(y_tmp,ny,gbl_min_step)

        do j = sy,ey
           iflow_urms(j) = DF_LundMatrix(1,1,j)
        enddo
        !$acc data copy(DF_Rnd2D,DF_ylen,DF_zlen,DF_uconv,DF_By,DF_Bz,DF_Fy,DF_LundMatrix) &
        !$acc copy(iflow_vf_old, iflow_vf_new, iflow_turb)

        call DFRandomField2D(ny,nz)

        call DFConvolution2D(sy,ey,sz,ez,iflow_vf_old)

        call DFRandomField2D(ny,nz)

        call DFConvolution2D(sy,ey,sz,ez,iflow_vf_new)

        call DFCastroTimeCorrelation(ik,c_rk,u_inf,dt,sy,ey,sz,ez,iflow_vf_old,iflow_vf_new)

        call DFEnforceReynoldsStresses2D(sy,ey,sz,ez,iflow_vf_new,iflow_turb)

        !$acc end data

#ifdef DEBUG
        if(nprocs == 1) then
        reyStressFile%name = 'rey_computed_inflow'
        reyStressFile%dir  = trim(data_dir)//'/INFLOW_REYNOLDS_STRESS'
        call OpenNewFile(ReyStressFile,it)
        do j=1,ny
         rmean   = 0._rp
         rmeansq = 0._rp
         rey12   = 0._rp
         rey13   = 0._rp
         rey23   = 0._rp
         do k=1,nz
           do m=1,3
            rmean(m)   = rmean(m)   + iflow_turb(m,j,k)
            rmeansq(m) = rmeansq(m) + iflow_turb(m,j,k)**2
           enddo
           rey12 = rey12+iflow_turb(1,j,k)*iflow_turb(2,j,k)
           rey13 = rey13+iflow_turb(1,j,k)*iflow_turb(3,j,k)
           rey23 = rey23+iflow_turb(2,j,k)*iflow_turb(3,j,k)
         enddo
         rmean   = rmean/(nz)
         rmeansq = rmeansq/(nz)
         rms     = sqrt(rmeansq-rmean*rmean)
         rey12   = rey12/(nz)
         rey13   = rey13/(nz)
         rey23   = rey23/(nz)
         rey12   = rey12-rmean(1)*rmean(2)
         rey13   = rey13-rmean(1)*rmean(3)
         rey23   = rey23-rmean(2)*rmean(3)
         write(reyStressFile%unit,*) &
         y(j),(rms(m)**2,m=1,3),rey12,rey13,rey23,(rmean(m),m=1,3)
        enddo
        call CloseFile(ReyStressFile)
        endif
#endif
        deallocate(y_tmp)

        return
end subroutine inflow_turbulent_init

subroutine inflow_turbulent_update()

        use df_module
        use FileModule
        use mesh_module, only: nx_xi, ny_xi

        implicit none
        real(rp) :: fs, fn
        integer  :: j,k,l

#ifdef DEBUG
        type(fileType)         :: turbInflow
        real(rp)               :: lcl_ekt, gbl_ekt
        integer                :: lcl_pts, gbl_pts
        real(rp)               :: lcl_ekt_max, gbl_ekt_max
        real(rp)               :: ek
        integer                :: err = 0
#endif

        call StartProfRange('inflow_turbulent_update')

        if(turb_inflow) then

          call DFRandomField2D(ny,nz)

          call DFConvolution2D(sy,ey,sz,ez,iflow_vf_new)

          call DFCastroTimeCorrelation(ik,c_rk,u_inf,dt,sy,ey,sz,ez,iflow_vf_old,iflow_vf_new)

          call DFEnforceReynoldsStresses2D(sy,ey,sz,ez,iflow_vf_new,iflow_turb)

          !$acc parallel default(present)
          !$acc loop gang, vector collapse(2)
          do    k = sz,ez
             do j = sy,ey
                fs = iflow_turb(1,j,k)
                fn = iflow_turb(2,j,k)
                iflow_turb(1,j,k) = nx_xi(sx,j)*fs - ny_xi(sx,j)*fn
                iflow_turb(2,j,k) = ny_xi(sx,j)*fs + nx_xi(sx,j)*fn
             enddo
          enddo
          !$acc end parallel

        else

          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = sz,ez
             do    j = sy,ey
                do l = 1,3
                   iflow_turb(l,j,k) = 0.0_rp
                enddo
             enddo
          enddo
          !$acc end parallel

        endif

#ifdef DEBUG
        turbInflow%name = 'inflowFluctuation'//trim(str(rank))
        turbInflow%dir  = trim(data_dir)//'/TURB_INFLOW/'
        call OpenNewFile(turbInflow,it)

        lcl_ekt     = 0.0_rp
        lcl_pts     = 0
        lcl_ekt_max = 0.0_rp
        do k    = sz,ez
           do j = sy,ey
              ek = (u_inf*iflow_turb(1,j,k))**2 + (u_inf*iflow_turb(2,j,k))**2 + (u_inf*iflow_turb(3,j,k))**2

              lcl_ekt     = lcl_ekt + ek
              lcl_ekt_max = max(ek,lcl_ekt_max)
              lcl_pts     = lcl_pts + 1
           enddo
        enddo
        if(mpi_flag) then
          call MPI_allreduce(lcl_ekt    ,gbl_ekt    ,1,MPI_RP, MPI_SUM, mpi_comm_cart, err)
          call MPI_allreduce(lcl_pts    ,gbl_pts    ,1,MPI_INTEGER         , MPI_SUM, mpi_comm_cart, err)
          call MPI_allreduce(lcl_ekt_max,gbl_ekt_max,1,MPI_RP, MPI_MAX, mpi_comm_cart, err)
        else
          gbl_ekt     = lcl_ekt
          gbl_pts     = lcl_pts
          gbl_ekt_max = lcl_ekt_max
        endif

        write(turbInflow%unit,*) '# turbulent mach nmber ', sqrt(gbl_ekt/real(gbl_pts,rp))
        write(turbInflow%unit,*) '# max turbulent energy ', gbl_ekt_max
        do k    = sz,ez
           do j = sy,ey
              write(turbInflow%unit,*) y(j), z(k), u_inf*iflow_turb(:,j,k)
           enddo
        enddo
        call CloseFile(turbInflow)
#endif

        call EndProfRange

        return
end subroutine inflow_turbulent_update

subroutine dump_inflow_mean_profile()
        implicit none

        character(len=512) :: fdir, fpath
        integer            :: j, iunit, ios
        logical            :: dir_exist

        if(rank /= root) return

        fdir = 'DATA/'//trim(data_dir)//'/INITIAL_CONDITION'
        inquire(file = trim(fdir), exist = dir_exist)
        if(.not. dir_exist) call execute_command_line('mkdir -p '//trim(fdir))

        fpath = trim(fdir)//'/inflow_mean_profile.dat'
        open(newunit = iunit, file = trim(fpath), status = 'replace', &
             action = 'write', iostat = ios)
        if(ios /= 0) then
           write(*,'(A)') ' dump_inflow_mean_profile: cannot open '//trim(fpath)
           return
        endif

        write(iunit, '(A)') '# inflow mean profile ('//trim(inflow_profile)//')'
        write(iunit, '(A)') '# y  rho  u/u_inf  v/u_inf  w/u_inf  T'
        do j = sy, ey
           write(iunit, *) y(j), iflow_mean(j,1),        &
                                 iflow_mean(j,2)/u_inf,  &
                                 iflow_mean(j,3)/u_inf,  &
                                 iflow_mean(j,4)/u_inf,  &
                                 iflow_mean(j,5)
        enddo
        close(iunit)

        return
end subroutine dump_inflow_mean_profile

subroutine df_state_filename(iter, fname)
        implicit none
        integer      , intent(in)  :: iter
        character(dl), intent(out) :: fname

        character(dl) :: iteration

        write(iteration,'(i0.7)') iter
        fname = 'DATA/'//trim(data_dir)//'/BINARY_DF/df'//trim(iteration)//'.bin'

        return
end subroutine df_state_filename

subroutine inflow_turbulent_dump_state
        use mpi
        use mpi_module, only: cart_coord, comm2Dyz, check_mpi_vec
        use df_module , only: DF_ndraw, DF_N

        implicit none

        character(dl) :: fname
        logical       :: is_writer
        integer       :: fh, int_size, int8_size, rp_size, header_count, ne
        integer       :: err, open_err_max, alr_err, hdr(6)
        integer       :: type_view
        integer(kind=MPI_OFFSET_KIND) :: disp
        integer, parameter :: array_rank = 3
        integer, dimension(array_rank) :: shp_view, shp_sub, start_view

        integer, parameter                :: NIO = 10
        integer, dimension(NIO)           :: errs
        character(len=28), dimension(NIO) :: lbls

        call StartProfRange('inflow_turbulent_dump_state')

        call df_state_filename(it, fname)

        is_writer = (cart_coord(1) == 0)
        errs = 0
        ne   = 0

        if(is_writer) then

        call MPI_FILE_OPEN(comm2Dyz, fname, MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                           MPI_INFO_NULL, fh, err)

        call MPI_Allreduce(err, open_err_max, 1, MPI_INTEGER, MPI_MAX, comm2Dyz, alr_err)
        if(open_err_max /= 0) then
          if(rank == root) write(*,'(1x,A,A)') ' WARNING: cannot write DF state ', trim(fname)
        else

        ne=ne+1; lbls(ne)='file_set_size'
        call MPI_FILE_SET_SIZE(fh, 0_MPI_OFFSET_KIND, errs(ne))

        ne=ne+1; lbls(ne)='type_size MPI_INTEGER'
        call MPI_TYPE_SIZE(MPI_INTEGER , int_size , errs(ne))
        ne=ne+1; lbls(ne)='type_size MPI_INTEGER8'
        call MPI_TYPE_SIZE(MPI_INTEGER8, int8_size, errs(ne))
        ne=ne+1; lbls(ne)='type_size MPI_RP'
        call MPI_TYPE_SIZE(MPI_RP      , rp_size  , errs(ne))

        hdr = (/DF_STATE_VERSION, ny, nz, DF_N, rp_size, it/)

        header_count = 0
        if(rank == root) header_count = 1

        ne=ne+1; lbls(ne)='write_at_all header'
        call MPI_FILE_WRITE_AT_ALL(fh, 0_MPI_OFFSET_KIND, hdr, 6*header_count, &
                                   MPI_INTEGER, MPI_STATUS_IGNORE, errs(ne))
        ne=ne+1; lbls(ne)='write_at_all ndraw'
        disp = 6*int_size
        call MPI_FILE_WRITE_AT_ALL(fh, disp, DF_ndraw, header_count, &
                                   MPI_INTEGER8, MPI_STATUS_IGNORE, errs(ne))

        shp_view   = (/3, ny, nz/)
        shp_sub    = (/3, ey-sy+1, ez-sz+1/)
        start_view = (/0, sy-1, sz-1/)

        ne=ne+1; lbls(ne)='type_create_subarray view'
        call MPI_TYPE_CREATE_SUBARRAY(array_rank, shp_view, shp_sub, start_view, &
                                      MPI_ORDER_FORTRAN, MPI_RP, type_view, errs(ne))
        call MPI_TYPE_COMMIT(type_view, err)

        disp = 6*int_size + int8_size
        ne=ne+1; lbls(ne)='file_set_view'
        call MPI_FILE_SET_VIEW(fh, disp, MPI_RP, type_view, "native", MPI_INFO_NULL, errs(ne))

        ne=ne+1; lbls(ne)='write_all plane'
        call MPI_FILE_WRITE_ALL(fh, iflow_vf_old, product(shp_sub), MPI_RP, &
                                MPI_STATUS_IGNORE, errs(ne))

        ne=ne+1; lbls(ne)='file_close'
        call MPI_FILE_CLOSE(fh, errs(ne))
        call MPI_TYPE_FREE(type_view, err)

        endif
        endif

        call check_mpi_vec('inflow_turbulent_dump_state', errs(1:ne), lbls(1:ne))

        call EndProfRange

        return
end subroutine inflow_turbulent_dump_state

subroutine inflow_turbulent_load_state
        use mpi
        use mpi_module, only: cart_coord, comm2Dyz, mpi_comm_cart, check_mpi_vec
        use df_module , only: DF_ndraw, DF_N

        implicit none

        character(dl) :: fname
        logical       :: is_writer, have_file
        integer       :: fh, int_size, int8_size, rp_size, ne
        integer       :: err, ok, ok_gbl, hdr(6)
        integer       :: type_view
        integer(kind=MPI_OFFSET_KIND) :: disp
        integer, parameter :: array_rank = 3
        integer, dimension(array_rank) :: shp_view, shp_sub, start_view

        integer, parameter                :: NIO = 8
        integer, dimension(NIO)           :: errs
        character(len=28), dimension(NIO) :: lbls

        if(len_trim(restart_file) == 0) return

        call df_state_filename(it, fname)

        have_file = .false.
        if(rank == root) inquire(file = trim(fname), exist = have_file)
        call MPI_BCAST(have_file, 1, MPI_LOGICAL, root, mpi_comm_cart, err)

        if(.not. have_file) then
          if(rank == root) then
            write(*,'(1x,A)')   ' WARNING: no digital-filter state next to this restart.'
            write(*,'(1x,A,A)') '          expected ', trim(fname)
            write(*,'(1x,A)')   '          the synthetic inflow restarts from a fresh slice.'
          endif
          return
        endif

        is_writer = (cart_coord(1) == 0)
        errs = 0
        ne   = 0
        ok   = 1

        call MPI_TYPE_SIZE(MPI_INTEGER , int_size , err)
        call MPI_TYPE_SIZE(MPI_INTEGER8, int8_size, err)
        call MPI_TYPE_SIZE(MPI_RP      , rp_size  , err)

        if(is_writer) then
          call MPI_FILE_OPEN(comm2Dyz, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, err)
          if(err /= 0) ok = 0
        endif

        hdr = 0
        if(rank == root .and. ok == 1) then
          call MPI_FILE_READ_AT(fh, 0_MPI_OFFSET_KIND, hdr, 6, MPI_INTEGER, MPI_STATUS_IGNORE, err)
          if(err /= 0) ok = 0
          disp = 6*int_size
          call MPI_FILE_READ_AT(fh, disp, DF_ndraw, 1, MPI_INTEGER8, MPI_STATUS_IGNORE, err)
          if(err /= 0) ok = 0
        endif
        call MPI_BCAST(hdr, 6, MPI_INTEGER, root, mpi_comm_cart, err)
        call MPI_ALLREDUCE(ok, ok_gbl, 1, MPI_INTEGER, MPI_MIN, mpi_comm_cart, err)
        ok = ok_gbl

        if(ok == 1) then
          if(hdr(1) /= DF_STATE_VERSION .or. hdr(2) /= ny .or. hdr(3) /= nz .or. &
             hdr(4) /= DF_N             .or. hdr(5) /= rp_size .or. hdr(6) /= it) ok = 0
        endif

        if(ok == 0) then
          if(rank == root) then
            write(*,'(1x,A,A)') ' WARNING: unusable digital-filter state ', trim(fname)
            write(*,'(1x,A)')   '          the synthetic inflow restarts from a fresh slice.'
          endif
          if(is_writer) call MPI_FILE_CLOSE(fh, err)
          return
        endif

        if(is_writer) then

        shp_view   = (/3, ny, nz/)
        shp_sub    = (/3, ey-sy+1, ez-sz+1/)
        start_view = (/0, sy-1, sz-1/)

        ne=ne+1; lbls(ne)='type_create_subarray view'
        call MPI_TYPE_CREATE_SUBARRAY(array_rank, shp_view, shp_sub, start_view, &
                                      MPI_ORDER_FORTRAN, MPI_RP, type_view, errs(ne))
        call MPI_TYPE_COMMIT(type_view, err)

        disp = 6*int_size + int8_size
        ne=ne+1; lbls(ne)='file_set_view'
        call MPI_FILE_SET_VIEW(fh, disp, MPI_RP, type_view, "native", MPI_INFO_NULL, errs(ne))

        ne=ne+1; lbls(ne)='read_all plane'
        call MPI_FILE_READ_ALL(fh, iflow_vf_old, product(shp_sub), MPI_RP, &
                               MPI_STATUS_IGNORE, errs(ne))

        ne=ne+1; lbls(ne)='file_close'
        call MPI_FILE_CLOSE(fh, errs(ne))
        call MPI_TYPE_FREE(type_view, err)

        endif

        call MPI_BCAST(DF_ndraw, 1, MPI_INTEGER8, root, mpi_comm_cart, err)

        if(rank == root) write(*,'(1x,A,A)') ' Digital-filter state resumed from ', trim(fname)

        call check_mpi_vec('inflow_turbulent_load_state', errs(1:ne), lbls(1:ne))

        return
end subroutine inflow_turbulent_load_state

end module inflow_module
