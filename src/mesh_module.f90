module mesh_module
use mpi_module
use math_tools_module
use allocate_module

implicit none
private
real(rp), public, allocatable, dimension(:) :: x, y, z
real(rp), public, allocatable, dimension(:) :: xstep, ystep, zstep
real(rp), public, allocatable, dimension(:) :: xstep_i, ystep_i, zstep_i

real(rp), public, allocatable, dimension(:) :: x_csi , y_eta , z_zta  !< Jacobian of the trasformated coordinates
real(rp), public, allocatable, dimension(:) :: ix_csi , iy_eta , iz_zta  !< Inverse of the Jacobian of the trasformated coordinates
real(rp), public, allocatable, dimension(:) :: x_csi2, y_eta2, z_zta2 !< Jacobian of the trasformated coordinates
real(rp), public, allocatable, dimension(:) :: ixsteph, iysteph, izsteph !< grid step at the grid interface location
real(rp), public, allocatable, dimension(:) :: csistep_i, etastep_i, ztastep_i

public init_grid, compute_grid_point

contains
subroutine init_grid
! -----------------------------------------------------------
!
!       This suboutine creates the grid. 
!       The physical domain start from 1/2 and ends in n+1/2._rp 
!       This two node are the physical boundary. 
!
! -----------------------------------------------------------
  use parameters_module, only: Lx,xmin,xmax,nx,Ly,ymin,ymax,ny,Lz,zmin,zmax,nz,gbl_min_step,dims,stretching_par
  
  implicit none
  real(rp), dimension(3) :: lcl_min_step
  integer                :: lbxG, ubxG
  integer                :: lbyG, ubyG
  integer                :: lbzG, ubzG
  integer                :: err = 0

  lbxG = lbx - 3; ubxG = ubx + 3
  lbyG = lby - 3; ubyG = uby + 3
  lbzG = lbz - 3; ubzG = ubz + 3

  ! ========================== x coordinate ======================== !
  call AllocateReal(x        ,lbxG,ubxG)
  call AllocateReal(xstep    ,lbx,ubx)
  call AllocateReal(xstep_i  ,lbx,ubx)
  call AllocateReal(csistep_i,lbx,ubx)
  call AllocateReal(x_csi    ,lbx,ubx)
  call AllocateReal(ix_csi   ,lbx,ubx)
  call AllocateReal(x_csi2   ,lbx,ubx)
  call AllocateReal(ixsteph  ,lbx,ubx)
  
  Lx = xmax - xmin
  call compute_grid_point(x     ,nx,lbxG,ubxG,xmin,xmax,stretching_par,gridpoint_x)
  call compute_grid_diff1(x,x_csi ,nx,lbx ,ubx ,xmin,xmax,stretching_par,gridpoint_x)
  call compute_grid_diff2(x,x_csi2,nx,lbx ,ubx ,xmin,xmax,stretching_par,gridpoint_x)

  call compute_grid_steps(lbx,ubx,nx,x_csi,xstep,xstep_i,csistep_i,ixsteph,ix_csi)

  ! ========================== y coordinate ======================== !
  call AllocateReal(y        ,lbyG,ubyG)
  call AllocateReal(ystep    ,lby ,uby )
  call AllocateReal(ystep_i  ,lby ,uby )
  call AllocateReal(etastep_i,lby ,uby )
  call AllocateReal(y_eta    ,lby ,uby )
  call AllocateReal(iy_eta   ,lby ,uby )
  call AllocateReal(y_eta2   ,lby ,uby )
  call AllocateReal(iysteph  ,lby ,uby )
        
  Ly = ymax - ymin
  call compute_grid_point(y     ,ny,lbyG,ubyG,ymin,ymax,stretching_par,gridpoint_y)
  call compute_grid_diff1(y,y_eta ,ny,lby ,uby ,ymin,ymax,stretching_par,gridpoint_y)
  call compute_grid_diff2(y,y_eta2,ny,lby ,uby ,ymin,ymax,stretching_par,gridpoint_y)

  call compute_grid_steps(lby,uby,ny,y_eta,ystep,ystep_i,etastep_i,iysteph,iy_eta)

  ! ========================== z coordinate ======================== !
  call AllocateReal(z        ,lbzG,ubzG)
  call AllocateReal(zstep    ,lbz ,ubz )
  call AllocateReal(zstep_i  ,lbz ,ubz )
  call AllocateReal(ztastep_i,lbz ,ubz )
  call AllocateReal(z_zta    ,lbz ,ubz )
  call AllocateReal(iz_zta   ,lbz ,ubz )
  call AllocateReal(z_zta2   ,lbz ,ubz )
  call AllocateReal(izsteph  ,lbz ,ubz )

  if    (dims == 2) then

    z         = 0.0_rp
    zstep     = 1.0_rp
    zstep_i   = 1.0_rp
    ztastep_i = 1.0_rp
    z_zta     = 0.0_rp
    iz_zta    = 0.0_rp
    izsteph   = 0.0_rp

  elseif(dims == 3) then

    Lz = zmax - zmin
    call compute_grid_point(z     ,nz,lbzG,ubzG,zmin,zmax,stretching_par,gridpoint_z)
    call compute_grid_diff1(z,z_zta,nz,lbz ,ubz ,zmin,zmax,stretching_par,gridpoint_z)
    call compute_grid_diff2(z,z_zta2,nz,lbz ,ubz ,zmin,zmax,stretching_par,gridpoint_z)
    call compute_grid_steps(lbz,ubz,nz,z_zta,zstep,zstep_i,ztastep_i,izsteph,iz_zta)

  endif
  !
  ! compute minimum grid step
  !
  lcl_min_step(1) = minval(xstep(sx:ex),1)
  lcl_min_step(2) = minval(ystep(sy:ey),1)
  lcl_min_step(3) = minval(zstep(sz:ez),1)
  
  if(mpi_flag) then
    call MPI_allreduce(lcl_min_step, gbl_min_step, 3, MPI_RP, MPI_MIN, mpi_comm_cart, err)
  else
    gbl_min_step = lcl_min_step
  endif

  ! 
  ! === mid point interpolation coefficients
  ! 
  allocate(mid_point_lele_x(-fd_order/2+1:fd_order/2,lbx:ubx),stat=err)
  if(err.ne.0) stop 'Allocation error' 
  allocate(mid_point_lele_y(-fd_order/2+1:fd_order/2,lby:uby),stat=err)
  if(err.ne.0) stop 'Allocation error' 
  allocate(mid_point_lele_z(-fd_order/2+1:fd_order/2,lbz:ubz),stat=err)
  if(err.ne.0) stop 'Allocation error' 

  mid_point_lele_x = 0.0_rp
  mid_point_lele_y = 0.0_rp
  mid_point_lele_z = 0.0_rp

  call ComputeMidPointCoefficients(gridpoint_x,fd_order,lbx,ubx,x,mid_point_lele_x)
  call ComputeMidPointCoefficients(gridpoint_y,fd_order,lby,uby,y,mid_point_lele_y)
  call ComputeMidPointCoefficients(gridpoint_z,fd_order,lbz,ubz,z,mid_point_lele_z)

  return
end subroutine init_grid


subroutine compute_grid_point(x,nx,lb,ub,xmin,xmax,alpha,gridType)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x
        integer                            , intent(in)    :: nx, lb, ub
        real(rp)                           , intent(in)    :: xmin, xmax
        real(rp)                           , intent(in)    :: alpha
        character(*)                       , intent(in)    :: gridType

        ! local declarations
        selectcase(trim(gridType))
          case('uniform')
          call uniformGrid_f0(x,nx,lb,ub,xmin,xmax)

          case('cluster_two_end')
          call clusterTwoEnd_f0(x,nx,lb,ub,xmin,xmax,alpha)

          case('cluster_two_end_erf')
          call clusterTwoEndErf_f0(x,nx,lb,ub,xmin,xmax,alpha)

          case('cluster_one_end')
          call clusterOneEnd_f0(x,nx,lb,ub,xmin,xmax,alpha)

          case('cluster_in_the_centre')
          call clusterInTheCentre_f0(x,nx,lb,ub,xmin,xmax,alpha)

          case('cluster_in_the_centre_uniform')
          call clusterInTheCentreUniform_f0(x,nx,lb,ub,xmin,xmax,alpha)

          case default
          call ReadMeshFile(x,nx,lb,ub,xmin,xmax,GridType)

        end select

        return
end subroutine compute_grid_point

subroutine compute_grid_diff1(x,x_csi,nx,lb,ub,xmin,xmax,alpha,gridType)
        implicit none
        real(rp), dimension(:), allocatable, intent(in)    :: x
        real(rp), dimension(:), allocatable, intent(inout) :: x_csi
        integer                            , intent(in)    :: nx, lb, ub
        real(rp)                           , intent(in)    :: xmin, xmax
        real(rp)                           , intent(in)    :: alpha
        character(*)                       , intent(in)    :: gridType

        selectcase(trim(gridType))
          case('uniform')
          call uniformGrid_f1(x_csi,lb,ub,xmin,xmax)

          case('cluster_two_end')
          call clusterTwoEnd_f1(x_csi,nx,lb,ub,xmin,xmax,alpha)

          case('cluster_two_end_erf')
          call clusterTwoEndErf_f1(x_csi,nx,lb,ub,xmin,xmax,alpha)

          case('cluster_one_end')
          call clusterOneEnd_f1(x_csi,nx,lb,ub,xmin,xmax,alpha)

          case('cluster_in_the_centre_uniform')
          call clusterInTheCentreUniform_f1(x_csi,nx,lb,ub,xmin,xmax,alpha)

          case default
          call computeGridF1(nx,lb,ub,x,x_csi)


        end select


        return
end subroutine compute_grid_diff1

subroutine compute_grid_diff2(x,x_csi2,nx,lb,ub,xmin,xmax,alpha,gridType)
        implicit none
        real(rp), dimension(:), allocatable, intent(in)    :: x
        real(rp), dimension(:), allocatable, intent(inout) :: x_csi2
        integer                            , intent(in)    :: nx, lb, ub
        real(rp)                           , intent(in)    :: xmin, xmax
        real(rp)                           , intent(in)    :: alpha
        character(*)                       , intent(in)    :: gridType

        selectcase(trim(gridType))
          case('uniform')
          call uniformGrid_f2(x_csi2,lb,ub)

          case('cluster_two_end')
          call clusterTwoEnd_f2(x_csi2,nx,lb,ub,xmin,xmax,alpha)

          case('cluster_two_end_erf')
          call clusterTwoEndErf_f2(x_csi2,nx,lb,ub,xmin,xmax,alpha)

          case('cluster_one_end')
          call clusterOneEnd_f2(x_csi2,nx,lb,ub,xmin,xmax,alpha)

          case('cluster_in_the_centre_uniform')
          call clusterInTheCentreUniform_f2(x_csi2,nx,lb,ub,xmin,xmax,alpha)

          case default
          call computeGridF2(nx,lb,ub,x,x_csi2)


        end select


        return
end subroutine compute_grid_diff2


subroutine compute_grid_steps(lb,ub,nx,x_csi,xstep,xstep_i,csistep_i,ixsteph,ix_csi)

        implicit none
        real(rp), dimension(:), allocatable, intent(in)    :: x_csi
        real(rp), dimension(:), allocatable, intent(inout) :: xstep
        real(rp), dimension(:), allocatable, intent(inout) :: xstep_i
        real(rp), dimension(:), allocatable, intent(inout) :: csistep_i
        real(rp), dimension(:), allocatable, intent(inout) :: ixsteph
        real(rp), dimension(:), allocatable, intent(inout) :: ix_csi
        integer                            , intent(in)    :: nx, lb,ub
       
        real(rp), dimension(:), allocatable :: midC
        real(rp) :: delta_csi, xsteph_
        integer  :: i, l, fdL, fdR

        do i = lb,ub
              
           ! space between computational coordinate
           delta_csi = 1.0_rp/real(nx,rp)

           ! physical grid step
           xstep  (i) = delta_csi * x_csi(i)

           ! inverse of physical step
           xstep_i(i) = 1.0_rp/(xstep(i))

           ! inverse of physical step
           csistep_i(i) = 1.0_rp/delta_csi

           ix_csi(i) = 1.0_rp/x_csi(i)
        enddo

        !
        ! === computing half grid spacings
        !
        selectcase(fd_order)

          case(2)
            fdL = 0
            fdR = 1
            allocate(midC(fdL:fdR))
            midC(0) = 0.5_rp
            midC(1) = midC(0)

          case(4)
            fdL = -1
            fdR =  2
            allocate(midC(fdL:fdR))
            midC(-1) = -1.0_rp/16.0_rp
            midC( 0) = +9.0_rp/16.0_rp
            midC( 1) = midC( 0)
            midC( 2) = midC(-1)

          case(6)
            fdL = -2
            fdR =  3
            allocate(midC(fdL:fdR))
            midC(-2) =  3.0_rp /256.0_rp
            midC(-1) = -25.0_rp/256.0_rp
            midC( 0) =  75.0_rp/128.0_rp
            midC( 1) = midC( 0)
            midC( 2) = midC(-1)
            midC( 3) = midC(-2)
        
          case default
          stop 'fd_order not implemented'

        end select

        do i = lb + 3, ub - 3
           xsteph_ = 0.0_rp
           do l = fdL, fdR
              xsteph_ = xsteph_ + midC(l)*xstep(i+l)
           enddo
           ixsteph(i) = 1.0_rp/xsteph_
        enddo

        deallocate(midC)


        return
end subroutine compute_grid_steps






subroutine ReadMeshFile(x,nx,lb,ub,xmin,xmax,GridType)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x
        real(rp)                           , intent(in)    :: xmin, xmax
        integer                            , intent(in)    :: nx, lb, ub
        character(*)                       , intent(in)    :: GridType

        ! local declarations
        real(rp), dimension(:), allocatable :: tmp
        real(rp), parameter                 :: toll = 10.0_rp
        character(dl) :: meshFile
        character(9)  :: gchr
        integer       :: meshUnit = 51, err = 0, ngGrid = 7, i, pts
        real(rp)      :: rdummy, strGrid, endGrid
        integer       :: idummy, nlen
        logical       :: GridExist
        !
        ! === open mesh file
        !
        meshFile = trim(GridType)
        open(unit = meshUnit, file = trim(meshFile), status = 'old', iostat = err)
        if(err .ne. 0) then
          print*, ' Mesh file ', trim(GridType), ' has not been opened!'
          stop
        endif
        
        !
        ! === check if the grid is consistent with the input
        !
        pts = 0
        do 
          read(meshUnit,*,iostat=err) idummy, rdummy
          if(err .ne. 0) exit
          pts = pts + 1
        enddo
        rewind(MeshUnit)
        if(mpi_flag) call MPI_BARRIER(mpi_comm_world,err)
        
        if(pts - 2*ngGrid.ne. nx) then
          if(rank == root) then
            print*, ' --------------------------------------------------------- '
            print*, ' Number of input points: ', nx
            print*, ' Number of grid  points: ', pts - 14
            print*, ' The grid file is not consistent with the input grid points'
            print*, ' --------------------------------------------------------- '
            stop
          endif
        endif
        !
        ! === allocate a temporary array containing the whole grid
        !
        allocate(tmp(1-ngGrid:nx+ngGrid), stat=err)
        if(err .ne. 0) stop ' Allocation error in ReadMeshFile'
        !
        ! === read the file and store it in the tmp array
        !
        do i = 1-ngGrid,nx+ngGrid
           read(meshUnit,*) idummy, tmp(i)
        enddo
        close(MeshUnit)
        !
        ! === save the grid to a file in the data directory
        !
        nlen = len_trim(GridType)
        gchr = GridType(nlen-8:nlen+1)
        inquire(file = 'DATA/'//trim(data_dir)//'/GRID/'//trim(gchr), exist= GridExist)
        if(.not.GridExist) then

          call execute_command_line('mkdir -p DATA/'//trim(data_dir)//'/GRID')
          open(unit = meshUnit, &
              file = 'DATA/'//trim(data_dir)//'/GRID/'//trim(gchr), iostat = err)
          if(err .ne. 0) print*,  'WARNING: Mesh file has not been saved!'
          do i = 1-ngGrid,nx+ngGrid
             write(meshUnit,*) i, tmp(i)
          enddo
          close(MeshUnit)

          if(rank == root .and. err == 0) &
          print*, ' The grid file ', trim(GridType), ' has been saved in DATA/',trim(data_dir),'/GRID/'
        endif
        if(mpi_flag) call MPI_BARRIER(mpi_comm_world,err)
        !
        ! === check dimension consistency
        !
        strGrid = 0.5_rp*(tmp(0) + tmp(1))
        endGrid = 0.5_rp*(tmp(nx)+ tmp(nx+1))

        if(abs(strGrid-xmin) > toll .or. abs(endGrid-xmax)> toll) then
          print*, '--------------------------------------- '
          print*, ' The grid   is :', strGrid, ' x ', endGrid
          print*, ' The domain is :', xmin   , ' x ', xmax
          print*, ' The grid is not consistent with inputs!'
          print*, '--------------------------------------- '
          stop
        endif
        !
        ! === distribute the grid in the procs
        !
        do i = lb,ub
           x(i) = tmp(i)
        enddo


        deallocate(tmp)

        if(mpi_flag) call MPI_BARRIER(mpi_comm_world,err)

        return
end subroutine ReadMeshFile



subroutine computeGridF1(nx,lb,ub,x,x_csi)
        
        use parameters_module, only: fd_order

        implicit none
        real(rp), allocatable, dimension(:), intent(in)    :: x
        real(rp), allocatable, dimension(:), intent(inout) :: x_csi
        integer                            , intent(in)    :: nx,lb, ub

        ! local declarations
        real(rp), allocatable, dimension(:) :: central_1
        real(rp)                            :: idxi, df
        integer                             :: i, s, fdO
        
        
        fdO = int(fd_order/2)
        allocate(central_1(-fdO:fdO)); central_1 = 0.0_rp

        !
        ! === select the derivative order
        !
        selectcase(fd_order)
        case(2)
          central_1(-1) = - 0.5_rp
          central_1( 0) =   0.0_rp
          central_1( 1) = - central_1(-1)
        case(4)
          central_1(-2) =   1.0_rp/12.0_rp
          central_1(-1) = - 2.0_rp/ 3.0_rp
          central_1( 0) =   0.0_rp
          central_1( 1) = - central_1(-1)
          central_1( 2) = - central_1(-2)
        case(6)
          central_1(-3) = - 1.0_rp/60.0_rp
          central_1(-2) =   3.0_rp/20.0_rp
          central_1(-1) = - 3.0_rp/ 4.0_rp
          central_1( 0) =   0.0_rp
          central_1( 1) = - central_1(-1)
          central_1( 2) = - central_1(-2)
          central_1( 3) = - central_1(-3)
        endselect
        if(sum(central_1) > 1.0E-14_rp) stop ' ERROR in computeGridF1'

        !
        ! === compute the derivative
        !
        idxi = real(nx,rp)
        do i = lb,ub
        
           df = 0.0_rp
           do s = -fdO, fdO
              df = df + central_1(s) * x(i+s)
           enddo

           x_csi(i) = idxi*df

        enddo

        deallocate(central_1)

        return
end subroutine computeGridF1


subroutine computeGridF2(nx,lb,ub,x,x_csi2)

        use parameters_module, only: fd_order

        implicit none
        real(rp), allocatable, dimension(:), intent(in)    :: x
        real(rp), allocatable, dimension(:), intent(inout) :: x_csi2
        integer                            , intent(in)    :: nx,lb, ub

        ! local declarations
        real(rp), allocatable, dimension(:) :: central_2
        integer                             :: i, s, fdO
        real(rp)                            :: idxi2, d2f
        
        fdO   = int(fd_order/2)
        allocate(central_2(-fdO:fdO)); central_2 = 0.0_rp
        
        !
        ! === select the derivative order
        !
        selectcase(fd_order)
        case(2)
          central_2(-1) =   1.0_rp
          central_2( 0) = - 2.0_rp
          central_2( 1) =   central_2(-1)
        case(4)
          central_2(-2) = - 1.0_rp/12.0_rp
          central_2(-1) =   4.0_rp/ 3.0_rp
          central_2( 0) = - 5.0_rp/ 2.0_rp
          central_2( 1) =   central_2(-1)
          central_2( 2) =   central_2(-2)
        case(6)
          central_2(-3) =    1.0_rp/90.0_rp
          central_2(-2) = -  3.0_rp/20.0_rp
          central_2(-1) =    3.0_rp/ 2.0_rp
          central_2( 0) = - 49.0_rp/18.0_rp
          central_2( 1) = central_2(-1)
          central_2( 2) = central_2(-2)
          central_2( 3) = central_2(-3)
        endselect
        if(sum(central_2) > 1.0E-14_rp) stop ' ERROR in computeGridF2'

        !
        ! === compute the derivative
        !
        idxi2 = (real(nx,rp))**2
        do i = lb,ub
        
           d2f = 0.0_rp
           do s = -fdO, fdO
              d2f = d2f + central_2(s) * x(i+s)
           enddo

           x_csi2(i) = idxi2*d2f

        enddo

        return
end subroutine computeGridF2




subroutine ComputeMidPointCoefficients(gridpoint,fd_order,lbx,ubx,x,mid_point_coef)
        
        use matrix_inversion_module
               
        implicit none
        real(rp), allocatable, dimension(:)  , intent(in)    :: x
        integer                              , intent(in)    :: fd_order
        integer                              , intent(in)    :: lbx,ubx
        character(*)                         , intent(in)    :: gridpoint
        real(rp), allocatable, dimension(:,:), intent(inout) :: mid_point_coef
        
        real(rp), parameter :: toll = 1.0E-13_rp
        integer :: mm , i, l, ll
        real(rp), dimension(1:fd_order, 1:fd_order)   :: amat
        real(rp), dimension(1:fd_order, 1:fd_order)   :: imat
        real(rp), dimension(-fd_order/2+1:fd_order/2) :: c
        real(rp) :: xx, xxx, expl

        mm = fd_order/2

        if(trim(gridpoint)=='uniform') then
          !
          ! === analytical solution
          !
          selectcase(fd_order)
            case(2)
              c(0) = 0.5_rp
              c(1) = c(0)

            case(4)
              c(-1) = -1.0_rp/16.0_rp
              c( 0) = +9.0_rp/16.0_rp
              c( 1) = c( 0)
              c( 2) = c(-1)

            case(6)
              c(-2) =  3.0_rp /256.0_rp
              c(-1) = -25.0_rp/256.0_rp
              c( 0) =  75.0_rp/128.0_rp
              c( 1) = c( 0)
              c( 2) = c(-1)
              c( 3) = c(-2)

          endselect
          do i = lbx,ubx

             do l = -fd_order/2+1,fd_order/2
                mid_point_coef(l,i) = c(l)
             enddo

             if(abs(sum(mid_point_coef(:,i))-1.0_rp)>1.0E-14_rp) then
               print*, abs(sum(mid_point_coef(:,i)) - 1.0_rp)
               stop 'Mid Point Coefficients are umbalanced'
             endif

          enddo

        else
          !
          ! === non uniform grids (numerical solution)
          !
          do i = lbx,ubx-1

             ! coordinata della faccia
             xx = 0.5_rp*(x(i) + x(i+1))     
             do l = -mm+1,mm
                ! coordinata dei nodi
                xxx = x(i+l) 
                do ll = 1,fd_order

                   expl = ll-1
                   amat(l+mm,ll) = (xxx-xx)**expl

                enddo
             enddo

             call inverse(amat,imat,fd_order)

             do l = -mm+1,mm
                c(l) = imat(1,mm+l)
             enddo
             mid_point_coef(:,i) = c(:)

             if(abs(sum(mid_point_coef(:,i))-1.0_rp)>toll) then
               print*, abs(sum(mid_point_coef(:,i)) - 1.0_rp)
               stop 'Mid Point Coefficients are umbalanced'
             endif

          enddo

        endif

        return
end subroutine ComputeMidPointCoefficients






! --------------------------------------------------------------------------------
! --------------------------------------------------------------------------------
! --------------------------------------------------------------------------------
subroutine uniformGrid_f0(x,nx,lb,ub,xmin,xmax)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x
        integer                            , intent(in)    :: nx,lb,ub
        real(rp)                           , intent(in)    :: xmin, xmax
        ! local declarations
        real(rp) :: Lx, xi
        integer  :: i
        
        Lx = xmax - xmin
        do i = lb, ub

           xi   = (i-0.5_rp)/real(nx,rp)
           x(i) = xmin + Lx*xi

        enddo

        return
end subroutine uniformGrid_f0
!
subroutine uniformGrid_f1(x_csi,lb,ub,xmin,xmax)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x_csi
        integer                            , intent(in)    :: lb,ub
        real(rp)                           , intent(in)    :: xmin, xmax
        ! local declarations
        real(rp) :: Lx

        Lx = xmax - xmin
        x_csi(lb:ub) = Lx
        
        return
end subroutine uniformGrid_f1
!
subroutine uniformGrid_f2(x_csi2,lb,ub)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x_csi2
        integer                            , intent(in)    :: lb,ub
        
        x_csi2(lb:ub) = 0.0_rp
        
        return
end subroutine uniformGrid_f2
! --------------------------------------------------------------------------------
! --------------------------------------------------------------------------------
! --------------------------------------------------------------------------------
subroutine clusterTwoEnd_f0(x,nx,lb,ub,xmin,xmax,a)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x
        integer                            , intent(in)    :: nx,lb,ub
        real(rp)                           , intent(in)    :: xmin, xmax, a
        
        ! local declarations
        integer  :: i
        real(rp) :: xi, f_xi, Lx
        
        Lx = xmax - xmin
        do i = lb, ub

           xi = (i-0.5_rp)/real(nx,rp)

           f_xi = Lx*0.5_rp*(tanh(a*(xi-0.5_rp)))/(tanh(0.5_rp*a))

           if(i <  1) f_xi = 2*xmin - Lx*cluster_two_end(    -i+1,nx,a)
           if(i > nx) f_xi = 2*xmax - Lx*cluster_two_end(2*nx-i+1,nx,a)

           x(i) = f_xi

        enddo

        return
end subroutine clusterTwoEnd_f0
!
subroutine clusterTwoEnd_f1(x_csi,nx,lb,ub,xmin,xmax,a)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x_csi
        integer                            , intent(in)    :: nx,lb,ub
        real(rp)                           , intent(in)    :: xmin, xmax, a
        ! local declarations
        integer  :: i
        real(rp) :: xi,  df_dxi, Lx
        
        Lx = xmax - xmin
        do i = lb,ub

           xi = (i-0.5_rp)/real(nx,rp)

           df_dxi = 0.5_rp*a*(cotanh(0.5_rp*a))*(sech(a*(xi-0.5_rp)))**2

           x_csi(i) = Lx * df_dxi

        enddo

        return
end subroutine clusterTwoEnd_f1
!
subroutine clusterTwoEnd_f2(x_csi2,nx,lb,ub,xmin,xmax,a)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x_csi2
        integer                            , intent(in)    :: nx,lb,ub
        real(rp)                           , intent(in)    :: xmin, xmax, a
        ! local declarations
        integer  :: i
        real(rp) :: xi,  d2f_dxi2, Lx

        Lx = xmax - xmin
        
        do i = lb,ub
           xi = (i-0.5_rp)/real(nx,rp)

           d2f_dxi2 = -a**2*(cotanh(0.5_rp*a))   *&
                      tanh(a*(xi-0.5_rp))        *&
                      (sech(a*(xi-0.5_rp)))**2

           x_csi2(i) = Lx * d2f_dxi2
        enddo

        return
end subroutine clusterTwoEnd_f2
! --------------------------------------------------------------------------------
! --------------------------------------------------------------------------------
! --------------------------------------------------------------------------------
subroutine clusterTwoEndErf_f0(x,nx,lb,ub,xmin,xmax,a)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x
        integer                            , intent(in)    :: nx,lb,ub
        real(rp)                           , intent(in)    :: xmin, xmax, a
        
        ! local declarations
        integer  :: i
        real(rp) :: xi, f_xi, Lx
        
        Lx = xmax - xmin
        do i = lb, ub

           xi = (i-0.5_rp)/real(nx,rp)

           f_xi = Lx*0.5_rp*erf(a*(xi-0.5_rp))/erf(0.5_rp*a)

           if(i <  1) f_xi = 2*xmin - Lx*cluster_two_end_erf(    -i+1,nx,a)
           if(i > nx) f_xi = 2*xmax - Lx*cluster_two_end_erf(2*nx-i+1,nx,a)

           x(i) = f_xi

        enddo

        return
end subroutine clusterTwoEndErf_f0
!
subroutine clusterTwoEndErf_f1(x_csi,nx,lb,ub,xmin,xmax,a)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x_csi
        integer                            , intent(in)    :: nx,lb,ub
        real(rp)                           , intent(in)    :: xmin, xmax, a
        ! local declarations
        integer             :: i
        real(rp)            :: xi,  df_dxi, Lx
        real(rp), parameter :: pi = 2.0_rp*asin(1.0_rp) 
        real(rp), parameter :: p  = 1.0_rp/sqrt(pi)
        
        Lx = xmax - xmin
        do i = lb,ub

           xi = (i-0.5_rp)/real(nx,rp)

           df_dxi = p*a*exp(-a**2*(xi-0.5_rp)**2)/erf(0.5_rp*a)

           x_csi(i) = Lx * df_dxi

        enddo

        return
end subroutine clusterTwoEndErf_f1
!
subroutine clusterTwoEndErf_f2(x_csi2,nx,lb,ub,xmin,xmax,a)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x_csi2
        integer                            , intent(in)    :: nx,lb,ub
        real(rp)                           , intent(in)    :: xmin, xmax, a
        ! local declarations
        integer             :: i
        real(rp)            :: xi,  d2f_dxi2, Lx
        real(rp), parameter :: pi = 2.0_rp*asin(1.0_rp) 
        real(rp), parameter :: p  = 2.0_rp/sqrt(pi)

        Lx = xmax - xmin
        
        do i = lb,ub
           xi = (i-0.5_rp)/real(nx,rp)

           d2f_dxi2 = - p*a**3*(xi-0.5_rp)*exp(-a**2*(xi-0.5_rp)**2)/erf(0.5_rp*a)

           x_csi2(i) = Lx * d2f_dxi2
        enddo

        return
end subroutine clusterTwoEndErf_f2
! --------------------------------------------------------------------------------
! --------------------------------------------------------------------------------
! --------------------------------------------------------------------------------
subroutine clusterOneEnd_f0(x,nx,lb,ub,xmin,xmax,a)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x
        integer                            , intent(in)    :: nx,lb,ub
        real(rp)                           , intent(in)    :: xmin, xmax, a
        
        ! local declarations
        integer  :: i
        real(rp) :: xi, f_xi, Lx
        
        Lx = xmax - xmin
        do i = lb, ub

           xi = (i-0.5_rp)/real(nx,rp)

           f_xi = Lx*(1.0_rp+tanh((xi-1.0_rp)*a)/tanh(a))      

           if(i <  1) f_xi = 2*xmin - Lx*cluster_one_end(    -i+1,nx,a)
           if(i > nx) f_xi = 2*xmax - Lx*cluster_one_end(2*nx-i+1,nx,a)

           x(i) = f_xi

        enddo

        return
end subroutine clusterOneEnd_f0
!
subroutine clusterOneEnd_f1(x_csi,nx,lb,ub,xmin,xmax,a)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x_csi
        integer                            , intent(in)    :: nx,lb,ub
        real(rp)                           , intent(in)    :: xmin, xmax, a
        ! local declarations
        integer  :: i
        real(rp) :: xi,  df_dxi, Lx
        
        Lx = xmax - xmin
        do i = lb,ub

           xi = (i-0.5_rp)/real(nx,rp)

           df_dxi = a*cotanh(a)*(sech(a*(xi-1.0_rp)))**2

           x_csi(i) = Lx * df_dxi

        enddo

        return
end subroutine clusterOneEnd_f1
!
subroutine clusterOneEnd_f2(x_csi2,nx,lb,ub,xmin,xmax,a)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x_csi2
        integer                            , intent(in)    :: nx,lb,ub
        real(rp)                           , intent(in)    :: xmin, xmax, a
        ! local declarations
        integer  :: i
        real(rp) :: xi,  d2f_dxi2, Lx

        Lx = xmax - xmin
        
        do i = lb,ub
           xi = (i-0.5_rp)/real(nx,rp)

           d2f_dxi2 = (-2*a**2)               *&
                      cotanh(a)               *&
                      tanh(a*(xi-1.0_rp))        *&
                      (sech(a*(xi-1.0_rp)))**2

           x_csi2(i) = Lx * d2f_dxi2
        enddo

        return
end subroutine clusterOneEnd_f2
! --------------------------------------------------------------------------------
! --------------------------------------------------------------------------------
! --------------------------------------------------------------------------------
subroutine clusterInTheCentre_f0(x,nx,lb,ub,xmin,xmax,a)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x
        integer                            , intent(in)    :: nx,lb,ub
        real(rp)                           , intent(in)    :: xmin, xmax, a
        
        ! local declarations
        integer  :: i
        real(rp) :: xi, f_xi, Lx
        
        Lx = xmax - xmin
        do i = lb, ub

           xi = (i-0.5_rp)/real(nx,rp)

           f_xi = 0.0_rp
           if    (xi <= 0.5_rp) then
             f_xi = 0.5_rp*(-1.0_rp + tanh(2*a*xi      )/(tanh(a)))

           elseif(xi > 0.5_rp) then
             f_xi = 0.5_rp*( 1.0_rp + tanh(2*a*(xi-1.0_rp))/(tanh(a)))

           endif

           x(i) = Lx*f_xi

        enddo

        return
end subroutine clusterInTheCentre_f0
!
!subroutine clusterInTheCentre_f1(x_csi,nx,lb,ub,xmin,xmax,a)
!        implicit none
!        real(rp), dimension(:), allocatable, intent(inout) :: x_csi
!        integer                            , intent(in)    :: nx,lb,ub
!        real(rp)                           , intent(in)    :: xmin, xmax, a
!        ! local declarations
!        integer  :: i
!        real(rp) :: xi,  df_dxi, Lx
!        
!        Lx = xmax - xmin
!        do i = lb,ub
!
!           xi = (i-0.5_rp)/real(nx,rp)
!
!           df_dxi = 0.0_rp
!           if    (xi <= 0.5_rp) then
!             df_dxi = a*cotanh(a)*(sech(2*a*xi))**2
!
!           elseif(xi >  0.5_rp) then
!             df_dxi = a*cotanh(a)*(sech(2*a*(xi-1.0_rp)))**2
!
!           endif
!
!           x_csi(i) = Lx * df_dxi
!
!        enddo
!
!        return
!end subroutine clusterInTheCentre_f1
!!
!subroutine clusterInTheCentre_f2(x_csi2,nx,lb,ub,xmin,xmax,a)
!        implicit none
!        real(rp), dimension(:), allocatable, intent(inout) :: x_csi2
!        integer                            , intent(in)    :: nx,lb,ub
!        real(rp)                           , intent(in)    :: xmin, xmax, a
!        ! local declarations
!        integer  :: i
!        real(rp) :: xi,  d2f_dxi2, Lx
!
!        Lx = xmax - xmin
!        
!        do i = lb,ub
!
!           xi = (i-0.5_rp)/real(nx,rp)
!
!           d2f_dxi2 = 0.0_rp
!           if    (xi <= 0.5_rp) then
!             d2f_dxi2 = -4*a**2*cotanh(a)*tanh(2*a*xi)*(sech(2*a*xi))**2
!
!           elseif(xi >  0.5_rp) then
!             d2f_dxi2 = -4*a**2*cotanh(a)*tanh(2*a*(xi-1.0_rp))*(sech(2*a*(xi-1.0_rp)))**2
!
!           endif
!
!           x_csi2(i) = Lx * d2f_dxi2
!
!        enddo
!
!        return
!end subroutine clusterInTheCentre_f2
! --------------------------------------------------------------------------------
! --------------------------------------------------------------------------------
! --------------------------------------------------------------------------------
subroutine clusterInTheCentreUniform_f0(x,nx,lb,ub,xmin,xmax,a)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x
        integer                            , intent(in)    :: nx,lb,ub
        real(rp)                           , intent(in)    :: xmin, xmax, a
        
        ! local declarations
        integer  :: i
        real(rp) :: xmin_i, xmax_i
        real(rp) :: xmin_o, xmax_o
        real(rp) :: xx, xi, f_xi
        real(rp) :: lenght, alpha, alpha_i, alpha_o, eta, x_c, eta_c
        
        ! inner domain extrema
        xmin_i = -abs(xmin)
        xmax_i =  abs(xmin)

        ! outer domain extrema
        xmin_o = -abs(xmax)
        xmax_o =  abs(xmax)

        lenght = xmax_o - xmin_o
        Lx     = xmax - xmin
        x_c    = (xmax_o-xmin_i)/15.0_rp ! size of the inner region  minimum grid step

        eta_c = a
        if(eta_c > 0.4_rp .or. eta_c < 0.2_rp) then
          stop 'FATAL ERROR: Stretching parameter must be in the range [0.2_rp : 0.4_rp]'
        endif
        
        call newton_raphson(opt_alpha,xmax_i,x_c,eta_c,1000, 1.0E-14_rp, 10.0_rp, alpha_i)
        call newton_raphson(opt_alpha,xmax_o,x_c,eta_c,1000, 1.0E-14_rp, 10.0_rp, alpha_o)
        if(abs(alpha_i)>1.0E+14_rp) stop ' Error in cluster_in_the_centre_uniform'
        if(abs(alpha_o)>1.0E+14_rp) stop ' Error in cluster_in_the_centre_uniform'

        do i = lb, ub

           xi  = (i-0.5_rp)/real(nx,rp)
           eta = xi - 0.5_rp
           xx  = xmin_o + lenght*xi

           f_xi = 0.0_rp
           if(abs(eta) < eta_c) then
             ! CENTRAL REGION
             f_xi = eta*x_c/eta_c
   
           elseif(abs(eta) > eta_c .and. xx < 0.0_rp) then
             ! LEFT REGION
             alpha = alpha_i
   
             f_xi = eta/abs(eta)*xmax_i + tanh(alpha*(eta-0.5_rp*eta/abs(eta)))/(tanh(alpha*(0.5_rp-eta_c)))*(xmax_i-x_c)
   
           elseif(abs(eta) > eta_c .and. xx > 0.0_rp) then
             ! RIGHT REGION
             alpha = alpha_o
   
             f_xi = eta/abs(eta)*xmax_o + tanh(alpha*(eta-0.5_rp*eta/abs(eta)))/(tanh(alpha*(0.5_rp-eta_c)))*(xmax_o-x_c)
   
           endif


           x(i) = f_xi

        enddo

        return
end subroutine clusterInTheCentreUniform_f0
!
subroutine clusterInTheCentreUniform_f1(x_csi,nx,lb,ub,xmin,xmax,a)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x_csi
        integer                            , intent(in)    :: nx,lb,ub
        real(rp)                           , intent(in)    :: xmin, xmax, a
        ! local declarations
        integer  :: i
        real(rp) :: xx, xi,  df_dxi, Lx
        real(rp) :: xmin_i, xmax_i
        real(rp) :: xmin_o, xmax_o
        real(rp) :: lenght, alpha, alpha_i, alpha_o, eta, x_c, eta_c
        
        ! inner domain extrema
        xmin_i = -abs(xmin)
        xmax_i =  abs(xmin)

        ! outer domain extrema
        xmin_o = -abs(xmax)
        xmax_o =  abs(xmax)

        lenght = xmax_o - xmin_o
        Lx     = xmax - xmin
        x_c    = (xmax_o-xmin_i)/15.0_rp ! size of the inner region  minimum grid step

        eta_c = a
        if(eta_c > 0.4_rp .or. eta_c < 0.2_rp) then
          stop 'FATAL ERROR: Stretching parameter must be in the range [0.2_rp : 0.4_rp]'
        endif
        
        call newton_raphson(opt_alpha,xmax_i,x_c,eta_c,1000, 1.0E-14_rp, 10.0_rp, alpha_i)
        call newton_raphson(opt_alpha,xmax_o,x_c,eta_c,1000, 1.0E-14_rp, 10.0_rp, alpha_o)
        if(abs(alpha_i)>1.0E+14_rp) stop ' Error in cluster_in_the_centre_uniform'
        if(abs(alpha_o)>1.0E+14_rp) stop ' Error in cluster_in_the_centre_uniform'

        do i = lb, ub

           xi  = (i-0.5_rp)/real(nx,rp)
           eta = xi - 0.5_rp
           xx  = xmin_o + lenght*xi

           df_dxi = 0.0_rp
           if(abs(eta) < eta_c) then
             ! CENTRAL REGION
             df_dxi = x_c/eta_c
   
           elseif(abs(eta) > eta_c .and. xx < 0.0_rp) then
             ! LEFT REGION
             alpha  = alpha_i
   
             df_dxi = alpha*(sech( alpha*eta*(abs(eta)-0.5_rp)/abs(eta) ))**2/(tanh(alpha*(0.5_rp-eta_c)))*(xmax_i-x_c)
   
           elseif(abs(eta) > eta_c .and. xx > 0.0_rp) then
             ! RIGHT REGION
             alpha  = alpha_o
   
             df_dxi = alpha*(sech( alpha*eta*(abs(eta)-0.5_rp)/abs(eta) ))**2/(tanh(alpha*(0.5_rp-eta_c)))*(xmax_o-x_c)
   
           endif

           x_csi(i) = df_dxi

        enddo

        return
end subroutine clusterInTheCentreUniform_f1
!
subroutine clusterInTheCentreUniform_f2(x_csi2,nx,lb,ub,xmin,xmax,a)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x_csi2
        integer                            , intent(in)    :: nx,lb,ub
        real(rp)                           , intent(in)    :: xmin, xmax, a
        ! local declarations
        integer  :: i
        real(rp) :: xx, xi, d2f_dxi2, Lx, arg, denom
        real(rp) :: xmin_i, xmax_i
        real(rp) :: xmin_o, xmax_o
        real(rp) :: lenght, alpha, alpha_i, alpha_o, eta, x_c, eta_c

        ! inner domain extrema
        xmin_i = -abs(xmin)
        xmax_i =  abs(xmin)

        ! outer domain extrema
        xmin_o = -abs(xmax)
        xmax_o =  abs(xmax)

        lenght = xmax_o - xmin_o
        Lx     = xmax - xmin
        x_c    = (xmax_o-xmin_i)/15.0_rp ! size of the inner region  minimum grid step

        eta_c = a
        if(eta_c > 0.4_rp .or. eta_c < 0.2_rp) then
          stop 'FATAL ERROR: Stretching parameter must be in the range [0.2_rp : 0.4_rp]'
        endif
        
        call newton_raphson(opt_alpha,xmax_i,x_c,eta_c,1000, 1.0E-14_rp, 10.0_rp, alpha_i)
        call newton_raphson(opt_alpha,xmax_o,x_c,eta_c,1000, 1.0E-14_rp, 10.0_rp, alpha_o)
        if(abs(alpha_i)>1.0E+14_rp) stop ' Error in cluster_in_the_centre_uniform'
        if(abs(alpha_o)>1.0E+14_rp) stop ' Error in cluster_in_the_centre_uniform'

        do i = lb, ub

           xi  = (i-0.5_rp)/real(nx,rp)
           eta = xi - 0.5_rp
           xx  = xmin_o + lenght*xi

           d2f_dxi2 = 0.0_rp
           if(abs(eta) < eta_c) then
             ! CENTRAL REGION
             d2f_dxi2 = 0.0_rp

           elseif(abs(eta) > eta_c .and. xx < 0.0_rp) then
             ! LEFT REGION
             alpha = alpha_i
             denom = 1.0_rp/(tanh(alpha*(0.5_rp-eta_c)))*(xmax_i-x_c)
             arg   = alpha*eta*(abs(eta)-0.5_rp)/abs(eta)

             d2f_dxi2 = -2*alpha**2*tanh(arg)*(sech(arg))**2*denom

           elseif(abs(eta) > eta_c .and. xx > 0.0_rp) then
             ! RIGHT REGION
             alpha = alpha_o
             denom = 1.0_rp/(tanh(alpha*(0.5_rp-eta_c)))*(xmax_o-x_c)
             arg   = alpha*eta*(abs(eta)-0.5_rp)/abs(eta)

             d2f_dxi2 = -2*alpha**2*tanh(arg)*(sech(arg))**2*denom

           endif

           x_csi2(i) = d2f_dxi2

        enddo

        return
end subroutine clusterInTheCentreUniform_f2
! --------------------------------------------------------------------------------
! --------------------------------------------------------------------------------
! --------------------------------------------------------------------------------





! 
! === cluster two end functions
!
function cluster_two_end(i,nx,a) result(f_xi)
        implicit none
        integer , intent(in) :: i
        integer , intent(in) :: nx
        real(rp), intent(in) :: a
        real(rp)             :: f_xi
        real(rp)             :: xi

        xi = (i-0.5_rp)/real(nx,rp)

        f_xi = 0.5_rp*(tanh(a*(xi-0.5_rp)))/(tanh(0.5_rp*a))
        
        return
end function cluster_two_end

! 
! === cluster two end ERF functions
!
function cluster_two_end_erf(i,nx,a) result(f_xi)
        implicit none
        integer , intent(in) :: i
        integer , intent(in) :: nx
        real(rp), intent(in) :: a
        real(rp)             :: f_xi
        real(rp)             :: xi

        xi = (i-0.5_rp)/real(nx,rp)

        f_xi = 0.5_rp*erf(a*(xi-0.5_rp))/erf(0.5_rp*a)
        
        return
end function cluster_two_end_erf

! 
! === cluster one end functions
!
function cluster_one_end(i,nx,a) result(f_xi)
        implicit none
        integer , intent(in) :: i
        integer , intent(in) :: nx
        real(rp), intent(in) :: a
        real(rp)             :: f_xi
        real(rp)             :: xi

        xi   = (i-0.5_rp)/real(nx,rp)

        f_xi = (1.0_rp+tanh((xi-1.0_rp)*a)/tanh(a))      

        return
end function cluster_one_end


function opt_alpha(alpha,xmax,x_c,eta_c) result(f_alpha)
! -------------------------------------------------------------------------------------------
!       This function  grants the continuity of the first derivative of the mapping function
!       for cluster in the centre uniform grid. 
!       REF: Orlandi, fluid flow phenomena, pag 13
! -------------------------------------------------------------------------------------------
        implicit none
        real(rp), intent(in) :: alpha
        real(rp), intent(in) :: xmax
        real(rp), intent(in) :: x_c
        real(rp), intent(in) :: eta_c
        real(rp)             :: f_alpha

        real(rp) :: zta_f 

        zta_f = 0.5_rp-eta_c

        f_alpha = alpha*(xmax-x_c)/(sinh(alpha*zta_f)*cosh(alpha*zta_f)) - x_c/eta_c

        return
end function opt_alpha





end module mesh_module
