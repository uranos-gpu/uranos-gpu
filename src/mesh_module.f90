module mesh_module
use mpi_module
use allocate_module

implicit none
private
real(rp), public, allocatable, dimension(:) :: x, y, z
real(rp), public, allocatable, dimension(:) :: xstep, ystep, zstep
real(rp), public, allocatable, dimension(:) :: xstep_i, ystep_i, zstep_i

real(rp), public, allocatable, dimension(:) :: x_csi , y_eta , z_zta
real(rp), public, allocatable, dimension(:) :: ix_csi , iy_eta , iz_zta
real(rp), public, allocatable, dimension(:) :: x_csi2, y_eta2, z_zta2
real(rp), public, allocatable, dimension(:) :: ixsteph, iysteph, izsteph
real(rp), public, allocatable, dimension(:) :: csistep_i, etastep_i, ztastep_i

real(rp), public, allocatable, dimension(:,:) :: xx, yy
real(rp), public, allocatable, dimension(:,:) :: dxi_dx , dxi_dy
real(rp), public, allocatable, dimension(:,:) :: deta_dx, deta_dy
real(rp), public, allocatable, dimension(:,:) :: Jac2D
real(rp), public, allocatable, dimension(:,:) :: dx_dxi , dx_deta
real(rp), public, allocatable, dimension(:,:) :: dy_dxi , dy_deta
real(rp), public, allocatable, dimension(:,:) :: m_xi   , m_eta
real(rp), public, allocatable, dimension(:,:) :: nx_xi  , ny_xi
real(rp), public, allocatable, dimension(:,:) :: nx_eta , ny_eta
real(rp), public, allocatable, dimension(:,:) :: Gamma1_2D, Gamma2_2D
real(rp), public, allocatable, dimension(:,:) :: dxi_dxx , dxi_dyy , dxi_dxy
real(rp), public, allocatable, dimension(:,:) :: deta_dxx, deta_dyy, deta_dxy
real(rp), public, allocatable, dimension(:,:) :: gJ11, gJ22, gJ12
real(rp), public, allocatable, dimension(:,:) :: fmx_csix, fmx_csiy, fmx_etax, fmx_etay
real(rp), public, allocatable, dimension(:,:) :: fmx_cof1, fmx_cof2
real(rp), public, allocatable, dimension(:,:) :: fmy_csix, fmy_csiy, fmy_etax, fmy_etay
real(rp), public, allocatable, dimension(:,:) :: fmy_cof1, fmy_cof2
real(rp), public, allocatable, dimension(:)   :: fmz_zeta

logical, public :: non_orthogonal = .false.
logical, public :: inplane_axis_aligned = .false.

integer, parameter :: mesh_gh = 3

public init_grid, compute_grid_point, destroy_grid, west_wall_normal_distance, &
       bottom_wall_streamwise_distance

contains

subroutine init_grid
        use parameters_module, only: Lx,xmin,xmax,nx,Ly,ymin,ymax,ny,Lz,zmin,zmax,nz,gbl_min_step,dims,mesh_file

        implicit none
        real(rp), dimension(3) :: lcl_min_step
        integer                :: lbxG, ubxG
        integer                :: lbyG, ubyG
        integer                :: lbzG, ubzG
        integer                :: i, j, err = 0

        if(fd_order/2 > mesh_gh) then
          if(rank == root) print*, ' ERROR: fd_order =', fd_order, &
                                   ' is too large for the mesh ghost width (', mesh_gh, ')'
          call secure_stop
        endif

        lbxG = lbx - mesh_gh; ubxG = ubx + mesh_gh
        lbyG = lby - mesh_gh; ubyG = uby + mesh_gh
        lbzG = lbz - mesh_gh; ubzG = ubz + mesh_gh

        call AllocateReal(x        ,lbxG,ubxG)
        call AllocateReal(xstep    ,lbx ,ubx )
        call AllocateReal(xstep_i  ,lbx ,ubx )
        call AllocateReal(csistep_i,lbx ,ubx )
        call AllocateReal(x_csi    ,lbx ,ubx )
        call AllocateReal(ix_csi   ,lbx ,ubx )
        call AllocateReal(x_csi2   ,lbx ,ubx )
        call AllocateReal(ixsteph  ,lbx ,ubx )

        Lx = xmax - xmin
        call compute_grid_point(x     ,nx,lbxG,ubxG,xmin,xmax,gridpoint_x,'x')
        call compute_grid_diff1(x,x_csi ,nx,lbx,ubx,xmin,xmax,gridpoint_x)
        call compute_grid_diff2(x,x_csi2,nx,lbx,ubx,gridpoint_x)

        call compute_grid_steps(lbx,ubx,nx,x_csi,xstep,xstep_i,csistep_i,ixsteph,ix_csi)

        call AllocateReal(y        ,lbyG,ubyG)
        call AllocateReal(ystep    ,lby ,uby )
        call AllocateReal(ystep_i  ,lby ,uby )
        call AllocateReal(etastep_i,lby ,uby )
        call AllocateReal(y_eta    ,lby ,uby )
        call AllocateReal(iy_eta   ,lby ,uby )
        call AllocateReal(y_eta2   ,lby ,uby )
        call AllocateReal(iysteph  ,lby ,uby )

        Ly = ymax - ymin
        call compute_grid_point(y     ,ny,lbyG,ubyG,ymin,ymax,gridpoint_y,'y')
        call compute_grid_diff1(y,y_eta ,ny,lby,uby,ymin,ymax,gridpoint_y)
        call compute_grid_diff2(y,y_eta2,ny,lby,uby,gridpoint_y)

        call compute_grid_steps(lby,uby,ny,y_eta,ystep,ystep_i,etastep_i,iysteph,iy_eta)

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
          call compute_grid_point(z     ,nz,lbzG,ubzG,zmin,zmax,gridpoint_z,'z')
          call compute_grid_diff1(z,z_zta ,nz,lbz,ubz,zmin,zmax,gridpoint_z)
          call compute_grid_diff2(z,z_zta2,nz,lbz,ubz,gridpoint_z)
          call compute_grid_steps(lbz,ubz,nz,z_zta,zstep,zstep_i,ztastep_i,izsteph,iz_zta)

        endif

        lcl_min_step(1) = minval(xstep(sx:ex),1)
        lcl_min_step(2) = minval(ystep(sy:ey),1)
        lcl_min_step(3) = minval(zstep(sz:ez),1)

        if(mpi_flag) then
          call MPI_allreduce(lcl_min_step, gbl_min_step, 3, MPI_RP, MPI_MIN, mpi_comm_cart, err)
        else
          gbl_min_step = lcl_min_step
        endif

        call AllocateReal(mid_point_lele_x,-2,3,lbx,ubx)
        call AllocateReal(mid_point_lele_y,-2,3,lby,uby)
        call AllocateReal(mid_point_lele_z,-2,3,lbz,ubz)

        call ComputeMidPointCoefficients(gridpoint_x,fd_order,lbx,ubx,x,mid_point_lele_x)
        call ComputeMidPointCoefficients(gridpoint_y,fd_order,lby,uby,y,mid_point_lele_y)
        call ComputeMidPointCoefficients(gridpoint_z,fd_order,lbz,ubz,z,mid_point_lele_z)

        call AllocateReal(xx       ,lbxG,ubxG,lbyG,ubyG)
        call AllocateReal(yy       ,lbxG,ubxG,lbyG,ubyG)
        call AllocateReal(dxi_dx   ,lbx ,ubx ,lby ,uby )
        call AllocateReal(dxi_dy   ,lbx ,ubx ,lby ,uby )
        call AllocateReal(deta_dx  ,lbx ,ubx ,lby ,uby )
        call AllocateReal(deta_dy  ,lbx ,ubx ,lby ,uby )
        call AllocateReal(Jac2D    ,lbx ,ubx ,lby ,uby )
        call AllocateReal(dx_dxi   ,lbx ,ubx ,lby ,uby )
        call AllocateReal(dx_deta  ,lbx ,ubx ,lby ,uby )
        call AllocateReal(dy_dxi   ,lbx ,ubx ,lby ,uby )
        call AllocateReal(dy_deta  ,lbx ,ubx ,lby ,uby )
        call AllocateReal(m_xi     ,lbx ,ubx ,lby ,uby )
        call AllocateReal(m_eta    ,lbx ,ubx ,lby ,uby )
        call AllocateReal(nx_xi    ,lbx ,ubx ,lby ,uby )
        call AllocateReal(ny_xi    ,lbx ,ubx ,lby ,uby )
        call AllocateReal(nx_eta   ,lbx ,ubx ,lby ,uby )
        call AllocateReal(ny_eta   ,lbx ,ubx ,lby ,uby )
        call AllocateReal(Gamma1_2D,lbx ,ubx ,lby ,uby )
        call AllocateReal(Gamma2_2D,lbx ,ubx ,lby ,uby )
        call AllocateReal(dxi_dxx  ,lbx ,ubx ,lby ,uby )
        call AllocateReal(dxi_dyy  ,lbx ,ubx ,lby ,uby )
        call AllocateReal(dxi_dxy  ,lbx ,ubx ,lby ,uby )
        call AllocateReal(deta_dxx ,lbx ,ubx ,lby ,uby )
        call AllocateReal(deta_dyy ,lbx ,ubx ,lby ,uby )
        call AllocateReal(deta_dxy ,lbx ,ubx ,lby ,uby )
        call AllocateReal(gJ11     ,lbx ,ubx ,lby ,uby )
        call AllocateReal(gJ22     ,lbx ,ubx ,lby ,uby )
        call AllocateReal(gJ12     ,lbx ,ubx ,lby ,uby )

        if(len_trim(mesh_file) > 0) then
          call read_curvilinear_mesh(xx,yy,lbxG,ubxG,lbyG,ubyG)
        else
          do j = lbyG,ubyG
             do i = lbxG,ubxG
                xx(i,j) = x(i)
                yy(i,j) = y(j)
             enddo
          enddo
        endif

        call compute_inplane_metrics(xx,yy,lbxG,ubxG,lbyG,ubyG)

        if(trim(diffusion_scheme) == 'conservative') call compute_face_metrics_conservative

        return
end subroutine init_grid

subroutine west_wall_normal_distance(wdist)
        use parameters_module, only: ny, ymin, cart_dims

        implicit none
        real(rp), dimension(0:ny), intent(out) :: wdist

        real(rp), dimension(0:ny) :: wloc
        integer , parameter       :: i0 = 1
        real(rp)                  :: dxj, dyj
        integer                   :: j, ierr

        if(cart_dims(2) /= 1) then
           if(rank == root) print*, ' ERROR: west_wall_normal_distance requires ', &
                                    'cart_dims(2)==1 (no wall-normal decomposition)'
           call secure_stop
        endif

        if(lbx <= i0 .and. i0 <= ubx) then
           wloc(0) = 0.0_rp
           wloc(1) = yy(i0,1) - ymin
           do j = 2, ny
              dxj     = xx(i0,j) - xx(i0,j-1)
              dyj     = yy(i0,j) - yy(i0,j-1)
              wloc(j) = wloc(j-1) + sqrt(dxj*dxj + dyj*dyj)
           enddo
        else
           wloc = -huge(1.0_rp)
        endif

        if(mpi_flag) then
           call MPI_allreduce(wloc, wdist, ny+1, MPI_RP, MPI_MAX, mpi_comm_cart, ierr)
        else
           wdist = wloc
        endif

        return
end subroutine west_wall_normal_distance

subroutine bottom_wall_streamwise_distance(sdist)
        use parameters_module, only: nx

        implicit none
        real(rp), dimension(-GN:nx+GN), intent(out) :: sdist

        real(rp), dimension(-GN:nx+GN) :: dloc, dseg
        integer , parameter            :: j0 = 1
        real(rp)                       :: dxi, dyi
        integer                        :: i, ierr

        dloc = -huge(1.0_rp)
        do i = max(1-GN, lbx-mesh_gh+1), min(nx+GN, ubx+mesh_gh)
           dxi     = xx(i,j0) - xx(i-1,j0)
           dyi     = yy(i,j0) - yy(i-1,j0)
           dloc(i) = sqrt(dxi*dxi + dyi*dyi)
        enddo

        if(mpi_flag) then
           call MPI_allreduce(dloc, dseg, nx+2*GN+1, MPI_RP, MPI_MAX, mpi_comm_cart, ierr)
        else
           dseg = dloc
        endif

        do i = 0, 1-GN, -1
           if(dseg(i) < -1.0_rp) dseg(i) = dseg(i+1)
        enddo
        do i = 1, nx+GN
           if(dseg(i) < -1.0_rp) dseg(i) = dseg(i-1)
        enddo

        sdist(0) = 0.0_rp
        do i = 1, nx+GN
           sdist(i) = sdist(i-1) + dseg(i)
        enddo
        do i = -1, -GN, -1
           sdist(i) = sdist(i+1) - dseg(i+1)
        enddo

        return
end subroutine bottom_wall_streamwise_distance

subroutine destroy_grid
        implicit none

        call DeallocateReal(x);         call DeallocateReal(y);         call DeallocateReal(z)

        call DeallocateReal(x_csi);     call DeallocateReal(y_eta);     call DeallocateReal(z_zta)
        call DeallocateReal(ix_csi);    call DeallocateReal(iy_eta);    call DeallocateReal(iz_zta)
        call DeallocateReal(x_csi2);    call DeallocateReal(y_eta2);    call DeallocateReal(z_zta2)

        call DeallocateReal(xstep);     call DeallocateReal(ystep);     call DeallocateReal(zstep)
        call DeallocateReal(xstep_i);   call DeallocateReal(ystep_i);   call DeallocateReal(zstep_i)
        call DeallocateReal(csistep_i); call DeallocateReal(etastep_i); call DeallocateReal(ztastep_i)
        call DeallocateReal(ixsteph);   call DeallocateReal(iysteph);   call DeallocateReal(izsteph)

        call DeallocateReal(mid_point_lele_x)
        call DeallocateReal(mid_point_lele_y)
        call DeallocateReal(mid_point_lele_z)

        call DeallocateReal(xx);        call DeallocateReal(yy)
        call DeallocateReal(dxi_dx);    call DeallocateReal(dxi_dy)
        call DeallocateReal(deta_dx);   call DeallocateReal(deta_dy)
        call DeallocateReal(Jac2D)
        call DeallocateReal(dx_dxi);    call DeallocateReal(dx_deta)
        call DeallocateReal(dy_dxi);    call DeallocateReal(dy_deta)
        call DeallocateReal(m_xi);      call DeallocateReal(m_eta)
        call DeallocateReal(nx_xi);     call DeallocateReal(ny_xi)
        call DeallocateReal(nx_eta);    call DeallocateReal(ny_eta)
        call DeallocateReal(Gamma1_2D); call DeallocateReal(Gamma2_2D)
        call DeallocateReal(dxi_dxx);   call DeallocateReal(dxi_dyy);  call DeallocateReal(dxi_dxy)
        call DeallocateReal(deta_dxx);  call DeallocateReal(deta_dyy); call DeallocateReal(deta_dxy)
        call DeallocateReal(gJ11);      call DeallocateReal(gJ22);      call DeallocateReal(gJ12)

        return
end subroutine destroy_grid

subroutine read_curvilinear_mesh(xx,yy,lbxG,ubxG,lbyG,ubyG)
        use parameters_module, only: nx,ny,mesh_file

        implicit none
        integer                                   , intent(in)    :: lbxG,ubxG,lbyG,ubyG
        real(rp), dimension(lbxG:ubxG,lbyG:ubyG)  , intent(inout) :: xx, yy

        real(rp), allocatable, dimension(:,:) :: tmp_x, tmp_y
        integer  :: meshUnit = 53, err = 0, errFlag, ng, i, j
        integer  :: nxh, nyh, nzh

        ng      = GN + mesh_gh
        errFlag = 0

        call AllocateReal(tmp_x,1-ng,nx+ng,1-ng,ny+ng)
        call AllocateReal(tmp_y,1-ng,nx+ng,1-ng,ny+ng)

        if(rank == root) then

          open(unit = meshUnit, file = trim(mesh_file), status = 'old', iostat = err)
          if(err /= 0) then
            print*, ' Curvilinear mesh file ', trim(mesh_file), ' has not been opened!'
            errFlag = 1
          endif

          if(errFlag == 0) then
            read(meshUnit,*,iostat=err) nxh, nyh, nzh
            if(err /= 0 .or. nxh /= nx .or. nyh /= ny) then
              print*, ' --------------------------------------------------------- '
              print*, ' Mesh header     : ', nxh, nyh
              print*, ' Input grid (nx,ny): ', nx, ny
              print*, ' The 2D mesh file is not consistent with the input grid'
              print*, ' --------------------------------------------------------- '
              errFlag = 2
            endif
          endif

          if(errFlag == 0) then
            do j = 1-ng,ny+ng
               do i = 1-ng,nx+ng
                  read(meshUnit,*) tmp_x(i,j)
               enddo
            enddo
            do j = 1-ng,ny+ng
               do i = 1-ng,nx+ng
                  read(meshUnit,*) tmp_y(i,j)
               enddo
            enddo
            close(meshUnit)
          endif

        endif

        if(mpi_flag) call MPI_BCAST(errFlag, 1, MPI_INTEGER, root, mpi_comm_world, err)
        if(errFlag /= 0) call secure_stop

        if(mpi_flag) then
          call MPI_BCAST(tmp_x, size(tmp_x), MPI_RP, root, mpi_comm_world, err)
          call MPI_BCAST(tmp_y, size(tmp_y), MPI_RP, root, mpi_comm_world, err)
        endif

        do j = lbyG,ubyG
           do i = lbxG,ubxG
              xx(i,j) = tmp_x(i,j)
              yy(i,j) = tmp_y(i,j)
           enddo
        enddo

        call write_inplane_mesh(tmp_x,tmp_y,ng)

        call DeallocateReal(tmp_x)
        call DeallocateReal(tmp_y)

        return
end subroutine read_curvilinear_mesh

subroutine compute_inplane_metrics(xx,yy,lbxG,ubxG,lbyG,ubyG)
        use parameters_module, only: fd_order, central_1, central_2

        implicit none
        integer                                  , intent(in) :: lbxG,ubxG,lbyG,ubyG
        real(rp), dimension(lbxG:ubxG,lbyG:ubyG) , intent(in) :: xx, yy

        real(rp), parameter :: toll = 1.0E-12_rp
        integer  :: i, j, l, m, errFlag, errFlagG, err = 0
        logical  :: lcl_skew, lcl_axis
        real(rp) :: x_xi, x_eta, y_xi, y_eta, detJ2D, idetJ2D, m_nrm
        real(rp) :: x_xixi, x_etaeta, x_xieta, y_xixi, y_etaeta, y_xieta
        real(rp) :: g11p, g22p, g12p, lapl_x, lapl_y
        real(rp) :: csi_x, csi_y, eta_x, eta_y
        real(rp) :: Bx_xx, By_xx, Bx_yy, By_yy, Bx_xy, By_xy

        errFlag = 0

        do j = lby,uby
           do i = lbx,ubx

              x_xi = 0.0_rp; x_eta = 0.0_rp; y_xi = 0.0_rp; y_eta = 0.0_rp
              do l = 1,3
                 x_xi  = x_xi  + central_1(l)*(xx(i+l,j) - xx(i-l,j))
                 y_xi  = y_xi  + central_1(l)*(yy(i+l,j) - yy(i-l,j))
                 x_eta = x_eta + central_1(l)*(xx(i,j+l) - xx(i,j-l))
                 y_eta = y_eta + central_1(l)*(yy(i,j+l) - yy(i,j-l))
              enddo

              detJ2D = x_xi*y_eta - x_eta*y_xi
              if(detJ2D <= toll) errFlag = 1
              idetJ2D = 1.0_rp/detJ2D

              dxi_dx (i,j) =  y_eta*idetJ2D
              dxi_dy (i,j) = -x_eta*idetJ2D
              deta_dx(i,j) = -y_xi *idetJ2D
              deta_dy(i,j) =  x_xi *idetJ2D
              Jac2D  (i,j) =  idetJ2D

              dx_dxi (i,j) =  x_xi
              dx_deta(i,j) =  x_eta
              dy_dxi (i,j) =  y_xi
              dy_deta(i,j) =  y_eta

              m_nrm      = sqrt(dxi_dx(i,j)**2 + dxi_dy(i,j)**2)
              m_xi (i,j) = m_nrm
              nx_xi(i,j) = dxi_dx(i,j)/m_nrm
              ny_xi(i,j) = dxi_dy(i,j)/m_nrm

              m_nrm       = sqrt(deta_dx(i,j)**2 + deta_dy(i,j)**2)
              m_eta (i,j) = m_nrm
              nx_eta(i,j) = deta_dx(i,j)/m_nrm
              ny_eta(i,j) = deta_dy(i,j)/m_nrm

              x_xixi = central_2(0)*xx(i,j); x_etaeta = central_2(0)*xx(i,j)
              y_xixi = central_2(0)*yy(i,j); y_etaeta = central_2(0)*yy(i,j)
              do l = 1,3
                 x_xixi   = x_xixi   + central_2(l)*(xx(i+l,j) + xx(i-l,j))
                 x_etaeta = x_etaeta + central_2(l)*(xx(i,j+l) + xx(i,j-l))
                 y_xixi   = y_xixi   + central_2(l)*(yy(i+l,j) + yy(i-l,j))
                 y_etaeta = y_etaeta + central_2(l)*(yy(i,j+l) + yy(i,j-l))
              enddo

              x_xieta = 0.0_rp; y_xieta = 0.0_rp
              do m = 1,3
                 do l = 1,3
                    x_xieta = x_xieta + central_1(l)*central_1(m)* &
                              (xx(i+l,j+m) - xx(i-l,j+m) - xx(i+l,j-m) + xx(i-l,j-m))
                    y_xieta = y_xieta + central_1(l)*central_1(m)* &
                              (yy(i+l,j+m) - yy(i-l,j+m) - yy(i+l,j-m) + yy(i-l,j-m))
                 enddo
              enddo

              g11p = dxi_dx(i,j)**2  + dxi_dy(i,j)**2
              g22p = deta_dx(i,j)**2 + deta_dy(i,j)**2
              g12p = dxi_dx(i,j)*deta_dx(i,j) + dxi_dy(i,j)*deta_dy(i,j)

              lapl_x = g11p*x_xixi + 2.0_rp*g12p*x_xieta + g22p*x_etaeta
              lapl_y = g11p*y_xixi + 2.0_rp*g12p*y_xieta + g22p*y_etaeta

              Gamma1_2D(i,j) = -(dxi_dx (i,j)*lapl_x + dxi_dy (i,j)*lapl_y)
              Gamma2_2D(i,j) = -(deta_dx(i,j)*lapl_x + deta_dy(i,j)*lapl_y)

              csi_x = dxi_dx (i,j);  csi_y = dxi_dy (i,j)
              eta_x = deta_dx(i,j);  eta_y = deta_dy(i,j)

              Bx_xx = x_xixi*csi_x*csi_x + 2.0_rp*x_xieta*csi_x*eta_x + x_etaeta*eta_x*eta_x
              By_xx = y_xixi*csi_x*csi_x + 2.0_rp*y_xieta*csi_x*eta_x + y_etaeta*eta_x*eta_x
              Bx_yy = x_xixi*csi_y*csi_y + 2.0_rp*x_xieta*csi_y*eta_y + x_etaeta*eta_y*eta_y
              By_yy = y_xixi*csi_y*csi_y + 2.0_rp*y_xieta*csi_y*eta_y + y_etaeta*eta_y*eta_y
              Bx_xy = x_xixi*csi_x*csi_y + x_xieta*(csi_x*eta_y + eta_x*csi_y) + x_etaeta*eta_x*eta_y
              By_xy = y_xixi*csi_x*csi_y + y_xieta*(csi_x*eta_y + eta_x*csi_y) + y_etaeta*eta_x*eta_y

              dxi_dxx (i,j) = -(csi_x*Bx_xx + csi_y*By_xx)
              dxi_dyy (i,j) = -(csi_x*Bx_yy + csi_y*By_yy)
              dxi_dxy (i,j) = -(csi_x*Bx_xy + csi_y*By_xy)
              deta_dxx(i,j) = -(eta_x*Bx_xx + eta_y*By_xx)
              deta_dyy(i,j) = -(eta_x*Bx_yy + eta_y*By_yy)
              deta_dxy(i,j) = -(eta_x*Bx_xy + eta_y*By_xy)

              gJ11(i,j) = g11p*detJ2D
              gJ22(i,j) = g22p*detJ2D
              gJ12(i,j) = g12p*detJ2D

           enddo
        enddo

        if(mpi_flag) then
          call MPI_allreduce(errFlag, errFlagG, 1, MPI_INTEGER, MPI_MAX, mpi_comm_cart, err)
        else
          errFlagG = errFlag
        endif
        if(errFlagG /= 0) then
          if(rank == root) print*, ' ERROR: singular or left-handed in-plane mesh (detJ2D <=', toll, ')'
          call secure_stop
        endif

        lcl_skew = any(gJ12 /= 0.0_rp)
        if(mpi_flag) then
          call MPI_allreduce(lcl_skew, non_orthogonal, 1, MPI_LOGICAL, MPI_LOR, mpi_comm_cart, err)
        else
          non_orthogonal = lcl_skew
        endif
        if(rank == root) print '(a,l2)', '  Grid non-orthogonal (g12 cross viscous flux active) : ', non_orthogonal

        lcl_axis = all(dxi_dy == 0.0_rp) .and. all(deta_dx == 0.0_rp)
        if(mpi_flag) then
          call MPI_allreduce(lcl_axis, inplane_axis_aligned, 1, MPI_LOGICAL, MPI_LAND, mpi_comm_cart, err)
        else
          inplane_axis_aligned = lcl_axis
        endif
        if(rank == root) print '(a,l2)', '  Grid in-plane axis-aligned (grad(div v) Cartesian fast-path) : ', inplane_axis_aligned

        return
end subroutine compute_inplane_metrics

subroutine compute_face_metrics_conservative
        implicit none

        real(rp) :: cm
        integer  :: i, j, p, m, pm, pp, plo, phi

        call AllocateReal(fmx_csix, lbx,ubx, lby,uby)
        call AllocateReal(fmx_csiy, lbx,ubx, lby,uby)
        call AllocateReal(fmx_etax, lbx,ubx, lby,uby)
        call AllocateReal(fmx_etay, lbx,ubx, lby,uby)
        call AllocateReal(fmx_cof1, lbx,ubx, lby,uby)
        call AllocateReal(fmx_cof2, lbx,ubx, lby,uby)
        call AllocateReal(fmy_csix, lbx,ubx, lby,uby)
        call AllocateReal(fmy_csiy, lbx,ubx, lby,uby)
        call AllocateReal(fmy_etax, lbx,ubx, lby,uby)
        call AllocateReal(fmy_etay, lbx,ubx, lby,uby)
        call AllocateReal(fmy_cof1, lbx,ubx, lby,uby)
        call AllocateReal(fmy_cof2, lbx,ubx, lby,uby)
        call AllocateReal(fmz_zeta, lbz,ubz)

        plo = lbx + visc_face_reach - 1
        phi = ubx - visc_face_reach
        do   j = lby,uby
          do p = plo,phi
            do m = 1,visc_face_reach
              cm = visc_intp(m)
              pm = p+1-m
              pp = p+m
              fmx_csix(p,j) = fmx_csix(p,j) + cm*(dxi_dx (pp,j) + dxi_dx (pm,j))
              fmx_csiy(p,j) = fmx_csiy(p,j) + cm*(dxi_dy (pp,j) + dxi_dy (pm,j))
              fmx_etax(p,j) = fmx_etax(p,j) + cm*(deta_dx(pp,j) + deta_dx(pm,j))
              fmx_etay(p,j) = fmx_etay(p,j) + cm*(deta_dy(pp,j) + deta_dy(pm,j))
              fmx_cof1(p,j) = fmx_cof1(p,j) + cm*(dy_deta(pp,j) + dy_deta(pm,j))
              fmx_cof2(p,j) = fmx_cof2(p,j) - cm*(dx_deta(pp,j) + dx_deta(pm,j))
            enddo
          enddo
        enddo

        plo = lby + visc_face_reach - 1
        phi = uby - visc_face_reach
        do   p = plo,phi
          do i = lbx,ubx
            do m = 1,visc_face_reach
              cm = visc_intp(m)
              pm = p+1-m
              pp = p+m
              fmy_csix(i,p) = fmy_csix(i,p) + cm*(dxi_dx (i,pp) + dxi_dx (i,pm))
              fmy_csiy(i,p) = fmy_csiy(i,p) + cm*(dxi_dy (i,pp) + dxi_dy (i,pm))
              fmy_etax(i,p) = fmy_etax(i,p) + cm*(deta_dx(i,pp) + deta_dx(i,pm))
              fmy_etay(i,p) = fmy_etay(i,p) + cm*(deta_dy(i,pp) + deta_dy(i,pm))
              fmy_cof1(i,p) = fmy_cof1(i,p) - cm*(dy_dxi (i,pp) + dy_dxi (i,pm))
              fmy_cof2(i,p) = fmy_cof2(i,p) + cm*(dx_dxi (i,pp) + dx_dxi (i,pm))
            enddo
          enddo
        enddo

        plo = lbz + visc_face_reach - 1
        phi = ubz - visc_face_reach
        do p = plo,phi
          do m = 1,visc_face_reach
            cm = visc_intp(m)
            pm = p+1-m
            pp = p+m
            fmz_zeta(p) = fmz_zeta(p) + cm*(ztastep_i(pp)*iz_zta(pp) + ztastep_i(pm)*iz_zta(pm))
          enddo
        enddo

        return
end subroutine compute_face_metrics_conservative

subroutine compute_grid_point(x,nx,lb,ub,xmin,xmax,gridType,dir)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x
        integer                            , intent(in)    :: nx, lb, ub
        real(rp)                           , intent(in)    :: xmin, xmax
        character(*)                       , intent(in)    :: gridType
        character(1)             , optional, intent(in)    :: dir

        selectcase(trim(gridType))
          case('uniform')
          call uniformGrid_f0(x,nx,lb,ub,xmin,xmax)
          if(present(dir)) call write_the_grid(nx,xmin,xmax,dir)

          case default
          call ReadMeshFile(x,nx,lb,ub,xmin,xmax,gridType,dir)

        end select

        return
end subroutine compute_grid_point

subroutine compute_grid_diff1(x,x_csi,nx,lb,ub,xmin,xmax,gridType)
        implicit none
        real(rp), dimension(:), allocatable, intent(in)    :: x
        real(rp), dimension(:), allocatable, intent(inout) :: x_csi
        integer                            , intent(in)    :: nx, lb, ub
        real(rp)                           , intent(in)    :: xmin, xmax
        character(*)                       , intent(in)    :: gridType

        selectcase(trim(gridType))
          case('uniform')
          call uniformGrid_f1(x_csi,lb,ub,xmin,xmax)

          case default
          call computeGridF1(nx,lb,ub,x,x_csi)

        end select

        return
end subroutine compute_grid_diff1

subroutine compute_grid_diff2(x,x_csi2,nx,lb,ub,gridType)
        implicit none
        real(rp), dimension(:), allocatable, intent(in)    :: x
        real(rp), dimension(:), allocatable, intent(inout) :: x_csi2
        integer                            , intent(in)    :: nx, lb, ub
        character(*)                       , intent(in)    :: gridType

        selectcase(trim(gridType))
          case('uniform')
          call uniformGrid_f2(x_csi2,lb,ub)

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
        integer                            , intent(in)    :: nx, lb, ub

        real(rp), dimension(:), allocatable :: midC
        real(rp) :: delta_csi, xsteph_
        integer  :: i, l, fdL, fdR

        do i = lb,ub

           delta_csi = 1.0_rp/real(nx,rp)

           xstep(i) = delta_csi * x_csi(i)

           xstep_i(i) = 1.0_rp/xstep(i)

           csistep_i(i) = 1.0_rp/delta_csi

           ix_csi(i) = 1.0_rp/x_csi(i)
        enddo

        selectcase(fd_order)

          case(2)
            fdL = 0
            fdR = 1
            call AllocateReal(midC,fdL,fdR)
            midC(0) = 0.5_rp
            midC(1) = midC(0)

          case(4)
            fdL = -1
            fdR =  2
            call AllocateReal(midC,fdL,fdR)
            midC(-1) = -1.0_rp/16.0_rp
            midC( 0) = +9.0_rp/16.0_rp
            midC( 1) = midC( 0)
            midC( 2) = midC(-1)

          case(6)
            fdL = -2
            fdR =  3
            call AllocateReal(midC,fdL,fdR)
            midC(-2) =  3.0_rp /256.0_rp
            midC(-1) = -25.0_rp/256.0_rp
            midC( 0) =  75.0_rp/128.0_rp
            midC( 1) = midC( 0)
            midC( 2) = midC(-1)
            midC( 3) = midC(-2)

          case default
          if(rank == root) print*, ' ERROR in compute_grid_steps: fd_order =', fd_order, ' not implemented'
          call secure_stop

        end select

        do i = lb + 3, ub - 3
           xsteph_ = 0.0_rp
           do l = fdL, fdR
              xsteph_ = xsteph_ + midC(l)*xstep(i+l)
           enddo
           ixsteph(i) = 1.0_rp/xsteph_
        enddo

        call DeallocateReal(midC)

        return
end subroutine compute_grid_steps

subroutine ReadMeshFile(x,nx,lb,ub,xmin,xmax,gridType,dir)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x
        real(rp)                           , intent(in)    :: xmin, xmax
        integer                            , intent(in)    :: nx, lb, ub
        character(*)                       , intent(in)    :: gridType
        character(1)             , optional, intent(in)    :: dir

        real(rp), dimension(:), allocatable :: tmp
        real(rp), parameter                 :: toll = 10.0_rp
        integer       :: meshUnit = 51, err = 0, ngGrid, i, pts, errFlag
        real(rp)      :: rdummy, strGrid, endGrid
        integer       :: idummy

        ngGrid  = GN + mesh_gh
        errFlag = 0

        call AllocateReal(tmp,1-ngGrid,nx+ngGrid)

        if(rank == root) then

          open(unit = meshUnit, file = trim(gridType), status = 'old', iostat = err)
          if(err /= 0) then
            print*, ' Mesh file ', trim(gridType), ' has not been opened!'
            errFlag = 1
          endif

          if(errFlag == 0) then
            pts = 0
            do
              read(meshUnit,*,iostat=err) idummy, rdummy
              if(err /= 0) exit
              pts = pts + 1
            enddo
            rewind(meshUnit)

            if(pts - 2*ngGrid /= nx) then
              print*, ' --------------------------------------------------------- '
              print*, ' Number of input points: ', nx
              print*, ' Number of grid  points: ', pts - 2*ngGrid
              print*, ' The grid file is not consistent with the input grid points'
              print*, ' --------------------------------------------------------- '
              errFlag = 2
            endif
          endif

          if(errFlag == 0) then
            do i = 1-ngGrid,nx+ngGrid
               read(meshUnit,*) idummy, tmp(i)
            enddo
            close(meshUnit)

            strGrid = 0.5_rp*(tmp(0) + tmp(1))
            endGrid = 0.5_rp*(tmp(nx)+ tmp(nx+1))
            if(abs(strGrid-xmin) > toll .or. abs(endGrid-xmax) > toll) then
              print*, '--------------------------------------- '
              print*, ' The grid   is :', strGrid, ' x ', endGrid
              print*, ' The domain is :', xmin   , ' x ', xmax
              print*, ' The grid is not consistent with inputs!'
              print*, '--------------------------------------- '
              errFlag = 3
            endif
          endif

        endif

        if(mpi_flag) call MPI_BCAST(errFlag, 1, MPI_INTEGER, root, mpi_comm_world, err)
        if(errFlag /= 0) call secure_stop

        if(mpi_flag) call MPI_BCAST(tmp, size(tmp), MPI_RP, root, mpi_comm_world, err)

        do i = lb,ub
           x(i) = tmp(i)
        enddo

        if(present(dir)) call write_the_grid(nx,xmin,xmax,dir,tmp)

        call DeallocateReal(tmp)

        return
end subroutine ReadMeshFile

subroutine write_the_grid(npt,xmin,xmax,dir,grid_in)
        implicit none

        real(rp)    , intent(in)           :: xmax, xmin
        integer     , intent(in)           :: npt
        character(1), intent(in)           :: dir
        real(rp), dimension(:), optional, intent(in) :: grid_in

        real(rp), dimension(:), allocatable :: gbl_arr
        real(rp) :: xi, Lx
        logical  :: GridExist
        character(50) :: filename
        integer  :: ngGrid = 7, err
        integer  :: i, MeshUnit = 10

        allocate(gbl_arr(1-ngGrid:npt+ngGrid))

        if(present(grid_in)) then
          gbl_arr(1-ngGrid:npt+ngGrid) = grid_in
        else
          Lx = xmax - xmin
          do i = 1-ngGrid,npt+ngGrid
             xi   = (i-0.5_rp)/real(npt,rp)
             gbl_arr(i) = xmin + Lx*xi
          enddo
        endif

        selectcase(dir)
          case('x','X'); filename = 'Xgrid.txt'
          case('y','Y'); filename = 'Ygrid.txt'
          case('z','Z'); filename = 'Zgrid.txt'
          case default
            if(rank == root) print*, 'WARNING: write_the_grid called with unknown dir=', dir
            deallocate(gbl_arr); return
        end select

        if(rank == root) then
          inquire(file = 'DATA/'//trim(data_dir)//'/GRID/'//trim(filename), exist= GridExist)
          if(.not.GridExist) then

            call execute_command_line('mkdir -p DATA/'//trim(data_dir)//'/GRID', wait=.true.)

            open(unit = meshUnit, status='replace', action='write', &
                file = 'DATA/'//trim(data_dir)//'/GRID/'//trim(filename), iostat = err)
            if(err .ne. 0) then
              print*, 'WARNING: Mesh file ', trim(filename), ' has not been saved! (iostat=', err, ')'
            else
              do i = 1-ngGrid,npt+ngGrid
                 write(meshUnit,'(i8,1x,es25.16e3)',iostat=err) i, gbl_arr(i)
                 if(err /= 0) exit
              enddo
              close(MeshUnit)
              if(err == 0) then
                print*, ' The grid file ', trim(filename), ' has been saved in DATA/',trim(data_dir),'/GRID/'
              else
                print*, 'WARNING: write failure on ', trim(filename), ' (iostat=', err, ')'
              endif
            endif
          endif
        endif

        deallocate(gbl_arr)

        return
end subroutine write_the_grid

subroutine write_inplane_mesh(mx,my,ng)
        use parameters_module, only: nx, ny, nz, data_dir

        implicit none
        integer                                       , intent(in) :: ng
        real(rp), dimension(1-ng:nx+ng,1-ng:ny+ng)    , intent(in) :: mx, my

        logical :: GridExist
        integer :: i, j, err, MeshUnit = 11

        if(rank /= root) return

        inquire(file = 'DATA/'//trim(data_dir)//'/GRID/MeshXY.txt', exist = GridExist)
        if(GridExist) return

        call execute_command_line('mkdir -p DATA/'//trim(data_dir)//'/GRID', wait=.true.)

        open(unit = MeshUnit, status='replace', action='write', &
             file = 'DATA/'//trim(data_dir)//'/GRID/MeshXY.txt', iostat = err)
        if(err /= 0) then
          print*, 'WARNING: MeshXY.txt has not been saved! (iostat=', err, ')'
          return
        endif

        write(MeshUnit,'(3(i0,1x))') nx, ny, nz
        do j = 1-ng,ny+ng
           do i = 1-ng,nx+ng
              write(MeshUnit,'(es25.16e3)') mx(i,j)
           enddo
        enddo
        do j = 1-ng,ny+ng
           do i = 1-ng,nx+ng
              write(MeshUnit,'(es25.16e3)') my(i,j)
           enddo
        enddo
        close(MeshUnit)

        print*, ' The curvilinear mesh MeshXY.txt has been saved in DATA/',trim(data_dir),'/GRID/'

        return
end subroutine write_inplane_mesh

subroutine computeGridF1(nx,lb,ub,x,x_csi)
        use parameters_module, only: fd_order

        implicit none
        real(rp), allocatable, dimension(:), intent(in)    :: x
        real(rp), allocatable, dimension(:), intent(inout) :: x_csi
        integer                            , intent(in)    :: nx, lb, ub

        real(rp), allocatable, dimension(:) :: central_1
        real(rp)                            :: idxi, df
        integer                             :: i, s, fdO

        fdO = int(fd_order/2)
        call AllocateReal(central_1,-fdO,fdO)

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
        case default
          if(rank == root) print*, ' ERROR in computeGridF1: fd_order =', fd_order, ' not implemented'
          call secure_stop
        endselect
        if(abs(sum(central_1)) > 1.0E-14_rp) then
          if(rank == root) print*, ' ERROR in computeGridF1: unbalanced derivative coefficients'
          call secure_stop
        endif

        idxi = real(nx,rp)
        do i = lb,ub

           df = 0.0_rp
           do s = 1, fdO
              df = df + central_1(s) * (x(i+s) - x(i-s))
           enddo

           x_csi(i) = idxi*df

        enddo

        call DeallocateReal(central_1)

        return
end subroutine computeGridF1

subroutine computeGridF2(nx,lb,ub,x,x_csi2)
        use parameters_module, only: fd_order

        implicit none
        real(rp), allocatable, dimension(:), intent(in)    :: x
        real(rp), allocatable, dimension(:), intent(inout) :: x_csi2
        integer                            , intent(in)    :: nx, lb, ub

        real(rp), allocatable, dimension(:) :: central_2
        integer                             :: i, s, fdO
        real(rp)                            :: idxi2, d2f

        fdO = int(fd_order/2)
        call AllocateReal(central_2,-fdO,fdO)

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
        case default
          if(rank == root) print*, ' ERROR in computeGridF2: fd_order =', fd_order, ' not implemented'
          call secure_stop
        endselect
        if(abs(sum(central_2)) > 1.0E-14_rp) then
          if(rank == root) print*, ' ERROR in computeGridF2: unbalanced derivative coefficients'
          call secure_stop
        endif

        idxi2 = (real(nx,rp))**2
        do i = lb,ub

           d2f = central_2(0) * x(i)
           do s = 1, fdO
              d2f = d2f + central_2(s) * (x(i+s) + x(i-s))
           enddo

           x_csi2(i) = idxi2*d2f

        enddo

        call DeallocateReal(central_2)

        return
end subroutine computeGridF2

subroutine ComputeMidPointCoefficients(gridpoint,fd_order,lbx,ubx,x,mid_point_coef)
        use matrix_inversion_module

        implicit none
        real(rp), allocatable, dimension(:)  , intent(in)    :: x
        integer                              , intent(in)    :: fd_order
        integer                              , intent(in)    :: lbx, ubx
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
          do i = lbx,ubx-1

             xx = 0.5_rp*(x(i) + x(i+1))
             do l = -mm+1,mm
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
             do l = -mm+1,mm
                mid_point_coef(l,i) = c(l)
             enddo

             if(abs(sum(mid_point_coef(:,i))-1.0_rp)>toll) then
               print*, abs(sum(mid_point_coef(:,i)) - 1.0_rp)
               stop 'Mid Point Coefficients are umbalanced'
             endif

          enddo

        endif

        return
end subroutine ComputeMidPointCoefficients

subroutine uniformGrid_f0(x,nx,lb,ub,xmin,xmax)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x
        integer                            , intent(in)    :: nx, lb, ub
        real(rp)                           , intent(in)    :: xmin, xmax

        real(rp) :: Lx, xi
        integer  :: i

        Lx = xmax - xmin
        do i = lb, ub

           xi   = (i-0.5_rp)/real(nx,rp)
           x(i) = xmin + Lx*xi

        enddo

        return
end subroutine uniformGrid_f0

subroutine uniformGrid_f1(x_csi,lb,ub,xmin,xmax)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x_csi
        integer                            , intent(in)    :: lb, ub
        real(rp)                           , intent(in)    :: xmin, xmax

        real(rp) :: Lx

        Lx = xmax - xmin
        x_csi(lb:ub) = Lx

        return
end subroutine uniformGrid_f1

subroutine uniformGrid_f2(x_csi2,lb,ub)
        implicit none
        real(rp), dimension(:), allocatable, intent(inout) :: x_csi2
        integer                            , intent(in)    :: lb, ub

        x_csi2(lb:ub) = 0.0_rp

        return
end subroutine uniformGrid_f2

end module mesh_module
