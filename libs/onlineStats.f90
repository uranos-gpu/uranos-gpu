module onlineStats

use mpi
use statistics_module

implicit none

contains


subroutine compute_1DStats(phi,VIS,nx,nz,nVar,itStat,&
                          sx,ex,sy,ey,sz,ez,comm,vmean,mpi_flag)

        use parameters_module     , only: rp, gamma0
        use fluid_functions_module, only: Sutherland
        use storage_module        , only: central_1, central_fd_order
        use mesh_module           , only: ystep_i
        use mpi_module            , only: lby,uby

        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(in)    :: phi
        real(rp), dimension(:,:,:)  , allocatable, intent(in)    :: VIS
        real(rp), dimension(:,:)    , allocatable, intent(inout) :: vmean

        integer , intent(in) :: comm
        integer , intent(in) :: nx,nz,nVar,sx,ex,sy,ey,sz,ez, itStat
        logical , intent(in) :: mpi_flag

        ! local declarations
        integer, parameter :: prec = MPI_DOUBLE_PRECISION
        integer            :: err  = 0, npt

        real(rp), dimension(lby:uby,nVar) :: tmp_lcl
        real(rp), dimension(lby:uby,nVar) :: tmp_gbl
        real(rp) :: r_, u_, v_, w_, ek, p_, T_, re, ir, mL, mT
        real(rp) :: cl, r_l, rul, u_l, r_y, u_y, ruy
        real(rp) :: irNxNz, r_itStat, r_itStatm1
        integer  :: i,j,k, l

        tmp_lcl = 0.0_rp
        tmp_gbl = 0.0_rp
        do       k = sz,ez
           do    j = sy,ey
              do i = sx,ex 

                 ! compute primitives
                 r_ = phi(i,j,k,1)
                 ir = 1.0_rp/r_
                 u_ = phi(i,j,k,2)*ir
                 v_ = phi(i,j,k,3)*ir
                 w_ = phi(i,j,k,4)*ir
                 ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
                 re = phi(i,j,k,5)
                 p_ = (gamma0-1.0_rp)*(re - r_*ek)
                 T_ = p_*ir

                 mL = Sutherland(T_)
                 mT = VIS(i,j,k) - mL

                 ! conservative variables
                 tmp_lcl(j,1) = tmp_lcl(j,1) + r_
                 tmp_lcl(j,2) = tmp_lcl(j,2) + r_*u_
                 tmp_lcl(j,3) = tmp_lcl(j,3) + r_*v_
                 tmp_lcl(j,4) = tmp_lcl(j,4) + r_*w_
                 tmp_lcl(j,5) = tmp_lcl(j,5) + re

                 ! Reynolds stress pieces
                 tmp_lcl(j,6) = tmp_lcl(j,6) + (r_*u_)*(r_*u_)*ir
                 tmp_lcl(j,7) = tmp_lcl(j,7) + (r_*v_)*(r_*v_)*ir
                 tmp_lcl(j,8) = tmp_lcl(j,8) + (r_*w_)*(r_*w_)*ir
                 tmp_lcl(j,9) = tmp_lcl(j,9) + (r_*u_)*(r_*v_)*ir

                 ! density 
                 tmp_lcl(j,10) = tmp_lcl(j,10) + r_*r_

                 !pressure
                 tmp_lcl(j,11) = tmp_lcl(j,11) + p_
                 tmp_lcl(j,12) = tmp_lcl(j,12) + p_*p_

                 ! temperature
                 tmp_lcl(j,13) = tmp_lcl(j,13) + T_
                 tmp_lcl(j,14) = tmp_lcl(j,14) + r_*T_*T_

                 ! viscosity
                 tmp_lcl(j,15) = tmp_lcl(j,15) + mL
                 tmp_lcl(j,16) = tmp_lcl(j,16) + mT
        
                 ! derivatives
                 r_y = 0.0_rp
                 ruy = 0.0_rp
                 u_y = 0.0_rp
                 do l = - central_fd_order/2,central_fd_order/2

                    cl   = central_1(l)

                    r_l = phi(i,j+l,k,1)
                    rul = phi(i,j+l,k,2)
                    u_l = rul/r_l

                    r_y = r_y + cl * r_l
                    u_y = u_y + cl * u_l
                    ruy = ruy + cl * rul

                 enddo
                 r_y = r_y*ystep_i(j)
                 u_y = u_y*ystep_i(j)
                 ruy = ruy*ystep_i(j)

                 tmp_lcl(j,17) = tmp_lcl(j,17) + r_y
                 tmp_lcl(j,18) = tmp_lcl(j,18) + ruy

                 tmp_lcl(j,19) = tmp_lcl(j,19) + mL*u_y
                 tmp_lcl(j,20) = tmp_lcl(j,20) + mT*u_y

              enddo
           enddo
        enddo

        irNxNz = 1.0_rp/real(nx*nz,rp)

        if(mpi_flag) then
          npt = nVar*(uby-lby+1)
          call mpi_allreduce(tmp_lcl,tmp_gbl,npt,prec,mpi_sum,comm,err)
          tmp_gbl = tmp_gbl*irNxNz
        else
          tmp_gbl = tmp_lcl*irNxNz
        endif
        
        ! summing up statistics
        r_itStat   = real(itStat,rp)
        r_itStatm1 = real(itStat-1,rp)
        do    l = 1,nVar
           do j = sy,ey
              vmean(j,l) = (vmean(j,l)*r_itStatm1 + tmp_gbl(j,l))/r_itStat
           enddo
        enddo

        return
end subroutine compute_1DStats




subroutine compute_2DStats(phi,VIS,nz,nVar,itStat, &
                           sx,ex,sy,ey,sz,ez,comm,vmean,mpi_flag)
        
        use parameters_module     , only: rp, gamma0, central_fd_order, hybrid_weno
        use fluid_functions_module, only: Sutherland
        use mpi_module            , only: lbx,ubx, lby,uby
        use storage_module        , only: central_1, vmean2D_aux, SSENSOR, weno
        use mesh_module           , only: xstep_i, ystep_i, zstep_i

        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(in)    :: phi
        real(rp), dimension(:,:,:)  , allocatable, intent(in)    :: VIS
        real(rp), dimension(:,:,:)  , allocatable, intent(inout) :: vmean
        
        integer , intent(in) :: comm, itStat
        integer , intent(in) :: nz,nVar,sx,ex,sy,ey,sz,ez
        logical , intent(in) :: mpi_flag

        ! local declarations
        integer, parameter :: prec = MPI_DOUBLE_PRECISION
        integer            :: err  = 0, npt

        real(rp), dimension(lbx:ubx,lby:uby,nVar) :: tmp_lcl
        real(rp), dimension(lbx:ubx,lby:uby,nVar) :: tmp_gbl
        real(rp) :: r_, u_, v_, w_, ek, p_, T_, re, ir, mL, mT, ru, rv, rw
        real(rp) :: irNz, r_itStat, r_itStatm1, ve, c_
        integer  :: i,j,k, l
        
        real(rp) :: cl, idx, idy, idz
        real(rp) :: r_x, r_y, r_z
        real(rp) :: rux, ruy, ruz
        real(rp) :: rvx, rvy, rvz
        real(rp) :: rwx, rwy, rwz
        real(rp) :: omx, omy, omz
        real(rp) :: sensor, weno_flag

        tmp_lcl = 0.0_rp
        tmp_gbl = 0.0_rp
        do       k = sz,ez
           do    j = sy,ey
              do i = sx,ex
        
                 ! compute primitives
                 r_ = phi(i,j,k,1)
                 ir = 1.0_rp/r_
                 u_ = phi(i,j,k,2)*ir
                 v_ = phi(i,j,k,3)*ir
                 w_ = phi(i,j,k,4)*ir
                 ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
                 re = phi(i,j,k,5)
                 p_ = (gamma0-1.0_rp)*(re - r_*ek)
                 T_ = p_*ir
                 c_ = sqrt(gamma0*T_)
                 ve = sqrt(u_*u_ + v_*v_ + w_*w_)

                 mL = Sutherland(T_)
                 mT = VIS(i,j,k) - mL

                 ! conservative variables
                 tmp_lcl(i,j,1) = tmp_lcl(i,j,1) + r_
                 tmp_lcl(i,j,2) = tmp_lcl(i,j,2) + r_*u_
                 tmp_lcl(i,j,3) = tmp_lcl(i,j,3) + r_*v_
                 tmp_lcl(i,j,4) = tmp_lcl(i,j,4) + r_*w_
                 tmp_lcl(i,j,5) = tmp_lcl(i,j,5) + re
        
                 ! Reynolds stress pieces
                 tmp_lcl(i,j,6) = tmp_lcl(i,j,6) + (r_*u_)*(r_*u_)*ir
                 tmp_lcl(i,j,7) = tmp_lcl(i,j,7) + (r_*v_)*(r_*v_)*ir
                 tmp_lcl(i,j,8) = tmp_lcl(i,j,8) + (r_*w_)*(r_*w_)*ir
                 tmp_lcl(i,j,9) = tmp_lcl(i,j,9) + (r_*u_)*(r_*v_)*ir
        
                 ! density 
                 tmp_lcl(i,j,10) = tmp_lcl(i,j,10) + r_*r_

                 !pressure
                 tmp_lcl(i,j,11) = tmp_lcl(i,j,11) + p_
                 tmp_lcl(i,j,12) = tmp_lcl(i,j,12) + p_*p_

                 ! temperature
                 tmp_lcl(i,j,13) = tmp_lcl(i,j,13) + T_
                 tmp_lcl(i,j,14) = tmp_lcl(i,j,14) + r_*T_*T_

                 ! viscosity
                 tmp_lcl(i,j,15) = tmp_lcl(i,j,15) + mL
                 tmp_lcl(i,j,16) = tmp_lcl(i,j,16) + mT
        
                 ! energies
                 tmp_lcl(i,j,17) = tmp_lcl(i,j,17) + r_*T_

                 ! Mach numbers
                 tmp_lcl(i,j,18) = tmp_lcl(i,j,18) + ve/c_
                 tmp_lcl(i,j,19) = tmp_lcl(i,j,19) + r_*ve/c_

                 ! derivatives
                 idx = xstep_i(i)
                 idy = ystep_i(j)
                 idz = zstep_i(k)
                 !
                 r_x = 0.0_rp
                 rux = 0.0_rp
                 rvx = 0.0_rp
                 rwx = 0.0_rp
                 !
                 r_y = 0.0_rp
                 ruy = 0.0_rp
                 rvy = 0.0_rp
                 rwy = 0.0_rp
                 !
                 r_z = 0.0_rp
                 ruz = 0.0_rp
                 rvz = 0.0_rp
                 rwz = 0.0_rp
                 do l = -central_fd_order/2, central_fd_order/2

                    cl = central_1(l)

                    r_x = r_x + cl*phi(i+l,j,k,1)
                    rux = rux + cl*phi(i+l,j,k,2)
                    rvx = rvx + cl*phi(i+l,j,k,3)
                    rwx = rwx + cl*phi(i+l,j,k,4)

                    r_y = r_y + cl*phi(i,j+l,k,1)
                    ruy = ruy + cl*phi(i,j+l,k,2)
                    rvy = rvy + cl*phi(i,j+l,k,3)
                    rwy = rwy + cl*phi(i,j+l,k,4)

                    r_z = r_z + cl*phi(i,j,k+l,1)
                    ruz = ruz + cl*phi(i,j,k+l,2)
                    rvz = rvz + cl*phi(i,j,k+l,3)
                    rwz = rwz + cl*phi(i,j,k+l,4)

                 enddo
                 r_x = r_x*idx
                 rux = rux*idx
                 rvx = rvx*idx
                 rwx = rwx*idx
                 !
                 r_y = r_y*idy
                 ruy = ruy*idy
                 rvy = rvy*idy
                 rwy = rwy*idy
                 !
                 r_z = r_z*idz
                 ruz = ruz*idz
                 rvz = rvz*idz
                 rwz = rwz*idz

                 tmp_lcl(i,j,20) = tmp_lcl(i,j,20) + r_x
                 tmp_lcl(i,j,21) = tmp_lcl(i,j,21) + r_y
                 tmp_lcl(i,j,22) = tmp_lcl(i,j,22) + r_z
                 !
                 tmp_lcl(i,j,23) = tmp_lcl(i,j,23) + rux
                 tmp_lcl(i,j,24) = tmp_lcl(i,j,24) + ruy
                 tmp_lcl(i,j,25) = tmp_lcl(i,j,25) + ruz
                 !
                 tmp_lcl(i,j,26) = tmp_lcl(i,j,26) + rvx
                 tmp_lcl(i,j,27) = tmp_lcl(i,j,27) + rvy
                 tmp_lcl(i,j,28) = tmp_lcl(i,j,28) + rvz
                 !
                 tmp_lcl(i,j,29) = tmp_lcl(i,j,29) + rwx
                 tmp_lcl(i,j,30) = tmp_lcl(i,j,30) + rwy
                 tmp_lcl(i,j,31) = tmp_lcl(i,j,31) + rwz
        
                 ! weno sensor
                 if(hybrid_weno) then
                   sensor = SSENSOR(i,j,k)
                   tmp_lcl(i,j,32) = tmp_lcl(i,j,32) + sensor
                   tmp_lcl(i,j,33) = tmp_lcl(i,j,33) + sensor*sensor

                   weno_flag = real(weno%flag(i,j,k),rp)
                   tmp_lcl(i,j,34) = tmp_lcl(i,j,34) + weno_flag
                   tmp_lcl(i,j,35) = tmp_lcl(i,j,35) + weno_flag*weno_flag
                 endif

              enddo
           enddo
        enddo
        
        irNz = 1.0_rp/real(nz,rp)
        if(mpi_flag) then
          npt = nVar*(ubx-lbx+1)*(uby-lby+1)
          call mpi_allreduce(tmp_lcl,tmp_gbl,npt,prec,mpi_sum,comm,err)
          tmp_gbl = tmp_gbl*irNz
        else
          tmp_gbl = tmp_lcl*irNz
        endif

        ! summing up statistics
        r_itStat   = real(itStat,rp)
        r_itStatm1 = real(itStat-1,rp)
        do       l = 1,nVar
           do    j = sy,ey
              do i = sx,ex
                 vmean(i,j,l) = (vmean(i,j,l)*r_itStatm1 + tmp_gbl(i,j,l))/r_itStat
              enddo
           enddo
        enddo

        do j    = sy,ey
           do i = sx,ex
              
              r_ = vmean(i,j,1)
              ir = 1.0_rp/r_
              ru = vmean(i,j,2)
              rv = vmean(i,j,3)
              rw = vmean(i,j,4)

              ! u,v Favre
              vmean2D_aux(i,j,1) = ru*ir
              vmean2D_aux(i,j,2) = rv*ir

              ! u"u", v"v", w"w", u"v"
              vmean2D_aux(i,j,3) = (vmean(i,j,6) - ru*ru*ir)*ir
              vmean2D_aux(i,j,4) = (vmean(i,j,7) - rv*rv*ir)*ir
              vmean2D_aux(i,j,5) = (vmean(i,j,8) - rw*rw*ir)*ir
              vmean2D_aux(i,j,6) = (vmean(i,j,9) - ru*rv*ir)*ir

              ! r_x, r_y, r_z
              r_x = vmean(i,j,20)
              r_y = vmean(i,j,21)
              r_z = vmean(i,j,22)

              ! u_x, u_y, u_z
              vmean2D_aux(i,j, 7) = (vmean(i,j,23)*r_ - r_x*ru)*ir**2
              vmean2D_aux(i,j, 8) = (vmean(i,j,24)*r_ - r_y*ru)*ir**2
              vmean2D_aux(i,j, 9) = (vmean(i,j,25)*r_ - r_z*ru)*ir**2

              ! v_x, v_y, v_z
              vmean2D_aux(i,j,10) = (vmean(i,j,26)*r_ - r_x*rv)*ir**2
              vmean2D_aux(i,j,11) = (vmean(i,j,27)*r_ - r_y*rv)*ir**2
              vmean2D_aux(i,j,12) = (vmean(i,j,28)*r_ - r_z*rv)*ir**2

              ! w_x, w_y, w_z
              vmean2D_aux(i,j,13) = (vmean(i,j,29)*r_ - r_x*rw)*ir**2
              vmean2D_aux(i,j,14) = (vmean(i,j,30)*r_ - r_y*rw)*ir**2
              vmean2D_aux(i,j,15) = (vmean(i,j,31)*r_ - r_z*rw)*ir**2

              ! omega_x, omega_y, omega_z
              omx = vmean2D_aux(i,j,14) - vmean2D_aux(i,j,12)
              omy = vmean2D_aux(i,j, 9) - vmean2D_aux(i,j,13)
              omz = vmean2D_aux(i,j,10) - vmean2D_aux(i,j, 8)

              vmean2D_aux(i,j,16) = omx
              vmean2D_aux(i,j,17) = omy
              vmean2D_aux(i,j,18) = omz
              
              ! vorticity magnitude
              vmean2D_aux(i,j,19) = sqrt(omx*omx + omy*omy + omz*omz)

              ! T favre
              vmean2D_aux(i,j,20) = vmean(i,j,17)/vmean(i,j,1)

           enddo
        enddo

        return
end subroutine compute_2DStats


subroutine compute_xz_plane_stats(var,nz,nVar,itStat,comm,vmean,mpi_flag)
        
        use parameters_module, only: rp, MPI_RP
        use mpi_module       , only: sx,ex,sz,ez,lbx,ubx

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(in)    :: var
        real(rp), dimension(:,:)  , allocatable, intent(inout) :: vmean

        integer, intent(in) :: nz,nVar,itStat,comm
        logical, intent(in) :: mpi_flag

        ! local declarations
        real(rp), dimension(lbx:ubx,nVar) :: tmp_lcl, tmp_gbl
        real(rp) :: irNz, r_itStat, r_itStatm1
        integer  :: i,k, l, npt, err = 0

        tmp_lcl = 0.0_rp
        tmp_gbl = 0.0_rp
        do       l = 1,nVar
           do    k = sz,ez
              do i = sx,ex

                 tmp_lcl(i,l) = tmp_lcl(i,l) + var(i,k,l)

              enddo
           enddo
        enddo

        ! reduce the sum
        irNz = 1.0_rp/real(nz,rp)
        if(mpi_flag) then
          npt = nVar*(ubx-lbx+1)
          call mpi_allreduce(tmp_lcl,tmp_gbl,npt,MPI_RP,mpi_sum,comm,err)
          tmp_gbl = tmp_gbl*irNz
        else
          tmp_gbl = tmp_gbl*irNz
        endif

        ! summing up the statistics in time
        r_itStat   = real(itStat,rp)
        r_itStatm1 = real(itStat-1,rp)
        do    l = 1,nVar
           do i = sx,ex
              vmean(i,l) = (vmean(i,l)*r_itStatm1 + tmp_gbl(i,l))/r_itStat
           enddo
        enddo


        return
end subroutine compute_xz_plane_stats

subroutine compute_xz_plane_tch(var,nx,nz,nVar,itStat,comm,vmean,mpi_flag)
        
        use parameters_module, only: rp, MPI_RP
        use mpi_module       , only: sx,ex,sz,ez

        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(in)    :: var
        real(rp), dimension(:)    , allocatable, intent(inout) :: vmean

        integer, intent(in) :: nx,nz,nVar,itStat,comm
        logical, intent(in) :: mpi_flag

        ! local declarations
        real(rp), dimension(nVar) :: tmp_lcl, tmp_gbl
        real(rp) :: irNxNz, r_itStat, r_itStatm1
        integer  :: i,k, l, npt, err = 0

        tmp_lcl = 0.0_rp
        tmp_gbl = 0.0_rp
        do       l = 1,nVar
           do    k = sz,ez
              do i = sx,ex

                 tmp_lcl(l) = tmp_lcl(l) + var(i,k,l)

              enddo
           enddo
        enddo

        ! reduce the sum
        irNxNz = 1.0_rp/real(nx,rp)/real(nz,rp)
        if(mpi_flag) then
          npt = nVar
          call mpi_allreduce(tmp_lcl,tmp_gbl,npt,MPI_RP,mpi_sum,comm,err)
          tmp_gbl = tmp_gbl*irNxNz
        else
          tmp_gbl = tmp_gbl*irNxNz
        endif

        ! summing up the statistics in time
        r_itStat   = real(itStat,rp)
        r_itStatm1 = real(itStat-1,rp)
        do l = 1,nVar
           vmean(l) = (vmean(l)*r_itStatm1 + tmp_gbl(l))/r_itStat
        enddo


        return
end subroutine compute_xz_plane_tch








end module onlineStats
