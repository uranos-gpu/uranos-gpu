module output_tch
use parameters_module, only: rp

implicit none
private
public stat_turbulent_channel

type wall_type
  real(rp) :: rho
  real(rp) :: vis
  real(rp) :: tmp
  real(rp) :: u_t
  real(rp) :: tau
  real(rp) :: u_y
  real(rp) :: ReTau
  real(rp) :: MaTau
endtype

type bulk_type
  real(rp) :: vel
  real(rp) :: tmp
endtype

type plus_unit
  real(rp), dimension(:), allocatable :: y
  real(rp), dimension(:), allocatable :: yT
  real(rp), dimension(:), allocatable :: u
  real(rp), dimension(:), allocatable :: uVD
  real(rp), dimension(:), allocatable :: uT
end type plus_unit

type(wall_type) :: wl
type(bulk_type) :: bk
type(plus_unit) :: plus

contains
subroutine stat_turbulent_channel

        use parameters_module
        use mpi_module
        use storage_module
        use FileModule
        use onlineStats
        use integration_module    , only: mean_field
        use wmles_module          , only: compute_tau_wall_wmles1D

        implicit none

        if(cart_dims(2) > 1) &
        stop ' Turbulent channel online statistics supports only 2Ddecomp!'
        
        if(wmles) call compute_xz_plane_tch(WMLES_DATA,nx,nz,nVWMLESData,&
                              itStat,comm2Dxz,vmean0D_wmles,mpi_flag)

        call mean_field(U,sx,ex,xstep,sy,ey,ystep,sz,ez,zstep,bk%vel,mpi_flag)
        call mean_field(T,sx,ex,xstep,sy,ey,ystep,sz,ez,zstep,bk%tmp,mpi_flag)

        call symmetry(sy,ey,nvAve1D,vmean1D)


        if(rank == root) then
          allocate(plus%y(0:ny), plus%yT(0:ny), plus%u(0:ny), &
                   plus%uVD(0:ny), plus%uT(0:ny))
        
          call GetMeanWallQuantities(sy,vmean1D,wl)

          call ComputePlusUnits(sy,ey,y,vmean1D,wl,plus)
        
          call TransformVanDriest(sy,ey,vmean1D,wl,plus)

          call TransformTrettelLaarson(sy,ey,y,vmean1D,wl,plus)
        
          call write_wall_normal_stats(sy,ey,y,plus,vmean1D,wl,bk)
        
          call write_stress_budget(sy,ey,y,vmean1D)

          call write_monitor(time,wl,bk,force_turbulent_channel)

          deallocate(plus%y,plus%u,plus%uVD,plus%yT,plus%uT)

        endif
        

        return
end subroutine stat_turbulent_channel






subroutine GetMeanWallQuantities(jW,vmean1D,wl)
        
        use parameters_module     , only: gamma0, mu_inf, bc, wmles
        use fluid_functions_module, only: Sutherland
        use storage_module        , only: vmean0D_wmles

        implicit none
        
        real(rp), dimension(:,:), allocatable, intent(in)  :: vmean1D
        integer                              , intent(in)  :: jW
        type(wall_type)                      , intent(out) :: wl
        
        ! get wall quantities
        wl%tmp = 1.0_rp
        wl%rho = vmean1D(jW,11)/wl%tmp
        wl%vis = Sutherland(wl%tmp)
        wl%u_y = (vmean1D(jW,18)*vmean1D(jW,1) - vmean1D(jW,2)*vmean1D(jW,17))/(vmean1D(jW,1)**2)

        ! compute tau wall
        wl%tau = mu_inf*wl%u_y
        wl%u_t = sqrt(wl%tau/wl%rho)
        if(wmles) then
          wl%tau = vmean0D_wmles(3)
          wl%u_t = sqrt(wl%tau/wl%rho)
        endif
        
        ! get wall groups
        wl%ReTau = wl%rho*wl%u_t/(mu_inf*wl%vis)
        wl%MaTau = wl%u_t/sqrt(gamma0*wl%tmp)

        return
end subroutine GetMeanWallQuantities



subroutine ComputePlusUnits(sy,ey,y,vmean1D,wl,plus)
        implicit none
        real(rp), allocatable, dimension(:,:), intent(in)    :: vmean1D
        real(rp), allocatable, dimension(:)  , intent(in)    :: y
        integer                              , intent(in)    :: sy,ey
        type(wall_type)                      , intent(in)    :: wl
        type(plus_unit)                      , intent(inout) :: plus

        ! local declarations
        real(rp) :: yj
        integer  :: j

        plus%u = 0.0_rp
        plus%y = 0.0_rp
        do j = sy,ey
           yj        = y(j)+1.0_rp
           plus%y(j) = yj*wl%ReTau
           plus%u(j) = (vmean1D(j,2)/vmean1D(j,1))/wl%u_t
        enddo

        return
end subroutine ComputePlusUnits


subroutine TransformVanDriest(sy,ey,vmean1D,wl,plus)
        implicit none
        real(rp), allocatable, dimension(:,:), intent(in)    :: vmean1D
        integer                              , intent(in)    :: sy,ey
        type(wall_type)                      , intent(in)    :: wl
        type(plus_unit)                      , intent(inout) :: plus
        
        real(rp) :: uc, du
        integer  :: j

        plus%uVD = 0.0_rp
        do j = sy,ey
           du = plus%u(j) - plus%u(j-1)
           uc = sqrt(vmean1D(j,1)/wl%rho)
           plus%uVD(j) = plus%uVD(j-1) + uc*du
        enddo

        return
end subroutine TransformVanDriest



subroutine TransformTrettelLaarson(sy,ey,y,vmean1D,wl,plus)
        implicit none
        real(rp), allocatable, dimension(:,:), intent(in)    :: vmean1D
        real(rp), allocatable, dimension(:)  , intent(in)    :: y
        integer                              , intent(in)    :: sy,ey
        type(wall_type)                      , intent(in)    :: wl
        type(plus_unit)                      , intent(inout) :: plus

        ! local
        real(rp), dimension(sy-1:ey) :: ft, dft,gt
        real(rp)                     :: r_, ir, rr, nn
        real(rp)                     :: df, dg
        integer                      :: j

        do j = sy,ey
           r_ = vmean1D(j,1)
           ir = 1.0_rp/r_
           rr = sqrt(r_/wl%rho)
           nn = (vmean1D(j,15) + vmean1D(j,16))*ir*(wl%rho/wl%vis)

           ft(j) = (y(j) + 1.0_rp)/(rr*nn)
        enddo
        do j = sy+1,ey-1
           dft(j) = (ft(j+1) - ft(j-1))/(y(j+1) - y(j-1))
        enddo
        dft(1)  = (-1.5_rp*ft(1)+2.0_rp*ft(2)-0.5_rp*ft(3))/&
                  (-1.5_rp* y(1)+2.0_rp* y(2)-0.5_rp* y(3))
        dft(ey) = ( 0.5_rp*ft(ey-2)-2._rp*ft(ey-1)+1.5_rp*ft(ey))/&
                  ( 0.5_rp* y(ey-2)-2._rp* y(ey-1)+1.5_rp* y(ey))

        do j = sy,ey
           r_ = vmean1D(j,1)
           ir = 1.0_rp/r_
           rr = r_/wl%rho
           nn = (vmean1D(j,15) + vmean1D(j,16))*ir*(wl%rho/wl%vis)

           gt(j) = dft(j)*rr*nn
        enddo
              
        plus%yT(0) = 0.0_rp
        plus%uT(0) = 0.0_rp
        plus%yT(1) = plus%y(1)*dft(1)
        plus%uT(1) = plus%u(1)*gt (1)
        do j = sy+1,ey
           df = 0.5_rp*(dft(j)+dft(j-1))
           dg = 0.5_rp*(gt (j)+gt (j-1))
           plus%yT(j) = plus%yT(j-1) + (plus%y(j) - plus%y(j-1))*df
           plus%uT(j) = plus%uT(j-1) + (plus%u(j) - plus%u(j-1))*dg
        enddo

        return
end subroutine TransformTrettelLaarson










subroutine write_wall_normal_stats(sy,ey,y,plus,vmean1D,wl,bk)
        
        use parameters_module, only: rp, Mach, Reynolds, gamma0, data_dir, it
        use FileModule

        implicit none
        integer                              , intent(in) :: sy,ey
        real(rp), dimension(:,:), allocatable, intent(in) :: vmean1D
        real(rp), dimension(:)  , allocatable, intent(in) :: y
        type(plus_unit)                      , intent(in) :: plus
        type(wall_type)                      , intent(in) :: wl
        type(bulk_type)                      , intent(in) :: bk

        ! local
        type(FileType) :: WNF
        real(rp)       :: yj
        real(rp)       :: mean_rhuu, mean_rhvv, mean_rhww, mean_rhuv
        real(rp)       :: meanR_u_u, meanR_v_v, meanR_w_w, meanR_u_v
        real(rp)       :: r_, ir, uu, vv, ww, uv, R11, R22, R33, R12
        integer        :: j

        WNF%name = 'stat'
        WNF%dir  = trim(data_dir)//'/VEL_STATS'
        call OpenNewFile(WNF,it)


        write(WNF%unit,'(A)')       '# *************************************************'
        write(WNF%unit,'(A,f18.6)') '#  TURBULENT CHANNEL AT MACH ', Mach 
        write(WNF%unit,'(A)')       '#'
        write(WNF%unit,'(A)')       '#  author: Francesco De Vanna '
        write(WNF%unit,'(A)')       '#  e-mail: fra.devanna@gmail.com'
        write(WNF%unit,'(A)')       '#'
        write(WNF%unit,'(A,f18.6)') '#  FRICTION REYNOLDS NUMBER  ', wl%ReTau
        write(WNF%unit,'(A,f18.6)') '#  BULK     REYNOLDS NUMBER  ', 2*Reynolds
        write(WNF%unit,'(A,f18.6)') '#  FRICTION MACH     NUMBER  ', wl%MaTau
        write(WNF%unit,'(A,f18.6)') '#  BULK     MACH     NUMBER  ', bk%vel/sqrt(gamma0*bk%tmp)
        write(WNF%unit,'(A)')       '# *************************************************'
        write(WNF%unit,'(A)')       '# Column 1  : y '
        write(WNF%unit,'(A)')       '# Column 2  : y+'
        write(WNF%unit,'(A)')       '# Column 3  : yT+'
        write(WNF%unit,'(A)')       '# Column 4  : u+'
        write(WNF%unit,'(A)')       '# Column 5  : u_VD+'
        write(WNF%unit,'(A)')       '# Column 6  : u_Trettel+'
        write(WNF%unit,'(A)')       '# Column 7  : tau11'
        write(WNF%unit,'(A)')       '# Column 8  : tau22'
        write(WNF%unit,'(A)')       '# Column 9  : tau33'
        write(WNF%unit,'(A)')       '# Column 10 : tau12'
        write(WNF%unit,'(A)')       '# Column 11 : rho/rhoWall'
        write(WNF%unit,'(A)')       '# Column 12 : Temp/TempWall'
        write(WNF%unit,'(A)')       '# Column 13 : MuMol/MuMolWall'
        write(WNF%unit,'(A)')       '# Column 14 : MuTurb/MuMol'
        write(WNF%unit,'(A)')       '# *************************************************'
        
        ! write to the file
        do j = sy,ey

           yj = y(j)+1.0_rp

           r_ = vmean1D(j,1)
           ir = 1.0_rp/r_

           ! mean(rho*u_i*u_j)
           mean_rhuu = vmean1D(j,6)
           mean_rhvv = vmean1D(j,7)
           mean_rhww = vmean1D(j,8)
           mean_rhuv = vmean1D(j,9)

           ! mean(rho)*tilde(u_i)*tilde(u_j)
           meanR_u_u = vmean1D(j,2)*vmean1D(j,2)*ir
           meanR_v_v = vmean1D(j,3)*vmean1D(j,3)*ir
           meanR_w_w = vmean1D(j,4)*vmean1D(j,4)*ir
           meanR_u_v = vmean1D(j,2)*vmean1D(j,3)*ir

           ! u_i''u_j''
           uu = (mean_rhuu - meanR_u_u)*ir
           vv = (mean_rhvv - meanR_v_v)*ir
           ww = (mean_rhww - meanR_w_w)*ir
           uv = (mean_rhuv - meanR_u_v)*ir

           R11 = (r_/wl%rho)*uu/wl%u_t**2
           R22 = (r_/wl%rho)*vv/wl%u_t**2
           R33 = (r_/wl%rho)*ww/wl%u_t**2
           R12 = (r_/wl%rho)*uv/wl%u_t**2

           write(WNF%unit,'(20e18.9)') yj, plus%y(j), plus%yT(j), &
                     plus%u(j), plus%uVD(j), plus%uT(j), &
                     R11, R22, R33, R12, r_/wl%rho, &
                     vmean1D(j,13), vmean1D(j,15), vmean1D(j,16)/vmean1D(j,15)
        enddo
        call CloseFile(WNF)
        





        return
end subroutine write_wall_normal_stats





subroutine write_stress_budget(sy,ey,y,vmean1D)
        
        use FileModule
        use parameters_module, only: rp, data_dir,it, Mach, mu_inf

        implicit none
        real(rp), allocatable, dimension(:,:), intent(in) :: vmean1D
        real(rp), allocatable, dimension(:)  , intent(in) :: y
        integer                              , intent(in) :: sy,ey

        ! local declaration
        type(FileType) :: BF
        real(rp)       :: yj,r_, ir, uv, R12, u_yj, u_yjp
        real(rp)       :: mMol, mTrb, mLUY, mTUY
        real(rp)       :: mean_rhuv, meanR_u_v
        integer        :: j

        !
        ! === open the file and write the header
        !
        BF%name = 'stat'
        BF%dir  = trim(data_dir)//'/BUD_STATS'
        call OpenNewFile(BF,it)

        write(BF%unit,'(A)')       '# *************************************************'
        write(BF%unit,'(A,f18.6)') '#  BUDGET STRESS FOR TURB. CHANNEL AT MACH ', Mach 
        write(BF%unit,'(A)')       '#'
        write(BF%unit,'(A)')       '#  author: Francesco De Vanna '
        write(BF%unit,'(A)')       '#  e-mail: fra.devanna@gmail.com'
        write(BF%unit,'(A)')       '#'
        write(BF%unit,'(A)')       '# *************************************************'
        write(BF%unit,'(A)')       '# Column 1  : y '
        write(BF%unit,'(A)')       '# Column 2  : mu_mol du/dy'
        write(BF%unit,'(A)')       '# Column 3  : mu_trb du/dy'
        write(BF%unit,'(A)')       '# Column 4  : Reynolds Stress'
        write(BF%unit,'(A)')       '# Column 5  : Total Stress'
        write(BF%unit,'(A)')       '# *************************************************'
        !
        ! === write to file
        !
        write(BF%unit,'(20e18.9)') 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 1.0_rp
        do j = sy,ey-1
        
           yj = 0.5_rp*(y(j)+y(j+1))+1.0_rp
           r_ = 0.5_rp*(vmean1D(j,1)+vmean1D(j+1,1))
           ir = 1.0_rp/r_
              
           ! Reynolds stress
           mean_rhuv = 0.5_rp*(vmean1D(j,9)+vmean1D(j+1,9))
           meanR_u_v = 0.5_rp*(vmean1D(j,2)+vmean1D(j+1,2))*&
                       0.5_rp*(vmean1D(j,3)+vmean1D(j+1,3))*ir
           uv        = (mean_rhuv - meanR_u_v)*ir
           R12       = (r_/wl%rho)*uv/wl%u_t**2

           mMol = 0.5_rp*(vmean1D(j,15)+vmean1D(j+1,15))
           mTrb = 0.5_rp*(vmean1D(j,16)+vmean1D(j+1,16))

           ! Viscous stress
           u_yj  = (vmean1D(j  ,18)*vmean1D(j  ,1) - vmean1D(j  ,2)*vmean1D(j  ,17))/(vmean1D(j  ,1)**2)
           u_yjp = (vmean1D(j+1,18)*vmean1D(j+1,1) - vmean1D(j+1,2)*vmean1D(j+1,17))/(vmean1D(j+1,1)**2)

           mLUY = mu_inf * mMol * 0.5_rp*(u_yj+u_yjp)/(wl%rho*wl%u_t**2)
           mTUY = mu_inf * mTrb * 0.5_rp*(u_yj+u_yjp)/(wl%rho*wl%u_t**2)

           write(BF%unit,'(20e18.9)') yj, mLUY, mTUY, - R12, &
                                            mLUY + mTUY - R12
        enddo
        write(BF%unit,'(20e18.9)') 2.0_rp, -1.0_rp, 0.0_rp, 0.0_rp, -1.0_rp

        call CloseFile(BF)

        return
end subroutine write_stress_budget



subroutine write_monitor(time,wl,bk,force_turbulent_channel)

        use parameters_module, only: rp, data_dir, it, gamma0, mu_inf
        use FileModule

        implicit none
        real(rp), intent(in) :: time
        real(rp), intent(in) :: force_turbulent_channel

        type(wall_type), intent(in) :: wl
        type(bulk_type), intent(in) :: bk
        type(FileType)              :: forceFile

        forceFile%name = 'TCH_MONITOR'
        forceFile%dir  = trim(data_dir)
        call AppendToFile(forceFile,it)
        
        ! write force turbulent channel
        write(forceFile%unit,10)      &
             it                      , &
             time*wl%u_t             , &
             wl%ReTau                , &
             wl%MaTau                , &
             force_turbulent_channel , &
             bk%vel/sqrt(gamma0)     , &
             2*bk%vel/mu_inf         , &
             bk%tmp                  , &
             2*(wl%u_t/bk%vel)**2

        call CloseFile(forceFile)

        10 format(I7,20e18.9)

        return
end subroutine write_monitor




subroutine symmetry(sy,ey,nv,v)
        
        use parameters_module, only: rp

        implicit none
        real(rp), dimension(:,:), allocatable, intent(inout) :: v
        integer                              , intent(in)    :: sy, ey, nv

        ! local declarations
        logical, dimension(nv) :: sim
        integer                :: id, j, j_sym
        
        sim = .true.
        sim( 3) = .false.
        sim( 9) = .false.
        sim(17) = .false.
        sim(18) = .false.
        sim(19) = .false.
        sim(20) = .false.
        
        do id = 1,nv

           if(sim(id)) then
             do j = sy,(ey-sy+1)/2
                j_sym = ey+1-j + (sy-1)

                v(j,id) = 0.5_rp*(v(j,id) + v(j_sym,id))
                v(j_sym,id) = v(j,id)
             enddo
            else
             do j = sy,(ey-sy+1)/2
                j_sym = ey+1-j + (sy-1)

                v(j,id) = 0.5_rp*(v(j,id) - v(j_sym,id))
                v(j_sym,id) = -v(j,id)
             enddo
           endif

        enddo

        return
end subroutine symmetry























end module output_tch
