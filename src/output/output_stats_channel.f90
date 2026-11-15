module output_stats_channel
use parameters_module, only: rp, Tref, vis_flag

implicit none
private
public stat_turbulent_channel

type wall_type
  real(rp) :: rho
  real(rp) :: vis
  real(rp) :: tmp
  real(rp) :: u_t
  real(rp) :: tau
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

        implicit none

        if(cart_dims(2) > 1) &
        stop ' Turbulent channel online statistics supports only 2Ddecomp!'

        if(itStat <= 1) return

        call bulk_mean(U,sy,ey,bk%vel)
        call bulk_mean(T,sy,ey,bk%tmp)

        call symmetry(sy,ey,nvAve1D,vmean1D)
        !$acc update device(vmean1D)

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

        use parameters_module     , only: mu_inf, bc, Tref
        use thermodynamics_module , only: sound_speed
        use fluid_functions_module, only: laminar_viscosity

        implicit none

        real(rp), dimension(:,:), allocatable, intent(in)  :: vmean1D
        integer                              , intent(in)  :: jW
        type(wall_type)                      , intent(out) :: wl

        wl%tmp = 1.0_rp
        wl%rho = vmean1D(jW,11)/wl%tmp
        wl%vis = laminar_viscosity(wl%tmp,Tref,vis_flag)

        wl%tau = mu_inf*(vmean1D(jW,19) + vmean1D(jW,20))
        wl%u_t = sqrt(wl%tau/wl%rho)

        wl%ReTau = wl%rho*wl%u_t/(mu_inf*wl%vis)
        wl%MaTau = wl%u_t/sound_speed(wl%tmp)

        return
end subroutine GetMeanWallQuantities

subroutine ComputePlusUnits(sy,ey,y,vmean1D,wl,plus)
        implicit none
        real(rp), allocatable, dimension(:,:), intent(in)    :: vmean1D
        real(rp), allocatable, dimension(:)  , intent(in)    :: y
        integer                              , intent(in)    :: sy,ey
        type(wall_type)                      , intent(in)    :: wl
        type(plus_unit)                      , intent(inout) :: plus

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
        plus%uVD(sy) = sqrt(vmean1D(sy,1)/wl%rho)*plus%u(sy)
        do j = sy+1,ey
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

        use parameters_module, only: rp, Mach, Reynolds, data_dir, it
        use thermodynamics_module, only: sound_speed
        use FileModule

        implicit none
        integer                              , intent(in) :: sy,ey
        real(rp), dimension(:,:), allocatable, intent(in) :: vmean1D
        real(rp), dimension(:)  , allocatable, intent(in) :: y
        type(plus_unit)                      , intent(in) :: plus
        type(wall_type)                      , intent(in) :: wl
        type(bulk_type)                      , intent(in) :: bk

        type(FileType) :: WNF
        real(rp)       :: yj
        real(rp)       :: mean_rhuu, mean_rhvv, mean_rhww, mean_rhuv
        real(rp)       :: mean_rhuuu, mean_rhvvv, mean_rhwww, mean_rhuuv
        real(rp)       :: mean_rhuuuu, mean_rhvvvv, mean_rhwwww
        real(rp)       :: meanR_u_u, meanR_v_v, meanR_w_w, meanR_u_v
        real(rp)       :: tilde_u, tilde_v, tilde_w
        real(rp)       :: r_, ir, uu, vv, ww, uv, R11, R22, R33, R12
        real(rp)       :: uuu, vvv, www, uuv
        reaL(rp)       :: uuuu, vvvv, wwww
        real(rp)       :: pp
        real(rp)       :: tilde_T, uT, vT, TTr, TTf, MM
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
        write(WNF%unit,'(A,f18.6)') '#  BULK     MACH     NUMBER  ', bk%vel/sound_speed(bk%tmp)
        write(WNF%unit,'(A,f18.6)') '#  FRICTION VELOCITY         ', wl%u_t
        write(WNF%unit,'(A,f18.6)') '#  RHO WALL                  ', wl%rho
        write(WNF%unit,'(A)')       '# -------------------------------------------------'
        write(WNF%unit,'(A)')       '# All quantities are non-dimensional. Running means over x, z and time,'
        write(WNF%unit,'(A)')       '# with the two halves of the channel folded by symmetry. y spans the'
        write(WNF%unit,'(A)')       '# full height, from the lower wall (y = 0) to the upper one (y = 2).'
        write(WNF%unit,'(A)')       "# <f> = Reynolds mean, tilde(f) = <rho f>/<rho> = Favre mean,"
        write(WNF%unit,'(A)')       "# f' = f - <f>, f'' = f - tilde(f). Wall values: T_w = 1, mu_w = mu(T_w),"
        write(WNF%unit,'(A)')       '# rho_w = <p>/T_w at the first node, c_w = c(T_w).'
        write(WNF%unit,'(A)')       '# tau_w = mu_inf <(mu_mol + mu_trb) du/dy> at the first node (half a cell'
        write(WNF%unit,'(A)')       '# off the wall), u_tau = sqrt(tau_w/rho_w).'
        write(WNF%unit,'(A)')       '# FRICTION REYNOLDS NUMBER = rho_w u_tau h/(mu_inf mu_w), h = 1.'
        write(WNF%unit,'(A)')       '# FRICTION MACH NUMBER     = u_tau/c_w.'
        write(WNF%unit,'(A)')       '# TURBULENT CHANNEL AT MACH = Mach of the namelist = <rho u>/(rho_w c_w), the'
        write(WNF%unit,'(A)')       '#   mass-flux Mach number. The bulk one of the literature is U_b/c_w = Mach/<rho>.'
        write(WNF%unit,'(A)')       '# BULK REYNOLDS NUMBER     = <rho u> 2h/(mu_inf mu_w) = 2 Reynolds of the namelist,'
        write(WNF%unit,'(A)')       '#   exact because the forcing keeps <rho u> fixed.'
        write(WNF%unit,'(A)')       '# BULK MACH NUMBER         = unweighted volume mean of u over c at the volume'
        write(WNF%unit,'(A)')       '#   mean of T. NOT the U_b/c_w of the literature: rebuild U_b = <rho u>/<rho>'
        write(WNF%unit,'(A)')       '#   from columns 4 and 11.'
        write(WNF%unit,'(A)')       '# *************************************************'
        write(WNF%unit,'(A)')       '# Column 1  : y, distance from the lower wall'
        write(WNF%unit,'(A)')       '# Column 2  : y+ = y rho_w u_tau/(mu_inf mu_w)'
        write(WNF%unit,'(A)')       '# Column 3  : yT+, Trettel-Larsson semi-local wall distance'
        write(WNF%unit,'(A)')       '# Column 4  : u+ = tilde(u)/u_tau'
        write(WNF%unit,'(A)')       '# Column 5  : u_VD+, van Driest transformed velocity'
        write(WNF%unit,'(A)')       '# Column 6  : u_TL+, Trettel-Larsson transformed velocity'
        write(WNF%unit,'(A)')       "# Column 7  : <rho> tilde(u'' u'')/(rho_w u_tau^2), Reynolds stress (not viscous)"
        write(WNF%unit,'(A)')       "# Column 8  : <rho> tilde(v'' v'')/(rho_w u_tau^2)"
        write(WNF%unit,'(A)')       "# Column 9  : <rho> tilde(w'' w'')/(rho_w u_tau^2)"
        write(WNF%unit,'(A)')       "# Column 10 : <rho> tilde(u'' v'')/(rho_w u_tau^2), negative in the lower half"
        write(WNF%unit,'(A)')       '# Column 11 : <rho>/rho_w'
        write(WNF%unit,'(A)')       '# Column 12 : <T>/T_w'
        write(WNF%unit,'(A)')       '# Column 13 : <mu_mol>/mu_w'
        write(WNF%unit,'(A)')       '# Column 14 : <mu_trb>/<mu_mol> (SGS viscosity, 0 in DNS)'
        write(WNF%unit,'(A)')       "# Column 15 : tilde(u'' u'' u'')/u_tau^3"
        write(WNF%unit,'(A)')       "# Column 16 : tilde(v'' v'' v'')/u_tau^3"
        write(WNF%unit,'(A)')       "# Column 17 : tilde(w'' w'' w'')/u_tau^3"
        write(WNF%unit,'(A)')       "# Column 18 : tilde(u'' u'' v'')/u_tau^3"
        write(WNF%unit,'(A)')       "# Column 19 : tilde(u'' u'' u'' u'')/u_tau^4"
        write(WNF%unit,'(A)')       "# Column 20 : tilde(v'' v'' v'' v'')/u_tau^4"
        write(WNF%unit,'(A)')       "# Column 21 : tilde(w'' w'' w'' w'')/u_tau^4"
        write(WNF%unit,'(A)')       "# Column 22 : <p' p'>, variance (not rms), not normalised: divide by tau_w^2 for p+^2"
        write(WNF%unit,'(A)')       '# Column 23 : <cp>'
        write(WNF%unit,'(A)')       '# Column 24 : <gamma>'
        write(WNF%unit,'(A)')       '# Column 25 : <mu_mol/rho>/(mu_w/rho_w)'
        write(WNF%unit,'(A)')       "# Column 26 : tilde(u'' T'')/(u_tau T_w)"
        write(WNF%unit,'(A)')       "# Column 27 : tilde(v'' T'')/(u_tau T_w)"
        write(WNF%unit,'(A)')       "# Column 28 : <T' T'>/T_w^2, Reynolds variance"
        write(WNF%unit,'(A)')       "# Column 29 : tilde(T'' T'')/T_w^2, Favre variance"
        write(WNF%unit,'(A)')       '# Column 30 : <M>, local Mach number |u|/c(T)'
        write(WNF%unit,'(A)')       "# Column 31 : <M' M'>, Reynolds variance"
        write(WNF%unit,'(A)')       '# *************************************************'

        do j = sy,ey

           yj = y(j)+1.0_rp

           r_ = vmean1D(j,1)
           ir = 1.0_rp/r_

           tilde_u = vmean1D(j,2) * ir
           tilde_v = vmean1D(j,3) * ir
           tilde_w = vmean1D(j,4) * ir

           mean_rhuu = vmean1D(j,6)
           mean_rhvv = vmean1D(j,7)
           mean_rhww = vmean1D(j,8)
           mean_rhuv = vmean1D(j,9)

           meanR_u_u = vmean1D(j,2)*vmean1D(j,2)*ir
           meanR_v_v = vmean1D(j,3)*vmean1D(j,3)*ir
           meanR_w_w = vmean1D(j,4)*vmean1D(j,4)*ir
           meanR_u_v = vmean1D(j,2)*vmean1D(j,3)*ir

           mean_rhuuu = vmean1D(j,21)
           mean_rhvvv = vmean1D(j,22)
           mean_rhwww = vmean1D(j,23)
           mean_rhuuv = vmean1D(j,24)

           mean_rhuuuu = vmean1D(j,25)
           mean_rhvvvv = vmean1D(j,26)
           mean_rhwwww = vmean1D(j,27)

           uu = (mean_rhuu - meanR_u_u)*ir
           vv = (mean_rhvv - meanR_v_v)*ir
           ww = (mean_rhww - meanR_w_w)*ir
           uv = (mean_rhuv - meanR_u_v)*ir

           R11 = (r_/wl%rho)*uu/wl%u_t**2
           R22 = (r_/wl%rho)*vv/wl%u_t**2
           R33 = (r_/wl%rho)*ww/wl%u_t**2
           R12 = (r_/wl%rho)*uv/wl%u_t**2

           uuu = (mean_rhuuu - 3.0_rp * mean_rhuu * tilde_u &
                 + 2.0_rp * vmean1D(j,2) * tilde_u * tilde_u) * ir

           vvv = (mean_rhvvv - 3.0_rp * mean_rhvv * tilde_v &
                 + 2.0_rp * vmean1D(j,3) * tilde_v * tilde_v) * ir

           www = (mean_rhwww - 3.0_rp * mean_rhww * tilde_w &
                 + 2.0_rp * vmean1D(j,4) * tilde_w * tilde_w) * ir

           uuv = (mean_rhuuv - 2.0_rp * mean_rhuv * tilde_u - mean_rhuu * tilde_v &
                 + 2.0_rp * vmean1D(j,2) * tilde_u * tilde_v) * ir

           uuuu = (mean_rhuuuu  - 4.0_rp * mean_rhuuu * tilde_u &
                  + 6.0_rp * mean_rhuu * tilde_u * tilde_u &
                   - 3.0_rp * vmean1D(j,2) * tilde_u * tilde_u * tilde_u) * ir

            vvvv = (mean_rhvvvv  - 4.0_rp * mean_rhvvv * tilde_v &
                   + 6.0_rp * mean_rhvv * tilde_v * tilde_v &
                   - 3.0_rp * vmean1D(j,3) * tilde_v * tilde_v * tilde_v) * ir

            wwww = (mean_rhwwww  - 4.0_rp * mean_rhwww * tilde_w &
                   + 6.0_rp * mean_rhww * tilde_w * tilde_w &
                   - 3.0_rp * vmean1D(j,4) * tilde_w * tilde_w * tilde_w) * ir

            pp = vmean1D(j,12) - (vmean1D(j,11) * vmean1D(j,11))

            tilde_T = vmean1D(j,11) * ir
            uT  = vmean1D(j,31) * ir - tilde_u * tilde_T
            vT  = vmean1D(j,32) * ir - tilde_v * tilde_T
            TTr = vmean1D(j,33) - vmean1D(j,13) * vmean1D(j,13)
            TTf = vmean1D(j,14) * ir - tilde_T * tilde_T
            MM  = vmean1D(j,35) - vmean1D(j,34) * vmean1D(j,34)

           write(WNF%unit,'(31e18.9)') yj               , &
                                       plus%y(j)        , &
                                       plus%yT(j)       , &
                                       plus%u(j)        , &
                                       plus%uVD(j)      , &
                                       plus%uT(j)       , &
                                       R11              , &
                                       R22              , &
                                       R33              , &
                                       R12              , &
                                       r_/wl%rho        , &
                                       vmean1D(j,13)/wl%tmp , &
                                       vmean1D(j,15)/wl%vis , &
                                       vmean1D(j,16)/vmean1D(j,15), &
                                       uuu/(wl%u_t**3)  , &
                                       vvv/(wl%u_t**3)  , &
                                       www/(wl%u_t**3)  , &
                                       uuv/(wl%u_t**3)  , &
                                       uuuu/(wl%u_t**4), &
                                       vvvv/(wl%u_t**4), &
                                       wwww/(wl%u_t**4), &
                                       pp               , &
                                       vmean1D(j,28)    , &
                                       vmean1D(j,29)    , &
                                       vmean1D(j,30)*wl%rho/wl%vis, &
                                       uT/(wl%u_t*wl%tmp), &
                                       vT/(wl%u_t*wl%tmp), &
                                       TTr/wl%tmp**2    , &
                                       TTf/wl%tmp**2    , &
                                       vmean1D(j,34)    , &
                                       MM
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

        type(FileType) :: BF
        real(rp)       :: yj,r_, ir, uv, R12
        real(rp)       :: mLUY, mTUY
        real(rp)       :: mean_rhuv, meanR_u_v
        integer        :: j

        BF%name = 'stat'
        BF%dir  = trim(data_dir)//'/BUD_STATS'
        call OpenNewFile(BF,it)

        write(BF%unit,'(A)')       '# *************************************************'
        write(BF%unit,'(A,f18.6)') '#  BUDGET STRESS FOR TURB. CHANNEL AT MACH ', Mach
        write(BF%unit,'(A)')       '#'
        write(BF%unit,'(A)')       '#  author: Francesco De Vanna '
        write(BF%unit,'(A)')       '#  e-mail: fra.devanna@gmail.com'
        write(BF%unit,'(A)')       '#'
        write(BF%unit,'(A)')       '# -------------------------------------------------'
        write(BF%unit,'(A)')       '# Mean shear stress balance, every term divided by tau_w (see VEL_STATS).'
        write(BF%unit,'(A)')       '# y is the distance from the lower wall, at the midpoints between nodes.'
        write(BF%unit,'(A)')       '# The exact total is 1 - y. The first and the last rows are the wall values'
        write(BF%unit,'(A)')       '# set by definition (y = 0 and y = 2), not computed.'
        write(BF%unit,'(A)')       "# <f> = Reynolds mean, f'' = fluctuation about the Favre mean."
        write(BF%unit,'(A)')       '# *************************************************'
        write(BF%unit,'(A)')       '# Column 1  : y'
        write(BF%unit,'(A)')       '# Column 2  : mu_inf <mu_mol du/dy>/tau_w, viscous stress'
        write(BF%unit,'(A)')       '# Column 3  : mu_inf <mu_trb du/dy>/tau_w, SGS stress (0 in DNS)'
        write(BF%unit,'(A)')       "# Column 4  : -<rho> tilde(u'' v'')/tau_w, Reynolds shear stress"
        write(BF%unit,'(A)')       '# Column 5  : total, columns 2 + 3 + 4'
        write(BF%unit,'(A)')       '# *************************************************'
        write(BF%unit,'(20e18.9)') 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 1.0_rp
        do j = sy,ey-1

           yj = 0.5_rp*(y(j)+y(j+1))+1.0_rp
           r_ = 0.5_rp*(vmean1D(j,1)+vmean1D(j+1,1))
           ir = 1.0_rp/r_

           mean_rhuv = 0.5_rp*(vmean1D(j,9)+vmean1D(j+1,9))
           meanR_u_v = 0.5_rp*(vmean1D(j,2)+vmean1D(j+1,2))*&
                       0.5_rp*(vmean1D(j,3)+vmean1D(j+1,3))*ir
           uv        = (mean_rhuv - meanR_u_v)*ir
           R12       = (r_/wl%rho)*uv/wl%u_t**2

           mLUY = mu_inf * 0.5_rp*(vmean1D(j,19)+vmean1D(j+1,19))/(wl%rho*wl%u_t**2)
           mTUY = mu_inf * 0.5_rp*(vmean1D(j,20)+vmean1D(j+1,20))/(wl%rho*wl%u_t**2)

           write(BF%unit,'(20e18.9)') yj, mLUY, mTUY, - R12, &
                                            mLUY + mTUY - R12
        enddo
        write(BF%unit,'(20e18.9)') 2.0_rp, -1.0_rp, 0.0_rp, 0.0_rp, -1.0_rp

        call CloseFile(BF)

        return
end subroutine write_stress_budget

subroutine write_monitor(time,wl,bk,force_turbulent_channel)

        use parameters_module, only: rp, data_dir, it, mu_inf
        use thermodynamics_module, only: sound_speed
        use FileModule

        implicit none
        real(rp), intent(in) :: time
        real(rp), intent(in) :: force_turbulent_channel

        type(wall_type), intent(in) :: wl
        type(bulk_type), intent(in) :: bk
        type(FileType)              :: forceFile
        logical                     :: new_file

        inquire(file = 'DATA/'//trim(data_dir)//'/TCH_MONITOR.txt', exist = new_file)
        new_file = (.not. new_file) .or. it == 0

        forceFile%name = 'TCH_MONITOR'
        forceFile%dir  = trim(data_dir)
        call AppendToFile(forceFile,it)

        if(new_file) then
          write(forceFile%unit,'(A)') '# TURBULENT CHANNEL MONITOR: one row per statistics sample (every StOut iterations).'
          write(forceFile%unit,'(A)') '# CUMULATIVE = from the running means since the statistics (re)start,'
          write(forceFile%unit,'(A)') '# INSTANTANEOUS = from the field at that iteration. Non-dimensional; see the'
          write(forceFile%unit,'(A)') '# VEL_STATS header for tau_w, u_tau and the wall values.'
          write(forceFile%unit,'(A)') '# <u>_V = unweighted volume mean of u, NOT the bulk velocity U_b = <rho u>/<rho>.'
          write(forceFile%unit,'(A)') '# Column 1 : iteration'
          write(forceFile%unit,'(A)') '# Column 2 : t u_tau/h, time in units of h/u_tau (u_tau cumulative)'
          write(forceFile%unit,'(A)') '# Column 3 : Re_tau, CUMULATIVE'
          write(forceFile%unit,'(A)') '# Column 4 : Ma_tau = u_tau/c_w, CUMULATIVE'
          write(forceFile%unit,'(A)') '# Column 5 : streamwise forcing that drives the flow, INSTANTANEOUS'
          write(forceFile%unit,'(A)') '# Column 6 : <u>_V/c_w, INSTANTANEOUS'
          write(forceFile%unit,'(A)') '# Column 7 : 2 <u>_V/mu_inf, Reynolds number on <u>_V, 2h and mu_w, INSTANTANEOUS'
          write(forceFile%unit,'(A)') '# Column 8 : volume mean of T, INSTANTANEOUS'
          write(forceFile%unit,'(A)') '# Column 9 : 2 (u_tau/<u>_V)^2 = tau_w/(rho_w <u>_V^2/2), skin friction on rho_w'
          write(forceFile%unit,'(A)') '#            and <u>_V (u_tau CUMULATIVE, <u>_V INSTANTANEOUS)'
        endif

        write(forceFile%unit,10)      &
             it                      , &
             time*wl%u_t             , &
             wl%ReTau                , &
             wl%MaTau                , &
             force_turbulent_channel , &
             bk%vel/sound_speed(1.0_rp), &
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

        logical, dimension(nv) :: sim
        integer                :: id, j, j_sym

        sim = .true.
        sim( 3) = .false.
        sim( 9) = .false.
        sim(17) = .false.
        sim(18) = .false.
        sim(19) = .false.
        sim(20) = .false.
        sim(22) = .false.
        sim(24) = .false.
        sim(32) = .false.

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

subroutine bulk_mean(f,js,je,mean)
        use parameters_module, only: MPI_RP
        use mpi_module       , only: sx,ex,sz,ez, mpi_comm_world, mpi_sum
        use mesh_module      , only: xstep, ystep, zstep
        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(in)  :: f
        integer                                , intent(in)  :: js, je
        real(rp)                               , intent(out) :: mean
        real(rp) :: lcl(2), gbl(2), s_f, s_v, dv
        integer  :: i, j, k, err

        s_f = 0.0_rp
        s_v = 0.0_rp
        !$acc parallel default(present) copy(s_f,s_v)
        !$acc loop gang, vector collapse(3) reduction(+:s_f,s_v)
        do       k = sz,ez
           do    j = js,je
              do i = sx,ex
                 dv  = xstep(i) * ystep(j) * zstep(k)
                 s_f = s_f + f(i,j,k) * dv
                 s_v = s_v + dv
              enddo
           enddo
        enddo
        !$acc end parallel

        lcl = (/s_f, s_v/)
        call MPI_allreduce(lcl, gbl, 2, MPI_RP, MPI_SUM, mpi_comm_world, err)
        mean = gbl(1)/gbl(2)

        return
end subroutine bulk_mean

end module output_stats_channel
