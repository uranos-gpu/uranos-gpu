module post_stat_output_bl_module
use parameters_module
use mpi_module
use storage_module
use post_stat_storage_module
use mesh_module
implicit none

private
public write_bl_streamwise_stats, write_bl_wall_normal_stats, &
       write_bl_wall_normal_stats_wmles


contains
subroutine write_bl_streamwise_stats

        use Fluid_functions_module, only: Sutherland
        use FileModule

        implicit none
        
        real(rp) :: rcF, T_w, r_w, m_w, ReyTu, delta_nu
        real(rp) :: tauW, utau, udel, uu, ufav099, ufav100, dely, unum, uden, d99
        integer  :: i, j, j99, l, id

        real(rp) :: dudyw, dy

        type(FileType) :: swFile

        swFile%name = 'sw'
        swFile%dir  = trim(data_dir)//'/POST_FRICTION'
        call OpenNewFile(swFile,itStat)

        do i = sx,ex

           rcF = Prandtl**(1.0_rp/3.0_rp)
           T_w = Trat*(1.0_rp + rcF*0.5_rp*(gamma0-1.0_rp)*Mach**2)
           r_w = vmean2D(i,sy,11)/T_w
           m_w = Sutherland(T_w)

           dudyw = (-22._rp*FavreVelU2D(i,1)+36._rp*FavreVelU2D(i,2)-18._rp*FavreVelU2D(i,3)+ 4._rp*FavreVelU2D(i,4))/12._rp
           dy    = (-22._rp*            y(1)+36._rp*            y(2)-18._rp*            y(3)+ 4._rp*            y(4))/12._rp
           dudyw = dudyw/dy

           tauW = mu_inf*m_w*dudyw
           if(wmles) tauW = vmean1D_wmles(i,3)
           utau = sqrt(abs(tauW)/r_w)
        
           !
           ! get the boundary layer thickness
           !
           do id = 1,ny
              if(y(id) > 7.0) exit
           enddo

           udel = 0.99_rp*FavreVelU2D(i,id)
           j99  = 1
           do j=1,ny-1
            uu = FavreVelU2D(i,j)
            if (uu>udel) then
             j99 = j-1
             exit
            endif
           enddo
           ufav099 = FavreVelU2D(i,j99)
           ufav100 = FavreVelU2D(i,j99+1)
           dely = y(j99+1)-y(j99)
           unum = udel-ufav099
           uden = ufav100-ufav099
           d99 = y(j99)+dely*(unum/uden)

           delta_nu = mu_inf*m_w/(r_w*utau)
           ReyTu    = d99/delta_nu

           StreamWiseStats(i,1) = x(i)
           StreamWiseStats(i,2) = tauW/q_inf
           StreamWiseStats(i,3) = delta_nu
           StreamWiseStats(i,4) = ReyTu
           StreamWiseStats(i,5) = utau
           StreamWiseStats(i,6) = vmean2D(i,sy,11)
           StreamWiseStats(i,7) = r_w
           StreamWiseStats(i,8) = d99

           if(wmles) then
             do l = 1,nvWmlesData
                StreamWiseStats(i,21+l) = vmean1D_wmles(i,l)
             enddo
           endif

           write(swFile%unit, '(100e18.9)') StreamWiseStats(i,:)

        enddo

        call CloseFile(swFile)
        
        return
end subroutine write_bl_streamwise_stats




subroutine write_bl_wall_normal_stats

        use FileModule

        implicit none
        
        type(FileType) :: statFile

        real(rp), parameter :: x_spacing = 1.0_rp
        real(rp)            :: xi, r_, ir, iq_inf, d99
        real(rp)            :: utau, delta_nu, r_w, T_, uf, T_w, rcF
        real(rp)            :: R11, R22, R33, R12, ekt
        real(rp)            :: ru, r_y, u_y, muL, muT
        integer             :: i, j, step_pt

        iq_inf = 1.0_rp/(0.5_rp*u_inf**2)
        rcF = Prandtl**(1.0_rp/3.0_rp)
        T_w = Trat*(1.0_rp + rcF*0.5_rp*(gamma0-1.0_rp)*Mach**2)

        xi  = xmin
        step_pt = nint(Nx/(Lx/x_spacing))
        do i = sx,ex, step_pt
        
           xi = x(i)

           statFile%name = 'stat'
           statFile%dir  = trim(data_dir)//'/POST_STATS_'//trim(str(nint(xi)))
           call OpenNewFile(statFile,itStat)

           write(statFile%unit,'(A)')       '# ******************************'
           write(statFile%unit,'(A,f18.6)') '#  TURBULENT BOUNDARY LAYER ', Mach 
           write(statFile%unit,'(A)')       '#'
           write(statFile%unit,'(A)')       '#  author: Francesco De Vanna '
           write(statFile%unit,'(A)')       '#  e-mail: fra.devanna@gmail.com'
           write(statFile%unit,'(A)')       '#'
           write(statFile%unit,'(A,f18.6)') '#  LOCATION%    ', x(i)/Lx*100
           write(statFile%unit,'(A,f18.6)') '#  FRICT. COEFF ', StreamWiseStats(i,2)
           write(statFile%unit,'(A,f18.6)') '#  FRICT. RENUM ', StreamWiseStats(i,4)
           write(statFile%unit,'(A)')       '# ******************************'
           write(statFile%unit,'(A)')       '# Column  1  : y '
           write(statFile%unit,'(A)')       '# Column  2  : ufavre/uinf'
           write(statFile%unit,'(A)')       '# Column  3  : tau11/ekinf'
           write(statFile%unit,'(A)')       '# Column  4  : tau22/ekinf'
           write(statFile%unit,'(A)')       '# Column  5  : tau33/ekinf'
           write(statFile%unit,'(A)')       '# Column  6  : tau12/ekinf'
           write(statFile%unit,'(A)')       '# Column  7  : ek/ekinf'
           write(statFile%unit,'(A)')       '# Column  8  : yplus'
           write(statFile%unit,'(A)')       '# Column  9  : uplus'
           write(statFile%unit,'(A)')       '# Column 10  : tau11plus'
           write(statFile%unit,'(A)')       '# Column 11  : tau22plus'
           write(statFile%unit,'(A)')       '# Column 12  : tau33plus'
           write(statFile%unit,'(A)')       '# Column 13  : tau12plus'
           write(statFile%unit,'(A)')       '# Column 14  : T/T_w'
           write(statFile%unit,'(A)')       '# Column 14  : T/T_w'
           write(statFile%unit,'(A)')       '# Column 15  : rho/rho_w'
           write(statFile%unit,'(A)')       '# Column 16  : u_rms+'
           write(statFile%unit,'(A)')       '# Column 17  : v_rms+'
           write(statFile%unit,'(A)')       '# Column 18  : w_rms+'
           write(statFile%unit,'(A)')       '# Column 18  : uvrms+'
           write(statFile%unit,'(A)')       '# ******************************'

           do j = sy,ey

              r_ = vmean2D(i,j,1)
              ir = 1.0_rp/vmean2D(i,j,1)
              uf = FavreVelU2D(i,j)

              R11 = FavreRUU_2D(i,j)
              R22 = FavreRVV_2D(i,j)
              R33 = FavreRWW_2D(i,j)
              R12 = FavreRUV_2D(i,j)
              ekt = FavreEkT_2D(i,j)
              T_  = vmean2D(i,j,13)
        
              delta_nu = StreamWiseStats(i,3)
              utau     = StreamWiseStats(i,5)
              r_w      = StreamWiseStats(i,7)
              d99      = StreamWiseStats(i,8)

              ! u_y derivative
              r_y = vmean2D(i,j,21)
              ru  = vmean2D(i,j,2)
              u_y = (vmean2D(i,j,24)*r_ - r_y*ru)*ir**2

              muL = vmean2D(i,j,15) ! laminar   viscosity
              muT = vmean2D(i,j,16) ! turbulent viscosity

              write(statFile%unit,'(100e18.9)') &
              y(j)/d99                        , & !1
              uf/u_inf                        , & !2
              r_*R11*iq_inf                   , & !3 
              r_*R22*iq_inf                   , & !4
              r_*R33*iq_inf                   , & !5
              r_*R12*iq_inf                   , & !6
              ekt*iq_inf                      , & !7
              y(j)/delta_nu                   , & !8
              FavreVelU2D(i,j)/utau           , & !9
              r_*R11/(r_w*utau**2)            , & !10
              r_*R22/(r_w*utau**2)            , & !11
              r_*R33/(r_w*utau**2)            , & !12
              r_*R12/(r_w*utau**2)            , & !13
              T_/T_w                          , & !14
              r_/r_w                          , & !15
              FavreUrms2D(i,j)/utau           , & !16
              FavreVrms2D(i,j)/utau           , & !17             
              FavreWrms2D(i,j)/utau           , & !18
              R12/utau**2                     , & !19
              !!! momentum balance variables
              r_*R12                          , & !20
              muL*u_y                         , & !21
              muT*u_y                         


              
        



           enddo

           call CloseFile(statFile)


        enddo

        return
end subroutine write_bl_wall_normal_stats



subroutine write_bl_wall_normal_stats_wmles

        use FileModule
        use wmles_module

        implicit none

        ! local declarations
        real(rp), parameter               :: x_spacing = 1.0_rp
        real(rp), dimension(1:wmles_npts) :: yc
        real(rp), dimension(1:wmles_npts) :: u_wm, T_wm
        real(rp)                          :: xi, iq_inf, T_w, rcF, tauW
        real(rp)                          :: delta_nu, utau
        integer                           :: i
        type(FileType)                    :: stF

        iq_inf = 1.0_rp/(0.5_rp*u_inf**2)
        rcF = Prandtl**(1.0_rp/3.0_rp)
        T_w = Trat*(1.0_rp + rcF*0.5_rp*(gamma0-1.0_rp)*Mach**2)

        xi  = xmin
        do
        
           i  = nint((xi - xmin)*nx/Lx)
           xi = xi + x_spacing
           if(xi > xmax) exit

           stF%name = 'stat'
           stF%dir  = trim(data_dir)//'/POST_WMLES_STATS_'//trim(str(nint(xi)))
           call OpenNewFile(stF,itStat)

           call compute_tau_wall_improved_wmles2D(i,vmean2D,T_w,yc,u_wm,T_wm,tauW)

           delta_nu = StreamWiseStats(i,3)
           utau     = StreamWiseStats(i,5)

           do j = 1,wmles_npts
              write(stF%unit,'(20e18.9)')  &
              yc(j),                       &
              u_wm(j)/u_inf,               &
              T_wm(j)/T_w,                 &
              yc(j)/delta_nu,              &
              u_wm(j)/utau

           enddo


           call CloseFile(stF)


        enddo




        return
end subroutine write_bl_wall_normal_stats_wmles




end module post_stat_output_bl_module
