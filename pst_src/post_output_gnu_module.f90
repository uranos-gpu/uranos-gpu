module post_output_gnu_module
! -----------------------------------------------------------
!
!       This module contains same routine to write inputs
!       for GNUPLOT.
!       
! -----------------------------------------------------------
use storage_module
use post_storage_module
use post_solutions_module
use norm_module
use real_to_integer_module

implicit none
private
public write_gnuplot, write_vs_time

contains
subroutine write_gnuplot
        implicit none

        select case(ic)

          case('periodic_euler_x')
          call write_periodic_euler('x')

          case('periodic_euler_y')
          call write_periodic_euler('y')

          case('periodic_euler_z')
          call write_periodic_euler('z')
                  
          case('shock_tube_x','lax_problem_x', 'shock_wave_interaction_x'                 , &
               'nscbc_perturbation_x', 'shock_vortex', 'shock_wave_x'       , 'shock_inflow')
          call print_primitives('x')

          case('shock_tube_y', 'lax_problem_y', 'shock_wave_interaction_y', 'shock_wave_y', &
               'nscbc_perturbation_y', 'advection_isentropic_vortex_y', 'isentropic_perturbation_y')
          call print_primitives('y')
                          
          case('shock_tube_z', 'lax_problem_z', 'shock_wave_z', 'shock_wave_interaction_z', &
               'advection_isentropic_vortex_z', 'isentropic_perturbation_z')
          call print_primitives('z')

          case('isentropic_vortex_x') 
          call print_isentropic_vortex

          case('couette_x','I_stokes_problem_x', 'II_stokes_problem_x', 'steady_couette_x', 'inflow_poiseuille')
          call print_stokes_problem('x')

          case('couette_y', 'I_stokes_problem_y', 'II_stokes_problem_y', 'steady_couette_y')
          call print_stokes_problem('y')

          case('couette_z', 'I_stokes_problem_z', 'II_stokes_problem_z', 'steady_couette_z')
          call print_stokes_problem('z')

          case('poiseuille_x')
          call print_poiseuille('x')
                  
          case('poiseuille_y')
          call print_poiseuille('y')
                  
          case('poiseuille_z')
          call print_poiseuille('z')

          case('turbulent_channel')
          if(bc(3) == 'aws_isothermal' .or. bc(3) == 'dws_isothermal') then
            call write_turbulent_channel_WallModelled 
          else
            call write_turbulent_channel_statistics
          endif
          call channel_flow_velocity_correlation

          case('swbli')
          call print_swbli_ascii

          case('hTurb')
          call print_hTurb_ascii

          case('turbulent_BL')
              if(hybrid_weno) call WritePDF(SSENSOR)
          !if(bc(3) == 'aws_adiabatic') then
          !  call write_turbulent_bdlayer_WallModelled
          !else
          !  call print_turbulent_boundary_layer_statistics
          !endif

          case('supersonic_intake')
          call write_sic_quantities
          call write_sic_centerline

          case('KolmogorovFlow')
          call print_KolmogorovFlow

        end select

        return
end subroutine write_gnuplot


subroutine write_periodic_euler(dir)
        use FileModule 

        implicit none
        character(1), intent(in) :: dir
        type(FileType)           :: PerEul
        integer                  :: i, j, k
        integer                  :: i0 = 1, j0 = 1, k0 = 1
        
        PerEul%name = trim(output_file_name)
        PerEul%dir  = trim(data_dir)//'/GNUPLOT'

        call OpenNewFile(PerEul,it)
        call gnuplot_error_header(PerEul%unit,error%rho)

        write(PerEul%unit,'(5A,1x)') ' #    coord ', '                density ', '          exact sol. ', &
                '    error', '                 hybrid weno flag'
               
        selectcase(dir)
          !
          ! === X CASE
          !
          case('x')
          if(hybrid_weno) then
            do i = sx,ex
               write(PerEul%unit,fmt=10) &
               x(i), phi(i,j0,k0,1), exact%rho(i,j0,k0), abs(error%rho(i,j0,k0)), weno_flag(i,j0,k0)
            enddo
          else
            do i = sx,ex
               write(PerEul%unit,fmt=11) x(i), phi(i,j0,k0,1), exact%rho(i,j0,k0), abs(error%rho(i,j0,k0))
            enddo
          endif
          !
          ! === Y CASE
          !
          case('y')
          do j = sy,ey
             write(PerEul%unit,fmt=10) y(j), phi(i0,j,k0,1), exact%rho(i0,j,k0), abs(error%rho(i0,j,k0))
          enddo
          !
          ! === Z CASE
          !
          case('z')
          do k = sz,ez
             write(PerEul%unit,fmt=10) z(k), phi(i0,j0,k,1), exact%rho(i0,j0,k), abs(error%rho(i0,j0,k))
          enddo

        endselect
                
        10 format(4e18.9,1i8)
        11 format(4e18.9)

        call CloseFile(PerEul)

        return
end subroutine write_periodic_euler


subroutine print_isentropic_vortex

        use FileModule

        implicit none
        type(FileType) :: VortexFile
        real(rp)       :: rho_, rho_exact
        integer        :: i, j0, k0

        VortexFile%name = trim(output_file_name)
        VortexFile%dir  = trim(data_dir)//'/GNUPLOT'

        call OpenNewFile(VortexFile,it)
        call gnuplot_error_header(VortexFile%unit,error%rho)

        j0 = ey/2
        k0 = ez/2

        write(VortexFile%unit,'(7A,1x)') ' #    x ', '                numerical sol. ','   exact sol.', &
                '        error', '             hybrid weno flag'
        
        if(hybrid_weno) then
          do i = sx,ex

             rho_      = phi(i,j0,k0,1)
             rho_exact = exact%rho(i,j0,k0)

             write(VortexFile%unit,10) x(i), rho_, rho_exact, abs(rho_-rho_exact), weno_flag(i,j0,k0)
          enddo
        else
          do i = sx,ex

             rho_      = phi(i,j0,k0,1)
             rho_exact = exact%rho(i,j0,k0)

             write(VortexFile%unit,11) x(i), rho_, rho_exact, abs(rho_-rho_exact)
          enddo
        endif

        call CloseFile(VortexFile)

        10 format(4e18.6,1i7)
        11 format(4e18.6)

        return
end subroutine print_isentropic_vortex

subroutine print_primitives(dir)
        use FileModule

        implicit none
        character(1), intent(in) :: dir
        type(FileType)           :: PrimFile
        integer                  :: i0 = 1, j0 = 1, k0 = 1
        integer                  :: i,j,k
        real(rp)                 :: rho_, i_rho, u_, v_, w_, p_, T_, c

        PrimFile%name = trim(output_file_name)
        PrimFile%dir  = trim(data_dir)//'/GNUPLOT'
        
        ! open the file
        call OpenNewFile(PrimFile,it)
        ! set the Gnuplot Header
        call GnuplotHeader(Primfile,sx,ex,sy,ey,sz,ez,time,ic)

        write(PrimFile%unit,'(8A,1x)') ' #    coord ', '                density ', '          pressure ', &
        '         temperature ', '      u speed ', '          speed of sound ', '   Mach', &
        '              Hybrid weno flag'
        
        selectcase(dir)
          !
          ! === X CASE
          !
          case('x')
          if(hybrid_weno) then
            do i=sx,ex

                rho_  = phi(i,j0,k0,1)
                i_rho = 1.0_rp/rho_

                u_    = phi(i,j0,k0,2)*i_rho
                v_    = phi(i,j0,k0,3)*i_rho
                w_    = phi(i,j0,k0,4)*i_rho
                p_    = (gamma0-1._rp)*(phi(i,j0,k0,5)-0.5_rp*rho_*(u_**2 + v_**2 + w_**2))
                T_    = p_*i_rho
                c     = sqrt(gamma0*T_)

                write(PrimFile%unit,10) x(i), rho_, p_, t_, u_, c, u_/c, weno_flag(i,j0,k0)

            enddo
          else
            do i=sx,ex

                rho_  = phi(i,j0,k0,1)
                i_rho = 1.0_rp/rho_

                u_    = phi(i,j0,k0,2)*i_rho
                v_    = phi(i,j0,k0,3)*i_rho
                w_    = phi(i,j0,k0,4)*i_rho
                p_    = (gamma0-1._rp)*(phi(i,j0,k0,5)-0.5_rp*rho_*(u_**2 + v_**2 + w_**2))
                T_    = p_*i_rho
                c     = sqrt(gamma0*T_)
                   
                write(PrimFile%unit,11) x(i), rho_, p_, t_, u_, c, u_/c

            enddo
          endif
          !
          ! === Y CASE
          !
          case('y')
          do j=sy,ey

              rho_  = phi(i0,j,k0,1)
              i_rho = 1.0_rp/rho_

              u_    = phi(i0,j,k0,2)*i_rho
              v_    = phi(i0,j,k0,3)*i_rho
              w_    = phi(i0,j,k0,4)*i_rho
              p_    = (gamma0-1._rp)*(phi(i0,j,k0,5)-0.5_rp*rho_*(u_**2 + v_**2 + w_**2))
              T_    = p_*i_rho
              c     = sqrt(gamma0*T_)
                  
              write(PrimFile%unit,11) y(j), rho_, p_, t_, v_, c, v_/c
          enddo
          !
          ! === Z CASE
          !
          case('z')
          do k=sz,ez

              rho_  = phi(i0,j0,k,1)
              i_rho = 1.0_rp/rho_

              u_    = phi(i0,j0,k,2)*i_rho
              v_    = phi(i0,j0,k,3)*i_rho
              w_    = phi(i0,j0,k,4)*i_rho
              p_    = (gamma0-1._rp)*(phi(i0,j0,k,5)-0.5_rp*rho_*(u_**2 + v_**2 + w_**2))
              T_    = p_*i_rho
              c     = sqrt(gamma0*T_)
                  
              write(PrimFile%unit,11) z(k), rho_, p_, t_, w_, c, w_/c
          enddo

        endselect

        call CloseFile(PrimFile)

        10 format(7e18.6,1i7)
        11 format(8e18.6)


        return
end subroutine print_primitives


subroutine print_poiseuille(dir)

        use FileModule

        implicit none
        character(1), intent(in) :: dir
        type(FileType)           :: PoisFile
        integer                  :: i0 = 1, j0 = 1, k0 = 1
        integer                  :: i,j,k,l
        real(rp)                 :: u_, u_exact
        real(rp)                 :: v_, v_exact
        real(rp)                 :: t_, t_exact, yjh, tauW = 0.0_rp
        real(rp)                 :: u_Rey, u_y, mu_, mVis_l, u_l, cl1, cl2, idyR
        
        PoisFile%name = trim(output_file_name)
        PoisFile%dir  = trim(data_dir)//'/GNUPLOT'

        ! open the file
        call OpenNewFile(PoisFile,it)

        write(PoisFile%unit,*) '# Velocity error'
        call gnuplot_error_header(PoisFile%unit,error%vel)
        write(PoisFile%unit,*) '# Temperatur error'
        call gnuplot_error_header(PoisFile%unit,error%tmp)

        write(PoisFile%unit,'(7A,1x)') ' #    coord ', '                u ','                exact u', &
                 '           T', '                 exact T'
        
        selectcase(dir)
          !
          ! === X CASE
          !
          case('x')
          if(Rey_average) then
            do j = sy,ey

               u_      = phi(i0,j,k0,2)/phi(i0,j,k0,1)
               u_exact = exact%vel(i0,j,k0)
               u_rey   = Re_av%rhu(i0,j,k0)/Re_av%rho(i0,j,k0)

               write(PoisFile%unit,10) y(j), u_/u_inf, u_exact/u_inf, u_Rey/u_inf

            enddo
          else
            do j = sy-1,ey+1

               u_      = phi(i0,j,k0,2)/phi(i0,j,k0,1)
               u_exact = exact%vel(i0,j,k0)

               t_      = T(i0,j,k0)
               t_exact = exact%tmp(i0,j,k0)

               idyR = 1.0_rp/(0.5_rp*(ystep(j) + ystep(j+1)))
               yjh  = 0.5_rp*(y(j) + y(j+1))

               u_y = 0.0_rp
               mu_ = 0.0_rp
               do l = fd_L,fd_R

                  cl1  = central_1_one_half(l)
                  cl2  = central_2_one_half(l)
        
                  ! u_y derivative
                  u_l = U(i0,j+l,k0)
                  u_y = u_y + cl2 * u_l * idyR

                  ! viscous components
                  mVis_l = 1.0_rp + 0.5_rp*cos(pi*y(j+l))
                  !mVis_l = 1.0_rp
                  !if(y(j+l) > -0.5_rp .and. y(j+l) < 0.5_rp) mVis_l = 4.0_rp

                  mu_ = mu_ + cl1 * mVis_l

               enddo
               if(j == sy-1) tauW = mu_*u_y

               write(PoisFile%unit,*) y(j), u_/u_inf, u_exact/u_inf, t_, t_exact, &
                       u_y/(2*u_inf), yjh, mu_*u_y/tauW, 0.5_rp*(y(j) + y(j+1))

            enddo
          endif
          !
          ! === Y CASE
          !
          case('y')
          do i = sx-1,ex+1

             v_      = phi(i,j0,k0,3)/phi(i,j0,k0,1)
             v_exact = exact%vel(i,j0,k0)

             t_      = T(i,j0,k0)
             t_exact = exact%tmp(i,j0,k0)

             write(PoisFile%unit,10) x(i), v_/u_inf, v_exact/u_inf, t_, t_exact

          enddo
          !
          ! === Z case
          !
          case('z')
          do k = sz-1,ez+1

             u_      = phi(i0,j0,k,2)/phi(i0,j0,k,1)
             u_exact = exact%vel(i0,j0,k)

             t_      = T(i0,j0,k)
             t_exact = exact%tmp(i0,j0,k)

             write(PoisFile%unit,10) z(k), u_/u_inf, u_exact/u_inf, t_, t_exact

          enddo

        endselect

        call CloseFile(PoisFile)

        10 format(7e18.9)

        return
end subroutine print_poiseuille


subroutine print_stokes_problem(dir)

        use FileModule

        implicit none
        character(1), intent(in) :: dir
        type(FileType)           :: StokesFile
        real(rp)                 :: u_, u_exact
        real(rp)                 :: v_, v_exact
        integer                  :: i0 = 1, j0 = 1, k0 = 1
        integer                  :: i,j,k,l 

        real(rp) :: u_y, mu_, cl1, cl2, u_l, mVis_l, idyR,tauMax
        real(rp), parameter :: visPar = 4.0_rp


        StokesFile%name = trim(output_file_name)
        StokesFile%dir  = trim(data_dir)//'/GNUPLOT'

        call OpenNewFile(StokesFile,it)
        call gnuplot_error_header(StokesFile%unit,error%vel)

        write(StokesFile%unit,'(7A,1x)') ' #    coord ', '                numerical sol. ','   exact sol.', &
                '        error'
        
        selectcase(dir)
          !
          ! === X CASE
          !
          case('x')
          tauMax = 0.0_rp
          do j = sy-1,ey+1

             u_      = phi(i0,j,k0,2)/phi(i0,j,k0,1)
             u_exact = exact%vel(i0,j,k0)

             idyR = 1.0_rp/(0.5_rp*(ystep(j) + ystep(j+1)))
             u_y = 0.0_rp
             mu_ = 0.0_rp
             do l = fd_L,fd_R

                cl1  = central_1_one_half(l)
                cl2  = central_2_one_half(l)
        
                ! u_y derivative
                u_l = U(i0,j+l,k0)
                u_y = u_y + cl2 * u_l * idyR

                ! viscous components
                mVis_l = 1.0_rp
                if(y(j+l) < 0.0_rp) mVis_l = visPar

                mu_ = mu_ + cl1 * mVis_l

             enddo

             tauMax = max(abs(mu_*u_y),tauMax)
          enddo

          do j = sy-1,ey+1

             u_      = phi(i0,j,k0,2)/phi(i0,j,k0,1)
             u_exact = exact%vel(i0,j,k0)

             idyR = 1.0_rp/(0.5_rp*(ystep(j) + ystep(j+1)))
             u_y = 0.0_rp
             mu_ = 0.0_rp
             do l = fd_L,fd_R

                cl1  = central_1_one_half(l)
                cl2  = central_2_one_half(l)
        
                ! u_y derivative
                u_l = U(i0,j+l,k0)
                u_y = u_y + cl2 * u_l * idyR

                ! viscous components
                mVis_l = 1.0_rp
                if(y(j+l) < 0.0_rp) mVis_l = visPar

                mu_ = mu_ + cl1 * mVis_l

             enddo

             write(StokesFile%unit,10) y(j), u_/u_inf, u_exact/u_inf, &
              abs(u_exact-u_)/u_inf, T(i0,j,k0), P(i0,j,k0)         , &
              0.5_rp*(y(j) + y(j+1)), mu_*u_y/tauMax

          enddo
          !
          ! === Y CASE
          !
          case('y')
          do i = sx-1,ex+1

             v_      = phi(i,j0,k0,3)/phi(i,j0,k0,1)
             v_exact = exact%vel(i,j0,k0)

             write(StokesFile%unit,10) x(i), v_/u_inf, v_exact/u_inf, abs(v_exact-v_)/u_inf

          enddo
          !
          ! === Z CASE
          !
          case('z')
          do k = sz-1,ez+1

             u_      = phi(i0,j0,k,2)/phi(i0,j0,k,1)
             u_exact = exact%vel(i0,j0,k)

             write(StokesFile%unit,10) z(k), u_/u_inf, u_exact/u_inf, abs(u_exact-u_)/u_inf

          enddo

        end select

        call CloseFile(StokesFile)

        10 format(8e18.6)

        return
end subroutine print_stokes_problem


subroutine print_KolmogorovFlow

        use FileModule

        implicit none
        type(FileType) :: kFile
        real(rp)       :: u_, ue
        integer        :: j, i0 = 1, k0 = 1

        kFile%name = trim(output_file_name)
        kFile%dir  = trim(data_dir)//'/GNUPLOT'

        call OpenNewFile(kFile,it)
        call gnuplot_error_header(kFile%unit,error%vel)

        write(kFile%unit,'(7A,1x)') ' #    coord ', '                numerical sol. ','   exact sol.', &
                '        error'

        do j = sy-1,ey+1

           u_ = phi(i0,j,k0,2)/phi(i0,j,k0,1)
           ue = exact%vel(i0,j,k0)

           write(kFile%unit,'(4e18.9)') y(j), u_/u_inf, ue/u_inf, abs(u_-ue)

        enddo


        return
end subroutine print_KolmogorovFlow


subroutine write_turbulent_channel_statistics

        use fluid_functions_module, only: iLogLaw, laminar_viscosity
        use math_tools_module     , only: newton_raphson
        use FileModule

        implicit none
        type(FileType)            :: PostVelStats
        integer , parameter       :: imax = 100
        real(rp), parameter       :: toll = 1.0E-14_rp
        real(rp), dimension(0:ny) :: yPlus, uPlus, uPlusVD
        real(rp)                  :: delta, rWall, uWall, wWall
        real(rp)                  :: u_yWl, tmpWl, mWall, nWall, tauWl
        real(rp)                  :: du, uc, r_, ir, d2, u2, w2
        real(rp)                  :: mean_rhuu, mean_rhvv, mean_rhww, mean_rhuv
        real(rp)                  :: meanR_u_u, meanR_v_v, meanR_w_w, meanR_u_v
        real(rp)                  :: uu, vv, ww, uv, R11, R22, R33, R12
        real(rp)                  :: u_tau, ReyTu, utModel, uP, yj
        integer                   :: i0, k0

        PostVelStats%dir  = trim(data_dir)//'/POST_VEL_STATS'
        PostVelStats%name = 'stat'
        call OpenNewFile(PostVelStats,it)

        i0 = ex/2
        k0 = ez/2

        ! compute shear quantities
        delta = 0.5_rp*(ymax-ymin)
        rWall = 1.0_rp
        uWall = Re_av%rhu(i0,sy,k0)/rWall
        wWall = Re_av%rhw(i0,sy,k0)/rWall
        u_yWl = 2*uWall/(y(sy) - y(sy-1))
        tmpWl = 1.0_rp
        mWall = laminar_viscosity(tmpWl,Tref,vis_flag)
        nWall = mu_inf*mWall/rWall
        tauWl = mu_inf*mWall*u_yWl
        u_tau = sqrt(tauWl/rWall)
        
        if(bc(3) == 'aws_isothermal') then
          d2 = y(sy+1) + 1.0_rp
          u2 = Re_av%rhu(i0,sy+1,k0)/Re_av%rho(i0,sy+1,k0)
          w2 = Re_av%rhw(i0,sy+1,k0)/Re_av%rho(i0,sy+1,k0)
          uP = sqrt(u2**2 + w2**2)
          call newton_raphson(iLogLaw,d2,uP,nWall,imax,toll,u_tau,utModel)
          u_tau = utModel
        endif
        ReyTu = rWall*u_tau/(mu_inf*mWall)

        ! compute uPlus and yPlus
        uPlus(0) = 0.0_rp
        yPlus(0) = 0.0_rp
        do j = sy, ey
           yj = y(j)+delta
           yPlus(j) = yj*ReyTu
           uPlus(j) = Re_av%rhu(i0,j,k0)/Re_av%rho(i0,j,k0)/u_tau
        enddo
        
        ! transform with VAN DRIEST
        uPlusVD = 0.0_rp
        do j = sy,ey
           du = uPlus(j)-uPlus(j-1)
           uc = sqrt(Re_av%rho(i0,j,k0)/rWall)
           uPlusVD(j) = uPlusVD(j-1)+uc*du
        enddo

        do j = sy, ey

           yj = y(j)+delta
           r_ = Re_av%rho(i0,j,k0)
           ir = 1.0_rp/r_

           ! mean(rho*ui*uj)
           mean_rhuu = Re_av%ruu(i0,j,k0)
           mean_rhvv = Re_av%rvv(i0,j,k0)
           mean_rhww = Re_av%rww(i0,j,k0)
           mean_rhuv = Re_av%ruv(i0,j,k0)
        
           ! mean(rhui)*mean(rhuj)
           meanR_u_u = Re_av%rhu(i0,j,k0)*Re_av%rhu(i0,j,k0)*ir
           meanR_v_v = Re_av%rhv(i0,j,k0)*Re_av%rhv(i0,j,k0)*ir
           meanR_w_w = Re_av%rhw(i0,j,k0)*Re_av%rhw(i0,j,k0)*ir
           meanR_u_v = Re_av%rhu(i0,j,k0)*Re_av%rhv(i0,j,k0)*ir
        
           ! u_i'*u_j'
           uu = (mean_rhuu - meanR_u_u)*ir
           vv = (mean_rhvv - meanR_v_v)*ir
           ww = (mean_rhww - meanR_w_w)*ir
           uv = (mean_rhuv - meanR_u_v)*ir

           R11 = (r_/rWall)*uu/u_tau**2
           R22 = (r_/rWall)*vv/u_tau**2
           R33 = (r_/rWall)*ww/u_tau**2
           R12 = (r_/rWall)*uv/u_tau**2

           write(PostvelStats%unit,'(20e18.9)') &
                 yj, yPlus(j), uPlus(j), uPlusVD(j), R11, R22, R33, R12, r_/rWall

        enddo
        call CloseFile(postVelStats)

        return
end subroutine write_turbulent_channel_statistics

subroutine write_vs_time

        use input_module, only: get_unit

        implicit none
        integer       :: f_unit
        character(dl) :: filename
        real(rp)      :: du_max, dv_max, dp_max

        call get_unit(f_unit)
        selectcase(ic)
          case('pirozzoli_vortex', 'isentropic_vortex_x', 'pirozzoli_vortex_y')

          if(.not.velocity) stop 'velocity is a presequite in write_vs_time: Pirozzoli_vortex'
          if(.not.pressure) stop 'pressure is a presequite in write_vs_time: Pirozzoli_vortex'

          filename = "DATA/"//trim(data_dir)//"/GNUPLOT/"//"quantities_vs_time.txt"
          if(f == file_ini) then

            open(unit=f_unit,file=filename,status='replace')
            write(f_unit,'(4A)') 'time', 'max|u-u_inf|/u_inf', 'max|v|/u_inf', 'max|p-1|'

          else
            open(unit=f_unit,file=filename,status='unknown',position='append')
          endif

          du_max = maxval(abs(U(:,:,:)-u_inf))/u_inf
          dv_max = maxval(abs(V(:,:,:)-0.0_rp))  /u_inf
          dp_max = maxval(abs(P(:,:,:)-1.0_rp))

          write(f_unit,'(4e18.6)') time, du_max, dv_max, dp_max

        end select

        close(f_unit)
        return
endsubroutine write_vs_time


subroutine gnuplot_error_header(funit,error)
        implicit none
        real(rp), allocatable, dimension(:,:,:), intent(in) :: error
        integer                                , intent(in) :: funit

        ! local declarations
        integer, dimension(3) :: lb,ub
        real(rp)              :: err_norm1, err_norm2, err_normi

        write(funit,'(A)'      ) " # ================================= #"
        write(funit,'(A,A)'    ) " # test     = ", trim(ic)
        write(funit,'(A,e18.6)') " # time     = ", time
        write(funit,'(A,I3)'   ) " # x_points = ", ex-sx+1
        write(funit,'(A,I3)'   ) " # y_points = ", ey-sy+1
        write(funit,'(A,I3)'   ) " # z_points = ", ez-sz+1
        write(funit,'(A)'      ) " # ================================= #"
        

        lb = (/sx,sy,sz/)
        ub = (/ex,ey,ez/)

        ! compute error in norm 1,2, infinit
        err_norm1 = p_norm  (error,1._rp,lb,ub)
        err_norm2 = p_norm  (error,2._rp,lb,ub)
        err_normi = inf_norm(error   ,lb,ub)

        write(funit, '(A,e16.6)') ' # dx   = ', x(sx+1) - x(sx)
        write(funit, '(A,e16.6)') ' # dy   = ', y(sy+1) - x(sy)
        write(funit, '(A,e16.6)') ' # dz   = ', z(sz+1) - x(sz)
        write(funit, '(A,e16.6)') ' # error norm - 1   = ', err_norm1
        write(funit, '(A,e16.6)') ' # error norm - 2   = ', err_norm2
        write(funit, '(A,e16.6)') ' # error norm - inf = ', err_normi
        write(funit, '(A)'      ) " # ================================= #"

        return
end subroutine gnuplot_error_header


subroutine print_swbli_ascii

        use FileModule
        use real_to_integer_module, only: nearest_integer_opt
        use storage_module        , only: U, P
        use fluid_functions_module, only: laminar_viscosity

        implicit none
        type(FileType)      :: prsFile, velFile, CfFile, gnuFile
        real(rp), parameter :: xsh_min = 0.20_rp
        real(rp), parameter :: xsh_max = 2.00_rp
        real(rp), parameter :: x1      = 0.60_rp
        real(rp), parameter :: x2      = 0.95_rp
        real(rp), parameter :: x3      = 1.60_rp
        real(rp)            :: xi, yj, p_, p0, i_u_inf
        real(rp)            :: i_dy, u_y, cf, mu_w
        integer             :: i, j, i1, i2, i3, fd, s
        !
        ! === pressure and density in some positions
        !
        gnuFile%name = 'wallNormal'
        gnuFile%dir  = trim(data_dir)//'/GNUPLOT' 
        call OpenNewFile(gnuFile,it)

        i1 = nearest_integer_opt(x,sx,ex,1.0_rp/3.0_rp*Lx)
        i2 = nearest_integer_opt(x,sx,ex,1.0_rp/2.0_rp*Lx)
        i3 = nearest_integer_opt(x,sx,ex,2.0_rp/3.0_rp*Lx)

        do j = sy-1,ey+1
           write(gnuFile%unit,*) y(j), P(i1,j,1), phi(i1,j,1,1), T(i1,j,1), &
                                       P(i2,j,1), phi(i2,j,1,1), T(i2,j,1), &
                                       P(i3,j,1), phi(i3,j,1,1), T(i3,j,1)
        enddo

        call CloseFile(gnuFile)


        !
        ! === stream wise pressure file
        !
        prsFile%name = 'pressure'
        prsFile%dir  = trim(data_dir)//'/PRESSURE'

        call OpenNewFile(prsFile,it)
        
        write(prsFile%unit,*) '# x ', ' pressure(x)'
        
        ! compute the minimum pressure along the plate
        p0 = huge(1.0_rp)
        do i = sx,ex
           xi = x(i)/abs(xmax)+1.0_rp
           if(xi > xsh_min .and. xi < xsh_max) then
           p0 = min(P(i,1,1),p0)
           endif
        enddo

        ! write to file
        do i = sx, ex
           xi = x(i)/abs(xmax)+1.0_rp
           if(xi > xsh_min .and. xi < xsh_max) then
              p_ = p(i,1,1)
              write(prsFile%unit, '(2e18.6)') xi, p_/p0
           endif
        enddo
        call CloseFile(prsFile) 
        !
        ! === stream-wise velocity in x1 = 0.6_rp, x2 = 0.95_rp, x3 = 1.6_rp
        !
        velFile%name = 'velocity'
        velFile%dir  = trim(data_dir)//'/VELOCITY'

        call OpenNewFile(velFile,it)

        write(velFile%unit,*) ' # y ', ' u(x1,y) ', ' U(x2,y) ', ' U(x3,y)'  

        i1 = nearest_integer_opt(x,sx,ex,x1)
        i2 = nearest_integer_opt(x,sx,ex,x2)
        i3 = nearest_integer_opt(x,sx,ex,x3)
        
        i_u_inf = 1.0_rp/u_inf
        do j = lby, uby
           yj = y(j)

           if(yj < 0.05_rp) then
             write(velFile%unit,'(4e18.6)') yj, U(i1,j,1)*i_u_inf, U(i2,j,1)*i_u_inf, U(i3,j,1)*i_u_inf
           endif

        enddo
        call CloseFile(velFile) 
        !
        ! === compute friction coefficient at wall location
        !
        CfFile%name = 'friction'
        CfFile%dir  = trim(data_dir)//'/FRICTION'

        call OpenNewFile(CfFile,it)
        
        i_dy = ystep_i(1)
        fd   = central_fd_order/2
        do i = sx,ex
           xi = x(i)/abs(xmax)+1.0_rp
            
             u_y = 0.0_rp
             do s = -central_fd_order/2, central_fd_order/2
                u_y = u_y + i_dy*central_1(s)*(U(i,1+s,1))
             enddo
                
             mu_w = 0.5_rp*(laminar_viscosity(T(i,1,1),Tref,vis_flag) + laminar_viscosity(T(i,0,1),Tref,vis_flag))
             cf = mu_inf * mu_w * u_y/q_inf

             write(CfFile%unit,'(3e18.6)') xi, cf, 0.644_rp/sqrt(Reynolds*xi)

        enddo
        call CloseFile(CfFile) 


        return
end subroutine print_swbli_ascii





subroutine print_hTurb_ascii

        use integration_module
        use post_statistic_time_module
        use FileModule
        use fluid_functions_module

        implicit none
        type(FileType) :: TrbMachFile
        type(FileType) :: TrbEpsilonFile
        type(FileType) :: TrbLengthsFile

        real(rp), dimension(3) :: u_rms
        real(rp), dimension(1) :: rh_mean
        real(rp)               :: urms, Mt, Re_lambda
        real(rp)               :: tp_mean, cs_mean, mu_mean, nu_mean
        real(rp)               :: eta_scale, lambda_scale
        real(rp)               :: eps_sol, eps_dil, eps_tot, eps_str
        
        if(.not.velocity)    stop ' Turbulent Mach number requires velocity = .true.'
        if(.not.temperature) stop ' Turbulent Mach number requires temperature = .true.'

        !
        ! === open output files
        !
        TrbMachFile%name    = 'TrbMach'
        TrbEpsilonFile%name = 'TrbEpsilon'
        TrbLengthsFile%name = 'TrbLengths'
        TrbMachFile%dir     = trim(data_dir)//'/TURBULENT_PARAMETERS'
        TrbEpsilonFile%dir  = trim(data_dir)//'/TURBULENT_PARAMETERS'
        TrbLengthsFile%dir  = trim(data_dir)//'/TURBULENT_PARAMETERS'

        call AppendToFile(TrbMachFile   ,it)
        call AppendToFile(TrbEpsilonFile,it)
        call AppendToFile(TrbLengthsFile,it)
        !
        ! === compute statistical quantities
        !
        call rms_field(U,sx,ex,xstep,sy,ey,ystep,sz,ez,zstep,u_rms(1))
        call rms_field(V,sx,ex,xstep,sy,ey,ystep,sz,ez,zstep,u_rms(2))
        call rms_field(W,sx,ex,xstep,sy,ey,ystep,sz,ez,zstep,u_rms(3))

        call mean_field(T      ,sx,ex,xstep,sy,ey,ystep,sz,ez,zstep,tp_mean)
        call mean_field(PHI    ,sx,ex,xstep,sy,ey,ystep,sz,ez,zstep,1,1,rh_mean)
        call mean_field(VIS,sx,ex,xstep,sy,ey,ystep,sz,ez,zstep,mu_mean)
        !
        ! === summing statistical parameters
        !
        urms    = sqrt((1.0_rp/3.0_rp)*(u_rms(1)**2 + u_rms(2)**2 + u_rms(3)**2))
        cs_mean = sqrt(gamma0*Tp_mean)
        nu_mean = mu_inf*mu_mean/rh_mean(1)
        Mt      = urms/cs_mean
        !
        ! === compute turbulent dissipation
        !
        call hTurbEpsilon(xstep,ystep,zstep,phi,U,V,W,VIS,rh_mean(1),&
                          eps_sol,eps_dil,eps_tot,eps_str,mpi_flag)
        
        ! Kolmogorov lenght scale
        eta_scale = (nu_mean**3/eps_tot)**(0.25_rp)
        ! Taylor lenght scale
        lambda_scale = urms*sqrt(15.0_rp*nu_mean/eps_tot)
        ! Reynolds Taylor
        Re_lambda = urms * lambda_scale/(nu_mean)

        ! print to files
        write(TrbMachFile%unit   ,'(3e18.6)') time, Mt, Re_lambda, urms
        write(TrbEpsilonFile%unit,'(4e18.6)') time, eps_dil/eps_sol, eps_tot, eps_str
        write(TrbLengthsFile%unit,'(4e18.6)') time, eta_scale, lambda_scale, lambda_scale/eta_scale

        ! close files
        call CloseFile(TrbMachFile)
        call CloseFile(TrbEpsilonFile)
        call CloseFile(TrbLengthsFile)

        return
end subroutine print_hTurb_ascii









subroutine print_turbulent_boundary_layer_statistics

        use FileModule
        use fluid_functions_module, only: laminar_viscosity
        use real_to_integer_module
        use interpolation_module

!        implicit none
!        type(fileType)      :: WallQuantities
!        type(fileType)      :: velFile, thdFile, rssFile
!        real(rp), parameter :: ivk = 1.0_rp/0.41_rp, C = 5.2_rp
!
!        real(rp), dimension(3)     :: xs
!        real(rp), dimension(0:ny)  :: uPlus, uPlusVD, yPlus, utemp
!        real(rp)                   :: tmpWl, mWall, rWall, u_yWl, tauWl, ReTWl, idvu
!        real(rp)                   :: u_tau, cf_Wl, r
!        real(rp)                   :: rhoY, tmpY
!        real(rp)                   :: du, uc, log_law, d99, u99, d99err
!        real(rp)                   :: r_, u_, v_, w_, uu, vv, ww, uv
!        real(rp)                   :: R11, R22, R33, R12
!        integer                    :: is, s, k0, j99, jj99, m = 4
!        
!        k0 = ez/2
!
!        ! Mean spanwise Friction
!        WallQuantities%name = 'WallQuantities'
!        WallQuantities%dir  = trim(data_dir)//'/POST_WALL_QUANTITIES'
!        call OpenNewFile(WallQuantities,it)
!
!        tmpWl = 1.0_rp + 0.5_rp*(gamma0-1.0_rp)*Prandtl**(1.0_rp/3.0_rp)*Mach**2
!        mWall = laminar_viscosity(tmpWl)
!        uTemp = 0.0_rp
!
!        do i = lbx,ubx
!        
!           rWall = Re_av%r(i,sy,k0)
!
!           u_yWl = Re_av%u(i,sy,k0)/y(sy)
!
!           tauWl = mu_inf*mWall*u_yWl
!        
!           u_tau = sqrt(tauWl/rWall)
!
!           cf_Wl = tauWl/q_inf
!
!           ! local delta99
!           uTemp(sy:ey) = Re_av%u(i,sy:ey,k0)
!           u99  = 0.99_rp*u_inf
!           call locate(uTemp,0,ny,u99,j99)
!           jj99 = min(max(j99-(m-1)/2,1),ny+1-m)
!           call polint(uTemp(jj99),y(jj99),m,u99,d99,d99err)
!           ReTWl = rWall*u_tau*d99/(mu_inf*mWall)
!
!           write(wallQuantities%unit,*) x(i), cf_wl, u_tau, ReTWl, rWall
!
!        enddo
!        call CloseFile(wallQuantities)
!
!
!        xs    = Lx*(/0.3_rp, 0.5_rp, 0.7_rp/)
!        uTemp = 0.0_rp
!
!        StLoop: do s = 1,size(xs)
!        
!           ! open Files
!           velFile%name = 'vel_'
!           velFile%dir  = trim(data_dir)//'/POST_VEL_ST'//trim(str(s))
!           call OpenNewFile(velFile,it)
!
!           thdFile%name = 'thd_'
!           thdFile%dir  = trim(data_dir)//'/POST_THD_ST'//trim(str(s))
!           call OpenNewFile(thdFile,it)
!
!           rssFile%name = 'rss_'
!           rssFile%dir  = trim(data_dir)//'/POST_RSS_ST'//trim(str(s))
!           call OpenNewFile(rssFile,it)
!        
!           ! find nearest location of the xs
!           is = nearest_integer_opt(x,sx,ex,xs(s))
!        
!           ! compute wall properties
!           r     = Prandtl**(1.0_rp/3.0_rp)
!           tmpWl = 1.0_rp + r*0.5_rp*(gamma0-1.0_rp)*Mach**2
!           mWall = laminar_viscosity(tmpWl)
!           rWall = Re_av%r(is,sy,k0)
!           u_yWl = Re_av%u(is,sy,k0)/y(sy)
!           tauWl = mu_inf*mWall*u_yWl
!           u_tau = sqrt(tauWl/rWall)
!
!           ! compute the local shear Reynolds number
!           idvu = rWall*u_tau/(mu_inf*mWall)
!
!           ! compute uPlu and yPlus
!           uPlus = 0.0_rp
!           yPlus = 0.0_rp
!           do j = sy,ey
!              yPlus(j) = y(j)*idvu
!              uPlus(j) = Re_av%u(is,j,k0)/u_tau
!           enddo
!        
!           ! transform with VAN DRIEST
!           uPlusVD = 0.0_rp
!           do j = sy,ey
!              du = uPlus(j)-uPlus(j-1)
!              uc = sqrt(Re_av%r(is,j,k0)/rWall)
!              uPlusVD(j) = uPlusVD(j-1)+uc*du
!           enddo
!
!           ! write velFile
!           do j = sy,ey
!              log_law = ivk * log(yPlus(j)) + C
!              write(velFile%unit,*) y(j), yPlus(j), uPlusVD(j), uPlus(j), log_law
!           enddo
!
!           ! write thdFile
!           do j = sy,ey
!              rhoY = Re_av%r(is,j,k0)/rWall
!              tmpY = Re_av%T(is,j,k0)/tmpWl
!              write(thdFile%unit,*) y(j), yPlus(j), tmpY, rhoY
!           enddo
!
!           ! write Reynolds Stress
!           do j = sy,ey
!
!              u_ = Re_av%u(is,j,k0)
!              v_ = Re_av%v(is,j,k0)
!              w_ = Re_av%w(is,j,k0)
!              r_ = Re_av%r(is,j,k0)
!
!              uu = Re_av%uu(is,j,k0)
!              vv = Re_av%vv(is,j,k0)
!              ww = Re_av%ww(is,j,k0)
!              uv = Re_av%uv(is,j,k0)
!
!              R11 = r_/rWall*(uu-u_*u_)/(u_tau**2)
!              R22 = r_/rWall*(vv-v_*v_)/(u_tau**2)
!              R33 = r_/rWall*(ww-w_*w_)/(u_tau**2)
!              R12 = r_/rWall*(uv-u_*v_)/(u_tau**2)
!              
!              write(rssFile%unit,*) y(j), yPlus(j), R11, R22, R33, R12
!
!           enddo
!
!           call CloseFile(velFile)
!           call CloseFile(thdFile)
!           call CloseFile(rssFile)
!
!
!        enddo StLoop

        return
end subroutine print_turbulent_boundary_layer_statistics



subroutine write_turbulent_channel_WallModelled
        
        use FileModule
        use fluid_functions_module, only: laminar_viscosity
        implicit none
        
        type(FileType) :: VelFile, BudFile
                
        real(rp), dimension(0:ny+1) :: molVis, trbVis, u_y
        real(rp), dimension(0:ny+1) :: yPlus , uPlus, uPlusVD
        integer                     :: i0 = 1, k0 = 1, l, j

        real(rp) :: cl1, cl2, u_l, idyR, yHalf, tWl
        real(rp) :: T_l, mVis_l, tVis_l, mMolUY, mTrbUY
        real(rp) :: rHalf, irhlf, mean_rhuv_half, meanR_u_v_half

        real(rp) :: r_, ir, rWl, u_tau, ReyTu, yj, du, uc
        real(rp) :: mean_rhuu, mean_rhvv, mean_rhww, mean_rhuv
        real(rp) :: meanR_u_u, meanR_v_v, meanR_w_w, meanR_u_v
        real(rp) :: uu, vv, ww, uv, R11, R22, R33, R12

        VelFile%name = 'stat'
        velFile%dir  = trim(data_dir)//'/POST_VEL_STATS'
        call OpenNewFile(VelFile,it)

        BudFile%name = 'stat'
        BudFile%dir  = trim(data_dir)//'/POST_BUD_STATS'
        call OpenNewFile(BudFile,it)
        
        ! compute velocity y-derivative
        do j = sy-1,ey+1
        
           idyR = 1.0_rp/(0.5_rp*(ystep(j) + ystep(j+1)))

           u_y   (j) = 0.0_rp
           molVis(j) = 0.0_rp
           trbVis(j) = 0.0_rp
           do l = fd_L,fd_R

              cl1  = central_1_one_half(l)
              cl2  = central_2_one_half(l)
        
              ! u_y derivative
              u_l = Re_av%rhu(i0,j+l,k0)/Re_av%rho(i0,j+l,k0)
              u_y(j) = u_y(j) + cl2 * u_l * idyR

              ! viscous components
              T_l = Re_av%tmp(i0,j+l,k0)
              mVis_l = laminar_viscosity(T_l,tref,vis_flag)

              tVis_l = Re_av%VIS(i0,j+l,k0) - mVis_l

              molVis(j) = molVis(j) + cl1 * mVis_l
              trbVis(j) = trbVis(j) + cl1 * tVis_l

           enddo
        enddo

        tWl   = mu_inf*(molVis(0))*u_y(0)
        rWl   = 0.5_rp*(Re_av%rho(i0,0,k0) + Re_av%rho(i0,1,k0))
        u_tau = sqrt(tWl/rWl)
        ReyTu = rWl*u_tau/mu_inf
        !
        ! === compute yPlus and uPlus
        !
        yPlus = 0.0_rp
        uPlus = 0.0_rp
        do j = sy-1,ey
           yj       = 0.5_rp*(y(j) + y(j+1)) + 1.0_rp
           yPlus(j) = yj*ReyTu
           uPlus(j) = 0.5_rp*(Re_av%rhu(i0,j  ,k0)/Re_av%rho(i0,j  ,k0) &
                         + Re_av%rhu(i0,j+1,k0)/Re_av%rho(i0,j+1,k0))/u_tau
        enddo
        !
        ! === compute uPlus Van Driest
        !
        uPlusVD = 0.0_rp
        do j = sy,ey
           du = uPlus(j) - uPlus(j-1)
           r_ = 0.5_rp*(Re_av%rho(i0,j,k0) + Re_av%rho(i0,j-1,k0))
           uc = sqrt(r_/rWl)
           uPlusVD(j) = uPlusVD(j-1) + uc*du
        enddo


        do j = sy,ey-1

           yj = 0.5_rp*(y(j) + y(j+1)) + 1.0_rp
           r_ = 0.5_rp*(Re_av%rho(i0,j,k0) + Re_av%rho(i0,j+1,k0))
           ir = 1.0_rp/r_

           ! mean(rho*ui*uj)
           mean_rhuu = 0.5_rp*(Re_av%ruu(i0,j,k0) + Re_av%ruu(i0,j+1,k0))
           mean_rhvv = 0.5_rp*(Re_av%rvv(i0,j,k0) + Re_av%rvv(i0,j+1,k0))
           mean_rhww = 0.5_rp*(Re_av%rww(i0,j,k0) + Re_av%rww(i0,j+1,k0))
           mean_rhuv = 0.5_rp*(Re_av%ruv(i0,j,k0) + Re_av%ruv(i0,j+1,k0))

           ! mean(rhui)*mean(rhuj)
           meanR_u_u = (0.5_rp*(Re_av%rhu(i0,j,k0)+Re_av%rhu(i0,j+1,k0)))**2*ir
           meanR_v_v = (0.5_rp*(Re_av%rhv(i0,j,k0)+Re_av%rhv(i0,j+1,k0)))**2*ir
           meanR_w_w = (0.5_rp*(Re_av%rhw(i0,j,k0)+Re_av%rhw(i0,j+1,k0)))**2*ir
           meanR_u_v = (0.5_rp*(Re_av%rhu(i0,j,k0)+Re_av%rhu(i0,j+1,k0)))*&
                       (0.5_rp*(Re_av%rhv(i0,j,k0)+Re_av%rhv(i0,j+1,k0)))*ir

           ! u_i'*u_j'
           uu = (mean_rhuu - meanR_u_u)*ir
           vv = (mean_rhvv - meanR_v_v)*ir
           ww = (mean_rhww - meanR_w_w)*ir
           uv = (mean_rhuv - meanR_u_v)*ir

           R11 = (r_/rWl)*uu/u_tau**2
           R22 = (r_/rWl)*vv/u_tau**2
           R33 = (r_/rWl)*ww/u_tau**2
           R12 = (r_/rWl)*uv/u_tau**2

           write(velFile%unit,'(20e18.9)') yj, yPlus(j), uPlus(j), uPlusVD(j), &
                   R11, R22, R33, R12, r_/rWl

        enddo


        !
        ! === compute stress budget
        !
        do j = sy-1,ey

           yHalf = 0.5_rp*(y(j) + y(j+1)) + 1.0_rp
           rHalf = 0.5_rp*(Re_av%rho(i0,j,k0) + Re_av%rho(i0,j+1,k0))
           irhlf = 1.0_rp/rHalf
        
           ! viscous stresses
           mTrbUY = mu_inf*TrbVis(j)*u_y(j)/tWl
           mMolUY = mu_inf*molVis(j)*u_y(j)/tWl

           ! turbulent stresses
           mean_rhuv_half = 0.5_rp*(Re_av%ruv(i0,j,k0) + Re_av%ruv(i0,j+1,k0))
           meanR_u_v_half = 0.5_rp*(Re_av%rhu(i0,j,k0) + Re_av%rhu(i0,j+1,k0))* &
                            0.5_rp*(Re_av%rhv(i0,j,k0) + Re_av%rhv(i0,j+1,k0))*irhlf

           uv = (mean_rhuv_half - meanR_u_v_half)*irhlf

           R12 = rHalf*uv/tWl
           if(j == 0 .or. j == ey) R12 = 0.0_rp
           if(j == 0 .or. j == ey) mTrbUy = 0.0_rp

           write(BudFile%unit,*) yHalf, mMolUY, mtrbUY, - R12, &
                                        mMolUY + mTrbUY - R12
        enddo

        
        call CloseFile(VelFile)
        call CloseFile(BudFile)

        return
end subroutine write_turbulent_channel_WallModelled



subroutine write_turbulent_bdlayer_WallModelled

        use FileModule
        use fluid_functions_module, only: laminar_viscosity
        use real_to_integer_module, only: locate, nearest_integer_opt
        use interpolation_module  , only: polint


        implicit none
        type(FileType) :: swFile, statFile

        integer, parameter :: j0 = 0, j1 = 1, k0 = 1
        integer, parameter :: nSt = 2

        real(rp), dimension(:), allocatable :: ReTauVEC, frictVEC
        real(rp), dimension(0:ny) :: uTemp, uPlusHalf, yPlusHalf, uPlusVD
        real(rp), dimension(nSt)  :: ReTgT

        real(rp) :: mean_rhuu, mean_rhvv, mean_rhww, mean_rhuv
        real(rp) :: meanR_u_u, meanR_v_v, meanR_w_w, meanR_u_v
        real(rp) :: uu, vv, ww, uv, R11, R22, R33, R12

        real(rp) :: uyW, rhW, muW, mtW, r_, ir, du, uc
        real(rp) :: cl1, cl2, u_l, mtl, mul, idyR, cfrc, ReyTu, idvu
        real(rp) :: tauW, utau, u99,d99,d99err, tmpW, t_
        integer  :: i,l,j99,jj99, is, s, m = 2

        tmpW = 1.0_rp + Prandtl**(1.0_rp/3.0_rp)*0.5_rp*(gamma0-1.0_rp)*Mach**2

        allocate(ReTauVec(sx:ex), frictVEC(sx:ex))

        swFile%name = 'sw'
        swFile%dir  = trim(data_dir)//'/POST_STREAMWISE'
        call OpenNewFile(swFile,it)

        write(swFile%unit,*) '# Column 1: x '
        write(swFile%unit,*) '# Column 2: friction coefficient'
        write(swFile%unit,*) '# Column 3: Reynolds tau'

        do i = sx,ex

           idyR = 1.0_rp/(0.5_rp*(ystep(j0) + ystep(j1)))

           uyW = 0.0_rp 
           rhW = 0.0_rp
           muW = 0.0_rp
           mtW = 0.0_rp
           do l = fd_L,fd_R

              cl1 = central_1_one_half(l)
              cl2 = central_2_one_half(l)

              ! u_y derivative at wall location
              u_l = Re_av%rhu(i,j0+l,k0)/Re_av%rho(i,j0+l,k0)
              uyW = uyW + cl2 * u_l * idyR

              ! density at wall location
              rhW = rhW + cl1 * Re_av%rho(i,j0+l,k0)
        
              ! molecular viscosity
              mul = laminar_viscosity(Re_av%tmp(i,j0+l,k0),Tref,vis_flag)
              muW = muW + cl1 * mul

              ! turbulent viscosity
              mtl = Re_av%VIS(i,j0+l,k0) - mul
              mtW = mtW + cl1 * mtl

           enddo

           tauW = mu_inf*(muW + mtW)*uyW
           uTau = sqrt(tauW/rhW)
           cfrc = tauW/q_inf

           uTemp(sy:ey) = Re_av%rhu(i,sy:ey,2)/Re_av%rho(i,sy:ey,1)
           u99          = 0.985_rp*u_inf
           call locate(uTemp,0,ny,u99,j99)
           jj99 = min(max(j99-(m-1)/2,1),ny+1-m)
           call polint(uTemp(jj99),y(jj99),m,u99,d99,d99err)

           ! compute the local shear reynolds number
           ReyTu = rhW*utau*d99/(mu_inf*muW)

           ReTauVEC(i) = ReyTu
           frictVEC(i) = cfrc

           write(swFile%unit,*) x(i), cfrc, ReyTu

        enddo
        call CloseFile(swFile)

        ! determine the location of the control station
        ReTgt = (/840.0_rp, 900.0_rp/)

        StLoop: do is = 1,nSt

          statFile%name = 'stat'
          statFile%dir  = trim(data_dir)//'/POST_STAT_RE'//trim(str(int(ReTgt(is))))
          call OpenNewFile(statFile,it)

          ! find integer of the StLoc station
          s = nearest_integer_opt(retauVEC,sx,ex,ReTgt(is))

          idyR = 1.0_rp/(0.5_rp*(ystep(j0) + ystep(j1)))

          uyW = 0.0_rp 
          rhW = 0.0_rp
          muW = 0.0_rp
          mtW = 0.0_rp
          do l = fd_L,fd_R

             cl1 = central_1_one_half(l)
             cl2 = central_2_one_half(l)

             ! u_y derivative at wall location
             u_l = Re_av%rhu(s,j0+l,k0)/Re_av%rho(s,j0+l,k0)
             uyW = uyW + cl2 * u_l * idyR

             ! density at wall location
             rhW = rhW + cl1 * Re_av%rho(s,j0+l,k0)

             ! molecular viscosity
             mul = laminar_viscosity(Re_av%tmp(s,j0+l,k0),Tref,vis_flag)
             muW = muW + cl1 * mul

             ! turbulent viscosity
             mtl = Re_av%VIS(s,j0+l,k0) - mul
             mtW = mtW + cl1 * mtl

          enddo
          tauW = mu_inf*(muW + mtW)*uyW
          uTau = sqrt(tauW/rhW)
          idvu = rhW*utau/(mu_inf*muW)

          ! compute uPlus and yPlus 
          yPlusHalf = 0.0_rp
          uPlusHalf = 0.0_rp
          do j = sy,ey-1
             yPlusHalf(j) = 0.5_rp*(y(j)+y(j+1))*idvu
             uPlusHalf(j) = 0.5_rp*(Re_av%rhu(s,j  ,k0)/Re_av%rho(s,j  ,k0) + &
                                 Re_av%rhu(s,j+1,k0)/Re_av%rho(s,j+1,k0))/utau
          enddo

          ! compute uPlus Van Driest
          uPlusVD = 0.0_rp
          do j = sy,ey
             r_ = 0.0_rp
             do l = fd_L,fd_R
                cl1 = central_1_one_half(l)
                r_  = r_  + cl1 * Re_av%rho(s,j+l,k0)
             enddo

             du = uPlusHalf(j) - uPlusHalf(j-1)
             uc = sqrt(r_/rhW)
             uPlusVD(j) = uPlusVD(j-1) + uc*du
          enddo

          write(statFile%unit,'(A)')       '# *************************************************'
          write(statFile%unit,'(A,f18.6)') '#  TURBULENT BOUNDARY LAYER ', Mach 
          write(statFile%unit,'(A)')       '#'
          write(statFile%unit,'(A)')       '#  author: Francesco De Vanna '
          write(statFile%unit,'(A)')       '#  e-mail: fra.devanna@gmail.com'
          write(statFile%unit,'(A)')       '#'
          write(statFile%unit,'(A,f18.6)') '#  LOCATION%                 ', x(s)/Lx*100
          write(statFile%unit,'(A,f18.6)') '#  FRICTION REYNOLDS NUMBER  ', ReTauVEC(s)
          write(statFile%unit,'(A,f18.6)') '#  FRICTION COEFFICIENT      ', frictVEC(s)
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
          write(statFile%unit,'(A)')       '# Column 10 : prsrms+'
          write(statFile%unit,'(A)')       '# Column 11 : tmprms+'
          write(statFile%unit,'(A)')       '# Column 12 : rhorms+'
          write(statFile%unit,'(A)')       '# Column 13 : MuMol/MuMolWall'
          write(statFile%unit,'(A)')       '# Column 14 : MuTurb/MuMol'
          write(statFile%unit,'(A)')       '# *************************************************'
        
          do j = 0,ey-1

             r_ = 0.0_rp 
             t_ = 0.0_rp
             rhW= 0.0_rp
             do l = fd_L,fd_R
                cl1 = central_1_one_half(l)
                r_  = r_  + cl1 * Re_av%rho(s,j +l,k0)
                t_  = t_  + cl1 * Re_av%Tmp(s,j +l,k0)
                rhW = rhW + cl1 * Re_av%rho(s,j0+l,k0)
             enddo
             ir = 1.0_rp/r_

             ! mean(rho*ui*uj)
             mean_rhuu = 0.5_rp*(Re_av%ruu(s,j,k0) + Re_av%ruu(s,j+1,k0))
             mean_rhvv = 0.5_rp*(Re_av%rvv(s,j,k0) + Re_av%rvv(s,j+1,k0))
             mean_rhww = 0.5_rp*(Re_av%rww(s,j,k0) + Re_av%rww(s,j+1,k0))
             mean_rhuv = 0.5_rp*(Re_av%ruv(s,j,k0) + Re_av%ruv(s,j+1,k0))

             ! mean(rhui)*mean(rhuj)
             meanR_u_u = (0.5_rp*(Re_av%rhu(s,j,k0)+Re_av%rhu(s,j+1,k0)))**2*ir
             meanR_v_v = (0.5_rp*(Re_av%rhv(s,j,k0)+Re_av%rhv(s,j+1,k0)))**2*ir
             meanR_w_w = (0.5_rp*(Re_av%rhw(s,j,k0)+Re_av%rhw(s,j+1,k0)))**2*ir
             meanR_u_v = (0.5_rp*(Re_av%rhu(s,j,k0)+Re_av%rhu(s,j+1,k0)))*&
                         (0.5_rp*(Re_av%rhv(s,j,k0)+Re_av%rhv(s,j+1,k0)))*ir

             ! u_i'*u_j'
             uu = (mean_rhuu - meanR_u_u)*ir
             vv = (mean_rhvv - meanR_v_v)*ir
             ww = (mean_rhww - meanR_w_w)*ir
             uv = (mean_rhuv - meanR_u_v)*ir

             R11 = (r_/rhW)*uu/utau**2
             R22 = (r_/rhW)*vv/utau**2
             R33 = (r_/rhW)*ww/utau**2
             R12 = (r_/rhW)*uv/utau**2

                
             write(statFile%unit,'(20e18.9)') 0.5_rp*(y(j)+y(j+1)), &
                     yPlusHalf(j), uPlusHalf(j), uPlusVD(j), R11, R22, R33, R12, &
                     r_/rhW, t_/tmpW

          enddo

          call CloseFile(statFile)

        enddo Stloop

        
        deallocate(ReTauVEC,frictVEC)

        return
end subroutine write_turbulent_bdlayer_WallModelled



subroutine write_sic_quantities
        
        use fileModule
        use math_tools_module, only: degtorad
        use real_to_integer_module, only: nearest_integer_opt

        implicit none

        type(fileType)      :: massFile
        real(rp), parameter :: L0 = 150.0_rp   , hc = 0.24_rp
        real(rp), parameter :: t1 = 10.00_rp   , t2 = 22.00_rp
        real(rp), parameter :: r1 = 52.80_rp/L0, r2 = 32.36_rp/L0
        real(rp)            :: b1, b2

        real(rp) :: xst, yst, mass_flow, dA, A
        integer  :: i0, js, je, j,k

        massFile%name = 'mass_flow'
        massFile%dir  = trim(data_dir)//'/MASS_FLOW'
        call AppendToFile(massFile,it)

        b1 = r1*sin(degtorad(t1))
        b2 = r2*sin(degtorad(t2))

        xst = 0.7_rp
        yst = b1+b2

        i0 = nearest_integer_opt(x,sx,ex,xst)
        js = nearest_integer_opt(y,sy,ey,yst)
        je = nearest_integer_opt(y,sy,ey,hc)

        mass_flow = 0.0_rp
        A = 0.0_rp
        do k    = sz,ez
           do j = js,je
               dA = ystep(j)*zstep(k)
               mass_flow = mass_flow + phi(i,j,k,2)*dA
               A = A + dA
           enddo
        enddo
        mass_flow = mass_flow/A

        write(massFile%unit,'(1i7,2e18.9)') it, time, mass_flow


        call CloseFile(massFile)

        





        return
end subroutine write_sic_quantities



subroutine write_sic_centerline

        use fileModule
        use real_to_integer_module, only: nearest_integer_opt

        implicit none
        type(FileType) :: outFile
        real(rp), allocatable, dimension(:), save :: U_cmean
        real(rp)                                  :: u_i
        integer :: js
        
        if(.not.allocated(U_cmean)) then
          allocate(U_cmean(sx:ex))
          U_cmean = 0.0_rp
        endif

        outFile%name = 'vel_centerline'
        outFile%dir  = trim(data_dir)//'/VEL_CENTERLINE'
        call OpenNewFile(outFile,it)

        js = nearest_integer_opt(y,sy,ey,0.18_rp)

        do i = sx,ex

           u_i = U(i,js,ez/2)

           U_cmean(i) = (U_cmean(i)*real(f-file_ini,rp) + u_i)/real(f-file_ini+1,rp)

           write(outFile%unit,'(3e18.6)') x(i), u_i, U_cmean(i)

        enddo


        call CloseFile(outFile)

        return
end subroutine write_sic_centerline



subroutine channel_flow_velocity_correlation
        
        use FileModule
        use post_statistic_time_module, only: Rey_averaged

        implicit none
        type(fileType) :: corrFileInner
        type(fileType) :: corrFileOuter

        integer :: i,j,k,r,s
        real(rp), allocatable, dimension(:,:,:) :: tmp
        real(rp) :: um2, co_uu, co_vv, co_ww

        allocate(tmp(lby:uby,lbz:ubz,3))
        
        tmp = 0.0_rp
        do j          = sy,ey
           do s       = 0,nz/2
              do k    = sz,ez
                 do i = sx,ex

                    r = s+k
                    if(r > ez) r = r-ez

                    tmp(j,s,1) = tmp(j,s,1) + U(i,j,k)*U(i,j,r)
                    tmp(j,s,2) = tmp(j,s,2) + V(i,j,k)*V(i,j,r)
                    tmp(j,s,3) = tmp(j,s,3) + W(i,j,k)*W(i,j,r)

                 enddo
              enddo
           enddo
        enddo
        tmp = tmp/real(nx*nz,rp)
        
        call Rey_averaged(f,file_ini,file_dt,time_from_restart,tmp,corr_v)
        
        corrFileInner%name = 'corrv'
        corrFileInner%dir  = trim(data_dir)//'/POST_CORR_INNER'
        call OpenNewFile(corrFileInner,it)

        corrFileOuter%name = 'corrv'
        corrFileOuter%dir  = trim(data_dir)//'/POST_CORR_OUTER'
        call OpenNewFile(corrFileOuter,it)

        do j    = sy,ey/2
           do s = 0,nz/2

              um2 = (Re_av%rhu(sx,j,sz)/Re_av%rho(sx,j,sz))**2
                
              co_uu = (corr_v(j,s,1)-um2)/(corr_v(j,0,1)-um2)
              co_vv = corr_v(j,s,2)/corr_v(j,0,2)
              co_ww = corr_v(j,s,3)/corr_v(j,0,3)

              write(corrFileInner%unit,'(5e18.6)') (y(j)+1.0_rp)*ReTau, s*zstep(1)*ReTau,&
                      co_uu, co_vv, co_ww

              write(corrFileOuter%unit,'(5e18.6)') (y(j)+1.0_rp), s*zstep(1),&
                      co_uu, co_vv, co_ww
           enddo

        enddo
        call CloseFile(corrFileInner)
        call CloseFile(corrFileOuter)


        deallocate(tmp)

        return
end subroutine channel_flow_velocity_correlation





subroutine WritePDF(var)
        implicit none
        real(rp), dimension(:,:,:), allocatable, intent(in) :: var

        ! local declaration
        real(rp), dimension(:)    , allocatable :: zones

        integer, parameter :: Nbin = 100

        real(rp) :: maxvar, minvar, dvar, zone_, val
        integer  :: plane, l
        type(FileType) :: PDFFile

        call AllocateReal(zones,1,Nbin)
        
        do j = 1,10
           plane = j

           maxvar = maxval(var(sx:ex,plane,sz:ez))
           minvar = minval(var(sx:ex,plane,sz:ez))

           dvar = 1.0_rp/Nbin

           zones = 0.0_rp
           do    k = sz,ez
              do i = sx,ex
                 
                 val = var(i,j,k)
                 zone_ = 0.0_rp
                 do l = 1,Nbin-1

                    if(val > zone_ .and. val < zone_+dvar) then
                      zones(l) = zones(l)+1
                    endif
                    zone_ = zone_ + dvar

                 enddo

              enddo
           enddo
           !zones = zones/maxval(zones)

           PDFFile%name = 'ducros_pdf'
           PDFFile%dir  = trim(data_dir)//'/DUCROS_PDF_PLANE'//trim(str(plane))
           call OpenNewFile(PDFFile,it)
           zone_ = 0.0_rp
           do l = 1,Nbin
              write(pdfFile%unit,*) zone_,zones(l)
              zone_ = zone_+dvar
           enddo
           call CloseFile(pdfFile)
        enddo

        call DeallocateReal(zones)

        return
end subroutine WritePDF







end module post_output_gnu_module
