program uranos
        use parameters_module
        use input_module
        use mesh_module
        use mpi_module
        use output_module
        use storage_module
        use time_module
        use init_bc_module
        use bc_module
        use rhs_module
        use update_variable_module
        use inflow_module
        use random_module
        use sgs_module
        use GetRetau_module
        use df_module
        use thermodynamics_module
        use profiling_module, only: FinalizeProfiling
        use viscous_shared_module, only: sgs_isotropic_setup
        implicit none

        call init_mpi
        call parse_command_line_args
        call read_input_namelist
        call init_non_calorically_coefficients
        call set_reference_quantities
        call resolve_scheme_flags
        call make_saving_directories
        call write_thermodynamics_profiles
        call init_Reynolds

        call init_cartesian_topology
        call init_subcartesian_topology
        call init_subdomain_boundaries
        call init_mpi_neighbours

        call init_mpi_exchange_data

        call init_GL_variables
        call init_FD_coefficients
        if(LES .and. yoshizawa) call sgs_isotropic_setup
        call init_grid
        call init_statistics_fields
        call init_RK_coefficients
        call init_DefaultSeed

        call init_inflow_profile

        call init_boundary_conditions

        call output_init

        if(turb_inflow) call inflow_turbulent_load_state

        !$acc data copyin(mid_point_lele_x,mid_point_lele_y,mid_point_lele_z) &
        !$acc copyin(xstep_i,ystep_i,zstep_i,ixsteph,iysteph,izsteph) &
        !$acc copyin(x,y,z,xstep,ystep,zstep,ix_csi,iy_eta,iz_zta) &
        !$acc copyin(csistep_i,y_eta,y_eta2,x_csi,x_csi2,etastep_i) &
        !$acc copyin(csistep_i, etastep_i, ztastep_i) &
        !$acc copyin(x_csi2, y_eta2, z_zta2) &
        !$acc copyin(dxi_dx,dxi_dy,deta_dx,deta_dy,m_xi,m_eta) &
        !$acc copyin(nx_xi,ny_xi,nx_eta,ny_eta) &
        !$acc copyin(Jac2D,dx_dxi,dx_deta,dy_dxi,dy_deta) &
        !$acc copyin(Gamma1_2D,Gamma2_2D,xx,yy) &
        !$acc copyin(dxi_dxx,dxi_dyy,dxi_dxy,deta_dxx,deta_dyy,deta_dxy) &
        !$acc copyin(gJ11,gJ22,gJ12) &
        !$acc copyin(central_1_one_half,central_2_one_half,central_1,central_2) &
        !$acc copyin(visc_divs,visc_stag,visc_intp) &
        !$acc copyin(visc_stag_b,visc_intp_b) &
        !$acc copyin(fward_1, bward_1, c1_wbias, c1_ebias) &
        !$acc copyin(weno_num,my_neighbour,bc) &
        !$acc copyin(a_rk,b_rk,c_rk) &
        !$acc copyin(pri_1D_x, pri_1D_y, pri_1D_z) &
        !$acc copyin(phi_arr_x, phi_arr_y, phi_arr_z) &
        !$acc copyin(flx_arr_x, flx_arr_y, flx_arr_z) &
        !$acc copyin(flx_x,flx_y,flx_z,ishock_x,ishock_y,ishock_z) &
        !$acc copyin(phi_bfr_send_E, phi_bfr_send_W, phi_bfr_recv_E, phi_bfr_recv_W) &
        !$acc copyin(phi_bfr_send_N, phi_bfr_send_S, phi_bfr_recv_N, phi_bfr_recv_S) &
        !$acc copyin(phi_bfr_send_B, phi_bfr_send_F, phi_bfr_recv_B, phi_bfr_recv_F) &
        !$acc copyin(bfr_send_E, bfr_send_W, bfr_recv_E, bfr_recv_W) &
        !$acc copyin(bfr_send_N, bfr_send_S, bfr_recv_N, bfr_recv_S) &
        !$acc copyin(bfr_send_B, bfr_send_F, bfr_recv_B, bfr_recv_F) &
        !$acc copyin(i13D_bfr_send_E, i13D_bfr_send_W, i13D_bfr_recv_E, i13D_bfr_recv_W) &
        !$acc copyin(i13D_bfr_send_N, i13D_bfr_send_S, i13D_bfr_recv_N, i13D_bfr_recv_S) &
        !$acc copyin(i13D_bfr_send_B, i13D_bfr_send_F, i13D_bfr_recv_B, i13D_bfr_recv_F) &
        !$acc copyin(DF_Rnd2D,DF_ylen,DF_zlen,DF_uconv,DF_By,DF_Bz,DF_Fy,DF_LundMatrix) &
        !$acc copyin(iflow_mean,iflow_turb,iflow_vf_old, iflow_vf_new,iflow_dTdu,iflow_urms) &
        !$acc copyin(weno_flag_xyz) &
        !$acc copy(phi,U,V,W,T,P,VIS,LMD,DIV,RHS,SSENSOR,weno_flag,phi_n,isoSGS) &
        !$acc copy(GAMMA_T,ENTHALPY) &
#ifdef NVIDIA
        !$acc copyin(generator)
#endif

       !$acc data copyin(vfclo_x,vfchi_x,vfcsl_x,vfcsh_x) &
       !$acc      copyin(vfclo_y,vfchi_y,vfcsl_y,vfcsh_y) &
       !$acc      copyin(vfclo_z,vfchi_z,vfcsl_z,vfcsh_z) &
       !$acc      copyin(fmx_csix,fmx_csiy,fmx_etax,fmx_etay,fmx_cof1,fmx_cof2) &
       !$acc      copyin(fmy_csix,fmy_csiy,fmy_etax,fmy_etay,fmy_cof1,fmy_cof2) &
       !$acc      copyin(fmz_zeta) if(diffusion_scheme == 'conservative')

       !$acc data copyin(vmean2D,stat_lcl2D,stat_gbl2D) if(allocated(vmean2D))
       !$acc data copyin(vmean1D,stat_pln1D,stat_lcl1D,stat_gbl1D) if(allocated(vmean1D))

        call init_mpi_windows

        if(enable_limiter) call limit_positivity

        call mpi_bc_communications
        call set_bc_conditions
        call update_all

        if(les) call compute_subgrid_model

        s_cpu_time = MPI_WTIME()

        call output_start

        do while(time .le. tmax .and. istop == 0 .and. it < itmax)

                call last_iteration(s_cpu_time,istop)

                if(logical_CFL) call compute_dt

                call init_runge_kutta

                do ik = 1, n_step
                   call rhs_navier_stokes
                   call runge_kutta
                   if(enable_limiter) call limit_positivity
                   call mpi_bc_communications
                   call set_bc_conditions
                   call update_all
                   if(les) call compute_subgrid_model

                enddo

                if(mod(it,rsOut) == 0 .and. enable_limiter) call report_positivity_limiter

                it = it + 1

                call output_step
        enddo

        call output_finalize

        !$acc end data
        !$acc end data
        !$acc end data
        !$acc end data

        if(enable_limiter) call finalize_positivity_limiter

        call FinalizeProfiling
        call end_all_variables
        call end_mpi

end program uranos
