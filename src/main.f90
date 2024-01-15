!================================================================================================
! Starting development: 01/04/2017, Venice, Italy.
!
!       ██╗   ██╗██████╗  █████╗ ███╗   ██╗ ██████╗ ███████╗
!       ██║   ██║██╔══██╗██╔══██╗████╗  ██║██╔═══██╗██╔════╝
!       ██║   ██║██████╔╝███████║██╔██╗ ██║██║   ██║███████╗
!       ██║   ██║██╔══██╗██╔══██║██║╚██╗██║██║   ██║╚════██║
!       ╚██████╔╝██║  ██║██║  ██║██║ ╚████║╚██████╔╝███████║
!        ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝ ╚═════╝ ╚══════╝
!                                                    
! Uranos (Unsteady Robust All-around Navier-StOkes Solver) is a fully-compressible Navier-Stokes 
! solver expecially developed to treat complex geometries. 
! The code was entirely written in Fortran90 by Francesco De Vanna. 
!
! To run the code ask to: fra.devanna@gmail.com
!================================================================================================
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
        use inflow_module
        use random_module
        use sgs_module
        use GetRetau_module
        use df_module
        implicit none

        ! Initiale mpi environment end read inputs
        call init_mpi
        call read_input_data
        call make_saving_directories
        call init_Reynolds

        !! Cartesian domain 3D splitting
        call init_cartesian_topology
        call init_subcartesian_topology
        call init_subdomain_boundaries
        call init_mpi_neighbours

        ! build datatype
        call init_mpi_exchange_data

        ! Initializing variables
        call init_GL_variables
        call init_statistics_fields
        call init_FD_coefficients
        call init_grid
        call init_RK_coefficients
        call init_DefaultSeed
        call initProbes(x,y,z,Probe)
        call init_sponge
        
        ! Initialize inflow
        if(inflow) call init_inflow_profile

        ! Initialize the solution
        call init_boundary_conditions


        ! Time loop ------------------------------
        !$acc data copyin(mid_point_lele_x,mid_point_lele_y,mid_point_lele_z) &
        !$acc copyin(xstep_i,ystep_i,zstep_i,ixsteph,iysteph,izsteph) &
        !$acc copyin(x,y,z,xstep,ystep,zstep,ix_csi,iy_eta,iz_zta) &
        !$acc copyin(csistep_i,y_eta,y_eta2,x_csi,x_csi2,etastep_i) &
        !$acc copyin(csistep_i, etastep_i, ztastep_i) &
        !$acc copyin(x_csi2, y_eta2, z_zta2) &
        !$acc copyin(central_1_one_half,central_2_one_half,central_1,central_2) &
        !$acc copyin(fward_1, bward_1) &
        !$acc copyin(weno_num,my_neighbour,bc) &
        !$acc copyin(a_rk,b_rk,c_rk) &
        !$acc copyin(mask_field) &
        !$acc copyin(tilde_op_x, tilde_op_y, tilde_op_z) &
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
        !$acc copyin(DF_Rnd2D,DF_ylen,DF_zlen,DF_By,DF_Bz,DF_Fy,DF_LundMatrix) &
        !$acc copyin(iflow_mean,iflow_turb,iflow_vf_old, iflow_vf_new) &   !!!!! >>>> WARNING HERE
        !$acc copyin(wmles_data,wmles_data_uw,wmles_data_lw) &
        !$acc copyin(uFarField,vFarField,pFarField,TFarField) &
        !$acc copyin(weno_flag_xyz, sponge_x) &
        !$acc copy(phi,U,V,W,T,P,VIS,LMD,DIV,RHS,SSENSOR,weno_flag,phi_n) &
#ifdef NVIDIA
        !$acc copyin(generator,phirand,psirand)
#endif

#ifdef AMD
       !$acc copyin(rnd2Dnum,rnd2Dnumptr,rnd2Dnumsize,rnd2Dnumshape,ic) &
       !$acc copyin(phirand,phirandnum,phirandsize,phirandnumptr) &
       !$acc copyin(psirand,psirandnum,psirandsize,psirandnumptr)
#endif
       !$acc data copyin(mask_GhostId,mask_GhostUnstr,mask_nGhost,mask_field) &

        ! init subgrid stresses
        if(les) call compute_subgrid_model

        s_cpu_time = MPI_WTIME()

        str_it_time = MPI_WTIME()

        do while(time .le. tmax .and. istop == 0 .and. it < itmax)
        
                call last_iteration(s_cpu_time,istop)
                
                ! compute dt under CFL condition

                if(logical_CFL) call compute_dt

                call init_runge_kutta

                ! print results sometimes
                if(printres) then
                  if(mod(it,itOut) == 0) call compute_average_iteration_time
                  if(mod(it,itout) == 0) call write_restart
                  if(mod(it,StOut) == 0 .and. stFlg) call write_restart_stat
                  if(mod(it,StOut) == 0 .and. stFlg) call write_statistics
                endif

                ! Runge - Kutta substeps
                do ik = 1, n_step
                   call rhs_navier_stokes
                   call runge_kutta
                   call mpi_bc_communications 
                   call set_bc_conditions
                   call update_all
                   if(les) call compute_subgrid_model

                enddo
               
                if(mod(it,rsOut) == 0 .and. printres) call dump_residuals

                it = it + 1
        enddo
        if(printres) then
          if(stFlg) call write_statistics
          call write_restart
        endif
        e_cpu_time = MPI_WTIME()
        !$acc end data
        !$acc end data
        !$acc end data

        call screen_elapse_time

        call end_all_variables
        call end_mpi

end program uranos


