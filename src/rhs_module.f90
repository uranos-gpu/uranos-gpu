module rhs_module
use parameters_module
use storage_module
use nvtx

implicit none
private
public rhs_navier_stokes, update_all, update_all_ghosts, mpi_wait_conservatives


contains
subroutine rhs_navier_stokes
! -----------------------------------------------------------------------
!
!       Computation of the RHS of the NAVIER-STOKES equations
!
! -----------------------------------------------------------------------
        use advection_module, only: advection_fluxes
        use viscous_module  , only: viscous_fluxes

        implicit none
        integer :: i,j,k,l

#ifdef TIME
        call mpi_stime(s_rhs_time)
#endif
#ifdef NVTX
        call nvtxStartRange("rhs_navier_stokes") 
#endif
        !$acc parallel default(present)
        !$acc loop gang, vector collapse(4)
        do          l = 1,5
           do       k = sz,ez
              do    j = sy,ey
                 do i = sx,ex
                    RHS(i,j,k,l) = 0.0_rp
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel loop
       
        ! === Convective fluxes computation
        call advection_fluxes(scheme)

        ! === Viscous fluxes computation
        if(viscous) call viscous_fluxes(dims)

        ! === Forcing terms computation
        call forcing_terms


#ifdef NVTX
        call nvtxEndRange 
#endif
#ifdef TIME
        call mpi_etime(s_rhs_time,t_rhs_calls,t_rhs_time)
#endif
        return
end subroutine rhs_navier_stokes



subroutine update_all
! -----------------------------------------------------------------------------
!       This subroutine update all the primitives variables in the inner region
! -----------------------------------------------------------------------------
        implicit none

#ifdef TIME
        call mpi_stime(s_upd_time)
#endif
#ifdef NVTX
        call nvtxStartRange("update_all") 
#endif
        
        !call compute_primitives(lbx,sy,sz,ubx,ey,ez)
        call compute_primitives(lbx,lby,lbz,ubx,uby,ubz)

#ifdef NVTX
        call nvtxEndRange
#endif
#ifdef TIME
        call mpi_etime(s_upd_time,t_upd_calls,t_upd_time)
#endif
        return
end subroutine update_all



subroutine mpi_wait_conservatives
! -------------------------------------------------------------
!       waiting the end of conservative communications
! -------------------------------------------------------------
        implicit none

#ifdef TIME
        call mpi_stime(s_wcs_time)
#endif
#ifdef NVTX
        call nvtxStartRange("mpi_wait_conservatives") 
#endif
#ifdef MPIBFR
        ! MPI BUFFER IMPLEMENTATION
        call mpi_wait_procs(req_array_yz)

        phi(:,sy-GN:sy-1,:,:) = phi_bfr_recv_S
        phi(:,ey+1:ey+GN,:,:) = phi_bfr_recv_N

        phi(:,:,sz-GN:sz-1,:) = phi_bfr_recv_B
        phi(:,:,ez+1:ez+GN,:) = phi_bfr_recv_F
#else 

        ! MPI DERIVED DATA TYPE IMPLEMENTATION
        call mpi_wait_procs(req_array_yz)

#endif
#ifdef NVTX
        call nvtxEndRange
#endif
#ifdef TIME
        call mpi_etime(s_wcs_time,t_wcs_calls,t_wcs_time)
#endif
        return
end subroutine mpi_wait_conservatives



subroutine update_all_ghosts
! -------------------------------------------------------------
!       Update the primitives variables in the ghost regions 
!       ATTENTION: only after 'mpi_wait_conservatives'
! -------------------------------------------------------------
        implicit none
#ifdef TIME
        call mpi_stime(s_upg_time)
#endif
#ifdef NVTX
        call nvtxStartRange("update_all_ghosts")
#endif

        ! south 
        call compute_primitives(lbx,lby,lbz,ubx,sy-1,ubz)
        ! north
        call compute_primitives(lbx,ey+1,lbz,ubx,uby,ubz)

        if(dims == 3) then
          ! backward
          call compute_primitives(lbx,lby,lbz,ubx,uby,sz-1)
          ! farward
          call compute_primitives(lbx,lby,ez+1,ubx,uby,ubz)
        endif

#ifdef NVTX
        call nvtxEndRange
#endif
#ifdef TIME
        call mpi_etime(s_upg_time,t_upg_calls,t_upg_time)
#endif
        return
end subroutine update_all_ghosts



subroutine compute_primitives(lo_1,lo_2,lo_3,hi_1,hi_2,hi_3)
        implicit none
        integer, intent(in) :: lo_1,lo_2,lo_3
        integer, intent(in) :: hi_1,hi_2,hi_3
        real(rp)            :: r_,ir,u_,v_,w_,p_,T_,mu_
        integer             :: i,j,k

#ifdef NVTX
        call nvtxStartRange("compute_primitives")
#endif

        !$omp parallel do collapse(3) default(private), shared(lo,hi,phi) &
        !$omp shared(U,V,W,P,T)

        !$acc parallel default(present)
        !$acc loop gang,vector collapse(3)
        do k       = lo_3,hi_3
           do j    = lo_2,hi_2
              do i = lo_1,hi_1

                 r_ = phi(i,j,k,1)
                 ir = 1.0_rp/r_
                 u_ = phi(i,j,k,2)*ir
                 v_ = phi(i,j,k,3)*ir
                 w_ = phi(i,j,k,4)*ir
                 p_ = cv_i*(phi(i,j,k,5) - 0.5_rp*r_*(u_*u_+v_*v_+w_*w_))
                 T_ = p_*ir

                 U(i,j,k) = u_
                 V(i,j,k) = v_
                 W(i,j,k) = w_
                 P(i,j,k) = p_
                 T(i,j,k) = T_

              enddo
           enddo
        enddo
        !$acc end parallel

        !$omp end parallel do

        if(viscous) then
        
          !$omp parallel do collapse(3) default(private), shared(lo,hi,VIS,LMD,T,k_inf)

          !$acc parallel default(present)
          !$acc loop gang,vector collapse(3)
          do k       = lo_3,hi_3
             do j    = lo_2,hi_2
                do i = lo_1,hi_1

                   T_  = T(i,j,k)
                   mu_ = T_*sqrt(T_) * (1.0_rp+suthc)/(T_+suthc)

                   VIS(i,j,k) = mu_
                   LMD(i,j,k) = k_inf * mu_

                enddo
             enddo
          enddo
          !$acc end parallel

          !$omp end parallel do

        endif

#ifdef NVTX
        call nvtxEndRange
#endif
        return
end subroutine compute_primitives



subroutine forcing_terms
        implicit none
        real(rp) :: f
        integer  :: i,j,k

#ifdef NVTX
        call nvtxStartRange("forcing_terms") 
#endif

        selectcase(ic)
          case('poiseuille_x', 'poiseuille_z','inflow_poiseuille')

            f = 2._rp*mu_inf *Mach*sqrt(gamma0)

            !$acc parallel default(present)
            !$acc loop gang, vector collapse(3)
            do       k = sz,ez
               do    j = sy,ey
                  do i = sx,ex
                    RHS(i,j,k,2) = RHS(i,j,k,2) + f
                  enddo
               enddo
            enddo
            !$acc end parallel loop

          case('poiseuille_y')

            f = 2._rp*mu_inf *Mach*sqrt(gamma0)

            !$acc parallel default(present)
            !$acc loop gang, vector collapse(3)
            do       k = sz,ez
               do    j = sy,ey
                  do i = sx,ex
                     RHS(i,j,k,3) = RHS(i,j,k,3) + f
                  enddo
               enddo
            enddo
            !$acc end parallel loop

          case('turbulent_channel')

            call ComputeTCHPressureGradient(phi,RHS,force_turbulent_channel)

          case('hTurb')

            call hTurb_ForcingTerm(phi,RHS)

          case('KolmogorovFlow')

            call kolmogorov_flow_forcing(RHS)
                
        endselect

#ifdef NVTX
        call nvtxEndRange 
#endif
return
end subroutine forcing_terms



subroutine ComputeTCHPressureGradient(phi,RHS,prs_grad)
! --------------------------------------------------------
!       Computation of the pressure gradient dp_dx for a 
!       channel driven by pressure
! --------------------------------------------------------
        
        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(in)    :: phi
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: RHS
        real(rp)                                 , intent(out)   :: prs_grad
        
        real(rp)           :: lcl_r_rhs, lcl_rurhs
        real(rp)           :: gbl_r_rhs, gbl_rurhs
        real(rp)           :: dv, u_
        integer            :: i,j,k

        lcl_r_rhs = 0.0_rp
        lcl_rurhs = 0.0_rp

        !$acc parallel default(present) copy(lcl_r_rhs,lcl_rurhs)
        !$acc loop gang, vector collapse(3) reduction(+:lcl_r_rhs,lcl_rurhs)
        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex

                 if(abs(y(j)) < 1.0_rp) then

                 dv = xstep(i)*ystep(j)*zstep(k)

                 lcl_r_rhs = lcl_r_rhs + RHS(i,j,k,1)*dv
                 lcl_rurhs = lcl_rurhs + RHS(i,j,k,2)*dv

                 endif

              enddo
           enddo
        enddo
        !$acc end parallel loop

        call MPI_allreduce(lcl_r_rhs,gbl_r_rhs,1,MPI_RP,mpi_sum,mpi_comm_cart,err)
        call MPI_allreduce(lcl_rurhs,gbl_rurhs,1,MPI_RP,mpi_sum,mpi_comm_cart,err)

        gbl_r_rhs = gbl_r_rhs/(Lx*2.0_rp*Lz)
        gbl_rurhs = gbl_rurhs/(Lx*2.0_rp*Lz)

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do       k = sz,ez
           do    j = sy,ey
              do i = sx,ex
                 if(abs(y(j)) < 1.0_rp) then

                   u_ = phi(i,j,k,2)/phi(i,j,k,1)

                   RHS(i,j,k,1) = RHS(i,j,k,1) - gbl_r_rhs
                   RHS(i,j,k,2) = RHS(i,j,k,2) - gbl_rurhs
                   RHS(i,j,k,5) = RHS(i,j,k,5) - gbl_rurhs*u_

                 endif
              enddo
           enddo
        enddo
        !$acc end parallel loop


        prs_grad = - gbl_rurhs

        return
end subroutine ComputeTCHPressureGradient



subroutine hTurb_ForcingTerm(phi,RHS)
! -------------------------------------------------------------------------
!       This subroutine compute the forcing term for Homogeneous turbulence
! -------------------------------------------------------------------------

        use integration_module
        use FileModule
        use random_module, only: rnd

        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(in)    :: phi
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: RHS

        ! local declarations
        real(rp), parameter :: a1 =   19.201_rp
        real(rp), parameter :: a2 = - 113.6716667_rp
        real(rp), parameter :: a3 = + 355.5_rp
        real(rp), parameter :: a4 = - 387.833333_rp

        real(rp) :: arm11, arm12, arm13
        real(rp) :: arm21, arm22, arm23
        real(rp) :: ph11 , ph12 , ph13
        real(rp) :: ph21 , ph22 , ph23, ph24, ph25, ph26
        real(rp) :: xi1  , xi2  , xi3
        real(rp) :: r_, ru, rv, rw, fx, fy, fz
        real(rp) :: iLx, iLy, iLz, A
        integer  :: i,j,k

        iLx = 1.0_rp/Lx
        iLy = 1.0_rp/Ly
        iLz = 1.0_rp/Lz

        A = (a1 + a2 * Mach + a3 * Mach**2 + a4 * Mach**3)*Mach**2
        
        ! compute the random phase for first mode
        ph11 = 2*pi*rnd()
        ph12 = 2*pi*rnd()
        ph13 = 2*pi*rnd()

        ! compute the random phase for second mode
        ph21 = 2*pi*rnd()
        ph22 = 2*pi*rnd()
        ph23 = 2*pi*rnd()
        ph24 = 2*pi*rnd()
        ph25 = 2*pi*rnd()
        ph26 = 2*pi*rnd()

       !$acc parallel default(present)
       !$acc loop gang, vector collapse(3) 
        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex

                 ! read the grid from RAM
                 xi1 = 2.0_rp*pi*x(i)*iLx
                 xi2 = 2.0_rp*pi*y(j)*iLy
                 xi3 = 2.0_rp*pi*z(k)*iLz

                 ! read density
                 r_ = phi(i,j,k,1)
                 ru = phi(i,j,k,2)
                 rv = phi(i,j,k,3)
                 rw = phi(i,j,k,4)

                 ! first armonic component
                 arm11 = A * sin(xi1 + ph11)
                 arm12 = A * sin(xi2 + ph12)
                 arm13 = A * sin(xi3 + ph13)

                 ! second armonic component
                 arm21 = A * sin(xi2 + ph21) * sin(xi3 + ph22)
                 arm22 = A * sin(xi1 + ph23) * sin(xi3 + ph24)
                 arm23 = A * sin(xi1 + ph25) * sin(xi2 + ph26)

                 ! solenoidal field
                 fx = arm12 + arm13 + arm21
                 fy = arm11 + arm13 + arm22
                 fz = arm11 + arm12 + arm23

                 ! compute forcing terms
                 RHS(i,j,k,2) = RHS(i,j,k,2) + r_*fx
                 RHS(i,j,k,3) = RHS(i,j,k,3) + r_*fy
                 RHS(i,j,k,4) = RHS(i,j,k,4) + r_*fz
                 RHS(i,j,k,5) = RHS(i,j,k,5) + ru*fx + rv*fy + rw*fz

              enddo
           enddo
        enddo
        !$acc end parallel loop

        return
end subroutine hTurb_ForcingTerm



subroutine kolmogorov_flow_forcing(RHS)

        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: RHS

        real(rp) :: fx
        integer  :: i,j,k

       !$acc parallel default(present)
       !$acc loop gang, vector collapse(3) 
        do       k = sz,ez
           do    j = sy,ey
              do i = sx,ex

                 fx = 0.1_rp*u_inf*sin(8*pi*y(j)/Ly)
                 RHS(i,j,k,2) = RHS(i,j,k,2) + fx
                 RHS(i,j,k,3) = 0.0_rp
                 RHS(i,j,k,4) = 0.0_rp
                 RHS(i,j,k,5) = 0.0_rp

              enddo
           enddo
        enddo
        !$acc end parallel loop

        return
end subroutine kolmogorov_flow_forcing





















end module rhs_module
