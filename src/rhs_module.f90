module rhs_module
use parameters_module
use storage_module
use profiling_module

implicit none
private
public rhs_navier_stokes, update_all, update_all_ghosts, mpi_wait_conservatives, dump_residuals


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

        call StartProfRange("rhs_navier_stokes") 

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
        !$acc end parallel
       
        ! === Convective fluxes computation
        call advection_fluxes(scheme)

        ! === Viscous fluxes computation
        if(viscous) call viscous_fluxes(dims)

        ! === Forcing terms computation
        call forcing_terms

        call EndProfRange 

        return
end subroutine rhs_navier_stokes



subroutine update_all
! -----------------------------------------------------------------------------
!       This subroutine update all the primitives variables in the inner region
! -----------------------------------------------------------------------------
        implicit none

        call StartProfRange("update_all") 
        
        !call compute_primitives(lbx,sy,sz,ubx,ey,ez)
        call compute_primitives(lbx,lby,lbz,ubx,uby,ubz)

        call EndProfRange

        return
end subroutine update_all



subroutine mpi_wait_conservatives
! -------------------------------------------------------------
!       waiting the end of conservative communications
! -------------------------------------------------------------
        implicit none

        call StartProfRange("mpi_wait_conservatives") 

        ! MPI BUFFER IMPLEMENTATION
        call mpi_wait_procs(req_array_yz)

        phi(:,sy-GN:sy-1,:,:) = phi_bfr_recv_S
        phi(:,ey+1:ey+GN,:,:) = phi_bfr_recv_N

        phi(:,:,sz-GN:sz-1,:) = phi_bfr_recv_B
        phi(:,:,ez+1:ez+GN,:) = phi_bfr_recv_F

        ! MPI DERIVED DATA TYPE IMPLEMENTATION
        call mpi_wait_procs(req_array_yz)

        call EndProfRange

        return
end subroutine mpi_wait_conservatives



subroutine update_all_ghosts
! -------------------------------------------------------------
!       Update the primitives variables in the ghost regions 
!       ATTENTION: only after 'mpi_wait_conservatives'
! -------------------------------------------------------------
        implicit none

        call StartProfRange("update_all_ghosts")

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

        call EndProfRange

        return
end subroutine update_all_ghosts



subroutine compute_primitives(lo_1,lo_2,lo_3,hi_1,hi_2,hi_3)
        
        use fluid_functions_module, only: laminar_viscosity

        implicit none
        integer, intent(in) :: lo_1,lo_2,lo_3
        integer, intent(in) :: hi_1,hi_2,hi_3
        real(rp)            :: r_,ir,u_,v_,w_,p_,T_,mu_
        integer             :: i,j,k

        call StartProfRange("compute_primitives")

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
                   mu_ = laminar_viscosity(T_,tref,vis_flag)

                   VIS(i,j,k) = mu_
                   LMD(i,j,k) = k_inf * mu_

                enddo
             enddo
          enddo
          !$acc end parallel

          !$omp end parallel do

        endif

        call EndProfRange

        return
end subroutine compute_primitives



subroutine forcing_terms
        implicit none
        real(rp) :: f
        integer  :: i,j,k

        call StartProfRange("forcing_terms") 

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
            !$acc end parallel

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
            !$acc end parallel

          case('turbulent_channel')

            call ComputeTCHPressureGradient(phi,RHS,force_turbulent_channel)


          case('hTurb')

            call hTurb_ForcingTerm(phi,RHS)

          case('KolmogorovFlow')

            call kolmogorov_flow_forcing(RHS)

          case('blades')

            call turbine_blades_sponge(RHS)
                
        endselect

       ! if(bc(4) == 'farFieldImposed') call farFieldImposedSponge(RHS)

        call EndProfRange 
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
        !$acc end parallel

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
        !$acc end parallel


        prs_grad = - gbl_rurhs

        return
end subroutine ComputeTCHPressureGradient


subroutine dump_residuals

        use FileModule
        use, intrinsic :: IEEE_Arithmetic

        implicit none
        integer  :: i,j,k
        real(rp) :: iNp
        real(rp) :: rhs_1, rhs_2, rhs_3, rhs_4, rhs_5
        real(rp) :: r_, ru, rv, rw, re
        real(rp) :: lcl_res_1, lcl_res_2, lcl_res_3, lcl_res_4, lcl_res_5
        real(rp) :: gbl_res_1, gbl_res_2, gbl_res_3, gbl_res_4, gbl_res_5
        real(rp) :: lcl_max_r_, lcl_max_ru, lcl_max_rv, lcl_max_rw, lcl_max_re
        real(rp) :: gbl_max_r_, gbl_max_ru, gbl_max_rv, gbl_max_rw, gbl_max_re
        real(rp) :: lcl_min_r_, lcl_min_ru, lcl_min_rv, lcl_min_rw, lcl_min_re
        real(rp) :: gbl_min_r_, gbl_min_ru, gbl_min_rv, gbl_min_rw, gbl_min_re
        type(FileType) :: residualFile
        
        lcl_res_1 = 0.0_rp
        lcl_res_2 = 0.0_rp
        lcl_res_3 = 0.0_rp
        lcl_res_4 = 0.0_rp
        lcl_res_5 = 0.0_rp

        lcl_max_r_ = 0.0_rp
        lcl_max_ru = 0.0_rp
        lcl_max_rv = 0.0_rp
        lcl_max_rw = 0.0_rp
        lcl_max_re = 0.0_rp

        lcl_min_r_ = huge(0.0_rp)
        lcl_min_ru = huge(0.0_rp)
        lcl_min_rv = huge(0.0_rp)
        lcl_min_rw = huge(0.0_rp)
        lcl_min_re = huge(0.0_rp)

        !$acc parallel default(present) &
        !$acc copy(lcl_res_1,lcl_res_2,lcl_res_3,lcl_res_4,lcl_res_5) &
        !$acc copy(lcl_max_r_,lcl_max_ru,lcl_max_rv,lcl_max_rw,lcl_max_re) &
        !$acc copy(lcl_min_r_,lcl_min_ru,lcl_min_rv,lcl_min_rw,lcl_min_re)
        !$acc loop gang, vector collapse(3) &
        !$acc reduction(+:lcl_res_1,lcl_res_2,lcl_res_3,lcl_res_4,lcl_res_5) &
        !$acc reduction(max:lcl_max_r_,lcl_max_ru,lcl_max_rv,lcl_max_rw,lcl_max_re) &
        !$acc reduction(min:lcl_min_r_,lcl_min_ru,lcl_min_rv,lcl_min_rw,lcl_min_re)
        do       k = sz,ez
           do    j = sy,ey
              do i = sx,ex
                        
                 rhs_1 = RHS(i,j,k,1)
                 rhs_2 = RHS(i,j,k,2)
                 rhs_3 = RHS(i,j,k,3)
                 rhs_4 = RHS(i,j,k,4)
                 rhs_5 = RHS(i,j,k,5)
                 
                 ! summing up residuals
                 lcl_res_1 = lcl_res_1 + rhs_1*rhs_1
                 lcl_res_2 = lcl_res_2 + rhs_2*rhs_2
                 lcl_res_3 = lcl_res_3 + rhs_3*rhs_3
                 lcl_res_4 = lcl_res_4 + rhs_4*rhs_4
                 lcl_res_5 = lcl_res_5 + rhs_5*rhs_5

                 ! getting max and min values of conservative variables
                 r_ = abs(phi(i,j,k,1))
                 ru = abs(phi(i,j,k,2))
                 rv = abs(phi(i,j,k,3))
                 rw = abs(phi(i,j,k,4))
                 re = abs(phi(i,j,k,5))

                 lcl_min_r_ = min(r_,lcl_min_r_)
                 lcl_min_ru = min(ru,lcl_min_ru)
                 lcl_min_rv = min(rv,lcl_min_rv)
                 lcl_min_rw = min(rw,lcl_min_rw)
                 lcl_min_re = min(re,lcl_min_re)

                 lcl_max_r_ = max(r_,lcl_max_r_)
                 lcl_max_ru = max(ru,lcl_max_ru)
                 lcl_max_rv = max(rv,lcl_max_rv)
                 lcl_max_rw = max(rw,lcl_max_rw)
                 lcl_max_re = max(re,lcl_max_re)

              enddo
           enddo
        enddo
        !$acc end parallel

        ! summing via MPI
        call MPI_allreduce(lcl_res_1, gbl_res_1, 1, MPI_RP, MPI_SUM, mpi_comm_cart, err)
        call MPI_allreduce(lcl_res_2, gbl_res_2, 1, MPI_RP, MPI_SUM, mpi_comm_cart, err)
        call MPI_allreduce(lcl_res_3, gbl_res_3, 1, MPI_RP, MPI_SUM, mpi_comm_cart, err)
        call MPI_allreduce(lcl_res_4, gbl_res_4, 1, MPI_RP, MPI_SUM, mpi_comm_cart, err)
        call MPI_allreduce(lcl_res_5, gbl_res_5, 1, MPI_RP, MPI_SUM, mpi_comm_cart, err)
        
        ! scaling residuals
        iNp = 1.0_rp/(real(nx*ny*nz,rp))

        gbl_res_1 = sqrt(gbl_res_1*iNp)*Dt
        gbl_res_2 = sqrt(gbl_res_2*iNp)*Dt
        gbl_res_3 = sqrt(gbl_res_3*iNp)*Dt
        gbl_res_4 = sqrt(gbl_res_4*iNp)*Dt
        gbl_res_5 = sqrt(gbl_res_5*iNp)*Dt

        if(ieee_is_nan(gbl_res_1)) stop 'divergency has been detected'
        if(ieee_is_nan(gbl_res_2)) stop 'divergency has been detected'
        if(ieee_is_nan(gbl_res_3)) stop 'divergency has been detected'
        if(ieee_is_nan(gbl_res_4)) stop 'divergency has been detected'
        if(ieee_is_nan(gbl_res_5)) stop 'divergency has been detected'

        ! get global max
        call MPI_allreduce(lcl_max_r_, gbl_max_r_, 1, MPI_RP, MPI_MAX, mpi_comm_cart, err)
        call MPI_allreduce(lcl_max_ru, gbl_max_ru, 1, MPI_RP, MPI_MAX, mpi_comm_cart, err)
        call MPI_allreduce(lcl_max_rv, gbl_max_rv, 1, MPI_RP, MPI_MAX, mpi_comm_cart, err)
        call MPI_allreduce(lcl_max_rw, gbl_max_rw, 1, MPI_RP, MPI_MAX, mpi_comm_cart, err)
        call MPI_allreduce(lcl_max_re, gbl_max_re, 1, MPI_RP, MPI_MAX, mpi_comm_cart, err)

        if(ieee_is_nan(gbl_max_r_)) stop 'divergency has been detected'
        if(ieee_is_nan(gbl_max_ru)) stop 'divergency has been detected'
        if(ieee_is_nan(gbl_max_rv)) stop 'divergency has been detected'
        if(ieee_is_nan(gbl_max_rw)) stop 'divergency has been detected'
        if(ieee_is_nan(gbl_max_re)) stop 'divergency has been detected'

        ! get global min
        call MPI_allreduce(lcl_min_r_, gbl_min_r_, 1, MPI_RP, MPI_MIN, mpi_comm_cart, err)
        call MPI_allreduce(lcl_min_ru, gbl_min_ru, 1, MPI_RP, MPI_MIN, mpi_comm_cart, err)
        call MPI_allreduce(lcl_min_rv, gbl_min_rv, 1, MPI_RP, MPI_MIN, mpi_comm_cart, err)
        call MPI_allreduce(lcl_min_rw, gbl_min_rw, 1, MPI_RP, MPI_MIN, mpi_comm_cart, err)
        call MPI_allreduce(lcl_min_re, gbl_min_re, 1, MPI_RP, MPI_MIN, mpi_comm_cart, err)

        if(rank == root) then
          residualFile%name = 'RESIDUALS_MONITOR'
          residualFile%dir  = trim(data_dir)
          call AppendToFile(residualFile,it)
          write(residualFile%unit,10) &
               it,            &
               gbl_res_1,     &
               gbl_res_2,     &
               gbl_res_3,     &
               gbl_res_4,     &
               gbl_res_5,     &
               gbl_max_r_,    & 
               gbl_min_r_,    &
               gbl_max_ru,    &
               gbl_min_ru,    &
               gbl_max_rv,    &
               gbl_min_rv,    &
               gbl_max_rw,    &
               gbl_min_rw,    &
               gbl_max_re,    &
               gbl_min_re

          call CloseFile(residualFile)
        endif

        10 format(I7,15e18.9)

        return
end subroutine dump_residuals
















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
        !$acc end parallel

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
        !$acc end parallel

        return
end subroutine kolmogorov_flow_forcing


subroutine turbine_blades_sponge(RHS)
        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: RHS

        real(rp) :: r_0, u_0, v_0, p_0
        real(rp) :: ru0, rv0, re0
        real(rp) :: sponge
        integer  :: i,j,k

        ! outlet conditions (to be determined from 1D or RANS models)
        r_0 = static_rho_outlet
        u_0 = vel_u_outlet
        v_0 = vel_v_outlet
        p_0 = static_prs_outlet

        ! conservative variables
        ru0 = r_0*u_0
        rv0 = r_0*v_0
        re0 = p_0/(gamma0-1.0_rp) + 0.5*r_0*(u_0*u_0 + v_0*v_0)

        !$acc parallel default(present)
        !$acc loop gang, vector, collapse(3)
        do   k = sz,ez
         do  j = sy,ey
          do i = sx,ex

             sponge = sponge_x(i)

             RHS(i,j,k,1) = RHS(i,j,k,1) - sponge*(phi(i,j,k,1) - r_0)
             RHS(i,j,k,2) = RHS(i,j,k,2) - sponge*(phi(i,j,k,2) - ru0)
             RHS(i,j,k,3) = RHS(i,j,k,3) - sponge*(phi(i,j,k,3) - rv0)
             RHS(i,j,k,5) = RHS(i,j,k,5) - sponge*(phi(i,j,k,5) - re0)

          enddo
         enddo
        enddo
        !$acc end parallel



        return
end subroutine turbine_blades_sponge



















end module rhs_module
