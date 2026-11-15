module rhs_module
use parameters_module
use storage_module
use profiling_module
use thermodynamics_module, only: enthalpy_t, gamma_Tref

implicit none
private
public rhs_navier_stokes, mpi_wait_conservatives

contains
subroutine rhs_navier_stokes
        use advection_module, only: advection_fluxes
        use viscous_module, only: viscous_fluxes

        implicit none

        call StartProfRange("rhs_navier_stokes")

        call advection_fluxes(scheme, ep_id, thermally_perfect_gas)

        if(mu_jump_ratio /= 1.0_rp) call impose_mu_jump

        if(viscous) call viscous_fluxes(dims)

        call forcing_terms

        call EndProfRange

        return
end subroutine rhs_navier_stokes

subroutine impose_mu_jump
        use mesh_module, only: yy
        implicit none
        integer :: i,j,k

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do       k = lbz,ubz
           do    j = lby,uby
              do i = lbx,ubx
                 if(abs(yy(i,j)) < mu_jump_loc) then
                   VIS(i,j,k) = mu_jump_ratio
                 else
                   VIS(i,j,k) = 1.0_rp
                 endif
              enddo
           enddo
        enddo
        !$acc end parallel

        return
end subroutine impose_mu_jump

subroutine mpi_wait_conservatives
        implicit none

        call StartProfRange("mpi_wait_conservatives")

        call mpi_wait_procs(req_array_yz)

        phi(:,sy-GN:sy-1,:,:) = phi_bfr_recv_S
        phi(:,ey+1:ey+GN,:,:) = phi_bfr_recv_N

        phi(:,:,sz-GN:sz-1,:) = phi_bfr_recv_B
        phi(:,:,ez+1:ez+GN,:) = phi_bfr_recv_F

        call mpi_wait_procs(req_array_yz)

        call EndProfRange

        return
end subroutine mpi_wait_conservatives

subroutine forcing_terms
        implicit none

        call StartProfRange("forcing_terms")

        selectcase(ic)
          case('turbulent_channel')

            call ComputeTCHPressureGradient(phi,RHS,force_turbulent_channel)

          case('KolmogorovFlow')

            call kolmogorov_flow_forcing(RHS)

        endselect

        call EndProfRange
return
end subroutine forcing_terms

subroutine ComputeTCHPressureGradient(phi,RHS,prs_grad)

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

subroutine kolmogorov_flow_forcing(RHS)

        use mesh_module, only: yy
        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: RHS

        real(rp) :: fx
        integer  :: i,j,k

       !$acc parallel default(present)
       !$acc loop gang, vector collapse(3)
        do       k = sz,ez
           do    j = sy,ey
              do i = sx,ex

                 fx = kolmogorov_forcing*u_inf*sin(8*pi*yy(i,j)/Ly)
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

end module rhs_module
