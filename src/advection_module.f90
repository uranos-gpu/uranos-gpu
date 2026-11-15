module advection_module
use parameters_module
use mpi_module
use storage_module
use profiling_module
use nscbc_boundary_conditions
use energy_preserving_JP_module
use energy_preserving_KEEP_module
use energy_preserving_AEC2CP_module
use energy_preserving_AECTP_module
use energy_preserving_PEPTP_module
use hybrid_wenoEP_JP_module
use hybrid_wenoEP_KEEP_module
use hybrid_wenoEP_AEC2CP_module
use hybrid_wenoEP_AECTP_module
use hybrid_wenoEP_PEPTP_module
implicit none
private

public advection_fluxes

contains
subroutine advection_fluxes(scheme, ep_id, flag)
        implicit none
        character(*), intent(in) :: scheme
        integer(1), intent(in)   :: ep_id
        integer, intent(in)      :: flag

        integer            :: istr,iend,jstr,jend,kstr,kend, ltot
        integer, parameter :: W = 1, E = 2, S = 3, N = 4
        integer, parameter :: NoOne = MPI_PROC_NULL
        integer            :: shock_recon
        integer            :: i,j,k,l

        call StartProfRange("advection_fluxes")

        selectcase(trim(scheme))
        case('hybrid_wenoEP')
          shock_recon = 1
        case('hybrid_tenoEP')
          shock_recon = 2
        case('hybrid_tenoaEP')
          shock_recon = 3
        endselect

        istr = sx
        iend = ex

        jstr = sy
        jend = ey

        kstr = sz
        kend = ez

        ltot = fd_order/2

        call set_nscbc_flux_ranges(istr, iend, jstr, jend)

        selectcase(trim(scheme))

          case('hybrid_wenoEP', 'hybrid_tenoEP', 'hybrid_tenoaEP')
            selectcase(ep_id)
            case(1)
              call hybrid_wenoEP_JP_X(weno_num,lbx,ubx,istr,iend,shock_recon, flag)
              call hybrid_wenoEP_JP_Y(weno_num,lby,uby,jstr,jend,shock_recon)
              if(dims==3) call hybrid_wenoEP_JP_Z(weno_num,lbz,ubz,kstr,kend,shock_recon)
            case(2)
              call hybrid_wenoEP_KEEP_X(weno_num,lbx,ubx,istr,iend,shock_recon, flag)
              call hybrid_wenoEP_KEEP_Y(weno_num,lby,uby,jstr,jend,shock_recon)
              if(dims==3) call hybrid_wenoEP_KEEP_Z(weno_num,lbz,ubz,kstr,kend,shock_recon)
            case(3)
              call hybrid_wenoEP_AEC2CP_X(weno_num,lbx,ubx,istr,iend,shock_recon, flag)
              call hybrid_wenoEP_AEC2CP_Y(weno_num,lby,uby,jstr,jend,shock_recon)
              if(dims==3) call hybrid_wenoEP_AEC2CP_Z(weno_num,lbz,ubz,kstr,kend,shock_recon)
            case(4)
              call hybrid_wenoEP_AECTP_X(weno_num,lbx,ubx,istr,iend,shock_recon, flag)
              call hybrid_wenoEP_AECTP_Y(weno_num,lby,uby,jstr,jend,shock_recon)
              if(dims==3) call hybrid_wenoEP_AECTP_Z(weno_num,lbz,ubz,kstr,kend,shock_recon)
            case(5)
              call hybrid_wenoEP_PEPTP_X(weno_num,lbx,ubx,istr,iend,shock_recon, flag)
              call hybrid_wenoEP_PEPTP_Y(weno_num,lby,uby,jstr,jend,shock_recon)
              if(dims==3) call hybrid_wenoEP_PEPTP_Z(weno_num,lbz,ubz,kstr,kend,shock_recon)
            case default
              if(rank == root) write(*,'(A,I0,A)') ' Blended (hybrid) schemes are implemented for'// &
                ' ep_type JP (1), KEEP (2), AEC (3=CP/4=TP) and PEP (5); got unknown ep_id=', ep_id, '.'
              call secure_stop
            endselect

          case('energy_preserving')
            selectcase(ep_id)
            case(1)
              call energy_preserving_JP_X(lbx,ubx,istr,iend)
              call energy_preserving_JP_Y(lby,uby,jstr,jend)
              if(dims==3) call energy_preserving_JP_Z(lbz,ubz,kstr,kend)
            case(2)
              call energy_preserving_KEEP_X(lbx,ubx,istr,iend)
              call energy_preserving_KEEP_Y(lby,uby,jstr,jend)
              if(dims==3) call energy_preserving_KEEP_Z(lbz,ubz,kstr,kend)
            case(3)
              call energy_preserving_AEC2CP_X(lbx,ubx,istr,iend)
              call energy_preserving_AEC2CP_Y(lby,uby,jstr,jend)
              if(dims==3) call energy_preserving_AEC2CP_Z(lbz,ubz,kstr,kend)
            case(4)
              call energy_preserving_AECTP_X(lbx,ubx,istr,iend)
              call energy_preserving_AECTP_Y(lby,uby,jstr,jend)
              if(dims==3) call energy_preserving_AECTP_Z(lbz,ubz,kstr,kend)
            case(5)
              call energy_preserving_PEPTP_X(lbx,ubx,istr,iend)
              call energy_preserving_PEPTP_Y(lby,uby,jstr,jend)
              if(dims==3) call energy_preserving_PEPTP_Z(lbz,ubz,kstr,kend)
            endselect

          case('none')
            !$acc parallel default(present)
            !$acc loop gang, vector collapse(4)
            do       l = 1,5
              do     k = sz,ez
                do   j = sy,ey
                  do i = sx,ex
                    RHS(i,j,k,l) = 0.0_rp
                  enddo
                enddo
              enddo
            enddo
          !$acc end parallel

          case default
            if(rank == root) write(*,'(A,A,A)') ' Scheme ', '"'//trim(scheme)//'"',  ' is not implemented.'
            call secure_stop

        endselect

        call set_nscbc_boundary_conditions(RHS)

        call EndProfRange

        return
end subroutine advection_fluxes

subroutine set_nscbc_flux_ranges(istr, iend, jstr, jend)
        implicit none
        integer, intent(inout) :: istr, iend, jstr, jend
        integer, parameter     :: NoOne = MPI_PROC_NULL
        integer                :: f

        do f = 1, 4
          if(my_neighbour(f) /= NoOne)  cycle
          if(.not. is_nscbc_bc(bc(f)))  cycle
          selectcase(f)
          case(1); istr = sx+1; call set_RHS_BC_x(sx)
          case(2); iend = ex-1; call set_RHS_BC_x(ex)
          case(3); jstr = sy+1
          case(4); jend = ey-1
          endselect
        enddo

end subroutine set_nscbc_flux_ranges

logical function is_nscbc_bc(bc_name)
        implicit none
        character(*), intent(in) :: bc_name
        is_nscbc_bc = .false.
        selectcase(trim(bc_name))
        case('nscbc_outflow', 'nscbc_inflow_generalized')
          is_nscbc_bc = .true.
        endselect
end function is_nscbc_bc

subroutine set_RHS_BC_x(x_i)

   implicit none
   integer, intent(in) :: x_i

   integer :: j,k,l

   call StartProfRange("set_RHS_BC_x")

   !$acc parallel default(present)
   !$acc loop gang, vector collapse(3)
   do l = 1,5
     do k = sz,ez
        do j = sy,ey
            RHS(x_i,j,k,l) = 0.0_rp
        enddo
     enddo
   enddo
   !$acc end parallel

   call EndProfRange

   return
end subroutine set_RHS_BC_x

end module advection_module
