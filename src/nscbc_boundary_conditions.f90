module nscbc_boundary_conditions
use parameters_module
use mpi_module
use mesh_module
use storage_module
use inflow_module
use profiling_module
use thermodynamics_module, only:enthalpy_t, gamma_Tref
use euler_matrices_module, only: euler_system_n
implicit none

private
public set_nscbc_boundary_conditions

contains
subroutine set_nscbc_boundary_conditions(RHS)
        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: RHS

        integer, dimension(6)      :: bound_nodes
        integer, dimension(6)      :: bound_norms
        character(1), dimension(6) :: bound_faces
        integer                    :: bnode, bnorm, f
        character(1)               :: bface

        call StartProfRange("set_nscbc_boundary_conditions")

        bound_nodes = (/ sx , ex , sy , ey , sz ,  ez/)
        bound_norms = (/ -1 ,  1 , -1 ,  1 , -1 ,  1 /)
        bound_faces = (/ 'W', 'E', 'S', 'N', 'B', 'F'/)

        do f = 1,6
           if(my_neighbour(f) == MPI_PROC_NULL) then

             bface = bound_faces(f)
             bnode = bound_nodes(f)
             bnorm = bound_norms(f)

             select case(bc(f))

               case('nscbc_inflow_generalized')
                 call nscbc_inflow_generalized(RHS,bface,bnode,bnorm)

               case('nscbc_outflow')
                 call nscbc_outflow(RHS,bface,bnode,bnorm)

             end select

           endif
        enddo

        call EndProfRange

        return
end subroutine set_nscbc_boundary_conditions

subroutine nscbc_outflow(RHS,bface,bnode,bnorm)

        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: RHS
        character(1)                             , intent(in)    :: bface
        integer                                  , intent(in)    :: bnode, bnorm

        real(rp), dimension(5,5) :: Rmatrix, Lmatrix
        real(rp), dimension(5) :: uxi, uxi_char, lambda_i

        real(rp) :: idx, idy, FU
        real(rp) :: r_, ir, u_, v_, w_, ek, p_, T_, c, cc, g_, h_
        real(rp) :: nx_, ny_, un
        real(rp) :: ddye, ddxe, ddyx, ddxx, jc, ru, rv, rEp
        integer  :: i,j,k, mm, kk, jj, b1, b2

        call StartProfRange("nscbc_outflow")

        select case(bface)
          case('E','W')

          !$acc parallel default(present)
          !$acc loop gang, vector collapse(2) private(Rmatrix,Lmatrix,uxi,uxi_char,lambda_i)
          do k    = sz,ez
             do j = sy,ey

                b1 = bnode - bnorm
                b2 = bnode - 2*bnorm

                r_ = phi(bnode,j,k,1)
                ir = 1.0_rp/r_
                u_ = U(bnode,j,k)
                v_ = V(bnode,j,k)
                w_ = W(bnode,j,k)
                p_ = P(bnode,j,k)
                g_ = GAMMA_T(bnode,j,k)

                ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
                h_ = ENTHALPY(bnode,j,k) + ek
                T_ = p_*ir
                c  = sqrt(g_ * T_)
                cc = c*c

                nx_ = nx_xi(bnode,j)
                ny_ = ny_xi(bnode,j)
                idx = m_xi (bnode,j)
                un  = nx_*u_ + ny_*v_

                call euler_system_n(nx_, ny_, u_, v_, w_, T_, g_, h_, Rmatrix, Lmatrix)
                uxi(1) =   bnorm * (0.5_rp*phi(b2,j,k,1)-2.0_rp*phi(b1,j,k,1)+1.5_rp*phi(bnode,j,k,1))*idx
                uxi(2) =   bnorm * (0.5_rp*phi(b2,j,k,2)-2.0_rp*phi(b1,j,k,2)+1.5_rp*phi(bnode,j,k,2))*idx
                uxi(3) =   bnorm * (0.5_rp*phi(b2,j,k,3)-2.0_rp*phi(b1,j,k,3)+1.5_rp*phi(bnode,j,k,3))*idx
                uxi(4) =   bnorm * (0.5_rp*phi(b2,j,k,4)-2.0_rp*phi(b1,j,k,4)+1.5_rp*phi(bnode,j,k,4))*idx
                uxi(5) =   bnorm * (0.5_rp*phi(b2,j,k,5)-2.0_rp*phi(b1,j,k,5)+1.5_rp*phi(bnode,j,k,5))*idx
                do kk = 1, 5
                    uxi_char(kk) = 0.0_rp
                      do jj = 1, 5
                          uxi_char(kk) = uxi_char(kk) + Lmatrix(kk, jj) * uxi(jj)
                      end do
                end do
                lambda_i(1) = 0.5_rp * (un - c + bnorm * abs(un - c))
                lambda_i(2) = 0.5_rp * (un     + bnorm * abs(un)    )
                lambda_i(3) = lambda_i(2)
                lambda_i(4) = lambda_i(2)
                lambda_i(5) = 0.5_rp * (un + c + bnorm * abs(un + c))

                uxi_char = lambda_i * uxi_char

                do mm = 1, 5
                  FU = 0.0_rp
                    !$acc loop seq
                    do kk = 1, 5
                      FU = FU + Rmatrix(mm,kk) * uxi_char(kk)
                    end do
                  RHS(bnode,j,k,mm) = RHS(bnode,j,k,mm) - FU
                end do

                ddye = central_1(1)*(dy_deta(bnode+1,j)-dy_deta(bnode-1,j)) &
                     + central_1(2)*(dy_deta(bnode+2,j)-dy_deta(bnode-2,j)) &
                     + central_1(3)*(dy_deta(bnode+3,j)-dy_deta(bnode-3,j))
                ddxe = central_1(1)*(dx_deta(bnode+1,j)-dx_deta(bnode-1,j)) &
                     + central_1(2)*(dx_deta(bnode+2,j)-dx_deta(bnode-2,j)) &
                     + central_1(3)*(dx_deta(bnode+3,j)-dx_deta(bnode-3,j))
                jc   = Jac2D(bnode,j)
                ru   = phi(bnode,j,k,2)
                rv   = phi(bnode,j,k,3)
                rEp  = phi(bnode,j,k,5) + p_
                RHS(bnode,j,k,1) = RHS(bnode,j,k,1) - jc*( ru        *ddye - rv        *ddxe)
                RHS(bnode,j,k,2) = RHS(bnode,j,k,2) - jc*((ru*u_+p_) *ddye - (rv*u_)   *ddxe)
                RHS(bnode,j,k,3) = RHS(bnode,j,k,3) - jc*((ru*v_)    *ddye - (rv*v_+p_)*ddxe)
                RHS(bnode,j,k,4) = RHS(bnode,j,k,4) - jc*((ru*w_)    *ddye - (rv*w_)   *ddxe)
                RHS(bnode,j,k,5) = RHS(bnode,j,k,5) - jc*((u_*rEp)   *ddye - (v_*rEp)  *ddxe)

             enddo
          enddo
          !$acc end parallel

        case('N','S')

          !$acc parallel default(present)
          !$acc loop gang, vector collapse(2) private(Rmatrix,Lmatrix,uxi,uxi_char,lambda_i)
          do k    = sz,ez
             do i = sx,ex

                b1 = bnode - bnorm
                b2 = bnode - 2*bnorm

                r_ = phi(i,bnode,k,1)
                ir = 1.0_rp/r_
                u_ = U(i,bnode,k)
                v_ = V(i,bnode,k)
                w_ = W(i,bnode,k)
                p_ = P(i,bnode,k)
                g_ = GAMMA_T(i,bnode,k)

                ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
                h_ = ENTHALPY(i,bnode,k) + ek
                T_ = p_*ir
                c  = sqrt(g_* T_)
                cc = c*c

                nx_ = nx_eta(i,bnode)
                ny_ = ny_eta(i,bnode)
                idy = m_eta (i,bnode)
                un  = nx_*u_ + ny_*v_

                call euler_system_n(nx_, ny_, u_, v_, w_, T_, g_, h_, Rmatrix, Lmatrix)

                uxi(1) =   bnorm * (0.5_rp*phi(i,b2,k,1)-2.0_rp*phi(i,b1,k,1)+1.5_rp*phi(i,bnode,k,1))*idy
                uxi(2) =   bnorm * (0.5_rp*phi(i,b2,k,2)-2.0_rp*phi(i,b1,k,2)+1.5_rp*phi(i,bnode,k,2))*idy
                uxi(3) =   bnorm * (0.5_rp*phi(i,b2,k,3)-2.0_rp*phi(i,b1,k,3)+1.5_rp*phi(i,bnode,k,3))*idy
                uxi(4) =   bnorm * (0.5_rp*phi(i,b2,k,4)-2.0_rp*phi(i,b1,k,4)+1.5_rp*phi(i,bnode,k,4))*idy
                uxi(5) =   bnorm * (0.5_rp*phi(i,b2,k,5)-2.0_rp*phi(i,b1,k,5)+1.5_rp*phi(i,bnode,k,5))*idy
                do kk = 1, 5
                    uxi_char(kk) = 0.0_rp
                    do jj = 1, 5
                        uxi_char(kk) = uxi_char(kk) + Lmatrix(kk, jj) * uxi(jj)
                    end do
                end do
                lambda_i(1) = 0.5_rp * (un - c + bnorm * abs(un - c))
                lambda_i(2) = 0.5_rp * (un     + bnorm * abs(un)    )
                lambda_i(3) = lambda_i(2)
                lambda_i(4) = lambda_i(2)
                lambda_i(5) = 0.5_rp * (un + c + bnorm * abs(un + c))

                uxi_char = lambda_i * uxi_char

                do mm = 1, 5
                  FU = 0.0_rp
                    !$acc loop seq
                    do kk = 1, 5
                      FU = FU + Rmatrix(mm,kk) * uxi_char(kk)
                    end do
                  RHS(i,bnode,k,mm) = RHS(i,bnode,k,mm) - FU
                end do

                ddyx = central_1(1)*(dy_dxi(i,bnode+1)-dy_dxi(i,bnode-1)) &
                     + central_1(2)*(dy_dxi(i,bnode+2)-dy_dxi(i,bnode-2)) &
                     + central_1(3)*(dy_dxi(i,bnode+3)-dy_dxi(i,bnode-3))
                ddxx = central_1(1)*(dx_dxi(i,bnode+1)-dx_dxi(i,bnode-1)) &
                     + central_1(2)*(dx_dxi(i,bnode+2)-dx_dxi(i,bnode-2)) &
                     + central_1(3)*(dx_dxi(i,bnode+3)-dx_dxi(i,bnode-3))
                jc   = Jac2D(i,bnode)
                ru   = phi(i,bnode,k,2)
                rv   = phi(i,bnode,k,3)
                rEp  = phi(i,bnode,k,5) + p_
                RHS(i,bnode,k,1) = RHS(i,bnode,k,1) - jc*( rv        *ddxx - ru        *ddyx)
                RHS(i,bnode,k,2) = RHS(i,bnode,k,2) - jc*((rv*u_)    *ddxx - (ru*u_+p_)*ddyx)
                RHS(i,bnode,k,3) = RHS(i,bnode,k,3) - jc*((rv*v_+p_) *ddxx - (ru*v_)   *ddyx)
                RHS(i,bnode,k,4) = RHS(i,bnode,k,4) - jc*((rv*w_)    *ddxx - (ru*w_)   *ddyx)
                RHS(i,bnode,k,5) = RHS(i,bnode,k,5) - jc*((v_*rEp)   *ddxx - (u_*rEp)  *ddyx)

             enddo
          enddo
          !$acc end parallel

        case default
          print*, 'nscbc_outflow is not implemented for face ', trim(bface)
          stop
          continue

        endselect

        call EndProfRange

        return
end subroutine nscbc_outflow

subroutine nscbc_inflow_generalized(RHS,bface,bnode,bnorm)
        implicit none

        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: RHS
        character(1)                             , intent(in)    :: bface
        integer                                  , intent(in)    :: bnode, bnorm

        real(rp), dimension(5,5) :: Rmatrix, Lmatrix
        real(rp), dimension(5) :: uxi,uxi_char, uxo, uxo_char, lambda_o, lambda_i

        real(rp) :: idx, r0, u0, v0, w0, p0, T0, FU
        real(rp) :: rm, um, vm, wm, Tm
        real(rp) ::     uf, vf, wf, th, sg
        real(rp) :: q1, q2, q3, q4, q5
        real(rp) :: r_, ir, u_, v_, w_, ek, p_, T_, c, g_, h_, cc
        real(rp) :: igm1
        real(rp) :: nx_, ny_, un
        real(rp) :: ddye, ddxe, jc, ru, rv, rEp
        integer  :: j,k, mm, kk, jj, b1, b2

        call StartProfRange("nscbc_inflow_generalized")

        igm1  = 1.0_rp/(gamma_Tref - 1.0_rp)

        if(turb_inflow) call inflow_turbulent_update()

        select case(bface)
          case('W')

          !$acc parallel default(present)
          !$acc loop gang, vector collapse(2) private(Rmatrix,Lmatrix,uxi,uxi_char,uxo,uxo_char,lambda_o,lambda_i)
          do k    = sz,ez
             do j = sy,ey

                rm = iflow_mean(j,1)
                um = iflow_mean(j,2)
                vm = iflow_mean(j,3)
                wm = iflow_mean(j,4)
                Tm = iflow_mean(j,5)

                uf = 0.0_rp
                vf = 0.0_rp
                wf = 0.0_rp
                T0 = Tm
                if(j > sy) then
                  uf = u_inf*iflow_turb(1,j,k)
                  vf = u_inf*iflow_turb(2,j,k)
                  wf = u_inf*iflow_turb(3,j,k)
                  th = iflow_dTdu(j)*uf/Tm
                  sg = iflow_dTdu(j)*u_inf*iflow_urms(j)/Tm
                  T0 = Tm*exp(th - 0.5_rp*sg*sg)
                endif

                p0 = rm * Tm
                u0 = um + uf
                v0 = vm + vf
                w0 = wm + wf
                r0 = p0/T0

                q1 = r0
                q2 = r0*u0
                q3 = r0*v0
                q4 = r0*w0
                if (thermally_perfect_gas > 0 ) then
                  q5 = r0 * enthalpy_t(T0) - p0 + 0.5_rp*r0*(u0*u0 + v0*v0 + w0*w0)
                else
                  q5 = p0*igm1 + 0.5_rp*r0*(u0*u0 + v0*v0 + w0*w0)
                endif

                b1 = bnode - bnorm
                b2 = bnode - 2*bnorm

                r_ = phi(bnode,j,k,1)
                ir = 1.0_rp/r_
                u_ = U(bnode,j,k)
                v_ = V(bnode,j,k)
                w_ = W(bnode,j,k)
                p_ = P(bnode,j,k)

                g_ = GAMMA_T(bnode,j,k)

                ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
                h_ = ENTHALPY(bnode,j,k) + ek
                T_ = p_*ir
                c  = sqrt(g_ * T_)
                cc = c*c

                nx_ = nx_xi(bnode,j)
                ny_ = ny_xi(bnode,j)
                idx = m_xi (bnode,j)
                un  = nx_*u_ + ny_*v_

                call euler_system_n(nx_, ny_, u_, v_, w_, T_, g_, h_, Rmatrix, Lmatrix)

                uxi(1) =   bnorm * (0.5_rp*phi(b2,j,k,1)-2.0_rp*phi(b1,j,k,1)+1.5_rp*phi(bnode,j,k,1))*idx
                uxi(2) =   bnorm * (0.5_rp*phi(b2,j,k,2)-2.0_rp*phi(b1,j,k,2)+1.5_rp*phi(bnode,j,k,2))*idx
                uxi(3) =   bnorm * (0.5_rp*phi(b2,j,k,3)-2.0_rp*phi(b1,j,k,3)+1.5_rp*phi(bnode,j,k,3))*idx
                uxi(4) =   bnorm * (0.5_rp*phi(b2,j,k,4)-2.0_rp*phi(b1,j,k,4)+1.5_rp*phi(bnode,j,k,4))*idx
                uxi(5) =   bnorm * (0.5_rp*phi(b2,j,k,5)-2.0_rp*phi(b1,j,k,5)+1.5_rp*phi(bnode,j,k,5))*idx
                uxo(1) = - bnorm * (phi(bnode,j,k,1)-q1)*idx
                uxo(2) = - bnorm * (phi(bnode,j,k,2)-q2)*idx
                uxo(3) = - bnorm * (phi(bnode,j,k,3)-q3)*idx
                uxo(4) = - bnorm * (phi(bnode,j,k,4)-q4)*idx
                uxo(5) = - bnorm * (phi(bnode,j,k,5)-q5)*idx
                do kk = 1, 5
                    uxi_char(kk) = 0.0_rp
                    uxo_char(kk) = 0.0_rp
                    do jj = 1, 5
                      uxi_char(kk) = uxi_char(kk) + Lmatrix(kk, jj) * uxi(jj)
                      uxo_char(kk) = uxo_char(kk) + Lmatrix(kk, jj) * uxo(jj)
                    end do
                end do
                lambda_i(1) = 0.5_rp * (un - c + bnorm * abs(un - c))
                lambda_i(2) = 0.5_rp * (un     + bnorm * abs(un)    )
                lambda_i(3) = lambda_i(2)
                lambda_i(4) = lambda_i(2)
                lambda_i(5) = 0.5_rp * (un + c + bnorm * abs(un + c))

                lambda_o(1) = 0.5_rp * (un - c - bnorm * abs(un - c))
                lambda_o(2) = 0.5_rp * (un     - bnorm * abs(un)    )
                lambda_o(3) = lambda_o(2)
                lambda_o(4) = lambda_o(2)
                lambda_o(5) = 0.5_rp * (un + c - bnorm * abs(un + c))

                uxi_char = lambda_i * uxi_char
                uxo_char = lambda_o * uxo_char
                do mm = 1, 5
                  FU = 0.0_rp
                    !$acc loop seq
                    do kk = 1, 5
                      FU = FU + Rmatrix(mm,kk) * (uxi_char(kk) + uxo_char(kk))
                    end do
                  RHS(bnode,j,k,mm) = RHS(bnode,j,k,mm) - FU
                end do

                ddye = central_1(1)*(dy_deta(bnode+1,j)-dy_deta(bnode-1,j)) &
                     + central_1(2)*(dy_deta(bnode+2,j)-dy_deta(bnode-2,j)) &
                     + central_1(3)*(dy_deta(bnode+3,j)-dy_deta(bnode-3,j))
                ddxe = central_1(1)*(dx_deta(bnode+1,j)-dx_deta(bnode-1,j)) &
                     + central_1(2)*(dx_deta(bnode+2,j)-dx_deta(bnode-2,j)) &
                     + central_1(3)*(dx_deta(bnode+3,j)-dx_deta(bnode-3,j))
                jc   = Jac2D(bnode,j)
                ru   = phi(bnode,j,k,2)
                rv   = phi(bnode,j,k,3)
                rEp  = phi(bnode,j,k,5) + p_
                RHS(bnode,j,k,1) = RHS(bnode,j,k,1) - jc*( ru        *ddye - rv        *ddxe)
                RHS(bnode,j,k,2) = RHS(bnode,j,k,2) - jc*((ru*u_+p_) *ddye - (rv*u_)   *ddxe)
                RHS(bnode,j,k,3) = RHS(bnode,j,k,3) - jc*((ru*v_)    *ddye - (rv*v_+p_)*ddxe)
                RHS(bnode,j,k,4) = RHS(bnode,j,k,4) - jc*((ru*w_)    *ddye - (rv*w_)   *ddxe)
                RHS(bnode,j,k,5) = RHS(bnode,j,k,5) - jc*((u_*rEp)   *ddye - (v_*rEp)  *ddxe)

             enddo
          enddo
          !$acc end parallel

        case default
          print*, 'nscbc_inflow_generalized is not implemented for face ', trim(bface)
          stop
          continue

        endselect

        call EndProfRange

        return
end subroutine nscbc_inflow_generalized

end module nscbc_boundary_conditions
