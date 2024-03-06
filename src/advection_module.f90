module advection_module
use parameters_module
use mpi_module
use storage_module
use flux_module
use profiling_module
use nscbc_boundary_conditions
!$ use omp_lib
implicit none
private

public advection_fluxes, CharacteristicRHSx, CharacteristicRHSy, &
       rhs_linear_ode, &
       central_finite_differences, rhs_linear_advection

contains
subroutine advection_fluxes(scheme)
        implicit none
        character(*), intent(in) :: scheme

        ! local declarations
        integer            :: istr,iend,jstr,jend,kstr,kend, ltot
        integer, parameter :: W = 1, E = 2, S = 3, N = 4
        integer, parameter :: NoOne = MPI_PROC_NULL
        integer            :: shock_recon, f

        integer, dimension(6)      :: bound_nodes
        integer, dimension(6)      :: bound_norms
        character(1), dimension(6) :: bound_faces
        integer                    :: bnode, bnorm
        character(1)               :: bface
        real(rp)                   :: Twall

        call StartProfRange("advection_fluxes")

        ! set the bound nodes, norms and faces
        bound_nodes = (/ sx , ex , sy , ey , sz ,  ez/)
        bound_norms = (/ -1 ,  1 , -1 ,  1 , -1 ,  1 /)
        bound_faces = (/ 'W', 'E', 'S', 'N', 'B', 'F'/)

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
        
        ! nscbc outflow
        if(my_neighbour(E) == NoOne .and. bc(E) == 'nscbc_outflow') iend = ex-1
        if(my_neighbour(N) == NoOne .and. bc(N) == 'nscbc_outflow') jend = ey-1
        if(my_neighbour(W) == NoOne .and. bc(W) == 'nscbc_inflow') istr = sx+1

        ! nscbc wall
        if(my_neighbour(S) == NoOne .and. bc(S) == 'nscbc_isothermal_wall') jstr = sy+1
        if(my_neighbour(N) == NoOne .and. bc(N) == 'nscbc_isothermal_wall') jend = ey-1

        selectcase(trim(scheme))

          case('hybrid_wenoEP', 'hybrid_tenoEP', 'hybrid_tenoaEP')
            call hybrid_wenoEPx(weno_num,lbx,ubx,istr,iend,shock_recon)
            call hybrid_wenoEPy(weno_num,lby,uby,jstr,jend,shock_recon)
            if(dims==3) then
            call hybrid_wenoEPz(weno_num,lbz,ubz,kstr,kend,shock_recon)
            endif

          case('energy_preserving')
            call energy_preservingX(lbx,ubx,ltot,istr,iend,tilde_op_x,pri_1D_x,RHS)
            call energy_preservingY(lby,uby,ltot,jstr,jend,tilde_op_y,pri_1D_y,RHS)
            if(dims==3)then
            call energy_preservingZ(lbz,ubz,ltot,kstr,kend,tilde_op_z,pri_1D_z,RHS)
            endif

          case('none')
            continue !zeroing the RHS

          case default
            if(rank == root) write(*,'(A,A,A)') ' Scheme ', '"'//trim(scheme)//'"',  ' is not implemented.'
            call secure_stop

        endselect

        ! nscbc outflow
        if(my_neighbour(E) == NoOne .and. bc(E) == 'nscbc_outflow') call BCRelax_X(RHS)
        if(my_neighbour(N) == NoOne .and. bc(N) == 'nscbc_outflow') call BCRelax_Y(RHS)
        if(my_neighbour(W) == NoOne .and. bc(W) == 'nscbc_inflow') call BCRelax_inflowX(RHS)
        
        do f = 1,6
           if(my_neighbour(f) == MPI_PROC_NULL) then

             bnode = bound_nodes(f)
             bnorm = bound_norms(f)
             bface = bound_faces(f)

             select case(bc(f))

               case('nscbc_isothermal_wall') 
                 Twall = 1.0_rp
                 call refl_wall(phi,bnode,bnorm,bface,Twall)      

             end select

           endif
        enddo

        call EndProfRange

        return
end subroutine advection_fluxes




!function compute_gang_size(dir, ltot) result(size)
!  implicit none
!  integer, intent(in) :: dir, ltot
!
!  integer :: loop_bound, array_size
!  integer, parameter :: max_memory=2147483647
!  integer :: size
!
!  if (dir .eq. 1 ) then ! X
!     loop_bound = (ez-sz+1)*(ey-sy+1)
!     array_size = (ubx-lbx+1)*ltot*5
!  else if (dir .eq. 2 ) then ! Y
!     loop_bound = (ez-sz+1)*(ex-sx+1)
!     array_size = (uby-lby+1)*ltot*5
!  else if (dir .eq. 3 ) then ! Y
!     loop_bound = (ex-sx+1)*(ey-sy+1)
!     array_size = (ubz-lbz+1)*ltot*5
!  else
!     stop "Invalid input to compute_gang_size"
!  end if
!
!  size = max_memory/8/array_size
!  size = min(size, loop_bound)
!  size = min(size, 64*1024-1) ! Limit to 64k gangs
!
!end function compute_gang_size



subroutine hybrid_wenoEPx(weno_num,lbx,ubx,istr,iend,shock_recon)

        implicit none

        integer, parameter :: vec_size  = 32
        integer, parameter :: chunkSize = 2*vec_size - 2*GN

        integer, intent(in) :: lbx,ubx
        integer, intent(in) :: weno_num
        integer, intent(in) :: istr,iend
        integer, intent(in) :: shock_recon

        real(rp),   dimension(3,1-GN:chunkSize+GN,5) :: tilde_op
        real(rp),   dimension(1-GN:chunkSize+GN,6)   :: pri_1D
        real(rp),   dimension(1-GN:chunkSize+GN,5)   :: phi_arr, flx_arr, flx
        integer(1), dimension(1-GN:chunkSize+GN)     :: weno_flag_chunk

        ! local declarations
        real(rp), parameter    :: gm1 = gamma0-1.0_rp
        real(rp), parameter    :: hgm = gamma0/gm1

        real(rp), dimension(5,-3:4) :: fvp, fvm, fp, fm
        real(rp), dimension(5)      :: fhatp, fhatm
        real(rp), dimension(5,5)    :: right, left
        real(rp), dimension(-3:3)   :: central_1_p

        real(rp) :: phi1, phi2, phi3, phi4, phi5, inner_sum
        real(rp) :: sqrt_rho0, sqrt_rho1, ff, pp, hFn
        real(rp) :: uroe, vroe, wroe, hroe, ekroe, croe, den_roe

        real(rp), parameter :: one8 = 1.0_rp/8.0_rp
        real(rp) :: ir, u_, v_, w_, ek, p_, idx
        real(rp) :: weight, lambda_max
        real(rp) :: vnc, hgU, hgV, hgW, hgE, cc, un_c2
        integer  :: i,j,k,l,m,il, n, ii, s
        integer  :: str_x, end_x, igbl, istr_x, iend_x
        integer(1), parameter :: one = 1

        call StartProfRange("hybrid_wenoEPx")

        !$acc parallel vector_length(vec_size) default(present) &
        !$acc private(tilde_op,pri_1D,phi_arr) &
        !$acc private(flx_arr,flx,weno_flag_chunk) &
        !$acc create(central_1_p)
        !$acc cache(central_1_p,tilde_op,pri_1D,phi_arr,flx_arr,flx,weno_flag_chunk)
        !$acc loop gang collapse(2)
        do k   = sz,ez
          do j = sy,ey
             central_1_p = 2*central_1

             do str_x = lbx,ubx,chunkSize
                end_x = min(str_x+ChunkSize, ubx-2*GN+1) - str_x

                ! correction for NSCBC
                istr_x = max(str_x+GN       , istr     ) - str_x - GN + 1
                iend_x = min(str_x+ChunkSize, iend-GN+1) - str_x
                if(iend_x - istr_x > -1) then
                !$acc loop vector
                do i = 1-GN,end_x+GN
                   igbl = i + str_x-1+GN

                   phi1 = phi(igbl,j,k,1)
                   phi2 = phi(igbl,j,k,2)
                   phi3 = phi(igbl,j,k,3)
                   phi4 = phi(igbl,j,k,4)
                   phi5 = phi(igbl,j,k,5)

                   phi_arr(i,1) = phi1
                   phi_arr(i,2) = phi2
                   phi_arr(i,3) = phi3
                   phi_arr(i,4) = phi4
                   phi_arr(i,5) = phi5

                   ir = 1.0_rp/phi1
                   u_ = phi2*ir
                   v_ = phi3*ir
                   w_ = phi4*ir
                   ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
                   p_ = gm1 * (phi5 - phi1*ek)

                   flx_arr(i,1) =  phi2
                   flx_arr(i,2) =  phi2   * u_ + p_
                   flx_arr(i,3) =  phi2   * v_
                   flx_arr(i,4) =  phi2   * w_
                   flx_arr(i,5) = (phi5+p_)*u_

                   pri_1D(i,1) = p_
                   pri_1D(i,2) = u_
                   pri_1D(i,3) = v_
                   pri_1D(i,4) = w_
                   pri_1D(i,5) = hgm*p_*ir + ek
                   pri_1D(i,6) = phi1

                   weno_flag_chunk(i) = iand(weno_flag_xyz(igbl,j,k),one)

                enddo

             !$acc loop vector collapse(2)
             do i = 1-GN,end_x+GN-3 !lbx,ubx-3
                do l = 1,3

                   il = i+l

                   weight = one8 * (pri_1D(i,6) + pri_1D(il,6)) &
                                 * (pri_1D(i,2) + pri_1D(il,2))

                   tilde_op(l,i,1) = 2*weight
                   tilde_op(l,i,2) = weight*(pri_1D(i,2) + pri_1D(il,2)) + 0.5_rp*(pri_1D(i,1) + pri_1D(il,1))
                   tilde_op(l,i,3) = weight*(pri_1D(i,3) + pri_1D(il,3))
                   tilde_op(l,i,4) = weight*(pri_1D(i,4) + pri_1D(il,4))
                   tilde_op(l,i,5) = weight*(pri_1D(i,5) + pri_1D(il,5))

                enddo
             enddo

            ! !$acc loop vector
            ! do i = 0,end_x ! sx-1,ex
            !    ishock(i) = 1
            !    !$acc loop seq
            !    do ii = i-weno_num+1,i+weno_num
            !       ishock(i) = min(weno_flag_chunk(ii),ishock(i))
            !    enddo
            ! enddo

             !$acc loop vector
             do i = 0,end_x ! sx-1,ex
             if (weno_flag_chunk(i) == weno_smooth) then

                do n = 1,5
                   hFn = 0.0_rp
                   do l = 1, 3
                      inner_sum = 0.0_rp
                      do m = 0, l-1
                         inner_sum = inner_sum + tilde_op(l,i-m,n)
                      enddo

                      hFn = hFn + central_1_p(l)*inner_sum
                   enddo
                   flx(i,n) = hFn
                enddo
             endif
             enddo

             !$acc loop vector private(right,left,fvp,fvm,fp,fm,fhatp,fhatm)
             do i = 0,end_x ! sx-1,ex
             if (weno_flag_chunk(i) == weno_shock) then

                !
                ! ROE MEAN STATE
                !
                sqrt_rho0 = sqrt(pri_1D(i,6))
                sqrt_rho1 = sqrt(pri_1D(i+1,6))
                den_roe   = 1.0_rp/(sqrt_rho0 + sqrt_rho1)

                uroe = den_roe*(pri_1D(i,2)*sqrt_rho1+pri_1D(i+1,2)*sqrt_rho0)
                vroe = den_roe*(pri_1D(i,3)*sqrt_rho1+pri_1D(i+1,3)*sqrt_rho0)
                wroe = den_roe*(pri_1D(i,4)*sqrt_rho1+pri_1D(i+1,4)*sqrt_rho0)
                hroe = den_roe*(pri_1D(i,5)*sqrt_rho1+pri_1D(i+1,5)*sqrt_rho0)
                ekroe= 0.5_rp*(uroe*uroe + vroe*vroe + wroe*wroe)
                croe = sqrt(gm1*(hroe-ekroe))
                !
                ! FILL EIGEN MATRIXES
                !
                vnc = uroe*croe
                hgU = gm1*uroe
                hgV = gm1*vroe
                hgW = gm1*wroe
                hgE = gm1*ekroe

                right(1,1)=1.0_rp
                right(1,2)=1.0_rp
                right(1,3)=1.0_rp
                right(1,4)=0.0_rp
                right(1,5)=0.0_rp
                !
                right(2,1)=uroe-croe
                right(2,2)=uroe
                right(2,3)=uroe+croe
                right(2,4)=0.0_rp
                right(2,5)=0.0_rp
                !
                right(3,1)=vroe
                right(3,2)=vroe
                right(3,3)=vroe
                right(3,4)=-1.0_rp
                right(3,5)=0.0_rp
                !
                right(4,1)=wroe
                right(4,2)=wroe
                right(4,3)=wroe
                right(4,4)=0.0_rp
                right(4,5)=1.0_rp
                !
                right(5,1)=hroe-vnc
                right(5,2)=ekroe
                right(5,3)=hroe+vnc
                right(5,4)=-vroe
                right(5,5)=wroe

                cc    = croe*croe
                un_c2 = 1.0_rp/cc

                left(1,1)=0.5_rp*un_c2*(hgE + vnc)
                left(1,2)=0.5_rp*un_c2*(-hgU-croe)
                left(1,3)=0.5_rp*un_c2*(-hgV)
                left(1,4)=0.5_rp*un_c2*(-hgW)
                left(1,5)=0.5_rp*un_c2*gm1
                !
                left(2,1)=un_c2*(cc  - hgE)
                left(2,2)=un_c2*hgU
                left(2,3)=un_c2*hgV
                left(2,4)=un_c2*hgW
                left(2,5)=un_c2*(-gm1)
                !
                left(3,1)=0.5_rp*un_c2*(hgE - vnc)
                left(3,2)=0.5_rp*un_c2*(-hgU+croe)
                left(3,3)=0.5_rp*un_c2*(-hgV)
                left(3,4)=0.5_rp*un_c2*(-hgW)
                left(3,5)=0.5_rp*un_c2*gm1
                !
                left(4,1)=vroe
                left(4,2)=0.0_rp
                left(4,3)=-1.0_rp
                left(4,4)=0.0_rp
                left(4,5)=0.0_rp
                !
                left(5,1)=-wroe
                left(5,2)=0.0_rp
                left(5,3)=0.0_rp
                left(5,4)=1.0_rp
                left(5,5)=0.0_rp
                !
                ! LAX FRIEDRICHS
                !
                lambda_max = max(abs(uroe-croe),abs(uroe+croe))

                do n = 1,5
                   do ii = -weno_num+1,weno_num
                      ff = flx_arr(i+ii,n)
                      pp = phi_arr(i+ii,n) * lambda_max

                      fvp(n,ii) = 0.5_rp * (ff + pp)
                      fvm(n,ii) = 0.5_rp * (ff - pp)
                   enddo
                enddo
                !
                ! WENO - TENO RECONSTRUCTION
                !

                do s = -weno_num+1, weno_num
                   do m = 1,5
                      fp(m,s) = 0.0_rp
                      fm(m,s) = 0.0_rp
                      do l = 1,5
                      fp(m,s) = fp(m,s) + left(m,l)*fvp(l,s)
                      fm(m,s) = fm(m,s) + left(m,l)*fvm(l,s)
                      enddo
                   enddo
                enddo

                call weno_reconstruction(weno_num,shock_recon,fp,fm,fhatp,fhatm)

                ! return in the orinal system
                do l = 1,5
                   flx(i,l) = 0.0_rp
                   do m = 1,5
                      flx(i,l) = flx(i,l) + right(l,m) * (fhatp(m) + fhatm(m))
                   enddo
                enddo
             endif
             enddo

             !$acc loop vector
             !do i = istr,iend
             do i = istr_x,iend_x

                igbl = i + str_x-1+GN
                
                idx = xstep_i(igbl)
                RHS(igbl,j,k,1) = - idx * (flx(i,1) - flx(i-1,1))
                RHS(igbl,j,k,2) = - idx * (flx(i,2) - flx(i-1,2))
                RHS(igbl,j,k,3) = - idx * (flx(i,3) - flx(i-1,3))
                RHS(igbl,j,k,4) = - idx * (flx(i,4) - flx(i-1,4))
                RHS(igbl,j,k,5) = - idx * (flx(i,5) - flx(i-1,5))

             enddo
             endif
             enddo ! chuncks
          enddo ! j
        enddo ! k
        !$acc end parallel

        !$omp end parallel do

        call EndProfRange

        return
end subroutine hybrid_wenoEPx

subroutine hybrid_wenoEPy(weno_num,lby,uby,jstr,jend,shock_recon)

        implicit none
        
        integer, parameter :: vec_size  = 64
        integer, parameter :: chunkSize = vec_size - 2*GN

        integer, intent(in) :: lby,uby
        integer, intent(in) :: weno_num
        integer, intent(in) :: jstr,jend
        integer, intent(in) :: shock_recon

        real(rp),   dimension(3,1-GN:chunkSize+GN,5) :: tilde_op
        real(rp),   dimension(1-GN:chunkSize+GN,6)   :: pri_1D
        real(rp),   dimension(1-GN:chunkSize+GN,5)   :: flx_arr, flx
        real(rp),   dimension(0:3,1-GN:chunkSize+GN,5)   :: phi_arr
        integer(1), dimension(1-GN:chunkSize+GN)     :: weno_flag_chunk

        ! local declarations
        real(rp), parameter    :: gm1 = gamma0-1.0_rp
        real(rp), parameter    :: hgm = gamma0/gm1

        real(rp), dimension(5,-3:4) :: fvp, fvm, fp, fm
        real(rp), dimension(5)      :: fhatp, fhatm
        real(rp), dimension(5,5)    :: right, left
        real(rp), dimension(-3:3)   :: central_1_p

        real(rp) :: phi1, phi2, phi3, phi4, phi5, inner_sum
        real(rp) :: sqrt_rho0, sqrt_rho1, ff, pp, hGn
        real(rp) :: uroe, vroe, wroe, hroe, ekroe, croe, den_roe

        real(rp), parameter :: one8 = 1.0_rp/8.0_rp
        real(rp) :: ir, u_, v_, w_, ek, p_, idy
        real(rp) :: weight, lambda_max
        real(rp) :: vnc, hgU, hgV, hgW, hgE, cc, un_c2
        integer  :: i,ii,j,k,l,m,jl, n, jj, s,iimax
        integer  :: str_y, end_y, jgbl, jstr_y, jend_y
        integer(1), parameter :: two = 2
        integer(1), parameter :: two_weno_smooth = 2*weno_smooth
        integer(1), parameter :: two_weno_shock  = 2*weno_shock

        call StartProfRange("hybrid_wenoEPy")
        
        !$acc parallel vector_length(vec_size) default(present) &
        !$acc private(tilde_op,pri_1D,phi_arr) &
        !$acc private(flx_arr,flx,weno_flag_chunk) &
        !$acc create(central_1_p)
        !$acc cache(central_1_p,tilde_op,pri_1D,phi_arr,flx_arr,flx,weno_flag_chunk)
        !$acc loop gang collapse(2)
        do k   = sz,ez
          do i = sx,ex,4
             central_1_p = 2*central_1

             iimax = min(3,ex-i)

             do str_y = lby,uby,chunkSize
                end_y = min(str_y+ChunkSize, uby-2*GN+1) - str_y

                ! correction for NSCBC
                jstr_y = max(str_y+GN       , jstr     ) - str_y - GN + 1
                jend_y = min(str_y+ChunkSize, jend-GN+1) - str_y
                if(jend_y - jstr_y > -1) then
                !$acc loop vector collapse(2)
                do j = 1-GN,end_y+GN
                   do ii = 0,iimax
                      jgbl = j + str_y-1+GN
                      phi_arr(ii,j,1) = phi(i+ii,jgbl,k,1)
                      phi_arr(ii,j,2) = phi(i+ii,jgbl,k,2)
                      phi_arr(ii,j,3) = phi(i+ii,jgbl,k,3)
                      phi_arr(ii,j,4) = phi(i+ii,jgbl,k,4)
                      phi_arr(ii,j,5) = phi(i+ii,jgbl,k,5)
                   enddo
                enddo
                do ii=0,iimax
                   !$acc loop vector
                   do j = 1-GN,end_y+GN
                      jgbl = j + str_y-1+GN

                      phi1 = phi_arr(ii,j,1)
                      phi2 = phi_arr(ii,j,2)
                      phi3 = phi_arr(ii,j,3)
                      phi4 = phi_arr(ii,j,4)
                      phi5 = phi_arr(ii,j,5)

                      ir = 1.0_rp/phi1
                      u_ = phi2*ir
                      v_ = phi3*ir
                      w_ = phi4*ir
                      ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
                      p_ = gm1 * (phi5 - phi1*ek)

                      flx_arr(j,1) =  phi3
                      flx_arr(j,2) =  phi3    * u_
                      flx_arr(j,3) =  phi3    * v_ + p_
                      flx_arr(j,4) =  phi3    * w_
                      flx_arr(j,5) = (phi5+p_)* v_

                      pri_1D(j,1) = p_
                      pri_1D(j,2) = u_
                      pri_1D(j,3) = v_
                      pri_1D(j,4) = w_
                      pri_1D(j,5) = hgm*p_*ir + ek
                      pri_1D(j,6) = phi1

                      weno_flag_chunk(j) = iand(weno_flag_xyz(i+ii,jgbl,k),two)

                   enddo

                   !$acc loop vector collapse(2)
                   do j = 1-GN,end_y+GN-3
                      do l = 1,3

                         jl = j+l

                         weight = one8 * (pri_1D(j,6) + pri_1D(jl,6)) &
                                       * (pri_1D(j,3) + pri_1D(jl,3))

                         tilde_op(l,j,1) = 2*weight
                         tilde_op(l,j,2) = weight*(pri_1D(j,2) + pri_1D(jl,2))
                         tilde_op(l,j,3) = weight*(pri_1D(j,3) + pri_1D(jl,3)) + 0.5_rp*(pri_1D(j,1) + pri_1D(jl,1))
                         tilde_op(l,j,4) = weight*(pri_1D(j,4) + pri_1D(jl,4))
                         tilde_op(l,j,5) = weight*(pri_1D(j,5) + pri_1D(jl,5))

                      enddo
                   enddo

                   !!$acc loop vector
                   !do j = 0,end_y ! sy-1,ey
                   !   ishock(j) = 1
                   !   !$acc loop seq
                   !   do jj = j-weno_num+1,j+weno_num
                   !      ishock(j) = min(weno_flag_chunk(jj),ishock(j))
                   !   enddo
                   !enddo

                   !$acc loop vector
                   do j = 0,end_y ! sy-1,ey
                      if(weno_flag_chunk(j) == two_weno_smooth) then

                        do n = 1,5
                           hGn = 0.0_rp
                           do l = 1, 3
                              inner_sum = 0.0_rp
                              do m = 0, l-1
                                 inner_sum = inner_sum + tilde_op(l,j-m,n)
                              enddo

                              hGn = hGn + central_1_p(l)*inner_sum
                           enddo
                           flx(j,n) = hGn
                        enddo
                      endif
                   enddo

                   !$acc loop vector private(right,left,fvp,fvm,fp,fm,fhatp,fhatm)
                   do j = 0,end_y ! sy-1,ey
                      if(weno_flag_chunk(j) == two_weno_shock) then

                      !
                      ! ROE MEAN STATE
                      !
                      sqrt_rho0 = sqrt(pri_1D(j,6))
                      sqrt_rho1 = sqrt(pri_1D(j+1,6))
                      den_roe   = 1.0_rp/(sqrt_rho0 + sqrt_rho1)

                      uroe = den_roe*(pri_1D(j,2)*sqrt_rho1+pri_1D(j+1,2)*sqrt_rho0)
                      vroe = den_roe*(pri_1D(j,3)*sqrt_rho1+pri_1D(j+1,3)*sqrt_rho0)
                      wroe = den_roe*(pri_1D(j,4)*sqrt_rho1+pri_1D(j+1,4)*sqrt_rho0)
                      hroe = den_roe*(pri_1D(j,5)*sqrt_rho1+pri_1D(j+1,5)*sqrt_rho0)
                      ekroe= 0.5_rp*(uroe*uroe + vroe*vroe + wroe*wroe)
                      croe = sqrt(gm1*(hroe-ekroe))
                      !
                      ! FILL EIGEN MATRIXES
                      !
                      vnc = vroe*croe
                      hgU = gm1*uroe
                      hgV = gm1*vroe
                      hgW = gm1*wroe
                      hgE = gm1*ekroe

                      right(1,1)=1.0_rp
                      right(1,2)=1.0_rp
                      right(1,3)=1.0_rp
                      right(1,4)=0.0_rp
                      right(1,5)=0.0_rp
                      !
                      right(2,1)=uroe
                      right(2,2)=uroe
                      right(2,3)=uroe
                      right(2,4)=1.0_rp
                      right(2,5)=0.0_rp
                      !
                      right(3,1)=vroe-croe
                      right(3,2)=vroe
                      right(3,3)=vroe+croe
                      right(3,4)=0.0_rp
                      right(3,5)=0.0_rp
                      !
                      right(4,1)=wroe
                      right(4,2)=wroe
                      right(4,3)=wroe
                      right(4,4)=0.0_rp
                      right(4,5)=-1.0_rp
                      !
                      right(5,1)=hroe-vnc
                      right(5,2)=ekroe
                      right(5,3)=hroe+vnc
                      right(5,4)=uroe
                      right(5,5)=-wroe

                      ! Left eigen matrix
                      cc    = croe*croe
                      un_c2 = 1.0_rp/cc

                      left(1,1)=0.5_rp*un_c2*(hgE+vnc)
                      left(1,2)=0.5_rp*un_c2*(-hgU)
                      left(1,3)=0.5_rp*un_c2*(-hgV-croe)
                      left(1,4)=0.5_rp*un_c2*(-hgW)
                      left(1,5)=0.5_rp*un_c2*gm1
                      !
                      left(2,1)=un_c2*(cc -hgE)
                      left(2,2)=un_c2*hgU
                      left(2,3)=un_c2*hgV
                      left(2,4)=un_c2*hgW
                      left(2,5)=un_c2*(-gm1)
                      !
                      left(3,1)=0.5_rp*un_c2*(hgE-vnc)
                      left(3,2)=0.5_rp*un_c2*(-hgU)
                      left(3,3)=0.5_rp*un_c2*(-hgV+croe)
                      left(3,4)=0.5_rp*un_c2*(-hgW)
                      left(3,5)=0.5_rp*un_c2*gm1
                      !
                      left(4,1)=-uroe
                      left(4,2)=1.0_rp
                      left(4,3)=0.0_rp
                      left(4,4)=0.0_rp
                      left(4,5)=0.0_rp
                      !
                      left(5,1)=wroe
                      left(5,2)=0.0_rp
                      left(5,3)=0.0_rp
                      left(5,4)=-1.0_rp
                      left(5,5)=0.0_rp
                      !
                      ! LAX FRIEDRICHS
                      !
                      lambda_max = max(abs(vroe-croe),abs(vroe+croe))

                      do n = 1,5
                         do jj = -weno_num+1,weno_num
                            ff = flx_arr(j+jj,n)
                            pp = phi_arr(ii,j+jj,n) * lambda_max

                            fvp(n,jj) = 0.5_rp * (ff + pp)
                            fvm(n,jj) = 0.5_rp * (ff - pp)
                         enddo
                      enddo
                      !
                      ! WENO - TENO RECONSTRUCTION
                      !

                      do s = -weno_num+1, weno_num
                         do m = 1,5
                            fp(m,s) = 0.0_rp
                            fm(m,s) = 0.0_rp
                            do l = 1,5
                            fp(m,s) = fp(m,s) + left(m,l)*fvp(l,s)
                            fm(m,s) = fm(m,s) + left(m,l)*fvm(l,s)
                            enddo
                         enddo
                      enddo

                      call weno_reconstruction(weno_num,shock_recon,fp,fm,fhatp,fhatm)

                      ! return in the orinal system
                      do l = 1,5
                         flx(j,l) = 0.0_rp
                         do m = 1,5
                            flx(j,l) = flx(j,l) + right(l,m) * (fhatp(m) + fhatm(m))
                         enddo
                      enddo
                      endif
                   enddo ! j
                
                   !$acc loop vector
                   do j = jstr_y,jend_y
                      jgbl = j + str_y-1+GN
                      idy = ystep_i(jgbl)
                      phi_arr(ii,j,1) = -idy * (flx(j,1) - flx(j-1,1))
                      phi_arr(ii,j,2) = -idy * (flx(j,2) - flx(j-1,2))
                      phi_arr(ii,j,3) = -idy * (flx(j,3) - flx(j-1,3))
                      phi_arr(ii,j,4) = -idy * (flx(j,4) - flx(j-1,4))
                      phi_arr(ii,j,5) = -idy * (flx(j,5) - flx(j-1,5))
                   enddo
                enddo ! ii

                !$acc loop vector collapse(2)
                ! do j = jstr,jend
                do j = jstr_y,jend_y
                   do ii = 0,iimax
                      jgbl = j + str_y-1+GN
                      RHS(i+ii,jgbl,k,1) = RHS(i+ii,jgbl,k,1) + phi_arr(ii,j,1)
                      RHS(i+ii,jgbl,k,2) = RHS(i+ii,jgbl,k,2) + phi_arr(ii,j,2)
                      RHS(i+ii,jgbl,k,3) = RHS(i+ii,jgbl,k,3) + phi_arr(ii,j,3)
                      RHS(i+ii,jgbl,k,4) = RHS(i+ii,jgbl,k,4) + phi_arr(ii,j,4)
                      RHS(i+ii,jgbl,k,5) = RHS(i+ii,jgbl,k,5) + phi_arr(ii,j,5)
                   enddo
                enddo
                endif
             enddo ! chuncks
          enddo ! i
        enddo ! k
        !$acc end parallel

        !$omp end parallel do

        call EndProfRange

        return
end subroutine hybrid_wenoEPy


subroutine hybrid_wenoEPz(weno_num,lbz,ubz,kstr,kend,shock_recon)

        implicit none
        integer, parameter :: vec_size = 64
        integer, parameter :: chunkSize = vec_size - 2*GN

        integer, intent(in) :: lbz,ubz
        integer, intent(in) :: weno_num
        integer, intent(in) :: kstr,kend
        integer, intent(in) :: shock_recon

        real(rp),   dimension(3,1-GN:chunkSize+GN,5) :: tilde_op
        real(rp),   dimension(1-GN:chunkSize+GN,6)   :: pri_1D
        real(rp),   dimension(1-GN:chunkSize+GN,5)   :: flx_arr, flx
        real(rp),   dimension(0:3,1-GN:chunkSize+GN,5)   :: phi_arr
        
        integer(1), dimension(1-GN:chunkSize+GN)     :: weno_flag_chunk

        ! local declarations
        real(rp), parameter    :: gm1 = gamma0-1.0_rp
        real(rp), parameter    :: hgm = gamma0/gm1

        real(rp), dimension(5,-3:4) :: fvp, fvm, fp, fm
        real(rp), dimension(5)      :: fhatp, fhatm
        real(rp), dimension(5,5)    :: right, left
        real(rp), dimension(-3:3)   :: central_1_p

        real(rp) :: phi1, phi2, phi3, phi4, phi5, inner_sum
        real(rp) :: sqrt_rho0, sqrt_rho1, ff, pp, hHn
        real(rp) :: uroe, vroe, wroe, hroe, ekroe, croe, den_roe

        real(rp), parameter :: one8 = 1.0_rp/8.0_rp
        real(rp) :: ir, u_, v_, w_, ek, p_, idz
        real(rp) :: weight, lambda_max
        real(rp) :: vnc, hgU, hgV, hgW, hgE, cc, un_c2
        integer  :: i,ii,j,k,l,m,kl, n, kk, s, iimax
        integer  :: str_z, end_z, kgbl, kstr_z, kend_z
        integer(1), parameter :: four = 4
        integer(1), parameter :: four_weno_smooth = four*weno_smooth
        integer(1), parameter :: four_weno_shock  = four*weno_shock

        call StartProfRange("hybrid_wenoEPz")

        !$acc parallel vector_length(vec_size) default(present) &
        !$acc private(tilde_op,pri_1D,phi_arr) &
        !$acc private(flx_arr,flx,weno_flag_chunk) &
        !$acc create(central_1_p)
        !$acc cache(central_1_p,tilde_op,pri_1D,phi_arr,flx_arr,flx,weno_flag_chunk)
        !$acc loop gang collapse(2)
        do j    = sy,ey
           do i = sx,ex,4
              central_1_p = 2*central_1

              iimax = min(3,ex-i)

              do str_z = lbz,ubz,chunkSize
                 end_z = min(str_z+ChunkSize, ubz-2*GN+1) - str_z

                 ! correction for NSCBC
                 kstr_z = max(str_z+GN       , kstr     ) - str_z - GN + 1
                 kend_z = min(str_z+ChunkSize, kend-GN+1) - str_z
                 if(kend_z - kstr_z > -1) then
                 !$acc loop vector collapse(2)
                 do k = 1-GN,end_z+GN
                    do ii=0,iimax
                       kgbl = k + str_z-1+GN
                       phi_arr(ii,k,1) = phi(i+ii,j,kgbl,1)
                       phi_arr(ii,k,2) = phi(i+ii,j,kgbl,2)
                       phi_arr(ii,k,3) = phi(i+ii,j,kgbl,3)
                       phi_arr(ii,k,4) = phi(i+ii,j,kgbl,4)
                       phi_arr(ii,k,5) = phi(i+ii,j,kgbl,5)
                    enddo
                 enddo
                 do ii=0,iimax
                    !$acc loop vector
                    do k = 1-GN,end_z+GN
                       kgbl = k + str_z-1+GN

                       phi1 = phi_arr(ii,k,1)
                       phi2 = phi_arr(ii,k,2)
                       phi3 = phi_arr(ii,k,3)
                       phi4 = phi_arr(ii,k,4)
                       phi5 = phi_arr(ii,k,5)

                       ir = 1.0_rp/phi1
                       u_ = phi2*ir
                       v_ = phi3*ir
                       w_ = phi4*ir
                       ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
                       p_ = gm1 * (phi5 - phi1*ek)

                       flx_arr(k,1) =  phi4
                       flx_arr(k,2) =  phi4    * u_
                       flx_arr(k,3) =  phi4    * v_
                       flx_arr(k,4) =  phi4    * w_ + p_
                       flx_arr(k,5) = (phi5+p_)* w_

                       pri_1D(k,1) = p_
                       pri_1D(k,2) = u_
                       pri_1D(k,3) = v_
                       pri_1D(k,4) = w_
                       pri_1D(k,5) = hgm*p_*ir + ek
                       pri_1D(k,6) = phi1

                       weno_flag_chunk(k) = iand(weno_flag_xyz(i+ii,j,kgbl),four)

                    enddo

                    !$acc loop vector collapse(2)
                    do k = 1-GN,end_z+GN-3 !lbz, ubz-3
                       do l = 1,3

                          kl = k+l

                          weight = one8 * (pri_1D(k,6) + pri_1D(kl,6)) &
                               * (pri_1D(k,4) + pri_1D(kl,4))

                          tilde_op(l,k,1) = 2*weight
                          tilde_op(l,k,2) = weight*(pri_1D(k,2) + pri_1D(kl,2))
                          tilde_op(l,k,3) = weight*(pri_1D(k,3) + pri_1D(kl,3))
                          tilde_op(l,k,4) = weight*(pri_1D(k,4) + pri_1D(kl,4)) + 0.5_rp*(pri_1D(k,1) + pri_1D(kl,1))
                          tilde_op(l,k,5) = weight*(pri_1D(k,5) + pri_1D(kl,5))

                       enddo
                    enddo

                    !!$acc loop vector
                    !do k = 0,end_z !sz-1,ez
                    !   ishock(k) = 1
                    !   do kk = k-weno_num+1,k+weno_num
                    !      ishock(k) = min(weno_flag_chunk(kk),ishock(k))
                    !   enddo
                    !enddo

                    !$acc loop vector
                    do k = 0,end_z !sz-1,ez
                       if(weno_flag_chunk(k) == four_weno_smooth) then

                          do n = 1,5
                             hHn = 0.0_rp
                             do l = 1, 3
                                inner_sum = 0.0_rp
                                do m = 0, l-1
                                   inner_sum = inner_sum + tilde_op(l,k-m,n)
                                enddo

                                hHn = hHn + central_1_p(l)*inner_sum
                             enddo
                             flx(k,n) = hHn
                          enddo
                       endif
                    enddo

                    !$acc loop vector private(right,left,fvp,fvm,fp,fm,fhatp,fhatm)
                    do k = 0,end_z !sz-1,ez
                       if (weno_flag_chunk(k) == four_weno_shock) then
                          !
                          ! ROE MEAN STATE
                          !
                          sqrt_rho0 = sqrt(pri_1D(k,6))
                          sqrt_rho1 = sqrt(pri_1D(k+1,6))
                          den_roe   = 1.0_rp/(sqrt_rho0 + sqrt_rho1)

                          uroe = den_roe*(pri_1D(k,2)*sqrt_rho1+pri_1D(k+1,2)*sqrt_rho0)
                          vroe = den_roe*(pri_1D(k,3)*sqrt_rho1+pri_1D(k+1,3)*sqrt_rho0)
                          wroe = den_roe*(pri_1D(k,4)*sqrt_rho1+pri_1D(k+1,4)*sqrt_rho0)
                          hroe = den_roe*(pri_1D(k,5)*sqrt_rho1+pri_1D(k+1,5)*sqrt_rho0)
                          ekroe= 0.5_rp*(uroe*uroe + vroe*vroe + wroe*wroe)
                          croe = sqrt(gm1*(hroe-ekroe))
                          !
                          ! FILL EIGEN MATRIXES
                          !
                          vnc = wroe*croe
                          hgU = gm1*uroe
                          hgV = gm1*vroe
                          hgW = gm1*wroe
                          hgE = gm1*ekroe

                          right(1,1)=1.0_rp
                          right(1,2)=1.0_rp
                          right(1,3)=1.0_rp
                          right(1,4)=0.0_rp
                          right(1,5)=0.0_rp
                          !
                          right(2,1)=uroe
                          right(2,2)=uroe
                          right(2,3)=uroe
                          right(2,4)=-1.0_rp
                          right(2,5)=0.0_rp
                          !
                          right(3,1)=vroe
                          right(3,2)=vroe
                          right(3,3)=vroe
                          right(3,4)=0.0_rp
                          right(3,5)=1.0_rp
                          !
                          right(4,1)=wroe-croe
                          right(4,2)=wroe
                          right(4,3)=wroe+croe
                          right(4,4)=0.0_rp
                          right(4,5)=0.0_rp
                          !
                          right(5,1)=hroe-vnc
                          right(5,2)=ekroe
                          right(5,3)=hroe+vnc
                          right(5,4)=-uroe
                          right(5,5)=vroe

                          ! Left eigen matrix
                          cc    = croe*croe
                          un_c2 = 1.0_rp/cc

                          left(1,1)=0.5_rp*un_c2*(hgE + vnc)
                          left(1,2)=0.5_rp*un_c2*(-hgU)
                          left(1,3)=0.5_rp*un_c2*(-hgV)
                          left(1,4)=0.5_rp*un_c2*(-hgW-croe)
                          left(1,5)=0.5_rp*un_c2*gm1
                          !
                          left(2,1)=un_c2*(cc  - hgE)
                          left(2,2)=un_c2*hgU
                          left(2,3)=un_c2*hgV
                          left(2,4)=un_c2*hgW
                          left(2,5)=un_c2*(-gm1)
                          !
                          left(3,1)=0.5_rp*un_c2*(hgE - vnc)
                          left(3,2)=0.5_rp*un_c2*(-hgU)
                          left(3,3)=0.5_rp*un_c2*(-hgV)
                          left(3,4)=0.5_rp*un_c2*(-hgW+croe)
                          left(3,5)=0.5_rp*un_c2*gm1
                          !
                          left(4,1)=uroe
                          left(4,2)=-1.0_rp
                          left(4,3)=0.0_rp
                          left(4,4)=0.0_rp
                          left(4,5)=0.0_rp
                          !
                          left(5,1)=-vroe
                          left(5,2)=0.0_rp
                          left(5,3)=1.0_rp
                          left(5,4)=0.0_rp
                          left(5,5)=0.0_rp
                          !
                          ! LAX FRIEDRICHS
                          !
                          lambda_max = max(abs(wroe-croe),abs(wroe+croe))

                          do n = 1,5
                             do kk = -weno_num+1,weno_num
                                ff = flx_arr(k+kk,n)
                                pp = phi_arr(ii,k+kk,n) * lambda_max

                                fvp(n,kk) = 0.5_rp * (ff + pp)
                                fvm(n,kk) = 0.5_rp * (ff - pp)
                             enddo
                          enddo
                          !
                          ! WENO - TENO RECONSTRUCTION
                          !

                          do s = -weno_num+1, weno_num
                             do m = 1,5
                                fp(m,s) = 0.0_rp
                                fm(m,s) = 0.0_rp
                                do l = 1,5
                                   fp(m,s) = fp(m,s) + left(m,l)*fvp(l,s)
                                   fm(m,s) = fm(m,s) + left(m,l)*fvm(l,s)
                                enddo
                             enddo
                          enddo

                          call weno_reconstruction(weno_num,shock_recon,fp,fm,fhatp,fhatm)

                          ! return in the orinal system
                          do l = 1,5
                             flx(k,l) = 0.0_rp
                             do m = 1,5
                                flx(k,l) = flx(k,l) + right(l,m) * (fhatp(m) + fhatm(m))
                             enddo
                          enddo
                       endif
                    enddo ! k

                    !$acc loop vector 
                    do k = kstr_z,kend_z
                       kgbl = k + str_z-1+GN                       
                       idz = zstep_i(kgbl)
                       phi_arr(ii,k,1) =  - idz * (flx(k,1) - flx(k-1,1))
                       phi_arr(ii,k,2) =  - idz * (flx(k,2) - flx(k-1,2))
                       phi_arr(ii,k,3) =  - idz * (flx(k,3) - flx(k-1,3))
                       phi_arr(ii,k,4) =  - idz * (flx(k,4) - flx(k-1,4))
                       phi_arr(ii,k,5) =  - idz * (flx(k,5) - flx(k-1,5))
                    enddo
                 enddo !! ii loop

                 !$acc loop vector collapse(2)
                 !do k = kstr,kend
                 do k = kstr_z,kend_z
                    do ii=0,iimax
                       kgbl = k + str_z-1+GN                       
                       RHS(i+ii,j,kgbl,1) = RHS(i+ii,j,kgbl,1) + phi_arr(ii,k,1)
                       RHS(i+ii,j,kgbl,2) = RHS(i+ii,j,kgbl,2) + phi_arr(ii,k,2)
                       RHS(i+ii,j,kgbl,3) = RHS(i+ii,j,kgbl,3) + phi_arr(ii,k,3)
                       RHS(i+ii,j,kgbl,4) = RHS(i+ii,j,kgbl,4) + phi_arr(ii,k,4)
                       RHS(i+ii,j,kgbl,5) = RHS(i+ii,j,kgbl,5) + phi_arr(ii,k,5)
                    enddo
                 enddo
              endif
              enddo ! chuncks
           enddo ! i
        enddo ! j
        !$acc end parallel

        !$omp end parallel do

        call EndProfRange

        return
end subroutine hybrid_wenoEPz





subroutine weno_reconstruction(weno_num,shock_recon,fp,fm,fhatp,fhatm)
! ---------------------------------------------------------------------------
!
!       This subroutine performs the WENO reconstruction in a cell bound.
!       The cell bound is identified by the b parameter.
!
!       b = 0  -> left  cell bound
!       b = -1 -> right cell bound
!
! ---------------------------------------------------------------------------
        !$acc routine seq
        implicit none

        integer                                       , intent(in ) :: shock_recon, weno_num
        real(rp), dimension(5,-3:4)                   , intent(in ) :: fp, fm
        real(rp), dimension(5)                        , intent(out) :: fhatp, fhatm

        ! local declarations
        real(rp), parameter          :: varepsilon = 1.0E-12_rp  !< pay attention on this in parallel
        real(rp), parameter          :: is1par = 13.0_rp/12.0_rp, is2par = 0.25_rp

        real(rp) :: q0, q1, q2, q3, Is0, Is1, Is2, Is3
        real(rp) :: alpha0, alpha1, alpha2, alpha3
        real(rp) :: omega0, omega1, omega2, omega3
        real(rp) :: isumAlpha, absTerm0, absTerm1

        integer :: s

        real(rp), parameter :: a1 = 10.5_rp, a2 = 3.5_rp, cr = 0.25_rp, csi = 1.0E-3_rp
        real(rp) :: eta0, eta1, eta2, eps, rm, dsum
        real(rp) :: ct
        integer  :: d0, d1, d2, d3


        ct = 1.0E-7_rp

        if(shock_recon == 1) then

        if    (weno_num == 2) then
        do s = 1, 5
           ! --- WENO plus part --- !

           ! polinomia
           q0 = aweno3(0,0)*fp(s,-1) + aweno3(0,1)*fp(s,0)
           q1 = aweno3(1,0)*fp(s, 0) + aweno3(1,1)*fp(s,1)

           ! smoothness index
           Is0 = (fp(s,-1)-fp(s,0))**2
           Is1 = (fp(s, 0)-fp(s,1))**2

           ! alpha
           !do r = 0, 1
           !   alpha(r) = cweno3(r)/(IS(r) + varepsilon)**2   !<WENO JS
           !enddo

           absTerm0 = abs(Is0 - Is1)
           alpha0 = cweno3(0)*(1.0_rp + absTerm0/(Is0 + varepsilon))  !<WENO Z
           alpha1 = cweno3(1)*(1.0_rp + absTerm0/(Is1 + varepsilon))  !<WENO Z

           ! omega
           isumAlpha = 1.0_rp/(alpha0 + alpha1)
           omega0 = alpha0*isumAlpha
           omega1 = alpha1*isumAlpha

           ! WENO plus reconstruction
           fhatp(s) = omega0*q0 + omega1*q1

           ! --- WENO minus part --- !

           ! polinomia
           q0 = aweno3(0,0)*fm(s,2) + aweno3(0,1)*fm(s,1)
           q1 = aweno3(1,0)*fm(s,1) + aweno3(1,1)*fm(s,0)

           ! smoothness indes
           Is0 = (fm(s,2)-fm(s,1))**2
           Is1 = (fm(s,1)-fm(s,0))**2

           ! alpha
           !do r = 0, 1
           !   alpha(r) = cweno3(r)/(IS(r) + varepsilon)**2   !<WENO JS
           !enddo

           absTerm0 = abs(Is0 - Is1)
           alpha0 = cweno3(0)*(1.0_rp + absTerm0/(Is0 + varepsilon))  !<WENO Z
           alpha1 = cweno3(1)*(1.0_rp + absTerm0/(Is1 + varepsilon))  !<WENO Z

           ! omega
           isumAlpha = 1.0_rp/(alpha0 + alpha1)
           omega0 = alpha0*isumAlpha
           omega1 = alpha1*isumAlpha

           ! WENO minus riconstruction
           fhatm(s) = omega0*q0 + omega1*q1

        enddo


        elseif(weno_num == 3) then
        do s = 1, 5
           ! --- WENO plus part --- !

           ! polinomia
           q0 = aweno5(0,0)*fp(s,-2) + aweno5(0,1)*fp(s,-1) + aweno5(0,2)*fp(s,0)
           q1 = aweno5(1,0)*fp(s,-1) + aweno5(1,1)*fp(s, 0) + aweno5(1,2)*fp(s,1)
           q2 = aweno5(2,0)*fp(s, 0) + aweno5(2,1)*fp(s, 1) + aweno5(2,2)*fp(s,2)

           ! smoothness index
           IS0 = is1par*(fp(s,-2)-2.0_rp*fp(s,-1)+fp(s,0))**2 + is2par*(    fp(s,-2)-4.0_rp*fp(s,-1)+3.0_rp*fp(s,0))**2
           IS1 = is1par*(fp(s,-1)-2.0_rp*fp(s, 0)+fp(s,1))**2 + is2par*(    fp(s,-1)                       -fp(s,1))**2
           IS2 = is1par*(fp(s, 0)-2.0_rp*fp(s, 1)+fp(s,2))**2 + is2par*(3.0_rp*fp(s, 0)-4.0_rp*fp(s, 1)    +fp(s,2))**2

           ! alpha
           !alphap = cweno5/(varepsilon+ISp)**2   !< WENO JS
           absTerm0 = abs(IS0 - IS2)
           alpha0 = absTerm0/(IS0 + varepsilon)
           alpha0 = cweno5(0)*(1.0_rp + alpha0*alpha0)
           alpha1 = absTerm0/(IS1 + varepsilon)
           alpha1 = cweno5(1)*(1.0_rp + alpha1*alpha1)
           alpha2 = absTerm0/(IS2 + varepsilon)
           alpha2 = cweno5(2)*(1.0_rp + alpha2*alpha2)

           ! omega
           isumAlpha = 1.0_rp/(alpha0 + alpha1 + alpha2)
           omega0 = alpha0*isumAlpha
           omega1 = alpha1*isumAlpha
           omega2 = alpha2*isumAlpha

           ! WENO plus reconstruction
           fhatp(s) = omega0*q0 + omega1*q1 + omega2*q2

           ! --- WENO minus part --- !

           ! polinomia
           q0 = aweno5(0,0)*fm(s,3) + aweno5(0,1)*fm(s,2) + aweno5(0,2)*fm(s, 1)
           q1 = aweno5(1,0)*fm(s,2) + aweno5(1,1)*fm(s,1) + aweno5(1,2)*fm(s, 0)
           q2 = aweno5(2,0)*fm(s,1) + aweno5(2,1)*fm(s,0) + aweno5(2,2)*fm(s,-1)

           ! smoothness indes
           IS0 = is1par*(fm(s,3)-2.0_rp*fm(s,2)+fm(s, 1))**2 + is2par*(    fm(s,3)-4.0_rp*fm(s,2)+3.0_rp*fm(s, 1))**2
           IS1 = is1par*(fm(s,2)-2.0_rp*fm(s,1)+fm(s, 0))**2 + is2par*(    fm(s,2)                      -fm(s, 0))**2
           IS2 = is1par*(fm(s,1)-2.0_rp*fm(s,0)+fm(s,-1))**2 + is2par*(3.0_rp*fm(s,1)-4.0_rp*fm(s,0)    +fm(s,-1))**2

           ! alpha
           !alpham = cweno5/(varepsilon+ISm)**2  !< WENO JS
           absTerm0 = abs(IS0 - IS2)
           alpha0 = absTerm0/(IS0 + varepsilon)
           alpha0 = cweno5(0)*(1.0_rp + alpha0*alpha0)
           alpha1 = absTerm0/(IS1 + varepsilon)
           alpha1 = cweno5(1)*(1.0_rp + alpha1*alpha1)
           alpha2 = absTerm0/(IS2 + varepsilon)
           alpha2 = cweno5(2)*(1.0_rp + alpha2*alpha2)

           ! omega
           isumAlpha = 1.0_rp/(alpha0+alpha1+alpha2)
           omega0 = alpha0*isumAlpha
           omega1 = alpha1*isumAlpha
           omega2 = alpha2*isumAlpha

           ! WENO minus riconstruction
           fhatm(s) = omega0*q0 + omega1*q1 + omega2*q2

        enddo
        elseif(weno_num == 4) then

        do s = 1, 5
           ! --- WENO plus part --- !

           ! polinomia
           q0 = aweno7(0,0)*fp(s,-3) + aweno7(0,1)*fp(s,-2) + aweno7(0,2)*fp(s,-1) + aweno7(0,3)*fp(s,0)
           q1 = aweno7(1,0)*fp(s,-2) + aweno7(1,1)*fp(s,-1) + aweno7(1,2)*fp(s, 0) + aweno7(1,3)*fp(s,1)
           q2 = aweno7(2,0)*fp(s,-1) + aweno7(2,1)*fp(s, 0) + aweno7(2,2)*fp(s, 1) + aweno7(2,3)*fp(s,2)
           q3 = aweno7(3,0)*fp(s, 0) + aweno7(3,1)*fp(s, 1) + aweno7(3,2)*fp(s, 2) + aweno7(3,3)*fp(s,3)

           ! smoothness index
           Is0 = fp(s,-3)*( 547.0_rp*fp(s,-3)- 3882.0_rp*fp(s,-2)+ 4642.0_rp*fp(s,-1)-1854.0_rp*fp(s,0)) + &
                 fp(s,-2)*(                    7043.0_rp*fp(s,-2)-17246.0_rp*fp(s,-1)+7042.0_rp*fp(s,0)) + &
                 fp(s,-1)*(                                       11003.0_rp*fp(s,-1)-9402.0_rp*fp(s,0)) + &
                 fp(s, 0)*(                                                           2107.0_rp*fp(s,0))
           Is1 = fp(s,-2)*( 267.0_rp*fp(s,-2)- 1642.0_rp*fp(s,-1)+ 1602.0_rp*fp(s, 0)- 494.0_rp*fp(s,1)) + &
                 fp(s,-1)*(                    2843.0_rp*fp(s,-1)- 5966.0_rp*fp(s, 0)+1922.0_rp*fp(s,1)) + &
                 fp(s, 0)*(                                        3443.0_rp*fp(s, 0)-2522.0_rp*fp(s,1)) + &
                 fp(s, 1)*(                                                            547.0_rp*fp(s,1))
           Is2 = fp(s,-1)*( 547.0_rp*fp(s,-1)- 2522.0_rp*fp(s, 0)+ 1922.0_rp*fp(s, 1)- 494.0_rp*fp(s,2)) + &
                 fp(s, 0)*(                    3443.0_rp*fp(s, 0)- 5966.0_rp*fp(s, 1)+1602.0_rp*fp(s,2)) + &
                 fp(s, 1)*(                                        2843.0_rp*fp(s, 1)-1642.0_rp*fp(s,2)) + &
                 fp(s, 2)*(                                                            267.0_rp*fp(s,2))
           Is3 = fp(s, 0)*(2107.0_rp*fp(s, 0)- 9402.0_rp*fp(s, 1)+ 7042.0_rp*fp(s, 2)-1854.0_rp*fp(s,3)) + &
                 fp(s, 1)*(                   11003.0_rp*fp(s, 1)-17246.0_rp*fp(s, 2)+4642.0_rp*fp(s,3)) + &
                 fp(s, 2)*(                                        7043.0_rp*fp(s, 2)-3882.0_rp*fp(s,3)) + &
                 fp(s, 3)*(                                                            547.0_rp*fp(s,3))

           ! alpha
           !do r = 0, 3
           !   alpha(r) = cweno7(r)/(IS(r) + varepsilon)**2   !<WENO JS
           !enddo

           absTerm0 = abs(Is0             - Is3)
           absTerm1 = abs(Is0 - Is1 - Is2 + Is3)

           !WENO Z
           alpha0 = absTerm0/(IS0 + varepsilon)
           alpha0 = cweno7(0)*(1.0_rp + alpha0*alpha0)
           alpha1 = absTerm1/(IS1 + varepsilon)
           alpha1 = cweno7(1)*(1.0_rp + alpha1*alpha1)
           alpha2 = absTerm0/(IS2 + varepsilon)
           alpha2 = cweno7(2)*(1.0_rp + alpha2*alpha2)
           alpha3 = absTerm1/(IS3 + varepsilon)
           alpha3 = cweno7(3)*(1.0_rp + alpha3*alpha3)

           ! omega
           isumAlpha = 1.0_rp/(alpha0 + alpha1 + alpha2 + alpha3)
           omega0 = alpha0*isumAlpha
           omega1 = alpha1*isumAlpha
           omega2 = alpha2*isumAlpha
           omega3 = alpha3*isumAlpha

           ! WENO plus reconstruction
           fhatp(s) = omega0*q0 + omega1*q1 + omega2*q2 + omega3*q3

           ! --- WENO minus part --- !

           ! polinomia
           q0 = aweno7(0,0)*fm(s,4) + aweno7(0,1)*fm(s,3) + aweno7(0,2)*fm(s, 2) + aweno7(0,3)*fm(s, 1)
           q1 = aweno7(1,0)*fm(s,3) + aweno7(1,1)*fm(s,2) + aweno7(1,2)*fm(s, 1) + aweno7(1,3)*fm(s, 0)
           q2 = aweno7(2,0)*fm(s,2) + aweno7(2,1)*fm(s,1) + aweno7(2,2)*fm(s, 0) + aweno7(2,3)*fm(s,-1)
           q3 = aweno7(3,0)*fm(s,1) + aweno7(3,1)*fm(s,0) + aweno7(3,2)*fm(s,-1) + aweno7(3,3)*fm(s,-2)

           ! smoothness indes
           Is0 = fm(s, 4)*( 547.0_rp*fm(s,4)- 3882.0_rp*fm(s,3)+ 4642.0_rp*fm(s, 2)-1854.0_rp*fm(s, 1)) + &
                 fm(s, 3)*(                   7043.0_rp*fm(s,3)-17246.0_rp*fm(s, 2)+7042.0_rp*fm(s, 1)) + &
                 fm(s, 2)*(                                     11003.0_rp*fm(s, 2)-9402.0_rp*fm(s, 1)) + &
                 fm(s, 1)*(                                                         2107.0_rp*fm(s, 1))
           Is1 = fm(s, 3)*( 267.0_rp*fm(s,3)- 1642.0_rp*fm(s,2)+ 1602.0_rp*fm(s, 1)- 494.0_rp*fm(s, 0)) + &
                 fm(s, 2)*(                   2843.0_rp*fm(s,2)- 5966.0_rp*fm(s, 1)+1922.0_rp*fm(s, 0)) + &
                 fm(s, 1)*(                                      3443.0_rp*fm(s, 1)-2522.0_rp*fm(s, 0)) + &
                 fm(s, 0)*(                                                          547.0_rp*fm(s, 0))
           Is2 = fm(s, 2)*( 547.0_rp*fm(s,2)- 2522.0_rp*fm(s,1)+ 1922.0_rp*fm(s, 0)- 494.0_rp*fm(s,-1)) + &
                 fm(s, 1)*(                   3443.0_rp*fm(s,1)- 5966.0_rp*fm(s, 0)+1602.0_rp*fm(s,-1)) + &
                 fm(s, 0)*(                                      2843.0_rp*fm(s, 0)-1642.0_rp*fm(s,-1)) + &
                 fm(s,-1)*(                                                          267.0_rp*fm(s,-1))
           Is3 = fm(s, 1)*(2107.0_rp*fm(s,1)- 9402.0_rp*fm(s,0)+ 7042.0_rp*fm(s,-1)-1854.0_rp*fm(s,-2)) + &
                 fm(s, 0)*(                  11003.0_rp*fm(s,0)-17246.0_rp*fm(s,-1)+4642.0_rp*fm(s,-2)) + &
                 fm(s,-1)*(                                      7043.0_rp*fm(s,-1)-3882.0_rp*fm(s,-2)) + &
                 fm(s,-2)*(                                                          547.0_rp*fm(s,-2))

           ! alpha
           !do r = 0, 3
           !   alpha(r) = cweno7(r)/(IS(r) + varepsilon)**2   !<WENO JS
           !enddo

           absTerm0 = abs(Is0             - Is3)
           absTerm1 = abs(Is0 - Is1 - Is2 + Is3)

           !WENO Z
           alpha0 = absTerm0/(IS0 + varepsilon)
           alpha0 = cweno7(0)*(1.0_rp + alpha0*alpha0)
           alpha1 = absTerm1/(IS1 + varepsilon)
           alpha1 = cweno7(1)*(1.0_rp + alpha1*alpha1)
           alpha2 = absTerm0/(IS2 + varepsilon)
           alpha2 = cweno7(2)*(1.0_rp + alpha2*alpha2)
           alpha3 = absTerm1/(IS3 + varepsilon)
           alpha3 = cweno7(3)*(1.0_rp + alpha3*alpha3)

           ! omega
           isumAlpha = 1.0_rp/(alpha0 + alpha1 + alpha2 + alpha3)
           omega0 = alpha0*isumAlpha
           omega1 = alpha1*isumAlpha
           omega2 = alpha2*isumAlpha
           omega3 = alpha3*isumAlpha

           ! WENO minus reconstruction
           fhatm(s) = omega0*q0 + omega1*q1 + omega2*q2 + omega3*q3

        enddo


        endif

        elseif(shock_recon == 2) then !TENO

        if(weno_num == 3) then
        do s = 1, 5
           ! --- WENO plus part --- !

           ! polinomia
           q0 = aweno5(0,0)*fp(s,-2) + aweno5(0,1)*fp(s,-1) + aweno5(0,2)*fp(s,0)
           q1 = aweno5(1,0)*fp(s,-1) + aweno5(1,1)*fp(s, 0) + aweno5(1,2)*fp(s,1)
           q2 = aweno5(2,0)*fp(s, 0) + aweno5(2,1)*fp(s, 1) + aweno5(2,2)*fp(s,2)

           ! smoothness index
           IS0 = is1par*(fp(s,-2)-2.0_rp*fp(s,-1)+fp(s,0))**2 + is2par*(    fp(s,-2)-4.0_rp*fp(s,-1)+3.0_rp*fp(s,0))**2
           IS1 = is1par*(fp(s,-1)-2.0_rp*fp(s, 0)+fp(s,1))**2 + is2par*(    fp(s,-1)                       -fp(s,1))**2
           IS2 = is1par*(fp(s, 0)-2.0_rp*fp(s, 1)+fp(s,2))**2 + is2par*(3.0_rp*fp(s, 0)-4.0_rp*fp(s, 1)    +fp(s,2))**2

           ! alpha
           absTerm0 = abs(IS0 - IS2)
           alpha0 = 1.0_rp + (absTerm0/(IS0 + varepsilon))**6
           alpha1 = 1.0_rp + (absTerm0/(IS1 + varepsilon))**6
           alpha2 = 1.0_rp + (absTerm0/(IS2 + varepsilon))**6

           ! omega
           isumAlpha = 1.0_rp/(alpha0 + alpha1 + alpha2)
           omega0 = alpha0*isumAlpha
           omega1 = alpha1*isumAlpha
           omega2 = alpha2*isumAlpha

           ! delta
           if(omega0 < ct) then
              d0 = 0
           else
              d0 = 1
           endif

           if(omega1 < ct) then
              d1 = 0
           else
              d1 = 1
           endif

           if(omega2 < ct) then
              d2 = 0
           else
              d2 = 1
           endif

           dsum = 1.0_rp/(d0*cweno5(0) + d1*cweno5(1) + d2*cweno5(2))

           omega0 = d0*cweno5(0)*dsum
           omega1 = d1*cweno5(1)*dsum
           omega2 = d2*cweno5(2)*dsum

           ! WENO plus reconstruction
           fhatp(s) = omega0*q0 + omega1*q1 + omega2*q2

           ! --- WENO minus part --- !

           ! polinomia
           q0 = aweno5(0,0)*fm(s,3) + aweno5(0,1)*fm(s,2) + aweno5(0,2)*fm(s, 1)
           q1 = aweno5(1,0)*fm(s,2) + aweno5(1,1)*fm(s,1) + aweno5(1,2)*fm(s, 0)
           q2 = aweno5(2,0)*fm(s,1) + aweno5(2,1)*fm(s,0) + aweno5(2,2)*fm(s,-1)

           ! smoothness indes
           IS0 = is1par*(fm(s,3)-2.0_rp*fm(s,2)+fm(s, 1))**2 + is2par*(    fm(s,3)-4.0_rp*fm(s,2)+3.0_rp*fm(s, 1))**2
           IS1 = is1par*(fm(s,2)-2.0_rp*fm(s,1)+fm(s, 0))**2 + is2par*(    fm(s,2)                      -fm(s, 0))**2
           IS2 = is1par*(fm(s,1)-2.0_rp*fm(s,0)+fm(s,-1))**2 + is2par*(3.0_rp*fm(s,1)-4.0_rp*fm(s,0)    +fm(s,-1))**2

           ! alpha
           absTerm0 = abs(IS0 - IS2)
           alpha0 = 1.0_rp + (absTerm0/(IS0 + varepsilon))**6
           alpha1 = 1.0_rp + (absTerm0/(IS1 + varepsilon))**6
           alpha2 = 1.0_rp + (absTerm0/(IS2 + varepsilon))**6

           ! omega
           isumAlpha = 1.0_rp/(alpha0 + alpha1 + alpha2)
           omega0 = alpha0*isumAlpha
           omega1 = alpha1*isumAlpha
           omega2 = alpha2*isumAlpha

           ! delta
           if(omega0 < ct) then
              d0 = 0
           else
              d0 = 1
           endif

           if(omega1 < ct) then
              d1 = 0
           else
              d1 = 1
           endif

           if(omega2 < ct) then
              d2 = 0
           else
              d2 = 1
           endif

           dsum = 1.0_rp/(d0*cweno5(0) + d1*cweno5(1) + d2*cweno5(2))

           omega0 = d0*cweno5(0)*dsum
           omega1 = d1*cweno5(1)*dsum
           omega2 = d2*cweno5(2)*dsum

           ! WENO minus riconstruction
           fhatm(s) = omega0*q0 + omega1*q1 + omega2*q2

        enddo
        elseif(weno_num == 4) then
        do s = 1, 5
        ! --- WENO plus part --- !

        ! polinomia
        q0 = aweno7(0,0)*fp(s,-3) + aweno7(0,1)*fp(s,-2) + aweno7(0,2)*fp(s,-1) + aweno7(0,3)*fp(s,0)
        q1 = aweno7(1,0)*fp(s,-2) + aweno7(1,1)*fp(s,-1) + aweno7(1,2)*fp(s, 0) + aweno7(1,3)*fp(s,1)
        q2 = aweno7(2,0)*fp(s,-1) + aweno7(2,1)*fp(s, 0) + aweno7(2,2)*fp(s, 1) + aweno7(2,3)*fp(s,2)
        q3 = aweno7(3,0)*fp(s, 0) + aweno7(3,1)*fp(s, 1) + aweno7(3,2)*fp(s, 2) + aweno7(3,3)*fp(s,3)

        ! smoothness index
        Is0 = fp(s,-3)*( 547.0_rp*fp(s,-3)- 3882.0_rp*fp(s,-2)+ 4642.0_rp*fp(s,-1)-1854.0_rp*fp(s,0)) + &
              fp(s,-2)*(                    7043.0_rp*fp(s,-2)-17246.0_rp*fp(s,-1)+7042.0_rp*fp(s,0)) + &
              fp(s,-1)*(                                       11003.0_rp*fp(s,-1)-9402.0_rp*fp(s,0)) + &
              fp(s, 0)*(                                                           2107.0_rp*fp(s,0))
        Is1 = fp(s,-2)*( 267.0_rp*fp(s,-2)- 1642.0_rp*fp(s,-1)+ 1602.0_rp*fp(s, 0)- 494.0_rp*fp(s,1)) + &
              fp(s,-1)*(                    2843.0_rp*fp(s,-1)- 5966.0_rp*fp(s, 0)+1922.0_rp*fp(s,1)) + &
              fp(s, 0)*(                                        3443.0_rp*fp(s, 0)-2522.0_rp*fp(s,1)) + &
              fp(s, 1)*(                                                            547.0_rp*fp(s,1))
        Is2 = fp(s,-1)*( 547.0_rp*fp(s,-1)- 2522.0_rp*fp(s, 0)+ 1922.0_rp*fp(s, 1)- 494.0_rp*fp(s,2)) + &
              fp(s, 0)*(                    3443.0_rp*fp(s, 0)- 5966.0_rp*fp(s, 1)+1602.0_rp*fp(s,2)) + &
              fp(s, 1)*(                                        2843.0_rp*fp(s, 1)-1642.0_rp*fp(s,2)) + &
              fp(s, 2)*(                                                            267.0_rp*fp(s,2))
        Is3 = fp(s, 0)*(2107.0_rp*fp(s, 0)- 9402.0_rp*fp(s, 1)+ 7042.0_rp*fp(s, 2)-1854.0_rp*fp(s,3)) + &
              fp(s, 1)*(                   11003.0_rp*fp(s, 1)-17246.0_rp*fp(s, 2)+4642.0_rp*fp(s,3)) + &
              fp(s, 2)*(                                        7043.0_rp*fp(s, 2)-3882.0_rp*fp(s,3)) + &
              fp(s, 3)*(                                                            547.0_rp*fp(s,3))

        absTerm0 = abs(Is0             - Is3)
        absTerm1 = abs(Is0 - Is1 - Is2 + Is3)

        ! alpha
        alpha0 = 1.0_rp + (absTerm0/(Is0 + varepsilon))**6
        alpha1 = 1.0_rp + (absTerm1/(Is1 + varepsilon))**6
        alpha2 = 1.0_rp + (absTerm0/(Is2 + varepsilon))**6
        alpha3 = 1.0_rp + (absTerm1/(Is3 + varepsilon))**6

        ! omega
        isumAlpha = 1.0_rp/(alpha0 + alpha1 + alpha2 + alpha3)
        omega0 = alpha0*isumAlpha
        omega1 = alpha1*isumAlpha
        omega2 = alpha2*isumAlpha
        omega3 = alpha3*isumAlpha

        ! delta
        if(omega0 < ct) then
        d0 = 0
        else
        d0 = 1
        endif

        if(omega1 < ct) then
        d1 = 0
        else
        d1 = 1
        endif

        if(omega2 < ct) then
        d2 = 0
        else
        d2 = 1
        endif

        if(omega3 < ct) then
        d3 = 0
        else
        d3 = 1
        endif

        dsum = 1.0_rp/(d0*cweno7(0) + d1*cweno7(1) + d2*cweno7(2) + d3*cweno7(3))

        omega0 = d0*cweno7(0)*dsum
        omega1 = d1*cweno7(1)*dsum
        omega2 = d2*cweno7(2)*dsum
        omega3 = d3*cweno7(3)*dsum

        ! WENO plus reconstruction
        fhatp(s) = omega0*q0 + omega1*q1 + omega2*q2 + omega3*q3

        ! --- WENO minus part --- !

        ! polinomia
        q0 = aweno7(0,0)*fm(s,4) + aweno7(0,1)*fm(s,3) + aweno7(0,2)*fm(s, 2) + aweno7(0,3)*fm(s, 1)
        q1 = aweno7(1,0)*fm(s,3) + aweno7(1,1)*fm(s,2) + aweno7(1,2)*fm(s, 1) + aweno7(1,3)*fm(s, 0)
        q2 = aweno7(2,0)*fm(s,2) + aweno7(2,1)*fm(s,1) + aweno7(2,2)*fm(s, 0) + aweno7(2,3)*fm(s,-1)
        q3 = aweno7(3,0)*fm(s,1) + aweno7(3,1)*fm(s,0) + aweno7(3,2)*fm(s,-1) + aweno7(3,3)*fm(s,-2)

        ! smoothness indes
        Is0 = fm(s, 4)*( 547.0_rp*fm(s,4)- 3882.0_rp*fm(s,3)+ 4642.0_rp*fm(s, 2)-1854.0_rp*fm(s, 1)) + &
              fm(s, 3)*(                   7043.0_rp*fm(s,3)-17246.0_rp*fm(s, 2)+7042.0_rp*fm(s, 1)) + &
              fm(s, 2)*(                                     11003.0_rp*fm(s, 2)-9402.0_rp*fm(s, 1)) + &
              fm(s, 1)*(                                                         2107.0_rp*fm(s, 1))
        Is1 = fm(s, 3)*( 267.0_rp*fm(s,3)- 1642.0_rp*fm(s,2)+ 1602.0_rp*fm(s, 1)- 494.0_rp*fm(s, 0)) + &
              fm(s, 2)*(                   2843.0_rp*fm(s,2)- 5966.0_rp*fm(s, 1)+1922.0_rp*fm(s, 0)) + &
              fm(s, 1)*(                                      3443.0_rp*fm(s, 1)-2522.0_rp*fm(s, 0)) + &
              fm(s, 0)*(                                                          547.0_rp*fm(s, 0))
        Is2 = fm(s, 2)*( 547.0_rp*fm(s,2)- 2522.0_rp*fm(s,1)+ 1922.0_rp*fm(s, 0)- 494.0_rp*fm(s,-1)) + &
              fm(s, 1)*(                   3443.0_rp*fm(s,1)- 5966.0_rp*fm(s, 0)+1602.0_rp*fm(s,-1)) + &
              fm(s, 0)*(                                      2843.0_rp*fm(s, 0)-1642.0_rp*fm(s,-1)) + &
              fm(s,-1)*(                                                          267.0_rp*fm(s,-1))
        Is3 = fm(s, 1)*(2107.0_rp*fm(s,1)- 9402.0_rp*fm(s,0)+ 7042.0_rp*fm(s,-1)-1854.0_rp*fm(s,-2)) + &
              fm(s, 0)*(                  11003.0_rp*fm(s,0)-17246.0_rp*fm(s,-1)+4642.0_rp*fm(s,-2)) + &
              fm(s,-1)*(                                      7043.0_rp*fm(s,-1)-3882.0_rp*fm(s,-2)) + &
              fm(s,-2)*(                                                          547.0_rp*fm(s,-2))

        absTerm0 = abs(Is0             - Is3)
        absTerm1 = abs(Is0 - Is1 - Is2 + Is3)

        ! alpha
        alpha0 = 1.0_rp + (absTerm0/(Is0 + varepsilon))**6
        alpha1 = 1.0_rp + (absTerm1/(Is1 + varepsilon))**6
        alpha2 = 1.0_rp + (absTerm0/(Is2 + varepsilon))**6
        alpha3 = 1.0_rp + (absTerm1/(Is3 + varepsilon))**6

        ! omega
        isumAlpha = 1.0_rp/(alpha0 + alpha1 + alpha2 + alpha3)
        omega0 = alpha0*isumAlpha
        omega1 = alpha1*isumAlpha
        omega2 = alpha2*isumAlpha
        omega3 = alpha3*isumAlpha

        ! delta
        if(omega0 < ct) then
        d0 = 0
        else
        d0 = 1
        endif

        if(omega1 < ct) then
        d1 = 0
        else
        d1 = 1
        endif

        if(omega2 < ct) then
        d2 = 0
        else
        d2 = 1
        endif

        if(omega3 < ct) then
        d3 = 0
        else
        d3 = 1
        endif

        dsum = 1.0_rp/(d0*cweno7(0) + d1*cweno7(1) + d2*cweno7(2) + d3*cweno7(3))

        omega0 = d0*cweno7(0)*dsum
        omega1 = d1*cweno7(1)*dsum
        omega2 = d2*cweno7(2)*dsum
        omega3 = d3*cweno7(3)*dsum

        ! WENO minus reconstruction
        fhatm(s) = omega0*q0 + omega1*q1 + omega2*q2 + omega3*q3

        enddo

        endif

        elseif(shock_recon == 3) then ! TENO A

        do s = 1, 5
           ! --- WENO plus part --- !

           ! polinomia
           q0 = aweno5(0,0)*fp(s,-2) + aweno5(0,1)*fp(s,-1) + aweno5(0,2)*fp(s,0)
           q1 = aweno5(1,0)*fp(s,-1) + aweno5(1,1)*fp(s, 0) + aweno5(1,2)*fp(s,1)
           q2 = aweno5(2,0)*fp(s, 0) + aweno5(2,1)*fp(s, 1) + aweno5(2,2)*fp(s,2)

           ! smoothness index
           IS0 = is1par*(fp(s,-2)-2.0_rp*fp(s,-1)+fp(s,0))**2 + is2par*(    fp(s,-2)-4.0_rp*fp(s,-1)+3.0_rp*fp(s,0))**2
           IS1 = is1par*(fp(s,-1)-2.0_rp*fp(s, 0)+fp(s,1))**2 + is2par*(    fp(s,-1)                       -fp(s,1))**2
           IS2 = is1par*(fp(s, 0)-2.0_rp*fp(s, 1)+fp(s,2))**2 + is2par*(3.0_rp*fp(s, 0)-4.0_rp*fp(s, 1)    +fp(s,2))**2

           ! alpha
           absTerm0 = abs(IS0 - IS2)
           alpha0 = 1.0_rp + (absTerm0/(IS0 + varepsilon))**6
           alpha1 = 1.0_rp + (absTerm0/(IS1 + varepsilon))**6
           alpha2 = 1.0_rp + (absTerm0/(IS2 + varepsilon))**6

           ! omega
           isumAlpha = 1.0_rp/(alpha0 + alpha1 + alpha2)
           omega0 = alpha0*isumAlpha
           omega1 = alpha1*isumAlpha
           omega2 = alpha2*isumAlpha

           ! adaptive ct
           eps = 0.9_rp*cr/(1.0_rp-0.9_rp*cr)*csi**2
           eta0 = (2.0_rp*abs((fp(s, 0)-fp(s, -1))*(fp(s, -1)-fp(s, -2))) + eps)/ &
                         ((fp(s, 0)-fp(s, -1))**2+(fp(s, -1)-fp(s, -2))**2+ eps)
           eta1 = (2.0_rp*abs((fp(s, 1)-fp(s,  0))*(fp(s,  0)-fp(s, -1))) + eps)/ &
                         ((fp(s, 1)-fp(s,  0))**2+(fp(s,  0)-fp(s, -1))**2+ eps)
           eta2 = (2.0_rp*abs((fp(s, 2)-fp(s,  1))*(fp(s,  1)-fp(s,  0))) + eps)/ &
                         ((fp(s, 2)-fp(s,  1))**2+(fp(s,  1)-fp(s,  0))**2+ eps)

           rm = 1.0_rp-min(1.0_rp,min(eta0,eta1,eta2)/cr)

           ct = 1.0_rp*10.0_rp**(-floor(a1-a2*(1.0_rp-((1.0_rp-rm)**4*(1.0_rp+4.0_rp*rm)))))

           ! delta
           if(omega0 < ct) then
              d0 = 0
           else
              d0 = 1
           endif

           if(omega1 < ct) then
              d1 = 0
           else
              d1 = 1
           endif

           if(omega2 < ct) then
              d2 = 0
           else
              d2 = 1
           endif

           dsum = 1.0_rp/(d0*cweno5(0) + d1*cweno5(1) + d2*cweno5(2))

           omega0 = d0*cweno5(0)*dsum
           omega1 = d1*cweno5(1)*dsum
           omega2 = d2*cweno5(2)*dsum

           ! WENO plus reconstruction
           fhatp(s) = omega0*q0 + omega1*q1 + omega2*q2

           ! --- WENO minus part --- !

           ! polinomia
           q0 = aweno5(0,0)*fm(s,3) + aweno5(0,1)*fm(s,2) + aweno5(0,2)*fm(s, 1)
           q1 = aweno5(1,0)*fm(s,2) + aweno5(1,1)*fm(s,1) + aweno5(1,2)*fm(s, 0)
           q2 = aweno5(2,0)*fm(s,1) + aweno5(2,1)*fm(s,0) + aweno5(2,2)*fm(s,-1)

           ! smoothness indes
           IS0 = is1par*(fm(s,3)-2.0_rp*fm(s,2)+fm(s, 1))**2 + is2par*(    fm(s,3)-4.0_rp*fm(s,2)+3.0_rp*fm(s, 1))**2
           IS1 = is1par*(fm(s,2)-2.0_rp*fm(s,1)+fm(s, 0))**2 + is2par*(    fm(s,2)                      -fm(s, 0))**2
           IS2 = is1par*(fm(s,1)-2.0_rp*fm(s,0)+fm(s,-1))**2 + is2par*(3.0_rp*fm(s,1)-4.0_rp*fm(s,0)    +fm(s,-1))**2

           ! alpha
           absTerm0 = abs(IS0 - IS2)
           alpha0 = 1.0_rp + (absTerm0/(IS0 + varepsilon))**6
           alpha1 = 1.0_rp + (absTerm0/(IS1 + varepsilon))**6
           alpha2 = 1.0_rp + (absTerm0/(IS2 + varepsilon))**6

           ! omega
           isumAlpha = 1.0_rp/(alpha0 + alpha1 + alpha2)
           omega0 = alpha0*isumAlpha
           omega1 = alpha1*isumAlpha
           omega2 = alpha2*isumAlpha

           ! adaptive ct
           eps = 0.9_rp*cr/(1.0_rp-0.9_rp*cr)*csi**2
           eta0 = (2.0_rp*abs((fm(s,  1)-fm(s, 2))*(fm(s, 2)-fm(s, 3))) + eps)/ &
                         ((fm(s,  1)-fm(s, 2))**2+(fm(s, 2)-fm(s, 3))**2+ eps)
           eta1 = (2.0_rp*abs((fm(s,  0)-fm(s, 1))*(fm(s, 1)-fm(s, 2))) + eps)/ &
                         ((fm(s,  0)-fm(s, 1))**2+(fm(s, 1)-fm(s, 2))**2+ eps)
           eta2 = (2.0_rp*abs((fm(s, -1)-fm(s, 0))*(fm(s, 0)-fm(s, 1))) + eps)/ &
                         ((fm(s, -1)-fm(s, 0))**2+(fm(s, 0)-fm(s, 1))**2+ eps)

           rm = 1.0_rp-min(1.0_rp,min(eta0,eta1,eta2)/cr)

           ct = 1.0_rp*10.0_rp**(-floor(a1-a2*(1.0_rp-((1.0_rp-rm)**4*(1.0_rp+4.0_rp*rm)))))

           ! delta
           if(omega0 < ct) then
              d0 = 0
           else
              d0 = 1
           endif

           if(omega1 < ct) then
              d1 = 0
           else
              d1 = 1
           endif

           if(omega2 < ct) then
              d2 = 0
           else
              d2 = 1
           endif

           dsum = 1.0_rp/(d0*cweno5(0) + d1*cweno5(1) + d2*cweno5(2))

           omega0 = d0*cweno5(0)*dsum
           omega1 = d1*cweno5(1)*dsum
           omega2 = d2*cweno5(2)*dsum

           ! WENO minus riconstruction
           fhatm(s) = omega0*q0 + omega1*q1 + omega2*q2

        enddo


        endif

return
end subroutine weno_reconstruction



!subroutine teno_reconstruction(weno_num,a,cw,fvp,fvm,left,right,hF)
!! ---------------------------------------------------------------------------
!!
!!       This subroutine performs the TENO reconstruction in a cell bound.
!!       The cell bound is identified by the b parameter.
!!
!!       b = 0  -> left  cell bound
!!       b = -1 -> right cell bound
!!
!! ---------------------------------------------------------------------------
!        !$acc routine seq
!        implicit none
!
!        integer                                        , intent(in)   :: weno_num
!        real(rp), dimension(0:weno_num-1,0:weno_num-1) , intent(in )  :: a
!        real(rp), dimension(0:weno_num-1)              , intent(in )  :: cw
!        real(rp), dimension(5,-3:4)                  , intent(in )  :: fvp, fvm
!        real(rp), dimension(5,5)                   , intent(in )  :: left, right
!        real(rp), dimension(5)                       , intent(out)  :: hF
!
!        ! local declarations
!        real(rp), parameter          :: varepsilon = 1.0E-12_rp  !< pay attention on this in parallel
!        real(rp), parameter          :: a1 = 10.5_rp, a2 = 3.5_rp, cr = 0.25_rp, csi = 1.0E-3_rp
!        real(rp), parameter          :: is1par = 13.0_rp/12.0_rp, is2par = 0.25_rp
!
!        real(rp) :: q0, q1, q2, q3, Is0, Is1, Is2, Is3
!        real(rp) :: alpha0, alpha1, alpha2, alpha3
!        real(rp) :: omega0, omega1, omega2, omega3
!        real(rp) :: isumAlpha, absTerm0, absTerm1, dsum
!        real(rp) :: eta0, eta1, eta2, eps, rm
!        real(rp) :: ct = 1.0E-7_rp
!
!        real(rp), dimension(5,-3:4) :: fp, fm
!        real(rp), dimension(5)      :: fhatp, fhatm
!        real(rp), dimension(5)      :: fhatp_b, fhatm_b
!        integer :: d0, d1, d2
!        integer :: s, l, m
!
!        do s = -weno_num+1, weno_num
!           do m = 1,5
!              fp(m,s) = 0.0_rp
!              fm(m,s) = 0.0_rp
!              do l = 1,5
!              fp(m,s) = fp(m,s) + left(m,l)*fvp(l,s)
!              fm(m,s) = fm(m,s) + left(m,l)*fvm(l,s)
!              enddo
!           enddo
!        enddo
!
!        do s = 1, 5
!           ! --- WENO plus part --- !
!
!           ! polinomia
!           q0 = a(0,0)*fp(s,-2) + a(0,1)*fp(s,-1) + a(0,2)*fp(s,0)
!           q1 = a(1,0)*fp(s,-1) + a(1,1)*fp(s, 0) + a(1,2)*fp(s,1)
!           q2 = a(2,0)*fp(s, 0) + a(2,1)*fp(s, 1) + a(2,2)*fp(s,2)
!
!           ! smoothness index
!           IS0 = is1par*(fp(s,-2)-2.0_rp*fp(s,-1)+fp(s,0))**2 + is2par*(    fp(s,-2)-4.0_rp*fp(s,-1)+3.0_rp*fp(s,0))**2
!           IS1 = is1par*(fp(s,-1)-2.0_rp*fp(s, 0)+fp(s,1))**2 + is2par*(    fp(s,-1)                 -fp(s,1))**2
!           IS2 = is1par*(fp(s, 0)-2.0_rp*fp(s, 1)+fp(s,2))**2 + is2par*(3.0_rp*fp(s, 0)-4.0_rp*fp(s, 1)    +fp(s,2))**2
!
!           ! alpha
!           absTerm0 = abs(IS0 - IS2)
!           alpha0 = 1.0_rp + (absTerm0/(IS0 + varepsilon))**6
!           alpha1 = 1.0_rp + (absTerm0/(IS1 + varepsilon))**6
!           alpha2 = 1.0_rp + (absTerm0/(IS2 + varepsilon))**6
!
!           ! omega
!           isumAlpha = 1.0_rp/(alpha0 + alpha1 + alpha2)
!           omega0 = alpha0*isumAlpha
!           omega1 = alpha1*isumAlpha
!           omega2 = alpha2*isumAlpha
!
!           ! delta
!           if(omega0 < ct) then
!              d0 = 0
!           else
!              d0 = 1
!           endif
!
!           if(omega1 < ct) then
!              d1 = 0
!           else
!              d1 = 1
!           endif
!
!           if(omega2 < ct) then
!              d2 = 0
!           else
!              d2 = 1
!           endif
!
!           dsum = 1.0_rp/(d0*cw(0) + d1*cw(1) + d2*cw(2))
!
!           omega0 = d0*cw(0)*dsum
!           omega1 = d1*cw(1)*dsum
!           omega2 = d2*cw(2)*dsum
!
!           ! WENO plus reconstruction
!           fhatp(s) = omega0*q0 + omega1*q1 + omega2*q2
!
!           ! --- WENO minus part --- !
!
!           ! polinomia
!           q0 = a(0,0)*fm(s,3) + a(0,1)*fm(s,2) + a(0,2)*fm(s, 1)
!           q1 = a(1,0)*fm(s,2) + a(1,1)*fm(s,1) + a(1,2)*fm(s, 0)
!           q2 = a(2,0)*fm(s,1) + a(2,1)*fm(s,0) + a(2,2)*fm(s,-1)
!
!           ! smoothness indes
!           IS0 = is1par*(fm(s,3)-2.0_rp*fm(s,2)+fm(s, 1))**2 + is2par*(    fm(s,3)-4.0_rp*fm(s,2)+3.0_rp*fm(s, 1))**2
!           IS1 = is1par*(fm(s,2)-2.0_rp*fm(s,1)+fm(s, 0))**2 + is2par*(    fm(s,2)                -fm(s, 0))**2
!           IS2 = is1par*(fm(s,1)-2.0_rp*fm(s,0)+fm(s,-1))**2 + is2par*(3.0_rp*fm(s,1)-4.0_rp*fm(s,0)    +fm(s,-1))**2
!
!           ! alpha
!           absTerm0 = abs(IS0 - IS2)
!           alpha0 = 1.0_rp + (absTerm0/(IS0 + varepsilon))**6
!           alpha1 = 1.0_rp + (absTerm0/(IS1 + varepsilon))**6
!           alpha2 = 1.0_rp + (absTerm0/(IS2 + varepsilon))**6
!
!           ! omega
!           isumAlpha = 1.0_rp/(alpha0 + alpha1 + alpha2)
!           omega0 = alpha0*isumAlpha
!           omega1 = alpha1*isumAlpha
!           omega2 = alpha2*isumAlpha
!
!           ! delta
!           if(omega0 < ct) then
!              d0 = 0
!           else
!              d0 = 1
!           endif
!
!           if(omega1 < ct) then
!              d1 = 0
!           else
!              d1 = 1
!           endif
!
!           if(omega2 < ct) then
!              d2 = 0
!           else
!              d2 = 1
!           endif
!
!           dsum = 1.0_rp/(d0*cw(0) + d1*cw(1) + d2*cw(2))
!
!           omega0 = d0*cw(0)*dsum
!           omega1 = d1*cw(1)*dsum
!           omega2 = d2*cw(2)*dsum
!
!           ! WENO minus riconstruction
!           fhatm(s) = omega0*q0 + omega1*q1 + omega2*q2
!
!        enddo
!
!        ! return in the orinal system
!        do l = 1,5
!           fhatp_b(l) = 0.0_rp
!           fhatm_b(l) = 0.0_rp
!           do m = 1,5
!           fhatp_b(l) = fhatp_b(l) + right(l,m)*fhatp(m)
!           fhatm_b(l) = fhatm_b(l) + right(l,m)*fhatm(m)
!           enddo
!           hF(l) = fhatp_b(l) + fhatm_b(l)
!        enddo
!
!return
!end subroutine teno_reconstruction






!subroutine tenoA_reconstruction(weno_num,a,cw,fvp,fvm,left,right,hF)
!! ---------------------------------------------------------------------------
!!
!!       This subroutine performs the TENO reconstruction in a cell bound.
!!       The cell bound is identified by the b parameter.
!!
!!       b = 0  -> left  cell bound
!!       b = -1 -> right cell bound
!!
!! ---------------------------------------------------------------------------
!        !$acc routine seq
!        implicit none
!
!        integer                                        , intent(in)   :: weno_num
!        real(rp), dimension(0:weno_num-1,0:weno_num-1) , intent(in )  :: a
!        real(rp), dimension(0:weno_num-1)              , intent(in )  :: cw
!        real(rp), dimension(5,-3:4)                  , intent(in )  :: fvp, fvm
!        real(rp), dimension(5,5)                   , intent(in )  :: left, right
!        real(rp), dimension(5)                       , intent(out)  :: hF
!
!        ! local declarations
!        real(rp), parameter          :: varepsilon = 1.0E-12_rp  !< pay attention on this in parallel
!        real(rp), parameter          :: a1 = 10.5_rp, a2 = 3.5_rp, cr = 0.25_rp, csi = 1.0E-3_rp
!        real(rp), parameter          :: is1par = 13.0_rp/12.0_rp, is2par = 0.25_rp
!
!        real(rp) :: q0, q1, q2, q3, Is0, Is1, Is2, Is3
!        real(rp) :: alpha0, alpha1, alpha2, alpha3
!        real(rp) :: omega0, omega1, omega2, omega3
!        real(rp) :: isumAlpha, absTerm0, absTerm1, dsum
!        real(rp) :: eta0, eta1, eta2, eps, rm
!        real(rp) :: ct = 1.0E-7_rp
!
!        real(rp), dimension(5,-3:4) :: fp, fm
!        real(rp), dimension(5)      :: fhatp, fhatm
!        real(rp), dimension(5)      :: fhatp_b, fhatm_b
!        integer :: d0, d1, d2
!        integer :: s, l, m
!
!        do s = -weno_num+1, weno_num
!           do m = 1,5
!              fp(m,s) = 0.0_rp
!              fm(m,s) = 0.0_rp
!              do l = 1,5
!              fp(m,s) = fp(m,s) + left(m,l)*fvp(l,s)
!              fm(m,s) = fm(m,s) + left(m,l)*fvm(l,s)
!              enddo
!           enddo
!        enddo
!
!        do s = 1, 5
!           ! --- WENO plus part --- !
!
!           ! polinomia
!           q0 = a(0,0)*fp(s,-2) + a(0,1)*fp(s,-1) + a(0,2)*fp(s,0)
!           q1 = a(1,0)*fp(s,-1) + a(1,1)*fp(s, 0) + a(1,2)*fp(s,1)
!           q2 = a(2,0)*fp(s, 0) + a(2,1)*fp(s, 1) + a(2,2)*fp(s,2)
!
!           ! smoothness index
!           IS0 = is1par*(fp(s,-2)-2.0_rp*fp(s,-1)+fp(s,0))**2 + is2par*(    fp(s,-2)-4.0_rp*fp(s,-1)+3.0_rp*fp(s,0))**2
!           IS1 = is1par*(fp(s,-1)-2.0_rp*fp(s, 0)+fp(s,1))**2 + is2par*(    fp(s,-1)                 -fp(s,1))**2
!           IS2 = is1par*(fp(s, 0)-2.0_rp*fp(s, 1)+fp(s,2))**2 + is2par*(3.0_rp*fp(s, 0)-4.0_rp*fp(s, 1)    +fp(s,2))**2
!
!           ! alpha
!           absTerm0 = abs(IS0 - IS2)
!           alpha0 = 1.0_rp + (absTerm0/(IS0 + varepsilon))**6
!           alpha1 = 1.0_rp + (absTerm0/(IS1 + varepsilon))**6
!           alpha2 = 1.0_rp + (absTerm0/(IS2 + varepsilon))**6
!
!           ! omega
!           isumAlpha = 1.0_rp/(alpha0 + alpha1 + alpha2)
!           omega0 = alpha0*isumAlpha
!           omega1 = alpha1*isumAlpha
!           omega2 = alpha2*isumAlpha
!
!           ! adaptive ct
!           eps = 0.9_rp*cr/(1.0_rp-0.9_rp*cr)*csi**2
!           eta0 = (2.0_rp*abs((fp(s, 0)-fp(s, -1))*(fp(s, -1)-fp(s, -2))) + eps)/ &
!                         ((fp(s, 0)-fp(s, -1))**2+(fp(s, -1)-fp(s, -2))**2+ eps)
!           eta1 = (2.0_rp*abs((fp(s, 1)-fp(s,  0))*(fp(s,  0)-fp(s, -1))) + eps)/ &
!                         ((fp(s, 1)-fp(s,  0))**2+(fp(s,  0)-fp(s, -1))**2+ eps)
!           eta2 = (2.0_rp*abs((fp(s, 2)-fp(s,  1))*(fp(s,  1)-fp(s,  0))) + eps)/ &
!                         ((fp(s, 2)-fp(s,  1))**2+(fp(s,  1)-fp(s,  0))**2+ eps)
!
!           rm = 1.0_rp-min(1.0_rp,min(eta0,eta1,eta2)/cr)
!
!           ct = 1.0_rp*10.0_rp**(-floor(a1-a2*(1.0_rp-((1.0_rp-rm)**4*(1.0_rp+4.0_rp*rm)))))
!
!           ! delta
!           if(omega0 < ct) then
!              d0 = 0
!           else
!              d0 = 1
!           endif
!
!           if(omega1 < ct) then
!              d1 = 0
!           else
!              d1 = 1
!           endif
!
!           if(omega2 < ct) then
!              d2 = 0
!           else
!              d2 = 1
!           endif
!
!           dsum = 1.0_rp/(d0*cw(0) + d1*cw(1) + d2*cw(2))
!
!           omega0 = d0*cw(0)*dsum
!           omega1 = d1*cw(1)*dsum
!           omega2 = d2*cw(2)*dsum
!
!           ! WENO plus reconstruction
!           fhatp(s) = omega0*q0 + omega1*q1 + omega2*q2
!
!           ! --- WENO minus part --- !
!
!           ! polinomia
!           q0 = a(0,0)*fm(s,3) + a(0,1)*fm(s,2) + a(0,2)*fm(s, 1)
!           q1 = a(1,0)*fm(s,2) + a(1,1)*fm(s,1) + a(1,2)*fm(s, 0)
!           q2 = a(2,0)*fm(s,1) + a(2,1)*fm(s,0) + a(2,2)*fm(s,-1)
!
!           ! smoothness indes
!           IS0 = is1par*(fm(s,3)-2.0_rp*fm(s,2)+fm(s, 1))**2 + is2par*(    fm(s,3)-4.0_rp*fm(s,2)+3.0_rp*fm(s, 1))**2
!           IS1 = is1par*(fm(s,2)-2.0_rp*fm(s,1)+fm(s, 0))**2 + is2par*(    fm(s,2)                -fm(s, 0))**2
!           IS2 = is1par*(fm(s,1)-2.0_rp*fm(s,0)+fm(s,-1))**2 + is2par*(3.0_rp*fm(s,1)-4.0_rp*fm(s,0)    +fm(s,-1))**2
!
!           ! alpha
!           absTerm0 = abs(IS0 - IS2)
!           alpha0 = 1.0_rp + (absTerm0/(IS0 + varepsilon))**6
!           alpha1 = 1.0_rp + (absTerm0/(IS1 + varepsilon))**6
!           alpha2 = 1.0_rp + (absTerm0/(IS2 + varepsilon))**6
!
!           ! omega
!           isumAlpha = 1.0_rp/(alpha0 + alpha1 + alpha2)
!           omega0 = alpha0*isumAlpha
!           omega1 = alpha1*isumAlpha
!           omega2 = alpha2*isumAlpha
!
!           ! adaptive ct
!           eps = 0.9_rp*cr/(1.0_rp-0.9_rp*cr)*csi**2
!           eta0 = (2.0_rp*abs((fm(s,  1)-fm(s, 2))*(fm(s, 2)-fm(s, 3))) + eps)/ &
!                         ((fm(s,  1)-fm(s, 2))**2+(fm(s, 2)-fm(s, 3))**2+ eps)
!           eta1 = (2.0_rp*abs((fm(s,  0)-fm(s, 1))*(fm(s, 1)-fm(s, 2))) + eps)/ &
!                         ((fm(s,  0)-fm(s, 1))**2+(fm(s, 1)-fm(s, 2))**2+ eps)
!           eta2 = (2.0_rp*abs((fm(s, -1)-fm(s, 0))*(fm(s, 0)-fm(s, 1))) + eps)/ &
!                         ((fm(s, -1)-fm(s, 0))**2+(fm(s, 0)-fm(s, 1))**2+ eps)
!
!           rm = 1.0_rp-min(1.0_rp,min(eta0,eta1,eta2)/cr)
!
!           ct = 1.0_rp*10.0_rp**(-floor(a1-a2*(1.0_rp-((1.0_rp-rm)**4*(1.0_rp+4.0_rp*rm)))))
!
!           ! delta
!           if(omega0 < ct) then
!              d0 = 0
!           else
!              d0 = 1
!           endif
!
!           if(omega1 < ct) then
!              d1 = 0
!           else
!              d1 = 1
!           endif
!
!           if(omega2 < ct) then
!              d2 = 0
!           else
!              d2 = 1
!           endif
!
!           dsum = 1.0_rp/(d0*cw(0) + d1*cw(1) + d2*cw(2))
!
!           omega0 = d0*cw(0)*dsum
!           omega1 = d1*cw(1)*dsum
!           omega2 = d2*cw(2)*dsum
!
!           ! WENO minus riconstruction
!           fhatm(s) = omega0*q0 + omega1*q1 + omega2*q2
!
!        enddo
!
!        ! return in the orinal system
!        do l = 1,5
!           fhatp_b(l) = 0.0_rp
!           fhatm_b(l) = 0.0_rp
!           do m = 1,5
!           fhatp_b(l) = fhatp_b(l) + right(l,m)*fhatp(m)
!           fhatm_b(l) = fhatm_b(l) + right(l,m)*fhatm(m)
!           enddo
!           hF(l) = fhatp_b(l) + fhatm_b(l)
!        enddo
!
!return
!end subroutine tenoA_reconstruction















!subroutine teno_reconstruction(weno_num,a,cw,fvp,fvm,left,right,hF)
!! ---------------------------------------------------------------------------
!!
!!       This subroutine performs the TENO reconstruction in a cell bound.
!!       The cell bound is identified by the b parameter.
!!
!!       b = 0  -> left  cell bound
!!       b = -1 -> right cell bound
!!
!! ---------------------------------------------------------------------------
!        !$acc routine seq
!        implicit none
!
!        integer                                        , intent(in)   :: weno_num
!        real(rp), dimension(0:weno_num-1,0:weno_num-1) , intent(in )  :: a
!        real(rp), dimension(0:weno_num-1)              , intent(in )  :: cw
!        real(rp), dimension(5,-3:4)                  , intent(in )  :: fvp, fvm
!        real(rp), dimension(5,5)                   , intent(in )  :: left, right
!        real(rp), dimension(5)                       , intent(out)  :: hF
!
!        ! local declarations
!        real(rp), parameter          :: varepsilon = 1.0E-12_rp  !< pay attention on this in parallel
!        real(rp), parameter          :: a1 = 10.5_rp, a2 = 3.5_rp, cr = 0.25_rp, csi = 1.0E-3_rp
!        real(rp), parameter          :: is1par = 13.0_rp/12.0_rp, is2par = 0.25_rp
!
!        real(rp) :: q0, q1, q2, q3, Is0, Is1, Is2, Is3
!        real(rp) :: alpha0, alpha1, alpha2, alpha3
!        real(rp) :: omega0, omega1, omega2, omega3
!        real(rp) :: isumAlpha, absTerm0, absTerm1, dsum
!        real(rp) :: eta0, eta1, eta2, eps, rm
!        real(rp) :: ct = 1.0E-7_rp
!
!        real(rp), dimension(5,-3:4) :: fp, fm
!        real(rp), dimension(5)      :: fhatp, fhatm
!        real(rp), dimension(5)      :: fhatp_b, fhatm_b
!        integer :: d0, d1, d2
!        integer :: s, l, m
!
!        do s = -weno_num+1, weno_num
!           do m = 1,5
!              fp(m,s) = 0.0_rp
!              fm(m,s) = 0.0_rp
!              do l = 1,5
!              fp(m,s) = fp(m,s) + left(m,l)*fvp(l,s)
!              fm(m,s) = fm(m,s) + left(m,l)*fvm(l,s)
!              enddo
!           enddo
!        enddo
!
!        do s = 1, 5
!           ! --- WENO plus part --- !
!
!           ! polinomia
!           q0 = a(0,0)*fp(s,-2) + a(0,1)*fp(s,-1) + a(0,2)*fp(s,0)
!           q1 = a(1,0)*fp(s,-1) + a(1,1)*fp(s, 0) + a(1,2)*fp(s,1)
!           q2 = a(2,0)*fp(s, 0) + a(2,1)*fp(s, 1) + a(2,2)*fp(s,2)
!
!           ! smoothness index
!           IS0 = is1par*(fp(s,-2)-2.0_rp*fp(s,-1)+fp(s,0))**2 + is2par*(    fp(s,-2)-4.0_rp*fp(s,-1)+3.0_rp*fp(s,0))**2
!           IS1 = is1par*(fp(s,-1)-2.0_rp*fp(s, 0)+fp(s,1))**2 + is2par*(    fp(s,-1)                 -fp(s,1))**2
!           IS2 = is1par*(fp(s, 0)-2.0_rp*fp(s, 1)+fp(s,2))**2 + is2par*(3.0_rp*fp(s, 0)-4.0_rp*fp(s, 1)    +fp(s,2))**2
!
!           ! alpha
!           absTerm0 = abs(IS0 - IS2)
!           alpha0 = 1.0_rp + (absTerm0/(IS0 + varepsilon))**6
!           alpha1 = 1.0_rp + (absTerm0/(IS1 + varepsilon))**6
!           alpha2 = 1.0_rp + (absTerm0/(IS2 + varepsilon))**6
!
!           ! omega
!           isumAlpha = 1.0_rp/(alpha0 + alpha1 + alpha2)
!           omega0 = alpha0*isumAlpha
!           omega1 = alpha1*isumAlpha
!           omega2 = alpha2*isumAlpha
!
!           ! delta
!           if(omega0 < ct) then
!              d0 = 0
!           else
!              d0 = 1
!           endif
!
!           if(omega1 < ct) then
!              d1 = 0
!           else
!              d1 = 1
!           endif
!
!           if(omega2 < ct) then
!              d2 = 0
!           else
!              d2 = 1
!           endif
!
!           dsum = 1.0_rp/(d0*cw(0) + d1*cw(1) + d2*cw(2))
!
!           omega0 = d0*cw(0)*dsum
!           omega1 = d1*cw(1)*dsum
!           omega2 = d2*cw(2)*dsum
!
!           ! WENO plus reconstruction
!           fhatp(s) = omega0*q0 + omega1*q1 + omega2*q2
!
!           ! --- WENO minus part --- !
!
!           ! polinomia
!           q0 = a(0,0)*fm(s,3) + a(0,1)*fm(s,2) + a(0,2)*fm(s, 1)
!           q1 = a(1,0)*fm(s,2) + a(1,1)*fm(s,1) + a(1,2)*fm(s, 0)
!           q2 = a(2,0)*fm(s,1) + a(2,1)*fm(s,0) + a(2,2)*fm(s,-1)
!
!           ! smoothness indes
!           IS0 = is1par*(fm(s,3)-2.0_rp*fm(s,2)+fm(s, 1))**2 + is2par*(    fm(s,3)-4.0_rp*fm(s,2)+3.0_rp*fm(s, 1))**2
!           IS1 = is1par*(fm(s,2)-2.0_rp*fm(s,1)+fm(s, 0))**2 + is2par*(    fm(s,2)                -fm(s, 0))**2
!           IS2 = is1par*(fm(s,1)-2.0_rp*fm(s,0)+fm(s,-1))**2 + is2par*(3.0_rp*fm(s,1)-4.0_rp*fm(s,0)    +fm(s,-1))**2
!
!           ! alpha
!           absTerm0 = abs(IS0 - IS2)
!           alpha0 = 1.0_rp + (absTerm0/(IS0 + varepsilon))**6
!           alpha1 = 1.0_rp + (absTerm0/(IS1 + varepsilon))**6
!           alpha2 = 1.0_rp + (absTerm0/(IS2 + varepsilon))**6
!
!           ! omega
!           isumAlpha = 1.0_rp/(alpha0 + alpha1 + alpha2)
!           omega0 = alpha0*isumAlpha
!           omega1 = alpha1*isumAlpha
!           omega2 = alpha2*isumAlpha
!
!           ! delta
!           if(omega0 < ct) then
!              d0 = 0
!           else
!              d0 = 1
!           endif
!
!           if(omega1 < ct) then
!              d1 = 0
!           else
!              d1 = 1
!           endif
!
!           if(omega2 < ct) then
!              d2 = 0
!           else
!              d2 = 1
!           endif
!
!           dsum = 1.0_rp/(d0*cw(0) + d1*cw(1) + d2*cw(2))
!
!           omega0 = d0*cw(0)*dsum
!           omega1 = d1*cw(1)*dsum
!           omega2 = d2*cw(2)*dsum
!
!           ! WENO minus riconstruction
!           fhatm(s) = omega0*q0 + omega1*q1 + omega2*q2
!
!        enddo
!
!        ! return in the orinal system
!        do l = 1,5
!           fhatp_b(l) = 0.0_rp
!           fhatm_b(l) = 0.0_rp
!           do m = 1,5
!           fhatp_b(l) = fhatp_b(l) + right(l,m)*fhatp(m)
!           fhatm_b(l) = fhatm_b(l) + right(l,m)*fhatm(m)
!           enddo
!           hF(l) = fhatp_b(l) + fhatm_b(l)
!        enddo
!
!return
!end subroutine teno_reconstruction






!subroutine tenoA_reconstruction(weno_num,a,cw,fvp,fvm,left,right,hF)
!! ---------------------------------------------------------------------------
!!
!!       This subroutine performs the TENO reconstruction in a cell bound.
!!       The cell bound is identified by the b parameter.
!!
!!       b = 0  -> left  cell bound
!!       b = -1 -> right cell bound
!!
!! ---------------------------------------------------------------------------
!        !$acc routine seq
!        implicit none
!
!        integer                                        , intent(in)   :: weno_num
!        real(rp), dimension(0:weno_num-1,0:weno_num-1) , intent(in )  :: a
!        real(rp), dimension(0:weno_num-1)              , intent(in )  :: cw
!        real(rp), dimension(5,-3:4)                  , intent(in )  :: fvp, fvm
!        real(rp), dimension(5,5)                   , intent(in )  :: left, right
!        real(rp), dimension(5)                       , intent(out)  :: hF
!
!        ! local declarations
!        real(rp), parameter          :: varepsilon = 1.0E-12_rp  !< pay attention on this in parallel
!        real(rp), parameter          :: a1 = 10.5_rp, a2 = 3.5_rp, cr = 0.25_rp, csi = 1.0E-3_rp
!        real(rp), parameter          :: is1par = 13.0_rp/12.0_rp, is2par = 0.25_rp
!
!        real(rp) :: q0, q1, q2, q3, Is0, Is1, Is2, Is3
!        real(rp) :: alpha0, alpha1, alpha2, alpha3
!        real(rp) :: omega0, omega1, omega2, omega3
!        real(rp) :: isumAlpha, absTerm0, absTerm1, dsum
!        real(rp) :: eta0, eta1, eta2, eps, rm
!        real(rp) :: ct = 1.0E-7_rp
!
!        real(rp), dimension(5,-3:4) :: fp, fm
!        real(rp), dimension(5)      :: fhatp, fhatm
!        real(rp), dimension(5)      :: fhatp_b, fhatm_b
!        integer :: d0, d1, d2
!        integer :: s, l, m
!
!        do s = -weno_num+1, weno_num
!           do m = 1,5
!              fp(m,s) = 0.0_rp
!              fm(m,s) = 0.0_rp
!              do l = 1,5
!              fp(m,s) = fp(m,s) + left(m,l)*fvp(l,s)
!              fm(m,s) = fm(m,s) + left(m,l)*fvm(l,s)
!              enddo
!           enddo
!        enddo
!
!        do s = 1, 5
!           ! --- WENO plus part --- !
!
!           ! polinomia
!           q0 = a(0,0)*fp(s,-2) + a(0,1)*fp(s,-1) + a(0,2)*fp(s,0)
!           q1 = a(1,0)*fp(s,-1) + a(1,1)*fp(s, 0) + a(1,2)*fp(s,1)
!           q2 = a(2,0)*fp(s, 0) + a(2,1)*fp(s, 1) + a(2,2)*fp(s,2)
!
!           ! smoothness index
!           IS0 = is1par*(fp(s,-2)-2.0_rp*fp(s,-1)+fp(s,0))**2 + is2par*(    fp(s,-2)-4.0_rp*fp(s,-1)+3.0_rp*fp(s,0))**2
!           IS1 = is1par*(fp(s,-1)-2.0_rp*fp(s, 0)+fp(s,1))**2 + is2par*(    fp(s,-1)                 -fp(s,1))**2
!           IS2 = is1par*(fp(s, 0)-2.0_rp*fp(s, 1)+fp(s,2))**2 + is2par*(3.0_rp*fp(s, 0)-4.0_rp*fp(s, 1)    +fp(s,2))**2
!
!           ! alpha
!           absTerm0 = abs(IS0 - IS2)
!           alpha0 = 1.0_rp + (absTerm0/(IS0 + varepsilon))**6
!           alpha1 = 1.0_rp + (absTerm0/(IS1 + varepsilon))**6
!           alpha2 = 1.0_rp + (absTerm0/(IS2 + varepsilon))**6
!
!           ! omega
!           isumAlpha = 1.0_rp/(alpha0 + alpha1 + alpha2)
!           omega0 = alpha0*isumAlpha
!           omega1 = alpha1*isumAlpha
!           omega2 = alpha2*isumAlpha
!
!           ! adaptive ct
!           eps = 0.9_rp*cr/(1.0_rp-0.9_rp*cr)*csi**2
!           eta0 = (2.0_rp*abs((fp(s, 0)-fp(s, -1))*(fp(s, -1)-fp(s, -2))) + eps)/ &
!                         ((fp(s, 0)-fp(s, -1))**2+(fp(s, -1)-fp(s, -2))**2+ eps)
!           eta1 = (2.0_rp*abs((fp(s, 1)-fp(s,  0))*(fp(s,  0)-fp(s, -1))) + eps)/ &
!                         ((fp(s, 1)-fp(s,  0))**2+(fp(s,  0)-fp(s, -1))**2+ eps)
!           eta2 = (2.0_rp*abs((fp(s, 2)-fp(s,  1))*(fp(s,  1)-fp(s,  0))) + eps)/ &
!                         ((fp(s, 2)-fp(s,  1))**2+(fp(s,  1)-fp(s,  0))**2+ eps)
!
!           rm = 1.0_rp-min(1.0_rp,min(eta0,eta1,eta2)/cr)
!
!           ct = 1.0_rp*10.0_rp**(-floor(a1-a2*(1.0_rp-((1.0_rp-rm)**4*(1.0_rp+4.0_rp*rm)))))
!
!           ! delta
!           if(omega0 < ct) then
!              d0 = 0
!           else
!              d0 = 1
!           endif
!
!           if(omega1 < ct) then
!              d1 = 0
!           else
!              d1 = 1
!           endif
!
!           if(omega2 < ct) then
!              d2 = 0
!           else
!              d2 = 1
!           endif
!
!           dsum = 1.0_rp/(d0*cw(0) + d1*cw(1) + d2*cw(2))
!
!           omega0 = d0*cw(0)*dsum
!           omega1 = d1*cw(1)*dsum
!           omega2 = d2*cw(2)*dsum
!
!           ! WENO plus reconstruction
!           fhatp(s) = omega0*q0 + omega1*q1 + omega2*q2
!
!           ! --- WENO minus part --- !
!
!           ! polinomia
!           q0 = a(0,0)*fm(s,3) + a(0,1)*fm(s,2) + a(0,2)*fm(s, 1)
!           q1 = a(1,0)*fm(s,2) + a(1,1)*fm(s,1) + a(1,2)*fm(s, 0)
!           q2 = a(2,0)*fm(s,1) + a(2,1)*fm(s,0) + a(2,2)*fm(s,-1)
!
!           ! smoothness indes
!           IS0 = is1par*(fm(s,3)-2.0_rp*fm(s,2)+fm(s, 1))**2 + is2par*(    fm(s,3)-4.0_rp*fm(s,2)+3.0_rp*fm(s, 1))**2
!           IS1 = is1par*(fm(s,2)-2.0_rp*fm(s,1)+fm(s, 0))**2 + is2par*(    fm(s,2)                -fm(s, 0))**2
!           IS2 = is1par*(fm(s,1)-2.0_rp*fm(s,0)+fm(s,-1))**2 + is2par*(3.0_rp*fm(s,1)-4.0_rp*fm(s,0)    +fm(s,-1))**2
!
!           ! alpha
!           absTerm0 = abs(IS0 - IS2)
!           alpha0 = 1.0_rp + (absTerm0/(IS0 + varepsilon))**6
!           alpha1 = 1.0_rp + (absTerm0/(IS1 + varepsilon))**6
!           alpha2 = 1.0_rp + (absTerm0/(IS2 + varepsilon))**6
!
!           ! omega
!           isumAlpha = 1.0_rp/(alpha0 + alpha1 + alpha2)
!           omega0 = alpha0*isumAlpha
!           omega1 = alpha1*isumAlpha
!           omega2 = alpha2*isumAlpha
!
!           ! adaptive ct
!           eps = 0.9_rp*cr/(1.0_rp-0.9_rp*cr)*csi**2
!           eta0 = (2.0_rp*abs((fm(s,  1)-fm(s, 2))*(fm(s, 2)-fm(s, 3))) + eps)/ &
!                         ((fm(s,  1)-fm(s, 2))**2+(fm(s, 2)-fm(s, 3))**2+ eps)
!           eta1 = (2.0_rp*abs((fm(s,  0)-fm(s, 1))*(fm(s, 1)-fm(s, 2))) + eps)/ &
!                         ((fm(s,  0)-fm(s, 1))**2+(fm(s, 1)-fm(s, 2))**2+ eps)
!           eta2 = (2.0_rp*abs((fm(s, -1)-fm(s, 0))*(fm(s, 0)-fm(s, 1))) + eps)/ &
!                         ((fm(s, -1)-fm(s, 0))**2+(fm(s, 0)-fm(s, 1))**2+ eps)
!
!           rm = 1.0_rp-min(1.0_rp,min(eta0,eta1,eta2)/cr)
!
!           ct = 1.0_rp*10.0_rp**(-floor(a1-a2*(1.0_rp-((1.0_rp-rm)**4*(1.0_rp+4.0_rp*rm)))))
!
!           ! delta
!           if(omega0 < ct) then
!              d0 = 0
!           else
!              d0 = 1
!           endif
!
!           if(omega1 < ct) then
!              d1 = 0
!           else
!              d1 = 1
!           endif
!
!           if(omega2 < ct) then
!              d2 = 0
!           else
!              d2 = 1
!           endif
!
!           dsum = 1.0_rp/(d0*cw(0) + d1*cw(1) + d2*cw(2))
!
!           omega0 = d0*cw(0)*dsum
!           omega1 = d1*cw(1)*dsum
!           omega2 = d2*cw(2)*dsum
!
!           ! WENO minus riconstruction
!           fhatm(s) = omega0*q0 + omega1*q1 + omega2*q2
!
!        enddo
!
!        ! return in the orinal system
!        do l = 1,5
!           fhatp_b(l) = 0.0_rp
!           fhatm_b(l) = 0.0_rp
!           do m = 1,5
!           fhatp_b(l) = fhatp_b(l) + right(l,m)*fhatp(m)
!           fhatm_b(l) = fhatm_b(l) + right(l,m)*fhatm(m)
!           enddo
!           hF(l) = fhatp_b(l) + fhatm_b(l)
!        enddo
!
!return
!end subroutine tenoA_reconstruction



















subroutine central_finite_differences(lbx,ubx,lby,uby,lbz,ubz)
! -----------------------------------------------------------------------
!
!       This subroutine solves the RHS of Navier-Stokes equation with a
!       central finite difference scheme.
!
! -----------------------------------------------------------------------
        implicit none
        integer , intent(in) :: lbx, ubx, lby,uby, lbz,ubz
        real(rp), dimension(lbx:ubx,5) :: phi_arr_x,flx,flx_arr_x
        real(rp), dimension(lby:uby,5) :: phi_arr_y,fly,flx_arr_y
        real(rp), dimension(lbz:ubz,5) :: phi_arr_z,flz,flx_arr_z
        real(rp), dimension(lbx:ubx,6) :: pri_1D_x
        real(rp), dimension(lby:uby,6) :: pri_1D_y
        real(rp), dimension(lbz:ubz,6) :: pri_1D_z
        real(rp), dimension(5)         :: hF, hG, hH
        integer                        :: s

        !$omp parallel do collapse(2) default(private), shared(sx,ex,sy,ey,sz,ez,x,RHS,xstep_i),&
        !$omp shared(central_1_one_half,fd_L,fd_R,lbx,ubx,phi)
        do k = sz,ez
           do j = sy,ey

              phi_arr_x(lbx:ubx,:) = phi(lbx:ubx,j,k,:)
              call conv_flux_x(lbx,ubx,phi_arr_x,flx_arr_x,pri_1D_x)

              do i = sx-1, ex

                 hF = 0.0_rp
                 do s = fd_L, fd_R
                    hF(:) = hF(:) + central_1_one_half(s) * flx_arr_x(i+s,:)
                 enddo

                 flx(i,:) = hf(:)
                 RHS(i,j,k,:) = - xstep_i(i) * (flx(i,:) - flx(i-1,:))

              enddo
           enddo
        enddo
        !$omp end parallel do

        !$omp parallel do collapse(2) default(private), shared(sx,ex,sy,ey,sz,ez,y,RHS,ystep_i),&
        !$omp shared(central_1_one_half,fd_L,fd_R,lby,uby,phi)
        do k = sz,ez
           do i = sx,ex

              phi_arr_y(lby:uby,:) = phi(i,lby:uby,k,:)
              call conv_flux_y(lby,uby,phi_arr_y,flx_arr_y,pri_1D_y)

              do j = sy-1, ey

                 hG = 0.0_rp
                 do s = fd_L, fd_R

                    hG(:) = hG(:) + central_1_one_half(s) * flx_arr_y(j+s,:)

                 enddo

                 fly(j,:) = hG(:)
                 RHS(i,j,k,:) = RHS(i,j,k,:) - ystep_i(j) * (fly(j,:) - fly(j-1,:))

              enddo
           enddo
        enddo
        !$omp end parallel do

        if(dims == 3) then

          !$omp parallel do collapse(2) default(private), shared(sx,ex,sy,ey,sz,ez,z,RHS,zstep_i), &
          !$omp shared(central_1_one_half,fd_L,fd_R,lbz,ubz,phi)
          do j = sy,ey
             do i = sx,ex

                phi_arr_z(lbz:ubz,:) = phi(i,j,lbz:ubz,:)
                call conv_flux_z(lbz,ubz,phi_arr_z,flx_arr_z,pri_1D_z)

                do k = sz-1, ez

                   hH = 0.0_rp
                   do s = fd_L, fd_R

                      hH(:) = hH(:) + central_1_one_half(s) * flx_arr_z(k+s,:)

                   enddo

                   flz(k,:) = hH(:)
                   RHS(i,j,k,:) = RHS(i,j,k,:) - zstep_i(k) * (flz(k,:) - flz(k-1,:))

                enddo
             enddo
          enddo
          !$omp end parallel do

        endif

return
endsubroutine central_finite_differences


subroutine energy_preservingX(lbx,ubx,ltot,istr,iend,tilde_op_x,pri_1D_x,RHS)

        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: RHS
        real(rp), allocatable, dimension(:,:,:)  , intent(inout) :: tilde_op_x
        real(rp), allocatable, dimension(:,:)    , intent(inout) :: pri_1D_x

        integer, intent(in) :: lbx,ubx
        integer, intent(in) :: ltot
        integer, intent(in) :: istr,iend

        integer                             :: i,j,k,l, il,m, n

        real(rp) :: ir, u_, v_, w_, ek, p_, weight, inner_sum, hFn
        real(rp) :: phi1, phi2, phi3, phi4, phi5, idx
        real(rp), parameter :: gm1 = gamma0-1.0_rp
        real(rp), parameter :: hgm = gamma0/gm1
        real(rp), parameter :: one8= 1.0_rp/8.0_rp

        call StartProfRange("energy_preservingX")

        !$acc parallel default(present)
        !$acc loop gang collapse(2) private(pri_1D_x,tilde_op_x,flx_x)
        do k = sz,ez
           do j = sy,ey

              !$acc loop vector
              do i = lbx, ubx

                 phi1 = phi(i,j,k,1)
                 phi2 = phi(i,j,k,2)
                 phi3 = phi(i,j,k,3)
                 phi4 = phi(i,j,k,4)
                 phi5 = phi(i,j,k,5)

                 ir = 1.0_rp/phi1
                 u_ = phi2*ir
                 v_ = phi3*ir
                 w_ = phi4*ir
                 ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
                 p_ = gm1 * (phi5 - phi1*ek)

                 pri_1D_x(i,1) = p_
                 pri_1D_x(i,2) = u_
                 pri_1D_x(i,3) = v_
                 pri_1D_x(i,4) = w_
                 pri_1D_x(i,5) = hgm*p_*ir + ek
                 pri_1D_x(i,6) = phi1

              enddo

              !$acc loop vector collapse(2)
              do l = 1,ltot
                 do i = lbx, ubx-3

                    il = i+l

                    weight = one8 * (pri_1D_x(i,6) + pri_1D_x(il,6)) &
                                  * (pri_1D_x(i,2) + pri_1D_x(il,2))

                    tilde_op_x(i,l,1) = 2*weight
                    tilde_op_x(i,l,2) = weight*(pri_1D_x(i,2) + pri_1D_x(il,2)) + 0.5_rp*(pri_1D_x(i,1) + pri_1D_x(il,1))
                    tilde_op_x(i,l,3) = weight*(pri_1D_x(i,3) + pri_1D_x(il,3))
                    tilde_op_x(i,l,4) = weight*(pri_1D_x(i,4) + pri_1D_x(il,4))
                    tilde_op_x(i,l,5) = weight*(pri_1D_x(i,5) + pri_1D_x(il,5))

                 enddo
              enddo

              !$acc loop vector
              do i = sx-1,ex

                 do n = 1,5
                    hFn = 0.0_rp
                    do l = 1,Ltot
                       inner_sum = 0.0_rp
                       do m = 0,l-1
                          inner_sum = inner_sum + tilde_op_x(i-m,l,n)
                       enddo

                       hFn = hFn + 2*central_1(l)*inner_sum
                    enddo
                    flx_x(i,n) = hFn
                 enddo

               enddo

               !$acc loop vector
               do i = istr,iend
                  idx = xstep_i(i)
                  RHS(i,j,k,1) = - idx * (flx_x(i,1) - flx_x(i-1,1))
                  RHS(i,j,k,2) = - idx * (flx_x(i,2) - flx_x(i-1,2))
                  RHS(i,j,k,3) = - idx * (flx_x(i,3) - flx_x(i-1,3))
                  RHS(i,j,k,4) = - idx * (flx_x(i,4) - flx_x(i-1,4))
                  RHS(i,j,k,5) = - idx * (flx_x(i,5) - flx_x(i-1,5))
               enddo

           enddo
        enddo
        !$acc end parallel
        !$omp end parallel do

        call EndProfRange

        return
end subroutine energy_preservingX


subroutine energy_preservingY(lby,uby,ltot,jstr,jend,tilde_op_y,pri_1D_y,RHS)

        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: RHS
        real(rp), allocatable, dimension(:,:,:)  , intent(inout) :: tilde_op_y
        real(rp), allocatable, dimension(:,:)    , intent(inout) :: pri_1D_y

        integer, intent(in) :: lby, uby
        integer, intent(in) :: ltot
        integer, intent(in) :: jstr,jend

        integer                             :: i,j,k,l, m,jl, n

        real(rp) :: ir, u_, v_, w_, ek, p_, weight, inner_sum, hGn
        real(rp) :: phi1, phi2, phi3, phi4, phi5, idy
        real(rp), parameter :: gm1 = gamma0-1.0_rp
        real(rp), parameter :: hgm = gamma0/gm1
        real(rp), parameter :: one8= 1.0_rp/8.0_rp

        call StartProfRange("energy_preservingY")

        !$acc parallel default(present)
        !$acc loop gang collapse(2) private(pri_1D_y,tilde_op_y,flx_y)
        do k = sz,ez
           do i = sx,ex

              !$acc loop vector
              do j = lby,uby

                 phi1 = phi(i,j,k,1)
                 phi2 = phi(i,j,k,2)
                 phi3 = phi(i,j,k,3)
                 phi4 = phi(i,j,k,4)
                 phi5 = phi(i,j,k,5)

                 ir = 1.0_rp/phi1
                 u_ = phi2*ir
                 v_ = phi3*ir
                 w_ = phi4*ir
                 ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
                 p_ = gm1 * (phi5 - phi1*ek)

                 pri_1D_y(j,1) = p_
                 pri_1D_y(j,2) = u_
                 pri_1D_y(j,3) = v_
                 pri_1D_y(j,4) = w_
                 pri_1D_y(j,5) = hgm*p_*ir + ek
                 pri_1D_y(j,6) = phi1

              enddo

              !$acc loop vector collapse(2)
              do l = 1,ltot
                 do j = lby, uby-3

                    jl = j+l

                    weight = one8 * (pri_1D_y(j,6) + pri_1D_y(jl,6)) &
                                  * (pri_1D_y(j,3) + pri_1D_y(jl,3))

                    tilde_op_y(j,l,1) = 2*weight
                    tilde_op_y(j,l,2) = weight*(pri_1D_y(j,2) + pri_1D_y(jl,2))
                    tilde_op_y(j,l,3) = weight*(pri_1D_y(j,3) + pri_1D_y(jl,3)) + 0.5_rp*(pri_1D_y(j,1) + pri_1D_y(jl,1))
                    tilde_op_y(j,l,4) = weight*(pri_1D_y(j,4) + pri_1D_y(jl,4))
                    tilde_op_y(j,l,5) = weight*(pri_1D_y(j,5) + pri_1D_y(jl,5))

                 enddo
              enddo

              !$acc loop vector
              do j = sy-1,ey

                 do n = 1,5
                    hGn = 0.0_rp
                    do l = 1,Ltot
                       inner_sum = 0.0_rp
                       do m = 0,l-1
                          inner_sum = inner_sum + tilde_op_y(j-m,l,n)
                       enddo

                       hGn = hGn + 2*central_1(l)*inner_sum
                    enddo
                    flx_y(j,n) = hGn
                 enddo

              enddo

              !$acc loop vector
              do j = jstr,jend
                 idy = ystep_i(j)
                 RHS(i,j,k,1) = RHS(i,j,k,1) - idy * (flx_y(j,1) - flx_y(j-1,1))
                 RHS(i,j,k,2) = RHS(i,j,k,2) - idy * (flx_y(j,2) - flx_y(j-1,2))
                 RHS(i,j,k,3) = RHS(i,j,k,3) - idy * (flx_y(j,3) - flx_y(j-1,3))
                 RHS(i,j,k,4) = RHS(i,j,k,4) - idy * (flx_y(j,4) - flx_y(j-1,4))
                 RHS(i,j,k,5) = RHS(i,j,k,5) - idy * (flx_y(j,5) - flx_y(j-1,5))
              enddo

           enddo
        enddo
        !$acc end parallel
        !$omp end parallel do

        call EndProfRange

        return
end subroutine energy_preservingY

subroutine energy_preservingZ(lbz,ubz,ltot,kstr,kend,tilde_op_z,pri_1D_z,RHS)

        implicit none
        real(rp), allocatable, dimension(:,:,:,:), intent(inout) :: RHS
        real(rp), allocatable, dimension(:,:,:)  , intent(inout) :: tilde_op_z
        real(rp), allocatable, dimension(:,:)    , intent(inout) :: pri_1D_z

        integer, intent(in) :: lbz,ubz
        integer, intent(in) :: ltot
        integer, intent(in) :: kstr,kend

        integer                             :: i,j,k, l, m, kl, n

        real(rp) :: ir, u_, v_, w_, ek, p_, weight, inner_sum, hHn
        real(rp) :: phi1, phi2, phi3, phi4, phi5, idz
        real(rp), parameter :: gm1 = gamma0-1.0_rp
        real(rp), parameter :: hgm = gamma0/gm1
        real(rp), parameter :: one8= 1.0_rp/8.0_rp

        call StartProfRange("energy_preservingZ")

        !$acc parallel default(present)
        !$acc loop gang collapse(2) private(pri_1D_z,tilde_op_z,flx_z)
        do j = sy,ey
           do i = sx,ex

              !$acc loop vector
              do k = lbz,ubz

                 phi1 = phi(i,j,k,1)
                 phi2 = phi(i,j,k,2)
                 phi3 = phi(i,j,k,3)
                 phi4 = phi(i,j,k,4)
                 phi5 = phi(i,j,k,5)

                 ir = 1.0_rp/phi1
                 u_ = phi2*ir
                 v_ = phi3*ir
                 w_ = phi4*ir
                 ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
                 p_ = gm1 * (phi5 - phi1*ek)

                 pri_1D_z(k,1) = p_
                 pri_1D_z(k,2) = u_
                 pri_1D_z(k,3) = v_
                 pri_1D_z(k,4) = w_
                 pri_1D_z(k,5) = hgm*p_*ir + ek
                 pri_1D_z(k,6) = phi1

              enddo

              !$acc loop vector collapse(2)
              do l = 1,ltot
                 do k = lbz, ubz-3

                    kl = k+l

                    weight = one8 * (pri_1D_z(k,6) + pri_1D_z(kl,6)) &
                                  * (pri_1D_z(k,4) + pri_1D_z(kl,4))

                    tilde_op_z(k,l,1) = 2*weight
                    tilde_op_z(k,l,2) = weight*(pri_1D_z(k,2) + pri_1D_z(kl,2))
                    tilde_op_z(k,l,3) = weight*(pri_1D_z(k,3) + pri_1D_z(kl,3))
                    tilde_op_z(k,l,4) = weight*(pri_1D_z(k,4) + pri_1D_z(kl,4)) + 0.5_rp*(pri_1D_z(k,1) + pri_1D_z(kl,1))
                    tilde_op_z(k,l,5) = weight*(pri_1D_z(k,5) + pri_1D_z(kl,5))

                 enddo
              enddo

              !$acc loop vector
              do k = sz-1,ez

                 do n = 1,5
                    hHn = 0.0_rp
                    do l = 1,Ltot
                       inner_sum = 0.0_rp
                       do m = 0,l-1
                          inner_sum = inner_sum + tilde_op_z(k-m,l,n)
                       enddo

                       hHn = hHn + 2*central_1(l)*inner_sum
                    enddo
                    flx_z(k,n) = hHn
                 enddo
              enddo

              !$acc loop vector
              do k = kstr,kend
                 idz = zstep_i(k)
                 RHS(i,j,k,1) = RHS(i,j,k,1) - idz * (flx_z(k,1) - flx_z(k-1,1))
                 RHS(i,j,k,2) = RHS(i,j,k,2) - idz * (flx_z(k,2) - flx_z(k-1,2))
                 RHS(i,j,k,3) = RHS(i,j,k,3) - idz * (flx_z(k,3) - flx_z(k-1,3))
                 RHS(i,j,k,4) = RHS(i,j,k,4) - idz * (flx_z(k,4) - flx_z(k-1,4))
                 RHS(i,j,k,5) = RHS(i,j,k,5) - idz * (flx_z(k,5) - flx_z(k-1,5))
              enddo

           enddo
        enddo
        !$acc end parallel
        !$omp end parallel do

        call EndProfRange

        return
end subroutine energy_preservingZ



subroutine rhs_linear_ode
        implicit none

        do k = sz,ez
           do j = sy,ey
              do i = sx,ex

                 RHS(i,j,k,1) = - phi(i,j,k,1)
                 RHS(i,j,k,2) = - phi(i,j,k,2)
                 RHS(i,j,k,3) = - phi(i,j,k,3)
                 RHS(i,j,k,4) = - phi(i,j,k,4)
                 RHS(i,j,k,5) = - phi(i,j,k,5)

              enddo
           enddo
        enddo
        return
end subroutine rhs_linear_ode

subroutine rhs_linear_advection
        implicit none
        real(rp) :: df_dx, i_dx
        integer  :: s


        do k = sz,ez
           do j = sy,ey
              do i = sx,ex

                 i_dx = 1._rp/(0.5_rp*(x(i+1) - x(i-1)))

                 df_dx = 0._rp
                 do s = -central_fd_order/2, central_fd_order/2

                         df_dx = df_dx + i_dx * central_1(s) * phi(i+s,j,k,1)

                 enddo

                 RHS(i,j,k,1) = - df_dx
                 RHS(i,j,k,2) =   0._rp
                 RHS(i,j,k,3) =   0._rp
                 RHS(i,j,k,4) =   0._rp
                 RHS(i,j,k,5) =   0._rp

              enddo
           enddo
        enddo
return
end subroutine rhs_linear_advection







subroutine CharacteristicRHSX(RHS,beta)
        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: RHS
        real(rp)                                 , intent(in)    :: beta

        ! local declaration
        real(rp), dimension(5)    :: lambda, L, d
        real(rp), dimension(-2:0) :: cl
        real(rp) :: r_, ir, u_, v_, w_, p_, T_, c, ek
        real(rp) :: r_is, iris, u_is, v_is, w_is, p_is
        real(rp) :: r_x, p_x, u_x, v_x, w_x, bw_s
        real(rp) :: idx
        integer  :: i,j,k,s, is

        cl(-2) =  0.5_rp
        cl(-1) = -2.0_rp
        cl( 0) =  1.5_rp

        do k = sz,ez
           do j = sy,ey

              i = ex
              RHS(i,j,k,:) = beta * RHS(i,j,k,:)

              r_ = phi(i,j,k,1)
              ir = 1.0_rp/r_
              u_ = phi(i,j,k,2)*ir
              v_ = phi(i,j,k,3)*ir
              w_ = phi(i,j,k,4)*ir
              ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
              p_ = cv_i * (phi(i,j,k,5) - r_*ek)
              T_ = p_*ir
              c  = sqrt(gamma0*T_)

              idx = xstep_i(i)

              r_x = 0.0_rp
              p_x = 0.0_rp
              u_x = 0.0_rp
              v_x = 0.0_rp
              w_x = 0.0_rp
              do s = -2,0
               is = i + s
               bw_s = cl(s)

               r_is = phi(is,j,k,1)
               iris = 1.0_rp/r_is
               u_is = phi(is,j,k,2)*iris
               v_is = phi(is,j,k,3)*iris
               w_is = phi(is,j,k,4)*iris
               p_is = cv_i*(phi(is,j,k,5) - 0.5_rp*r_is*(u_is*u_is + v_is*v_is + w_is*w_is))
               r_x = r_x + bw_s * r_is
               p_x = p_x + bw_s * p_is
               u_x = u_x + bw_s * u_is
               v_x = v_x + bw_s * v_is
               w_x = w_x + bw_s * w_is
              enddo
              r_x = r_x * idx
              p_x = p_x * idx
              u_x = u_x * idx
              v_x = v_x * idx
              w_x = w_x * idx

              lambda(1) = 0.5_rp * (u_ - c + abs(u_ - c))
              lambda(2) = 0.5_rp * (u_     + abs(u_    ))
              lambda(3) = 0.5_rp * (u_     + abs(u_    ))
              lambda(4) = 0.5_rp * (u_     + abs(u_    ))
              lambda(5) = 0.5_rp * (u_ + c + abs(u_ + c))

              L(1) = lambda(1) * (p_x - r_*c*u_x)
              if(u_/c < 1.0_rp .and. L1_wave == 'poinsot_model') then
                 L(1) = 0.27_rp*(1.0_rp - (u_/c)**2) * c/Lx * (p_ - 1.0_rp)
              endif
              L(2) = lambda(2) * ((c**2)*r_x - p_x)
              L(3) = lambda(3) * (v_x)
              L(4) = lambda(4) * (w_x)
              L(5) = lambda(5) * (p_x + r_*c*u_x)

              d(1) = 1._rp/c**2 *(L(2)+0.5_rp*(L(5)+L(1)))
              d(2) = 1._rp/(2*r_*c)*(L(5)-L(1))
              d(3) = L(3)
              d(4) = L(4)
              d(5) = 0.5_rp*(L(5)+L(1))

              RHS(i,j,k,1) = RHS(i,j,k,1) -     d(1)
              RHS(i,j,k,2) = RHS(i,j,k,2) - (u_*d(1)+r_*d(2))
              RHS(i,j,k,3) = RHS(i,j,k,3) - (v_*d(1)+r_*d(3))
              RHS(i,j,k,4) = RHS(i,j,k,4) - (w_*d(1)+r_*d(4))
              RHS(i,j,k,5) = RHS(i,j,k,5) - (ek*d(1) + d(5)/(gamma0-1.0_rp) &
               + r_*u_*d(2) + r_*v_*d(3) + r_*w_*d(4))

          enddo
        enddo

        return
end subroutine CharacteristicRHSX


subroutine CharacteristicRHSy(RHS,beta)
        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: RHS
        real(rp)                                 , intent(in)    :: beta

        ! local declarations
        real(rp), dimension(5)    :: lambda, L, d
        real(rp), dimension(-2:0) :: cl
        real(rp) :: r_, ir, u_, v_, w_, p_, T_, c, ek
        real(rp) :: r_js, irjs, u_js, v_js, w_js, p_js
        real(rp) :: r_y, p_y, u_y, v_y, w_y, bw_s,idy
        integer  :: i,j,k,s,js

        cl(-2) =  0.5_rp
        cl(-1) = -2.0_rp
        cl( 0) =  1.5_rp

        do k = sz,ez
           do i = sx,ex

              j = ey
              RHS(i,j,k,:) = beta * RHS(i,j,k,:)

              r_ = phi(i,j,k,1)
              ir = 1.0_rp/r_
              u_ = phi(i,j,k,2)*ir
              v_ = phi(i,j,k,3)*ir
              w_ = phi(i,j,k,4)*ir
              ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
              p_ = cv_i*(phi(i,j,k,5) - r_*ek)
              T_ = p_*ir
              c  = sqrt(gamma0*T_)

              idy= ystep_i(j)

              r_y = 0.0_rp
              p_y = 0.0_rp
              u_y = 0.0_rp
              v_y = 0.0_rp
              w_y = 0.0_rp
              do s = -2,0
                 js   = j+s
                 bw_s = cl(s)

                 r_js = phi(i,js,k,1)
                 irjs = 1.0_rp/r_js
                 u_js = phi(i,js,k,2)*irjs
                 v_js = phi(i,js,k,3)*irjs
                 w_js = phi(i,js,k,4)*irjs
                 p_js = cv_i*(phi(i,js,k,5) &
                        - 0.5_rp*r_js*(u_js*u_js + v_js*v_js + w_js*w_js))

                 r_y = r_y + bw_s * r_js
                 p_y = p_y + bw_s * p_js
                 u_y = u_y + bw_s * u_js
                 v_y = v_y + bw_s * v_js
                 w_y = w_y + bw_s * w_js
              enddo
              r_y = r_y * idy
              p_y = p_y * idy
              u_y = u_y * idy
              v_y = v_y * idy
              w_y = w_y * idy

              lambda(1) = 0.5_rp * (v_ - c + abs(v_ - c))
              lambda(2) = 0.5_rp * (v_     + abs(v_    ))
              lambda(3) = 0.5_rp * (v_     + abs(v_    ))
              lambda(4) = 0.5_rp * (v_     + abs(v_    ))
              lambda(5) = 0.5_rp * (v_ + c + abs(v_ + c))

              L(1) = lambda(1) * (p_y - r_*c*v_y)
              if(v_/c < 1.0_rp .and. L1_wave == 'poinsot_model') then
                 L(1) = 0.27_rp*(1.0_rp - (v_/c)**2) * c/Ly * (p_ - 1.0_rp)
              endif
              L(2) = lambda(2) * (c**2*r_y - p_y)
              L(3) = lambda(3) * u_y
              L(4) = lambda(4) * w_y
              L(5) = lambda(5) * (p_y + r_*c*v_y)

              d(1) = 1.0_rp/c**2 * (L(2) + 0.5_rp*(L(5) + L(1)))
              d(2) = L(5) + L(1)
              d(3) = L(3)
              d(4) = 1.0_rp/(2*r_*c)*(L(5) - L(1))
              d(5) = L(4)

              RHS(i,j,k,1) = RHS(i,j,k,1) -     d(1)
              RHS(i,j,k,2) = RHS(i,j,k,2) - (u_*d(1) + r_*d(3))
              RHS(i,j,k,3) = RHS(i,j,k,3) - (v_*d(1) + r_*d(4))
              RHS(i,j,k,4) = RHS(i,j,k,4) - (w_*d(1) + r_*d(5))
              RHS(i,j,k,5) = RHS(i,j,k,5) - (ek*d(1) + d(2)/(gamma0-1.0_rp) &
                             + r_*u_*d(3) + r_*v_*d(4) + r_*w_*d(5))

           enddo
        enddo


        return
end subroutine CharacteristicRHSy




subroutine BCRelax_X(RHS)

        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: RHS

        ! local declarations
        real(rp), parameter      :: kk  =  gamma0-1.0_rp
        real(rp), parameter      :: ik  =  1.0_rp/kk
        real(rp), parameter      :: cl2 =  0.5_rp
        real(rp), parameter      :: cl1 = -2.0_rp
        real(rp), parameter      :: cl0 =  1.5_rp

        real(rp), dimension(5,5) :: PP, iP, LL, RR
        real(rp) :: uxi1,uxi_char1, uxo1, uxo_char1
        real(rp) :: uxi2,uxi_char2, uxo2, uxo_char2
        real(rp) :: uxi3,uxi_char3, uxo3, uxo_char3
        real(rp) :: uxi4,uxi_char4, uxo4, uxo_char4
        real(rp) :: uxi5,uxi_char5, uxo5, uxo_char5

        real(rp) :: idx, r0, u0, v0, w0, p0, alpha, FU
        real(rp) :: q1, q2, q3, q4, q5
        real(rp) :: e1, e2, e3, e4, e5
        real(rp) :: r_, ir, u_, v_, w_, ek, p_, T_, c, cc
        integer  :: i,j,k, mm, m

        call StartProfRange("BCRelax_X")

        r0 = 1.0_rp
        u0 = sqrt(gamma0)*Mach
        v0 = 0.0_rp
        w0 = 0.0_rp
        p0 = 1.0_rp

        q1 = r0
        q2 = r0*u0
        q3 = r0*v0
        q4 = r0*w0
        q5 = p0/(gamma0-1.0_rp) + 0.5_rp*r0*(u0**2 + v0**2 + w0**2)

        alpha = 0.0_rp
        if(L1_wave=='poinsot_model') alpha = 1.0_rp
        if(ic     =='smooth_body'  ) alpha = 0.1_rp

        i = ex
        idx = xstep_i(i)

        !$acc parallel default(present) &
        !$acc private(PP,iP,LL,RR)
        !$acc loop gang, vector collapse(2)
        do k    = sz,ez
           do j = sy,ey

              r_ = phi(i,j,k,1)
              ir = 1.0_rp/r_
              u_ = phi(i,j,k,2)*ir
              v_ = phi(i,j,k,3)*ir
              w_ = phi(i,j,k,4)*ir
              ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
              p_ = (gamma0-1._rp)*(phi(i,j,k,5)-r_*ek)
              T_ = p_*ir
              c  = sqrt(gamma0*T_)
              cc = c*c

              ! P matrix
              PP(1,1) = 1.0_rp
              PP(1,2) = 0.0_rp
              PP(1,3) = 0.0_rp
              PP(1,4) = 0.0_rp
              PP(1,5) = 0.0_rp
              !
              PP(2,1) = u_
              PP(2,2) = r_
              PP(2,3) = 0.0_rp
              PP(2,4) = 0.0_rp
              PP(2,5) = 0.0_rp
              !
              PP(3,1) = v_
              PP(3,2) = 0.0_rp
              PP(3,3) = r_
              PP(3,4) = 0.0_rp
              PP(3,5) = 0.0_rp
              !
              PP(4,1) = w_
              PP(4,2) = 0.0_rp
              PP(4,3) = 0.0_rp
              PP(4,4) = r_
              PP(4,5) = 0.0_rp
              !
              PP(5,1) = ek
              PP(5,2) = r_*u_
              PP(5,3) = r_*v_
              PP(5,4) = r_*w_
              PP(5,5) = ik

              ! inverse of P
              iP(1,1) = 1.0_rp
              iP(1,2) = 0.0_rp
              iP(1,3) = 0.0_rp
              iP(1,4) = 0.0_rp
              iP(1,5) = 0.0_rp
              !
              iP(2,1) = -u_*ir
              iP(2,2) = ir
              iP(2,3) = 0.0_rp
              iP(2,4) = 0.0_rp
              iP(2,5) = 0.0_rp
              !
              iP(3,1) = -v_*ir
              iP(3,2) = 0.0_rp
              iP(3,3) = ir
              iP(3,4) = 0.0_rp
              iP(3,5) = 0.0_rp
              !
              iP(4,1) = -w_*ir
              iP(4,2) = 0.0_rp
              iP(4,3) = 0.0_rp
              iP(4,4) = ir
              iP(4,5) = 0.0_rp
              !
              iP(5,1) = ek*kk
              iP(5,2) =-kk*u_
              iP(5,3) =-kk*v_
              iP(5,4) =-kk*w_
              iP(5,5) = kk

              ! left eigenvectors
              LL(1,1) = 0.0_rp
              LL(1,2) = -r_*c
              LL(1,3) = 0.0_rp
              LL(1,4) = 0.0_rp
              LL(1,5) = 1.0_rp
              !
              LL(2,1) = cc
              LL(2,2) = 0.0_rp
              LL(2,3) = 0.0_rp
              LL(2,4) = 0.0_rp
              LL(2,5) =-1.0_rp
              !
              LL(3,1) = 0.0_rp
              LL(3,2) = 0.0_rp
              LL(3,3) = 1.0_rp
              LL(3,4) = 0.0_rp
              LL(3,5) = 0.0_rp
              !
              LL(4,1) = 0.0_rp
              LL(4,2) = 0.0_rp
              LL(4,3) = 0.0_rp
              LL(4,4) = 1.0_rp
              LL(4,5) = 0.0_rp
              !
              LL(5,1) = 0.0_rp
              LL(5,2) = r_*c
              LL(5,3) = 0.0_rp
              LL(5,4) = 0.0_rp
              LL(5,5) = 1.0_rp

              !right eigenvector
              RR(1,1) = 0.5_rp/cc
              RR(1,2) = 1.0_rp/cc
              RR(1,3) = 0.0_rp
              RR(1,4) = 0.0_rp
              RR(1,5) = 0.5_rp/cc
              !
              RR(2,1) = -0.5_rp*ir/c
              RR(2,2) = 0.0_rp
              RR(2,3) = 0.0_rp
              RR(2,4) = 0.0_rp
              RR(2,5) = 0.5_rp*ir/c
              !
              RR(3,1) = 0.0_rp
              RR(3,2) = 0.0_rp
              RR(3,3) = 1.0_rp
              RR(3,4) = 0.0_rp
              RR(3,5) = 0.0_rp
              !
              RR(4,1) = 0.0_rp
              RR(4,2) = 0.0_rp
              RR(4,3) = 0.0_rp
              RR(4,4) = 1.0_rp
              RR(4,5) = 0.0_rp
              !
              RR(5,1) = 0.5_rp
              RR(5,2) = 0.0_rp
              RR(5,3) = 0.0_rp
              RR(5,4) = 0.0_rp
              RR(5,5) = 0.5_rp

              ! eigenvalues
              e1 = u_ - c
              e2 = u_
              e3 = u_
              e4 = u_
              e5 = u_ + c
              !
              ! compute derivative on physical space
              !
              ! inner
              uxi1 = (cl2*phi(i-2,j,k,1)+cl1*phi(i-1,j,k,1)+cl0*phi(i,j,k,1))*idx
              uxi2 = (cl2*phi(i-2,j,k,2)+cl1*phi(i-1,j,k,2)+cl0*phi(i,j,k,2))*idx
              uxi3 = (cl2*phi(i-2,j,k,3)+cl1*phi(i-1,j,k,3)+cl0*phi(i,j,k,3))*idx
              uxi4 = (cl2*phi(i-2,j,k,4)+cl1*phi(i-1,j,k,4)+cl0*phi(i,j,k,4))*idx
              uxi5 = (cl2*phi(i-2,j,k,5)+cl1*phi(i-1,j,k,5)+cl0*phi(i,j,k,5))*idx
              ! outer
              uxo1 = -(phi(i,j,k,1)-q1)*idx
              uxo2 = -(phi(i,j,k,2)-q2)*idx
              uxo3 = -(phi(i,j,k,3)-q3)*idx
              uxo4 = -(phi(i,j,k,4)-q4)*idx
              uxo5 = -(phi(i,j,k,5)-q5)*idx
              !
              ! derivative of characteristic variables
              !
              uxi_char1 = 0.0_rp
              uxi_char2 = 0.0_rp
              uxi_char3 = 0.0_rp
              uxi_char4 = 0.0_rp
              uxi_char5 = 0.0_rp

              uxo_char1 = 0.0_rp
              uxo_char2 = 0.0_rp
              uxo_char3 = 0.0_rp
              uxo_char4 = 0.0_rp
              uxo_char5 = 0.0_rp
              !$acc loop seq
              do mm=1,5
                 uxi_char1=uxi_char1+LL(1,mm)*iP(mm,1)*uxi1
                 uxi_char1=uxi_char1+LL(1,mm)*iP(mm,2)*uxi2
                 uxi_char1=uxi_char1+LL(1,mm)*iP(mm,3)*uxi3
                 uxi_char1=uxi_char1+LL(1,mm)*iP(mm,4)*uxi4
                 uxi_char1=uxi_char1+LL(1,mm)*iP(mm,5)*uxi5
                 !
                 uxi_char2=uxi_char2+LL(2,mm)*iP(mm,1)*uxi1
                 uxi_char2=uxi_char2+LL(2,mm)*iP(mm,2)*uxi2
                 uxi_char2=uxi_char2+LL(2,mm)*iP(mm,3)*uxi3
                 uxi_char2=uxi_char2+LL(2,mm)*iP(mm,4)*uxi4
                 uxi_char2=uxi_char2+LL(2,mm)*iP(mm,5)*uxi5
                 !
                 uxi_char3=uxi_char3+LL(3,mm)*iP(mm,1)*uxi1
                 uxi_char3=uxi_char3+LL(3,mm)*iP(mm,2)*uxi2
                 uxi_char3=uxi_char3+LL(3,mm)*iP(mm,3)*uxi3
                 uxi_char3=uxi_char3+LL(3,mm)*iP(mm,4)*uxi4
                 uxi_char3=uxi_char3+LL(3,mm)*iP(mm,5)*uxi5
                 !
                 uxi_char4=uxi_char4+LL(4,mm)*iP(mm,1)*uxi1
                 uxi_char4=uxi_char4+LL(4,mm)*iP(mm,2)*uxi2
                 uxi_char4=uxi_char4+LL(4,mm)*iP(mm,3)*uxi3
                 uxi_char4=uxi_char4+LL(4,mm)*iP(mm,4)*uxi4
                 uxi_char4=uxi_char4+LL(4,mm)*iP(mm,5)*uxi5
                 !
                 uxi_char5=uxi_char5+LL(5,mm)*iP(mm,1)*uxi1
                 uxi_char5=uxi_char5+LL(5,mm)*iP(mm,2)*uxi2
                 uxi_char5=uxi_char5+LL(5,mm)*iP(mm,3)*uxi3
                 uxi_char5=uxi_char5+LL(5,mm)*iP(mm,4)*uxi4
                 uxi_char5=uxi_char5+LL(5,mm)*iP(mm,5)*uxi5
                 !!!!!!
                 uxo_char1=uxo_char1+LL(1,mm)*iP(mm,1)*uxo1
                 uxo_char1=uxo_char1+LL(1,mm)*iP(mm,2)*uxo2
                 uxo_char1=uxo_char1+LL(1,mm)*iP(mm,3)*uxo3
                 uxo_char1=uxo_char1+LL(1,mm)*iP(mm,4)*uxo4
                 uxo_char1=uxo_char1+LL(1,mm)*iP(mm,5)*uxo5
                 !
                 uxo_char2=uxo_char2+LL(2,mm)*iP(mm,1)*uxo1
                 uxo_char2=uxo_char2+LL(2,mm)*iP(mm,2)*uxo2
                 uxo_char2=uxo_char2+LL(2,mm)*iP(mm,3)*uxo3
                 uxo_char2=uxo_char2+LL(2,mm)*iP(mm,4)*uxo4
                 uxo_char2=uxo_char2+LL(2,mm)*iP(mm,5)*uxo5
                 !
                 uxo_char3=uxo_char3+LL(3,mm)*iP(mm,1)*uxo1
                 uxo_char3=uxo_char3+LL(3,mm)*iP(mm,2)*uxo2
                 uxo_char3=uxo_char3+LL(3,mm)*iP(mm,3)*uxo3
                 uxo_char3=uxo_char3+LL(3,mm)*iP(mm,4)*uxo4
                 uxo_char3=uxo_char3+LL(3,mm)*iP(mm,5)*uxo5
                 !
                 uxo_char4=uxo_char4+LL(4,mm)*iP(mm,1)*uxo1
                 uxo_char4=uxo_char4+LL(4,mm)*iP(mm,2)*uxo2
                 uxo_char4=uxo_char4+LL(4,mm)*iP(mm,3)*uxo3
                 uxo_char4=uxo_char4+LL(4,mm)*iP(mm,4)*uxo4
                 uxo_char4=uxo_char4+LL(4,mm)*iP(mm,5)*uxo5
                 !
                 uxo_char5=uxo_char5+LL(5,mm)*iP(mm,1)*uxo1
                 uxo_char5=uxo_char5+LL(5,mm)*iP(mm,2)*uxo2
                 uxo_char5=uxo_char5+LL(5,mm)*iP(mm,3)*uxo3
                 uxo_char5=uxo_char5+LL(5,mm)*iP(mm,4)*uxo4
                 uxo_char5=uxo_char5+LL(5,mm)*iP(mm,5)*uxo5
              enddo
              !
              ! enforce LODI relations
              !
              if(e1 < 0.0_rp) uxi_char1 = alpha*uxo_char1
              if(e2 < 0.0_rp) uxi_char2 = alpha*uxo_char2
              if(e3 < 0.0_rp) uxi_char3 = alpha*uxo_char3
              if(e4 < 0.0_rp) uxi_char4 = alpha*uxo_char4
              if(e5 < 0.0_rp) uxi_char5 = alpha*uxo_char5

              uxi_char1 = e1 * uxi_char1
              uxi_char2 = e2 * uxi_char2
              uxi_char3 = e3 * uxi_char3
              uxi_char4 = e4 * uxi_char4
              uxi_char5 = e5 * uxi_char5
              !
              ! computing RHS
              !
              !$acc loop seq
              do m = 1,5
                 FU = 0.0_rp

                 FU = FU + PP(m,1)*RR(1,1)*uxi_char1
                 FU = FU + PP(m,1)*RR(1,2)*uxi_char2
                 FU = FU + PP(m,1)*RR(1,3)*uxi_char3
                 FU = FU + PP(m,1)*RR(1,4)*uxi_char4
                 FU = FU + PP(m,1)*RR(1,5)*uxi_char5

                 FU = FU + PP(m,2)*RR(2,1)*uxi_char1
                 FU = FU + PP(m,2)*RR(2,2)*uxi_char2
                 FU = FU + PP(m,2)*RR(2,3)*uxi_char3
                 FU = FU + PP(m,2)*RR(2,4)*uxi_char4
                 FU = FU + PP(m,2)*RR(2,5)*uxi_char5

                 FU = FU + PP(m,3)*RR(3,1)*uxi_char1
                 FU = FU + PP(m,3)*RR(3,2)*uxi_char2
                 FU = FU + PP(m,3)*RR(3,3)*uxi_char3
                 FU = FU + PP(m,3)*RR(3,4)*uxi_char4
                 FU = FU + PP(m,3)*RR(3,5)*uxi_char5

                 FU = FU + PP(m,4)*RR(4,1)*uxi_char1
                 FU = FU + PP(m,4)*RR(4,2)*uxi_char2
                 FU = FU + PP(m,4)*RR(4,3)*uxi_char3
                 FU = FU + PP(m,4)*RR(4,4)*uxi_char4
                 FU = FU + PP(m,4)*RR(4,5)*uxi_char5

                 FU = FU + PP(m,5)*RR(5,1)*uxi_char1
                 FU = FU + PP(m,5)*RR(5,2)*uxi_char2
                 FU = FU + PP(m,5)*RR(5,3)*uxi_char3
                 FU = FU + PP(m,5)*RR(5,4)*uxi_char4
                 FU = FU + PP(m,5)*RR(5,5)*uxi_char5

                 RHS(i,j,k,m) = RHS(i,j,k,m) - FU
              enddo

           enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine BCRelax_X




subroutine BCRelax_InflowX(RHS)

        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: RHS

        ! local declarations
        real(rp), parameter      :: kk = gamma0-1.0_rp
        real(rp), parameter      :: ik = 1.0_rp/kk
        real(rp), dimension(5,5) :: PP, iP, LL, RR
        real(rp), dimension(5)   :: uxi,uxi_char, uxo, uxo_char

        real(rp) :: idx, cl0, cl1, cl2, FU
        real(rp) :: q1, q2, q3, q4, q5
        real(rp) :: e1, e2, e3, e4, e5
        real(rp) :: r_, ir, u_, v_, w_, ek, p_, T_, c, cc
        integer  :: i,j,k,l, lll, mm, m, n

        call StartProfRange("BCRelax_InflowX")

        cl0 = -1.5_rp
        cl1 =  2.0_rp
        cl2 = -0.5_rp

        i = sx
        idx = xstep_i(i)

        !$acc parallel default(present) &
        !$acc private(PP,iP,LL,RR,uxi,uxo,uxi_char,uxo_char)
        !$acc loop gang, vector collapse(2)
        do k = sz,ez
           do j = sy,ey

              ! reference inflow quantities from the node i-1
              q1 = phi(i-1,j,k,1)
              q2 = phi(i-1,j,k,2)
              q3 = phi(i-1,j,k,3)
              q4 = phi(i-1,j,k,4)
              q5 = phi(i-1,j,k,5)

              r_ = phi(i,j,k,1)
              ir = 1.0_rp/r_
              u_ = phi(i,j,k,2)*ir
              v_ = phi(i,j,k,3)*ir
              w_ = phi(i,j,k,4)*ir
              ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
              p_ = (gamma0-1._rp)*(phi(i,j,k,5)-r_*ek)
              T_ = p_*ir
              c  = sqrt(gamma0*T_)
              cc = c*c

              ! P matrix
              PP(1,1) = 1.0_rp
              PP(1,2) = 0.0_rp
              PP(1,3) = 0.0_rp
              PP(1,4) = 0.0_rp
              PP(1,5) = 0.0_rp
              !
              PP(2,1) = u_
              PP(2,2) = r_
              PP(2,3) = 0.0_rp
              PP(2,4) = 0.0_rp
              PP(2,5) = 0.0_rp
              !
              PP(3,1) = v_
              PP(3,2) = 0.0_rp
              PP(3,3) = r_
              PP(3,4) = 0.0_rp
              PP(3,5) = 0.0_rp
              !
              PP(4,1) = w_
              PP(4,2) = 0.0_rp
              PP(4,3) = 0.0_rp
              PP(4,4) = r_
              PP(4,5) = 0.0_rp
              !
              PP(5,1) = ek
              PP(5,2) = r_*u_
              PP(5,3) = r_*v_
              PP(5,4) = r_*w_
              PP(5,5) = ik

              ! inverse of P
              iP(1,1) = 1.0_rp
              iP(1,2) = 0.0_rp
              iP(1,3) = 0.0_rp
              iP(1,4) = 0.0_rp
              iP(1,5) = 0.0_rp
              !
              iP(2,1) = -u_*ir
              iP(2,2) = ir
              iP(2,3) = 0.0_rp
              iP(2,4) = 0.0_rp
              iP(2,5) = 0.0_rp
              !
              iP(3,1) = -v_*ir
              iP(3,2) = 0.0_rp
              iP(3,3) = ir
              iP(3,4) = 0.0_rp
              iP(3,5) = 0.0_rp
              !
              iP(4,1) = -w_*ir
              iP(4,2) = 0.0_rp
              iP(4,3) = 0.0_rp
              iP(4,4) = ir
              iP(4,5) = 0.0_rp
              !
              iP(5,1) = ek*kk
              iP(5,2) = -kk*u_
              iP(5,3) = -kk*v_
              iP(5,4) = -kk*w_
              iP(5,5) = kk

              ! left eigenvectors
              LL(1,1) = 0.0_rp
              LL(1,2) = -r_*c
              LL(1,3) = 0.0_rp
              LL(1,4) = 0.0_rp
              LL(1,5) = 1.0_rp
              !
              LL(2,1) = cc
              LL(2,2) = 0.0_rp
              LL(2,3) = 0.0_rp
              LL(2,4) = 0.0_rp
              LL(2,5) = -1.0_rp
              !
              LL(3,1) = 0.0_rp
              LL(3,2) = 0.0_rp
              LL(3,3) = 1.0_rp
              LL(3,4) = 0.0_rp
              LL(3,5) = 0.0_rp
              !
              LL(4,1) = 0.0_rp
              LL(4,2) = 0.0_rp
              LL(4,3) = 0.0_rp
              LL(4,4) = 1.0_rp
              LL(4,5) = 0.0_rp
              !
              LL(5,1) = 0.0_rp
              LL(5,2) = r_*c
              LL(5,3) = 0.0_rp
              Ll(5,4) = 0.0_rp
              LL(5,5) = 1.0_rp

              !right eigenvector
              RR(1,1) = 0.5_rp/cc
              RR(1,2) = 1.0_rp/cc
              RR(1,3) = 0.0_rp
              RR(1,4) = 0.0_rp
              RR(1,5) = 0.5_rp/cc
              !
              RR(2,1) = -0.5_rp*ir/c
              RR(2,2) = 0.0_rp
              RR(2,3) = 0.0_rp
              RR(2,4) = 0.0_rp
              RR(2,5) = 0.5_rp*ir/c
              !
              RR(3,1) = 0.0_rp
              RR(3,2) = 0.0_rp
              RR(3,3) = 1.0_rp
              RR(3,4) = 0.0_rp
              RR(3,5) = 0.0_rp
              !
              RR(4,1) = 0.0_rp
              RR(4,2) = 0.0_rp
              RR(4,3) = 0.0_rp
              RR(4,4) = 1.0_rp
              RR(4,5) = 0.0_rp
              !
              RR(5,1) = 0.5_rp
              RR(5,2) = 0.0_rp
              RR(5,3) = 0.0_rp
              RR(5,4) = 0.0_rp
              RR(5,5) = 0.5_rp

              ! eigenvalues
              e1 = u_ - c
              e2 = u_
              e3 = u_
              e4 = u_
              e5 = u_ + c
              !
              ! compute derivative on physical space
              !
              ! inner
              uxi(1) = (cl0*phi(i,j,k,1)+cl1*phi(i+1,j,k,1)+cl2*phi(i+2,j,k,1))*idx
              uxi(2) = (cl0*phi(i,j,k,2)+cl1*phi(i+1,j,k,2)+cl2*phi(i+2,j,k,2))*idx
              uxi(3) = (cl0*phi(i,j,k,3)+cl1*phi(i+1,j,k,3)+cl2*phi(i+2,j,k,3))*idx
              uxi(4) = (cl0*phi(i,j,k,4)+cl1*phi(i+1,j,k,4)+cl2*phi(i+2,j,k,4))*idx
              uxi(5) = (cl0*phi(i,j,k,5)+cl1*phi(i+1,j,k,5)+cl2*phi(i+2,j,k,5))*idx
              ! outer
              uxo(1) = (phi(i,j,k,1)-q1)*idx
              uxo(2) = (phi(i,j,k,2)-q2)*idx
              uxo(3) = (phi(i,j,k,3)-q3)*idx
              uxo(4) = (phi(i,j,k,4)-q4)*idx
              uxo(5) = (phi(i,j,k,5)-q5)*idx
              !
              ! derivative of characteristic variables
              !
              !$acc loop seq
              do m=1,5
                 uxi_char(m) = 0.0_rp
                 uxo_char(m) = 0.0_rp
                 !$acc loop seq
                 do mm=1,5
                    !$acc loop seq
                    do lll=1,5
                       uxi_char(m)=uxi_char(m)+LL(m,mm)*iP(mm,lll)*uxi(lll)
                       uxo_char(m)=uxo_char(m)+LL(m,mm)*iP(mm,lll)*uxo(lll)
                    enddo
                 enddo
              enddo
              !
              ! enforce LODI relations
              !
              if(e1 > 0.0_rp) uxi_char(1) = uxo_char(1)
              if(e2 > 0.0_rp) uxi_char(2) = uxo_char(2)
              if(e3 > 0.0_rp) uxi_char(3) = uxo_char(3)
              if(e4 > 0.0_rp) uxi_char(4) = uxo_char(4)
              if(e5 > 0.0_rp) uxi_char(5) = uxo_char(5)

              uxi_char(1) = e1 * uxi_char(1)
              uxi_char(2) = e2 * uxi_char(2)
              uxi_char(3) = e3 * uxi_char(3)
              uxi_char(4) = e4 * uxi_char(4)
              uxi_char(5) = e5 * uxi_char(5)
              !
              ! computing RHS
              !
              !$acc loop seq
              do m = 1,5
                 FU = 0.0_rp
                 !$acc loop seq
                 do l = 1,5
                    !$acc loop seq
                    do n = 1,5
                       FU = FU + PP(m,n)*RR(n,l)*uxi_char(l)
                    enddo
                 enddo
                 RHS(i,j,k,m) = RHS(i,j,k,m) - FU
              enddo

           enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine BCRelax_InflowX






















subroutine BCRelax_Y(RHS)

        implicit none
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: RHS

        ! local declarations
        real(rp), parameter      :: kk = gamma0 - 1.0_rp
        real(rp), parameter      :: ik = 1.0_rp/kk
        real(rp), dimension(5,5) :: PP, iP, LL, RR
        real(rp)                 :: uyi1, uyi_char1, uyo_char1, uyo1
        real(rp)                 :: uyi2, uyi_char2, uyo_char2, uyo2
        real(rp)                 :: uyi3, uyi_char3, uyo_char3, uyo3
        real(rp)                 :: uyi4, uyi_char4, uyo_char4, uyo4
        real(rp)                 :: uyi5, uyi_char5, uyo_char5, uyo5

        real(rp) :: cl2, cl1, cl0, alpha, GU
        real(rp) :: q1, q2, q3, q4, q5
        real(rp) :: e1, e2, e3, e4, e5
        real(rp) :: r_, ir, u_, v_, w_, ek, p_, T_, c, cc
        real(rp) :: idy, r0, u0, v0, w0, p0
        integer  :: i,j,k,mm, m

        call StartProfRange("BCRelax_Y")

        cl2 =  0.5_rp
        cl1 = -2.0_rp
        cl0 =  1.5_rp

        r0 = 1.0_rp
        u0 = sqrt(gamma0)*Mach
        v0 = 0.0_rp
        w0 = 0.0_rp
        p0 = 1.0_rp

        q1 = r0
        q2 = r0*u0
        q3 = r0*v0
        q4 = r0*w0
        q5 = p0/(gamma0-1.0_rp) + 0.5_rp*r0*(u0**2 + v0**2 + w0**2)

        alpha = 0.0_rp
        if(Mach < 1.0_rp)     alpha = 1.0_rp
        if(ic=='smooth_body') alpha = 0.01_rp

        j = ey
        idy = ystep_i(j)

        !$acc parallel default(present) &
        !$acc private(PP,iP,LL,RR)
        !$acc loop gang, vector collapse(2)
        do k    = sz,ez
           do i = sx,ex

              r_ = phi(i,j,k,1)
              ir = 1.0_rp/r_
              u_ = phi(i,j,k,2)*ir
              v_ = phi(i,j,k,3)*ir
              w_ = phi(i,j,k,4)*ir
              ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
              p_ = (gamma0-1._rp)*(phi(i,j,k,5)-r_*ek)
              T_ = p_*ir
              c  = sqrt(gamma0*T_)
              cc = c*c

              ! P matrix
              PP(1,1) = 1.0_rp
              PP(1,2) = 0.0_rp
              PP(1,3) = 0.0_rp
              PP(1,4) = 0.0_rp
              PP(1,5) = 0.0_rp
              !
              PP(2,1) = u_
              PP(2,2) = r_
              PP(2,3) = 0.0_rp
              PP(2,4) = 0.0_rp
              PP(2,5) = 0.0_rp
              !
              PP(3,1) = v_
              PP(3,2) = 0.0_rp
              PP(3,3) = r_
              PP(3,4) = 0.0_rp
              PP(3,5) = 0.0_rp
              !
              PP(4,1) = w_
              PP(4,2) = 0.0_rp
              PP(4,3) = 0.0_rp
              PP(4,4) = r_
              PP(4,5) = 0.0_rp
              !
              PP(5,1) = ek
              PP(5,2) = r_*u_
              PP(5,3) = r_*v_
              PP(5,4) = r_*w_
              PP(5,5) = ik

              ! inverse of P
              iP(1,1) = 1.0_rp
              iP(1,2) = 0.0_rp
              iP(1,3) = 0.0_rp
              ip(1,4) = 0.0_rp
              iP(1,5) = 0.0_rp
              !
              iP(2,1) = -u_*ir
              iP(2,2) = ir
              iP(2,3) = 0.0_rp
              iP(2,4) = 0.0_rp
              iP(2,5) = 0.0_rp
              !
              iP(3,1) = -v_*ir
              iP(3,2) = 0.0_rp
              iP(3,3) = ir
              iP(3,4) = 0.0_rp
              iP(3,5) = 0.0_rp
              !
              iP(4,1) = -w_*ir
              iP(4,2) = 0.0_rp
              iP(4,3) = 0.0_rp
              iP(4,4) = ir
              iP(4,5) = 0.0_rp
              !
              iP(5,1) = ek*kk
              iP(5,2) = -kk*u_
              iP(5,3) = -kk*v_
              iP(5,4) = -kk*w_
              iP(5,5) = kk

              ! left eigenvectors
              LL(1,1) = 0.0_rp
              LL(1,2) = 0.0_rp
              LL(1,3) = -r_*c
              LL(1,4) = 0.0_rp
              LL(1,5) = 1.0_rp

              LL(2,1) = 0.0_rp
              LL(2,2) = 1.0_rp
              LL(2,3) = 0.0_rp
              LL(2,4) = 0.0_rp
              LL(2,5) = 0.0_rp

              LL(3,1) = cc
              LL(3,2) = 0.0_rp
              LL(3,3) = 0.0_rp
              LL(3,4) = 0.0_rp
              LL(3,5) = -1.0_rp

              LL(4,1) = 0.0_rp
              LL(4,2) = 0.0_rp
              LL(4,3) = 0.0_rp
              LL(4,4) = 1.0_rp
              LL(4,5) = 0.0_rp

              LL(5,1) = 0.0_rp
              LL(5,2) = 0.0_rp
              LL(5,3) = r_*c
              LL(5,4) = 0.0_rp
              LL(5,5) = 1.0_rp

              ! right eigenvectors
              RR(1,1) = 0.5_rp/cc
              RR(1,2) = 0.0_rp
              RR(1,3) = 1.0_rp/cc
              RR(1,4) = 0.0_rp
              RR(1,5) = 0.5_rp/cc
              !
              RR(2,1) = 0.0_rp
              RR(2,2) = 1.0_rp
              RR(2,3) = 0.0_rp
              RR(2,4) = 0.0_rp
              RR(2,5) = 0.0_rp
              !
              RR(3,1) = -0.5_rp*ir/c
              RR(3,2) = 0.0_rp
              RR(3,3) = 0.0_rp
              RR(3,4) = 0.0_rp
              RR(3,5) = 0.5_rp*ir/ c
              !
              RR(4,1) = 0.0_rp
              RR(4,2) = 0.0_rp
              RR(4,3) = 0.0_rp
              RR(4,4) = 1.0_rp
              RR(4,5) = 0.0_rp
              !
              RR(5,1) = 0.5_rp
              RR(5,2) = 0.0_rp
              RR(5,3) = 0.0_rp
              RR(5,4) = 0.0_rp
              RR(5,5) = 0.5_rp

              ! eigenvalues
              e1 = v_ - c
              e2 = v_
              e3 = v_
              e4 = v_
              e5 = v_ + c
              !
              ! compute derivative on physical space
              !
              ! inner
              uyi1 = (cl2*phi(i,j-2,k,1)+cl1*phi(i,j-1,k,1)+cl0*phi(i,j,k,1))*idy
              uyi2 = (cl2*phi(i,j-2,k,2)+cl1*phi(i,j-1,k,2)+cl0*phi(i,j,k,2))*idy
              uyi3 = (cl2*phi(i,j-2,k,3)+cl1*phi(i,j-1,k,3)+cl0*phi(i,j,k,3))*idy
              uyi4 = (cl2*phi(i,j-2,k,4)+cl1*phi(i,j-1,k,4)+cl0*phi(i,j,k,4))*idy
              uyi5 = (cl2*phi(i,j-2,k,5)+cl1*phi(i,j-1,k,5)+cl0*phi(i,j,k,5))*idy
              ! outer
              uyo1 = -(phi(i,j,k,1)-q1)*idy
              uyo2 = -(phi(i,j,k,2)-q2)*idy
              uyo3 = -(phi(i,j,k,3)-q3)*idy
              uyo4 = -(phi(i,j,k,4)-q4)*idy
              uyo5 = -(phi(i,j,k,5)-q5)*idy
              !
              ! derivative of characteristic variables
              !
              uyi_char1 = 0.0_rp
              uyi_char2 = 0.0_rp
              uyi_char3 = 0.0_rp
              uyi_char4 = 0.0_rp
              uyi_char5 = 0.0_rp

              uyo_char1 = 0.0_rp
              uyo_char2 = 0.0_rp
              uyo_char3 = 0.0_rp
              uyo_char4 = 0.0_rp
              uyo_char5 = 0.0_rp
              !$acc loop seq
              do mm = 1,5
                 uyi_char1=uyi_char1+LL(1,mm)*iP(mm,1)*uyi1
                 uyi_char2=uyi_char2+LL(2,mm)*iP(mm,1)*uyi1
                 uyi_char3=uyi_char3+LL(3,mm)*iP(mm,1)*uyi1
                 uyi_char4=uyi_char4+LL(4,mm)*iP(mm,1)*uyi1
                 uyi_char5=uyi_char5+LL(5,mm)*iP(mm,1)*uyi1
                 !
                 uyi_char1=uyi_char1+LL(1,mm)*iP(mm,2)*uyi2
                 uyi_char2=uyi_char2+LL(2,mm)*iP(mm,2)*uyi2
                 uyi_char3=uyi_char3+LL(3,mm)*iP(mm,2)*uyi2
                 uyi_char4=uyi_char4+LL(4,mm)*iP(mm,2)*uyi2
                 uyi_char5=uyi_char5+LL(5,mm)*iP(mm,2)*uyi2
                 !
                 uyi_char1=uyi_char1+LL(1,mm)*iP(mm,3)*uyi3
                 uyi_char2=uyi_char2+LL(2,mm)*iP(mm,3)*uyi3
                 uyi_char3=uyi_char3+LL(3,mm)*iP(mm,3)*uyi3
                 uyi_char4=uyi_char4+LL(4,mm)*iP(mm,3)*uyi3
                 uyi_char5=uyi_char5+LL(5,mm)*iP(mm,3)*uyi3
                 !
                 uyi_char1=uyi_char1+LL(1,mm)*iP(mm,4)*uyi4
                 uyi_char2=uyi_char2+LL(2,mm)*iP(mm,4)*uyi4
                 uyi_char3=uyi_char3+LL(3,mm)*iP(mm,4)*uyi4
                 uyi_char4=uyi_char4+LL(4,mm)*iP(mm,4)*uyi4
                 uyi_char5=uyi_char5+LL(5,mm)*iP(mm,4)*uyi4
                 !
                 uyi_char1=uyi_char1+LL(1,mm)*iP(mm,5)*uyi5
                 uyi_char2=uyi_char2+LL(2,mm)*iP(mm,5)*uyi5
                 uyi_char3=uyi_char3+LL(3,mm)*iP(mm,5)*uyi5
                 uyi_char4=uyi_char4+LL(4,mm)*iP(mm,5)*uyi5
                 uyi_char5=uyi_char5+LL(5,mm)*iP(mm,5)*uyi5
                 !!!!!
                 uyo_char1=uyo_char1+LL(1,mm)*iP(mm,1)*uyo1
                 uyo_char2=uyo_char2+LL(2,mm)*iP(mm,1)*uyo1
                 uyo_char3=uyo_char3+LL(3,mm)*iP(mm,1)*uyo1
                 uyo_char4=uyo_char4+LL(4,mm)*iP(mm,1)*uyo1
                 uyo_char5=uyo_char5+LL(5,mm)*iP(mm,1)*uyo1
                 !
                 uyo_char1=uyo_char1+LL(1,mm)*iP(mm,2)*uyo2
                 uyo_char2=uyo_char2+LL(2,mm)*iP(mm,2)*uyo2
                 uyo_char3=uyo_char3+LL(3,mm)*iP(mm,2)*uyo2
                 uyo_char4=uyo_char4+LL(4,mm)*iP(mm,2)*uyo2
                 uyo_char5=uyo_char5+LL(5,mm)*iP(mm,2)*uyo2
                 !
                 uyo_char1=uyo_char1+LL(1,mm)*iP(mm,3)*uyo3
                 uyo_char2=uyo_char2+LL(2,mm)*iP(mm,3)*uyo3
                 uyo_char3=uyo_char3+LL(3,mm)*iP(mm,3)*uyo3
                 uyo_char4=uyo_char4+LL(4,mm)*iP(mm,3)*uyo3
                 uyo_char5=uyo_char5+LL(5,mm)*iP(mm,3)*uyo3
                 !
                 uyo_char1=uyo_char1+LL(1,mm)*iP(mm,4)*uyo4
                 uyo_char2=uyo_char2+LL(2,mm)*iP(mm,4)*uyo4
                 uyo_char3=uyo_char3+LL(3,mm)*iP(mm,4)*uyo4
                 uyo_char4=uyo_char4+LL(4,mm)*iP(mm,4)*uyo4
                 uyo_char5=uyo_char5+LL(5,mm)*iP(mm,4)*uyo4
                 !
                 uyo_char1=uyo_char1+LL(1,mm)*iP(mm,5)*uyo5
                 uyo_char2=uyo_char2+LL(2,mm)*iP(mm,5)*uyo5
                 uyo_char3=uyo_char3+LL(3,mm)*iP(mm,5)*uyo5
                 uyo_char4=uyo_char4+LL(4,mm)*iP(mm,5)*uyo5
                 uyo_char5=uyo_char5+LL(5,mm)*iP(mm,5)*uyo5
              enddo
              !
              ! enforce LODI relations
              !
              if(e1 < 0.0_rp) uyi_char1 = alpha*uyo_char1
              if(e2 < 0.0_rp) uyi_char2 = alpha*uyo_char2
              if(e3 < 0.0_rp) uyi_char3 = alpha*uyo_char3
              if(e4 < 0.0_rp) uyi_char4 = alpha*uyo_char4
              if(e5 < 0.0_rp) uyi_char5 = alpha*uyo_char5

              uyi_char1 = e1 * uyi_char1
              uyi_char2 = e2 * uyi_char2
              uyi_char3 = e3 * uyi_char3
              uyi_char4 = e4 * uyi_char4
              uyi_char5 = e5 * uyi_char5
              !
              ! computing RHS
              !
              !$acc loop seq
              do m = 1,5
                 GU = 0.0_rp

                 GU = GU + PP(m,1)*RR(1,1)*uyi_char1
                 GU = GU + PP(m,1)*RR(1,2)*uyi_char2
                 GU = GU + PP(m,1)*RR(1,3)*uyi_char3
                 GU = GU + PP(m,1)*RR(1,4)*uyi_char4
                 GU = GU + PP(m,1)*RR(1,5)*uyi_char5

                 GU = GU + PP(m,2)*RR(2,1)*uyi_char1
                 GU = GU + PP(m,2)*RR(2,2)*uyi_char2
                 GU = GU + PP(m,2)*RR(2,3)*uyi_char3
                 GU = GU + PP(m,2)*RR(2,4)*uyi_char4
                 GU = GU + PP(m,2)*RR(2,5)*uyi_char5

                 GU = GU + PP(m,3)*RR(3,1)*uyi_char1
                 GU = GU + PP(m,3)*RR(3,2)*uyi_char2
                 GU = GU + PP(m,3)*RR(3,3)*uyi_char3
                 GU = GU + PP(m,3)*RR(3,4)*uyi_char4
                 GU = GU + PP(m,3)*RR(3,5)*uyi_char5

                 GU = GU + PP(m,4)*RR(4,1)*uyi_char1
                 GU = GU + PP(m,4)*RR(4,2)*uyi_char2
                 GU = GU + PP(m,4)*RR(4,3)*uyi_char3
                 GU = GU + PP(m,4)*RR(4,4)*uyi_char4
                 GU = GU + PP(m,4)*RR(4,5)*uyi_char5

                 GU = GU + PP(m,5)*RR(5,1)*uyi_char1
                 GU = GU + PP(m,5)*RR(5,2)*uyi_char2
                 GU = GU + PP(m,5)*RR(5,3)*uyi_char3
                 GU = GU + PP(m,5)*RR(5,4)*uyi_char4
                 GU = GU + PP(m,5)*RR(5,5)*uyi_char5

                 RHS(i,j,k,m) = RHS(i,j,k,m) - GU
              enddo


           enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine BCRelax_Y

























end module advection_module

