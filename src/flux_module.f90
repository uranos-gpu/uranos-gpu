module flux_module

implicit none
private
public conv_flux_x, conv_flux_y, conv_flux_z, energy_preserving_tilde_op!, compute_pressure_array
contains

subroutine conv_flux_x(lb,ub,phi_arr_x,flx_arr_x,pri_arr_x)

        use parameters_module, only: rp, gamma0

        implicit none
        integer                     , intent(in   ) :: lb,ub
        real(rp), dimension(lb:ub,5), intent(in   ) :: phi_arr_x
        real(rp), dimension(lb:ub,5), intent(inout) :: flx_arr_x
        real(rp), dimension(lb:ub,6), intent(inout) :: pri_arr_x

        real(rp), dimension(5) :: phi_
        real(rp), parameter    :: gm1 = gamma0-1.0_rp
        real(rp), parameter    :: hgm = gamma0/gm1
        real(rp)               :: r_, ir, u_, v_, w_, p_, ek
        integer                :: i

        do i = lb, ub

           phi_(:)  = phi_arr_x(i,:)
        
           r_ = phi_(1)
           ir = 1.0_rp/r_
           u_ = phi_(2)*ir
           v_ = phi_(3)*ir
           w_ = phi_(4)*ir
           ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
           p_ = gm1 * (phi_(5) - r_*ek)

           flx_arr_x(i,1) =  phi_(2)
           flx_arr_x(i,2) =  phi_(2)   * u_ + p_
           flx_arr_x(i,3) =  phi_(2)   * v_
           flx_arr_x(i,4) =  phi_(2)   * w_
           flx_arr_x(i,5) = (phi_(5)+p_)*u_

           pri_arr_x(i,1) = p_
           pri_arr_x(i,2) = u_
           pri_arr_x(i,3) = v_
           pri_arr_x(i,4) = w_
           pri_arr_x(i,5) = hgm*p_*ir + ek
           pri_arr_x(i,6) = r_

        enddo

        return
end subroutine conv_flux_x




subroutine conv_flux_y(lb,ub,phi_arr_y,flx_arr_y,pri_arr_y)

        use parameters_module, only: rp, gamma0

        implicit none
        integer                     , intent(in   ) :: lb,ub
        real(rp), dimension(lb:ub,5), intent(in   ) :: phi_arr_y
        real(rp), dimension(lb:ub,5), intent(inout) :: flx_arr_y
        real(rp), dimension(lb:ub,6), intent(inout) :: pri_arr_y

        real(rp), dimension(5) :: phi_
        real(rp), parameter    :: gm1 = gamma0-1.0_rp
        real(rp), parameter    :: hgm = gamma0/gm1
        real(rp)               :: r_, ir, u_, v_, w_, p_, ek
        integer                :: j

        do j = lb,ub

           phi_(:)  = phi_arr_y(j,:)

           r_ = phi_(1)
           ir = 1.0_rp/r_
           u_ = phi_(2)*ir
           v_ = phi_(3)*ir
           w_ = phi_(4)*ir
           ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
           p_ = gm1 * (phi_(5) - r_*ek)

           flx_arr_y(j,1) =  phi_(3)
           flx_arr_y(j,2) =  phi_(3)    * u_
           flx_arr_y(j,3) =  phi_(3)    * v_ + p_
           flx_arr_y(j,4) =  phi_(3)    * w_
           flx_arr_y(j,5) = (phi_(5)+p_)* v_

           pri_arr_y(j,1) = p_
           pri_arr_y(j,2) = u_
           pri_arr_y(j,3) = v_
           pri_arr_y(j,4) = w_
           pri_arr_y(j,5) = hgm*p_*ir + ek
           pri_arr_y(j,6) = r_

        enddo

        return
end subroutine conv_flux_y



subroutine conv_flux_z(lb,ub,phi_arr_z,flx_arr_z,pri_arr_z)

        use parameters_module, only: rp, gamma0

        implicit none
        integer                     , intent(in   ) :: lb,ub
        real(rp), dimension(lb:ub,5), intent(inout) :: phi_arr_z
        real(rp), dimension(lb:ub,5), intent(inout) :: flx_arr_z
        real(rp), dimension(lb:ub,6), intent(inout) :: pri_arr_z

        real(rp), dimension(5) :: phi_
        real(rp), parameter    :: gm1 = gamma0-1.0_rp
        real(rp), parameter    :: hgm = gamma0/gm1
        real(rp)               :: r_, ir, u_, v_, w_, p_, ek
        integer                :: k

        do k = lb,ub

           phi_(:) = phi_arr_z(k,:)

           r_ = phi_(1)
           ir = 1.0_rp/r_
           u_ = phi_(2)*ir
           v_ = phi_(3)*ir
           w_ = phi_(4)*ir
           ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
           p_ = gm1 * (phi_(5) - r_*ek)

           flx_arr_z(k,1) =  phi_(4)
           flx_arr_z(k,2) =  phi_(4)    * u_
           flx_arr_z(k,3) =  phi_(4)    * v_ 
           flx_arr_z(k,4) =  phi_(4)    * w_ + p_
           flx_arr_z(k,5) = (phi_(5)+p_)* w_

           pri_arr_z(k,1) = p_
           pri_arr_z(k,2) = u_
           pri_arr_z(k,3) = v_
           pri_arr_z(k,4) = w_
           pri_arr_z(k,5) = hgm*p_*ir + ek
           pri_arr_z(k,6) = r_

        enddo

        return
end subroutine conv_flux_z


subroutine energy_preserving_tilde_op(dir,ltot,lb,ub,pri_arr,tilde_op)
! ---------------------------------------------------------------------------
!       
!       Coomputation of the tilde operator for energy preserving scheme.
!
!       REF: Pirozzoli, "Generalized conservative approximations of split 
!             convective derivative operators, JCP 2010"
!
!       IN:  dir         !< phisycal direction (1,2,3)
!            lb,ub       !< lower ad upper bound of the arrays 
!            phi_array   !< array of conservative variable
!            prs_array   !< pressure array
!
!       OUT: tilde_op    !> tilde operator
!
! ---------------------------------------------------------------------------
        use parameters_module, only: rp

        implicit none
        integer                          , intent(in)    :: lb,ub
        integer                          , intent(in)    :: ltot
        real(rp), dimension(lb:ub,6)     , intent(in)    :: pri_arr
        real(rp), dimension(ltot,lb:ub,5), intent(inout) :: tilde_op
        integer                          , intent(in)    :: dir

        real(rp), dimension(5) :: fi_k0, fi_k1
        real(rp), parameter    :: one8 = 1.0_rp/8.0_rp
        real(rp)               :: r0, r1, weight
        integer                :: k, l, kl, id

        id = dir+1
        do k = lb, ub-3
        
           fi_k0(1) = 1.0_rp
           fi_k0(2) = pri_arr(k,2)
           fi_k0(3) = pri_arr(k,3)
           fi_k0(4) = pri_arr(k,4)
           fi_k0(5) = pri_arr(k,5)
           r0       = pri_arr(k,6)
        
           do l = 1,ltot

              kl = k+l

              fi_k1(1) = 1.0_rp
              fi_k1(2) = pri_arr(kl,2)
              fi_k1(3) = pri_arr(kl,3)
              fi_k1(4) = pri_arr(kl,4)
              fi_k1(5) = pri_arr(kl,5)
              r1       = pri_arr(kl,6)

              weight = one8 * (r0 + r1) * (fi_k0(id) + fi_k1(id))

              tilde_op(l,k,:) = weight * (fi_k0(:) + fi_k1(:))

           enddo
        enddo
        
        return
end subroutine energy_preserving_tilde_op












end module flux_module
