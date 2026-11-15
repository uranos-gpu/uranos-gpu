module euler_matrices_module
use parameters_module   , only: rp
use thermodynamics_module, only: gamma_Tref
implicit none

private
public euler_system_x, euler_system_y, euler_system_n
public euler_eigenmatrix_x, euler_eigenmatrix_y, euler_jacobian

contains

subroutine euler_system_x(r_, u_, v_, w_, p_, g_, h_, Rmatrix, Lmatrix)
        !$acc routine seq
        implicit none

        real(rp), dimension(1:5, 1:5), intent(out) :: Rmatrix, Lmatrix
        real(rp),                      intent(in)  :: r_, u_, v_, w_, p_, g_, h_

        real(rp) :: c, cc, ci, icc, V2
        real(rp) :: gm1, igm1
        real(rp) :: b1, b2, b3
        real(rp) :: b2u, b2v, b2w, hb2u, hb2v, hb2w
        real(rp) :: hb1, hb2, hci, huci, uc

        cc   = g_ * p_ / r_
        c    = sqrt(cc)
        ci   = 1.0_rp / c
        icc  = ci * ci
        gm1  = g_ - 1.0_rp
        igm1 = 1.0_rp / gm1
        V2   = u_*u_ + v_*v_ + w_*w_

        b2 = gm1 * icc
        b3 = h_ - cc * igm1
        b1 = 1.0_rp + b2 * (V2 - h_)

        b2u  = b2 * u_
        b2v  = b2 * v_
        b2w  = b2 * w_
        hb2u = 0.5_rp * b2u
        hb2v = 0.5_rp * b2v
        hb2w = 0.5_rp * b2w
        hb1  = 0.5_rp * b1
        hb2  = 0.5_rp * b2
        hci  = 0.5_rp * ci
        huci = u_ * hci
        uc   = u_ * c

        Rmatrix(1,1) = 1.0_rp
        Rmatrix(2,1) = u_ - c
        Rmatrix(3,1) = v_
        Rmatrix(4,1) = w_
        Rmatrix(5,1) = h_ - uc

        Rmatrix(1,2) = 1.0_rp
        Rmatrix(2,2) = u_
        Rmatrix(3,2) = v_
        Rmatrix(4,2) = w_
        Rmatrix(5,2) = b3

        Rmatrix(1,3) = 0.0_rp
        Rmatrix(2,3) = 0.0_rp
        Rmatrix(3,3) = 1.0_rp
        Rmatrix(4,3) = 0.0_rp
        Rmatrix(5,3) = v_

        Rmatrix(1,4) = 0.0_rp
        Rmatrix(2,4) = 0.0_rp
        Rmatrix(3,4) = 0.0_rp
        Rmatrix(4,4) = 1.0_rp
        Rmatrix(5,4) = w_

        Rmatrix(1,5) = 1.0_rp
        Rmatrix(2,5) = u_ + c
        Rmatrix(3,5) = v_
        Rmatrix(4,5) = w_
        Rmatrix(5,5) = h_ + uc

        Lmatrix(1,1) = hb1 + huci
        Lmatrix(2,1) = 1.0_rp - b1
        Lmatrix(3,1) = -v_
        Lmatrix(4,1) = -w_
        Lmatrix(5,1) = hb1 - huci

        Lmatrix(1,2) = -hb2u - hci
        Lmatrix(2,2) =  b2u
        Lmatrix(3,2) =  0.0_rp
        Lmatrix(4,2) =  0.0_rp
        Lmatrix(5,2) = -hb2u + hci

        Lmatrix(1,3) = -hb2v
        Lmatrix(2,3) =  b2v
        Lmatrix(3,3) =  1.0_rp
        Lmatrix(4,3) =  0.0_rp
        Lmatrix(5,3) = -hb2v

        Lmatrix(1,4) = -hb2w
        Lmatrix(2,4) =  b2w
        Lmatrix(3,4) =  0.0_rp
        Lmatrix(4,4) =  1.0_rp
        Lmatrix(5,4) = -hb2w

        Lmatrix(1,5) =  hb2
        Lmatrix(2,5) = -b2
        Lmatrix(3,5) =  0.0_rp
        Lmatrix(4,5) =  0.0_rp
        Lmatrix(5,5) =  hb2
end subroutine euler_system_x

subroutine euler_system_y(r_, u_, v_, w_, p_, g_, h_, Rmatrix, Lmatrix)
        !$acc routine seq
        implicit none

        real(rp), dimension(1:5, 1:5), intent(out) :: Rmatrix, Lmatrix
        real(rp),                      intent(in)  :: r_, u_, v_, w_, p_, g_, h_

        real(rp) :: c, cc, ci, icc, V2
        real(rp) :: gm1, igm1
        real(rp) :: b1, b2, b3
        real(rp) :: b2u, b2v, b2w, hb2u, hb2v, hb2w
        real(rp) :: hb1, hb2, hci, hvci, vc

        cc   = g_ * p_ / r_
        c    = sqrt(cc)
        ci   = 1.0_rp / c
        icc  = ci * ci
        gm1  = g_ - 1.0_rp
        igm1 = 1.0_rp / gm1
        V2   = u_*u_ + v_*v_ + w_*w_

        b2 = gm1 * icc
        b3 = h_ - cc * igm1
        b1 = 1.0_rp + b2 * (V2 - h_)

        b2u  = b2 * u_
        b2v  = b2 * v_
        b2w  = b2 * w_
        hb2u = 0.5_rp * b2u
        hb2v = 0.5_rp * b2v
        hb2w = 0.5_rp * b2w
        hb1  = 0.5_rp * b1
        hb2  = 0.5_rp * b2
        hci  = 0.5_rp * ci
        hvci = v_ * hci
        vc   = v_ * c

        Rmatrix(1,1) = 1.0_rp
        Rmatrix(2,1) = u_
        Rmatrix(3,1) = v_ - c
        Rmatrix(4,1) = w_
        Rmatrix(5,1) = h_ - vc

        Rmatrix(1,2) = 1.0_rp
        Rmatrix(2,2) = u_
        Rmatrix(3,2) = v_
        Rmatrix(4,2) = w_
        Rmatrix(5,2) = b3

        Rmatrix(1,3) = 0.0_rp
        Rmatrix(2,3) = 1.0_rp
        Rmatrix(3,3) = 0.0_rp
        Rmatrix(4,3) = 0.0_rp
        Rmatrix(5,3) = u_

        Rmatrix(1,4) = 0.0_rp
        Rmatrix(2,4) = 0.0_rp
        Rmatrix(3,4) = 0.0_rp
        Rmatrix(4,4) = 1.0_rp
        Rmatrix(5,4) = w_

        Rmatrix(1,5) = 1.0_rp
        Rmatrix(2,5) = u_
        Rmatrix(3,5) = v_ + c
        Rmatrix(4,5) = w_
        Rmatrix(5,5) = h_ + vc

        Lmatrix(1,1) = hb1 + hvci
        Lmatrix(2,1) = 1.0_rp - b1
        Lmatrix(3,1) = -u_
        Lmatrix(4,1) = -w_
        Lmatrix(5,1) = hb1 - hvci

        Lmatrix(1,2) = -hb2u
        Lmatrix(2,2) =  b2u
        Lmatrix(3,2) =  1.0_rp
        Lmatrix(4,2) =  0.0_rp
        Lmatrix(5,2) = -hb2u

        Lmatrix(1,3) = -hb2v - hci
        Lmatrix(2,3) =  b2v
        Lmatrix(3,3) =  0.0_rp
        Lmatrix(4,3) =  0.0_rp
        Lmatrix(5,3) = -hb2v + hci

        Lmatrix(1,4) = -hb2w
        Lmatrix(2,4) =  b2w
        Lmatrix(3,4) =  0.0_rp
        Lmatrix(4,4) =  1.0_rp
        Lmatrix(5,4) = -hb2w

        Lmatrix(1,5) =  hb2
        Lmatrix(2,5) = -b2
        Lmatrix(3,5) =  0.0_rp
        Lmatrix(4,5) =  0.0_rp
        Lmatrix(5,5) =  hb2
end subroutine euler_system_y

subroutine euler_system_n(nx_, ny_, u_, v_, w_, t_, g_, h_, Rmatrix, Lmatrix)
        !$acc routine seq
        implicit none

        real(rp), dimension(1:5, 1:5), intent(out) :: Rmatrix, Lmatrix
        real(rp),                      intent(in)  :: nx_, ny_, u_, v_, w_, t_, g_, h_

        real(rp) :: c, cc, ci, icc, V2
        real(rp) :: gm1, igm1
        real(rp) :: b1, b2, b3
        real(rp) :: un, ut
        real(rp) :: b2u, b2v, b2w
        real(rp) :: cnx, cny, cun

        cc   = g_ * t_
        c    = sqrt(cc)
        ci   = 1.0_rp / c
        icc  = ci * ci
        gm1  = g_ - 1.0_rp
        igm1 = 1.0_rp / gm1
        V2   = u_*u_ + v_*v_ + w_*w_

        un =  nx_*u_ + ny_*v_
        ut = -ny_*u_ + nx_*v_

        b2 = gm1 * icc
        b3 = h_ - cc * igm1
        b1 = 1.0_rp + b2 * (V2 - h_)

        b2u = b2 * u_
        b2v = b2 * v_
        b2w = b2 * w_
        cnx = c * nx_
        cny = c * ny_
        cun = c * un

        Rmatrix(1,1) = 1.0_rp
        Rmatrix(2,1) = u_ - cnx
        Rmatrix(3,1) = v_ - cny
        Rmatrix(4,1) = w_
        Rmatrix(5,1) = h_ - cun

        Rmatrix(1,2) = 1.0_rp
        Rmatrix(2,2) = u_
        Rmatrix(3,2) = v_
        Rmatrix(4,2) = w_
        Rmatrix(5,2) = b3

        Rmatrix(1,3) =  0.0_rp
        Rmatrix(2,3) = -ny_
        Rmatrix(3,3) =  nx_
        Rmatrix(4,3) =  0.0_rp
        Rmatrix(5,3) =  ut

        Rmatrix(1,4) = 0.0_rp
        Rmatrix(2,4) = 0.0_rp
        Rmatrix(3,4) = 0.0_rp
        Rmatrix(4,4) = 1.0_rp
        Rmatrix(5,4) = w_

        Rmatrix(1,5) = 1.0_rp
        Rmatrix(2,5) = u_ + cnx
        Rmatrix(3,5) = v_ + cny
        Rmatrix(4,5) = w_
        Rmatrix(5,5) = h_ + cun

        Lmatrix(1,1) = 0.5_rp*(b1 + un*ci)
        Lmatrix(1,2) = 0.5_rp*(-b2u - nx_*ci)
        Lmatrix(1,3) = 0.5_rp*(-b2v - ny_*ci)
        Lmatrix(1,4) = 0.5_rp*(-b2w)
        Lmatrix(1,5) = 0.5_rp*b2

        Lmatrix(2,1) = 1.0_rp - b1
        Lmatrix(2,2) = b2u
        Lmatrix(2,3) = b2v
        Lmatrix(2,4) = b2w
        Lmatrix(2,5) = -b2

        Lmatrix(3,1) = -ut
        Lmatrix(3,2) = -ny_
        Lmatrix(3,3) =  nx_
        Lmatrix(3,4) =  0.0_rp
        Lmatrix(3,5) =  0.0_rp

        Lmatrix(4,1) = -w_
        Lmatrix(4,2) =  0.0_rp
        Lmatrix(4,3) =  0.0_rp
        Lmatrix(4,4) =  1.0_rp
        Lmatrix(4,5) =  0.0_rp

        Lmatrix(5,1) = 0.5_rp*(b1 - un*ci)
        Lmatrix(5,2) = 0.5_rp*(-b2u + nx_*ci)
        Lmatrix(5,3) = 0.5_rp*(-b2v + ny_*ci)
        Lmatrix(5,4) = 0.5_rp*(-b2w)
        Lmatrix(5,5) = 0.5_rp*b2
end subroutine euler_system_n

subroutine euler_eigenmatrix_x(r_, u_, v_, w_, p_, RR, LL)
        !$acc routine seq
        implicit none

        real(rp), dimension(1:5, 1:5), intent(out) :: RR, LL
        real(rp), intent(in) :: r_, u_, v_, w_, p_

        real(rp) :: c, cc, ek, ir, T_

        ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
        ir = 1.0_rp/r_
        T_ = p_*ir
        c  = sqrt(gamma_Tref*T_)
        cc = c*c

        LL(1,1) = 0.0_rp
        LL(1,2) = -r_*c
        LL(1,3) = 0.0_rp
        LL(1,4) = 0.0_rp
        LL(1,5) = 1.0_rp
        LL(2,1) = cc
        LL(2,2) = 0.0_rp
        LL(2,3) = 0.0_rp
        LL(2,4) = 0.0_rp
        LL(2,5) = -1.0_rp
        LL(3,1) = 0.0_rp
        LL(3,2) = 0.0_rp
        LL(3,3) = 1.0_rp
        LL(3,4) = 0.0_rp
        LL(3,5) = 0.0_rp
        LL(4,1) = 0.0_rp
        LL(4,2) = 0.0_rp
        LL(4,3) = 0.0_rp
        LL(4,4) = 1.0_rp
        LL(4,5) = 0.0_rp
        LL(5,1) = 0.0_rp
        LL(5,2) = r_*c
        LL(5,3) = 0.0_rp
        Ll(5,4) = 0.0_rp
        LL(5,5) = 1.0_rp

        RR(1,1) = 0.5_rp/cc
        RR(1,2) = 1.0_rp/cc
        RR(1,3) = 0.0_rp
        RR(1,4) = 0.0_rp
        RR(1,5) = 0.5_rp/cc
        RR(2,1) = -0.5_rp*ir/c
        RR(2,2) = 0.0_rp
        RR(2,3) = 0.0_rp
        RR(2,4) = 0.0_rp
        RR(2,5) = 0.5_rp*ir/c
        RR(3,1) = 0.0_rp
        RR(3,2) = 0.0_rp
        RR(3,3) = 1.0_rp
        RR(3,4) = 0.0_rp
        RR(3,5) = 0.0_rp
        RR(4,1) = 0.0_rp
        RR(4,2) = 0.0_rp
        RR(4,3) = 0.0_rp
        RR(4,4) = 1.0_rp
        RR(4,5) = 0.0_rp
        RR(5,1) = 0.5_rp
        RR(5,2) = 0.0_rp
        RR(5,3) = 0.0_rp
        RR(5,4) = 0.0_rp
        RR(5,5) = 0.5_rp

        return
end subroutine euler_eigenmatrix_x

subroutine euler_eigenmatrix_y(r_, u_, v_, w_, p_, RR, LL)
        !$acc routine seq
        implicit none

        real(rp), dimension(1:5, 1:5), intent(out) :: RR, LL
        real(rp), intent(in) :: r_, u_, v_, w_, p_

        real(rp) :: c, cc, ek, ir, T_

        ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
        ir = 1.0_rp/r_
        T_ = p_*ir
        c  = sqrt(gamma_Tref*T_)
        cc = c*c

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

        RR(1,1) = 0.5_rp/cc
        RR(1,2) = 0.0_rp
        RR(1,3) = 1.0_rp/cc
        RR(1,4) = 0.0_rp
        RR(1,5) = 0.5_rp/cc
        RR(2,1) = 0.0_rp
        RR(2,2) = 1.0_rp
        RR(2,3) = 0.0_rp
        RR(2,4) = 0.0_rp
        RR(2,5) = 0.0_rp
        RR(3,1) = -0.5_rp*ir/c
        RR(3,2) = 0.0_rp
        RR(3,3) = 0.0_rp
        RR(3,4) = 0.0_rp
        RR(3,5) = 0.5_rp*ir/ c
        RR(4,1) = 0.0_rp
        RR(4,2) = 0.0_rp
        RR(4,3) = 0.0_rp
        RR(4,4) = 1.0_rp
        RR(4,5) = 0.0_rp
        RR(5,1) = 0.5_rp
        RR(5,2) = 0.0_rp
        RR(5,3) = 0.0_rp
        RR(5,4) = 0.0_rp
        RR(5,5) = 0.5_rp

        return
end subroutine euler_eigenmatrix_y

subroutine euler_jacobian(r_,u_,v_,w_,p_,PP,iP)
        !$acc routine seq

        implicit none
        real(rp), dimension(1:5,1:5), intent(out) :: PP, iP
        real(rp)                    , intent(in)  :: r_,u_,v_,w_,p_

        real(rp) :: kk, ik
        real(rp) :: c, cc, ek, ir, T_

        kk = gamma_Tref - 1.0_rp
        ik = 1.0_rp/kk

        ek = 0.5_rp*(u_*u_ + v_*v_ + w_*w_)
        ir = 1.0_rp/r_
        T_ = p_*ir
        c  = sqrt(gamma_Tref*T_)
        cc = c*c

        PP(1,1) = 1.0_rp
        PP(1,2) = 0.0_rp
        PP(1,3) = 0.0_rp
        PP(1,4) = 0.0_rp
        PP(1,5) = 0.0_rp
        PP(2,1) = u_
        PP(2,2) = r_
        PP(2,3) = 0.0_rp
        PP(2,4) = 0.0_rp
        PP(2,5) = 0.0_rp
        PP(3,1) = v_
        PP(3,2) = 0.0_rp
        PP(3,3) = r_
        PP(3,4) = 0.0_rp
        PP(3,5) = 0.0_rp
        PP(4,1) = w_
        PP(4,2) = 0.0_rp
        PP(4,3) = 0.0_rp
        PP(4,4) = r_
        PP(4,5) = 0.0_rp
        PP(5,1) = ek
        PP(5,2) = r_*u_
        PP(5,3) = r_*v_
        PP(5,4) = r_*w_
        PP(5,5) = ik

        iP(1,1) = 1.0_rp
        iP(1,2) = 0.0_rp
        iP(1,3) = 0.0_rp
        iP(1,4) = 0.0_rp
        iP(1,5) = 0.0_rp
        iP(2,1) = -u_*ir
        iP(2,2) = ir
        iP(2,3) = 0.0_rp
        iP(2,4) = 0.0_rp
        iP(2,5) = 0.0_rp
        iP(3,1) = -v_*ir
        iP(3,2) = 0.0_rp
        iP(3,3) = ir
        iP(3,4) = 0.0_rp
        iP(3,5) = 0.0_rp
        iP(4,1) = -w_*ir
        iP(4,2) = 0.0_rp
        iP(4,3) = 0.0_rp
        iP(4,4) = ir
        iP(4,5) = 0.0_rp
        iP(5,1) = ek*kk
        iP(5,2) =-kk*u_
        iP(5,3) =-kk*v_
        iP(5,4) =-kk*w_
        iP(5,5) = kk

        return
end subroutine euler_jacobian

end module euler_matrices_module
