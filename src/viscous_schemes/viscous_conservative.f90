module viscous_conservative_module
use profiling_module
use parameters_module
use mpi_module
use mpi_comm_module, only: mpi_share_planes
use mesh_module
use storage_module

implicit none
private
public viscous_faces_conservative
public viscous_flux_3D_conservative_xi_curv , viscous_flux_3D_conservative_xi_cart
public viscous_flux_3D_conservative_eta_curv, viscous_flux_3D_conservative_eta_cart
public viscous_flux_3D_conservative_zta_curv, viscous_flux_3D_conservative_zta_cart
public divergence_3D_conservative

contains

subroutine viscous_flux_3D_conservative_xi_curv
        implicit none

        real(rp), parameter :: two3 = 2.0_rp/3.0_rp

        integer , parameter :: vec_size  = 64
        integer , parameter :: chunkSize = vec_size - 2*visc_div_reach + 1
        integer , parameter :: nnd_max   = vec_size + 2*visc_face_reach - 1

        real(rp), dimension(0:vec_size-1,1:8) :: fc
        real(rp), dimension(0:nnd_max-1,1:8)  :: nd

        real(rp) :: csi_x, csi_y, eta_x, eta_y, cof1, cof2, zeta_z
        real(rp) :: muf, lmf, Jac, bl, cm, cn, fsg
        real(rp) :: duc, dvc, dwc, dtc
        real(rp) :: due, dve, dwe, dte
        real(rp) :: duz, dvz, dwz
        real(rp) :: uh, vh, wh
        real(rp) :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tx, ty
        real(rp) :: t11, t12, t13, t22, t23, div3
        real(rp) :: flx1, flx2, flx3, qflx, gflx, pflx, rti
        real(rp) :: dfl1, dfl2, dfl3, dhe, dph, dpr, cl, tnd
        real(rp) :: rhs_u, rhs_v, rhs_w, uu, vv, ww
        integer  :: i, j, k, l, m, n, s, p, pm, pp
        integer  :: vfp_lo, vfp_hi
        integer  :: nt, q, q0, nnd, ntm, ntp
        integer  :: ft, cs0, cs1, nfa

        vfp_lo = sx - GN + visc_face_reach - 1
        vfp_hi = ex + GN - visc_face_reach

        call StartProfRange("viscous_flux_3D_conservative_xi_curv")

        !$acc parallel vector_length(vec_size) default(present) private(fc,nd)
        !$acc cache(fc,nd)
        !$acc loop gang collapse(2)
        do   k = sz,ez
          do  j = sy,ey

            zeta_z = ztastep_i(k)*iz_zta(k)

            do cs0 = sx,ex,chunkSize
              cs1 = min(cs0+chunkSize-1, ex)
              nfa = (cs1-cs0+1) + 2*visc_div_reach - 1
              nnd = nfa + 2*visc_face_reach - 1
              q0  = cs0 - visc_div_reach - visc_face_reach + 1

              !$acc loop vector
              do nt = 0,nnd-1
                q = q0 + nt
                if(q >= lbx .and. q <= ubx) then
                  due = 0.0_rp
                  dve = 0.0_rp
                  dwe = 0.0_rp
                  dte = 0.0_rp
                  duz = 0.0_rp
                  dvz = 0.0_rp
                  dwz = 0.0_rp
                  !$acc loop seq
                  do n = 1,3
                    cn  = central_1(n)
                    due = due + cn*(U(q,j+n,k) - U(q,j-n,k))
                    dve = dve + cn*(V(q,j+n,k) - V(q,j-n,k))
                    dwe = dwe + cn*(W(q,j+n,k) - W(q,j-n,k))
                    dte = dte + cn*(T(q,j+n,k) - T(q,j-n,k))
                    duz = duz + cn*(U(q,j,k+n) - U(q,j,k-n))
                    dvz = dvz + cn*(V(q,j,k+n) - V(q,j,k-n))
                    dwz = dwz + cn*(W(q,j,k+n) - W(q,j,k-n))
                  enddo
                  nd(nt,1) = due
                  nd(nt,2) = dve
                  nd(nt,3) = dwe
                  nd(nt,4) = dte
                  nd(nt,5) = duz
                  nd(nt,6) = dvz
                  nd(nt,7) = dwz
                  nd(nt,8) = 1.0_rp/T(q,j,k)
                endif
              enddo

              !$acc loop vector
              do ft = 0,nfa-1
                p = cs0 - visc_div_reach + ft

                if(p < vfp_lo .or. p > vfp_hi) then

                  if(p < vfp_lo) then
                    flx1 = vfclo_x(1,j,k)
                    flx2 = vfclo_x(2,j,k)
                    flx3 = vfclo_x(3,j,k)
                    gflx = vfclo_x(4,j,k)
                    pflx = vfclo_x(8,j,k)
                    uh   = vfclo_x(5,j,k)
                    vh   = vfclo_x(6,j,k)
                    wh   = vfclo_x(7,j,k)
                  else
                    flx1 = vfchi_x(1,j,k)
                    flx2 = vfchi_x(2,j,k)
                    flx3 = vfchi_x(3,j,k)
                    gflx = vfchi_x(4,j,k)
                    pflx = vfchi_x(8,j,k)
                    uh   = vfchi_x(5,j,k)
                    vh   = vfchi_x(6,j,k)
                    wh   = vfchi_x(7,j,k)
                  endif

                else

                  csi_x = fmx_csix(p,j)
                  csi_y = fmx_csiy(p,j)
                  eta_x = fmx_etax(p,j)
                  eta_y = fmx_etay(p,j)
                  cof1  = fmx_cof1(p,j)
                  cof2  = fmx_cof2(p,j)
                  muf   = 0.0_rp
                  rti   = 0.0_rp
                  lmf   = 0.0_rp
                  uh    = 0.0_rp
                  vh    = 0.0_rp
                  wh    = 0.0_rp
                  duc   = 0.0_rp
                  dvc   = 0.0_rp
                  dwc   = 0.0_rp
                  dtc   = 0.0_rp
                  due   = 0.0_rp
                  dve   = 0.0_rp
                  dwe   = 0.0_rp
                  dte   = 0.0_rp
                  duz   = 0.0_rp
                  dvz   = 0.0_rp
                  dwz   = 0.0_rp

                  !$acc loop seq
                  do m = 1,visc_face_reach
                    cm  = visc_intp(m)
                    pm  = p+1-m
                    pp  = p+m
                    ntm = pm - q0
                    ntp = pp - q0

                    muf = muf + cm*(VIS(pp,j,k) + VIS(pm,j,k))
                    lmf = lmf + cm*(LMD(pp,j,k) + LMD(pm,j,k))

                    uh = uh + cm*(U(pp,j,k) + U(pm,j,k))
                    vh = vh + cm*(V(pp,j,k) + V(pm,j,k))
                    wh = wh + cm*(W(pp,j,k) + W(pm,j,k))

                    duc = duc + visc_stag(m)*(U(pp,j,k) - U(pm,j,k))
                    dvc = dvc + visc_stag(m)*(V(pp,j,k) - V(pm,j,k))
                    dwc = dwc + visc_stag(m)*(W(pp,j,k) - W(pm,j,k))
                    dtc = dtc + visc_stag(m)*(T(pp,j,k) - T(pm,j,k))

                    due = due + cm*(nd(ntp,1) + nd(ntm,1))
                    dve = dve + cm*(nd(ntp,2) + nd(ntm,2))
                    dwe = dwe + cm*(nd(ntp,3) + nd(ntm,3))
                    dte = dte + cm*(nd(ntp,4) + nd(ntm,4))
                    duz = duz + cm*(nd(ntp,5) + nd(ntm,5))
                    dvz = dvz + cm*(nd(ntp,6) + nd(ntm,6))
                    dwz = dwz + cm*(nd(ntp,7) + nd(ntm,7))
                  enddo

                  ux = csi_x*duc + eta_x*due
                  uy = csi_y*duc + eta_y*due
                  uz = zeta_z*duz
                  vx = csi_x*dvc + eta_x*dve
                  vy = csi_y*dvc + eta_y*dve
                  vz = zeta_z*dvz
                  wx = csi_x*dwc + eta_x*dwe
                  wy = csi_y*dwc + eta_y*dwe
                  wz = zeta_z*dwz
                  tx = csi_x*dtc + eta_x*dte
                  ty = csi_y*dtc + eta_y*dte

                  div3 = two3*(ux + vy + wz)
                  t11  = muf*(2.0_rp*ux - div3)
                  t12  = muf*(uy + vx)
                  t13  = muf*(uz + wx)
                  t22  = muf*(2.0_rp*vy - div3)
                  t23  = muf*(vz + wy)

                  flx1 = cof1*t11 + cof2*t12
                  flx2 = cof1*t12 + cof2*t22
                  flx3 = cof1*t13 + cof2*t23
                  qflx = lmf*(cof1*tx + cof2*ty)
                  !$acc loop seq
                  do m = 1,visc_face_reach
                    rti = rti + visc_intp(m)*(nd(p+m-q0,8) + nd(p+1-m-q0,8))
                  enddo
                  gflx = qflx*rti
                  pflx = gflx*dtc
                endif

                fc(ft,1) = flx1
                fc(ft,2) = flx2
                fc(ft,3) = flx3
                fc(ft,4) = gflx
                fc(ft,5) = uh
                fc(ft,6) = vh
                fc(ft,7) = wh
                fc(ft,8) = pflx
              enddo

              !$acc loop vector
              do i = cs0,cs1

                uu   = U(i,j,k)
                vv   = V(i,j,k)
                ww   = W(i,j,k)
                tnd  = T(i,j,k)
                dfl1 = 0.0_rp
                dfl2 = 0.0_rp
                dfl3 = 0.0_rp
                dhe  = 0.0_rp
                dpr  = 0.0_rp
                dph  = 0.0_rp

                !$acc loop seq
                do l = 1,visc_div_reach
                  bl = visc_divs(l)
                  cl = visc_intp(l)
                  !$acc loop seq
                  do s = 0,1
                    if(s == 1) then
                      p   =  i+l-1
                      fsg =  bl
                    else
                      p   =  i-l
                      fsg = -bl
                    endif
                    ft = p - cs0 + visc_div_reach

                    flx1 = fc(ft,1)
                    flx2 = fc(ft,2)
                    flx3 = fc(ft,3)
                    gflx = fc(ft,4)
                    pflx = fc(ft,8)
                    uh   = fc(ft,5)
                    vh   = fc(ft,6)
                    wh   = fc(ft,7)

                    dfl1 = dfl1 + fsg*flx1
                    dfl2 = dfl2 + fsg*flx2
                    dfl3 = dfl3 + fsg*flx3
                    dhe  = dhe  + fsg*gflx
                    dpr  = dpr  + cl*pflx
                    dph  = dph  + fsg*(flx1*(uh-uu) + flx2*(vh-vv) + flx3*(wh-ww))
                  enddo
                enddo

                Jac   = Jac2D(i,j)
                rhs_u = Jac*dfl1
                rhs_v = Jac*dfl2
                rhs_w = Jac*dfl3

                RHS(i,j,k,2) = RHS(i,j,k,2) + mu_inf*rhs_u
                RHS(i,j,k,3) = RHS(i,j,k,3) + mu_inf*rhs_v
                RHS(i,j,k,4) = RHS(i,j,k,4) + mu_inf*rhs_w
                RHS(i,j,k,5) = RHS(i,j,k,5) &
                             + mu_inf*(uu*rhs_u + vv*rhs_v + ww*rhs_w + Jac*(dph + tnd*dhe + dpr))

              enddo
            enddo
          enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine viscous_flux_3D_conservative_xi_curv

subroutine viscous_flux_3D_conservative_xi_cart
        implicit none

        real(rp), parameter :: two3 = 2.0_rp/3.0_rp

        integer , parameter :: vec_size  = 64
        integer , parameter :: chunkSize = vec_size - 2*visc_div_reach + 1
        integer , parameter :: nnd_max   = vec_size + 2*visc_face_reach - 1

        real(rp), dimension(0:vec_size-1,1:8) :: fc
        real(rp), dimension(0:nnd_max-1,1:5)  :: nd

        real(rp) :: csi_x, eta_y, cof1, zeta_z
        real(rp) :: muf, lmf, Jac, bl, cm, cn, fsg
        real(rp) :: duc, dvc, dwc, dtc
        real(rp) :: due, dve, dwe, dte
        real(rp) :: duz, dvz, dwz
        real(rp) :: uh, vh, wh
        real(rp) :: ux, uy, uz, vx, vy, wx, wz, tx
        real(rp) :: t11, t12, t13, div3
        real(rp) :: flx1, flx2, flx3, qflx, gflx, pflx, rti
        real(rp) :: dfl1, dfl2, dfl3, dhe, dph, dpr, cl, tnd
        real(rp) :: rhs_u, rhs_v, rhs_w, uu, vv, ww
        integer  :: i, j, k, l, m, n, s, p, pm, pp
        integer  :: vfp_lo, vfp_hi
        integer  :: nt, q, q0, nnd, ntm, ntp
        integer  :: ft, cs0, cs1, nfa

        vfp_lo = sx - GN + visc_face_reach - 1
        vfp_hi = ex + GN - visc_face_reach

        call StartProfRange("viscous_flux_3D_conservative_xi_cart")

        !$acc parallel vector_length(vec_size) default(present) private(fc,nd)
        !$acc cache(fc,nd)
        !$acc loop gang collapse(2)
        do   k = sz,ez
          do  j = sy,ey

            zeta_z = ztastep_i(k)*iz_zta(k)

            do cs0 = sx,ex,chunkSize
              cs1 = min(cs0+chunkSize-1, ex)
              nfa = (cs1-cs0+1) + 2*visc_div_reach - 1
              nnd = nfa + 2*visc_face_reach - 1
              q0  = cs0 - visc_div_reach - visc_face_reach + 1

              !$acc loop vector
              do nt = 0,nnd-1
                q = q0 + nt
                if(q >= lbx .and. q <= ubx) then
                  due = 0.0_rp
                  dve = 0.0_rp
                  duz = 0.0_rp
                  dwz = 0.0_rp
                  !$acc loop seq
                  do n = 1,3
                    cn  = central_1(n)
                    due = due + cn*(U(q,j+n,k) - U(q,j-n,k))
                    dve = dve + cn*(V(q,j+n,k) - V(q,j-n,k))
                    duz = duz + cn*(U(q,j,k+n) - U(q,j,k-n))
                    dwz = dwz + cn*(W(q,j,k+n) - W(q,j,k-n))
                  enddo
                  nd(nt,1) = due
                  nd(nt,2) = dve
                  nd(nt,3) = duz
                  nd(nt,4) = dwz
                  nd(nt,5) = 1.0_rp/T(q,j,k)
                endif
              enddo

              !$acc loop vector
              do ft = 0,nfa-1
                p = cs0 - visc_div_reach + ft

                if(p < vfp_lo .or. p > vfp_hi) then

                  if(p < vfp_lo) then
                    flx1 = vfclo_x(1,j,k)
                    flx2 = vfclo_x(2,j,k)
                    flx3 = vfclo_x(3,j,k)
                    gflx = vfclo_x(4,j,k)
                    pflx = vfclo_x(8,j,k)
                    uh   = vfclo_x(5,j,k)
                    vh   = vfclo_x(6,j,k)
                    wh   = vfclo_x(7,j,k)
                  else
                    flx1 = vfchi_x(1,j,k)
                    flx2 = vfchi_x(2,j,k)
                    flx3 = vfchi_x(3,j,k)
                    gflx = vfchi_x(4,j,k)
                    pflx = vfchi_x(8,j,k)
                    uh   = vfchi_x(5,j,k)
                    vh   = vfchi_x(6,j,k)
                    wh   = vfchi_x(7,j,k)
                  endif

                else

                  csi_x = fmx_csix(p,j)
                  eta_y = fmx_etay(p,j)
                  cof1  = fmx_cof1(p,j)
                  muf   = 0.0_rp
                  rti   = 0.0_rp
                  lmf   = 0.0_rp
                  uh    = 0.0_rp
                  vh    = 0.0_rp
                  wh    = 0.0_rp
                  duc   = 0.0_rp
                  dvc   = 0.0_rp
                  dwc   = 0.0_rp
                  dtc   = 0.0_rp
                  due   = 0.0_rp
                  dve   = 0.0_rp
                  dwe   = 0.0_rp
                  dte   = 0.0_rp
                  duz   = 0.0_rp
                  dvz   = 0.0_rp
                  dwz   = 0.0_rp

                  !$acc loop seq
                  do m = 1,visc_face_reach
                    cm  = visc_intp(m)
                    pm  = p+1-m
                    pp  = p+m
                    ntm = pm - q0
                    ntp = pp - q0

                    muf = muf + cm*(VIS(pp,j,k) + VIS(pm,j,k))
                    lmf = lmf + cm*(LMD(pp,j,k) + LMD(pm,j,k))
                    rti = rti + cm*(nd(ntp,5) + nd(ntm,5))

                    uh = uh + cm*(U(pp,j,k) + U(pm,j,k))
                    vh = vh + cm*(V(pp,j,k) + V(pm,j,k))
                    wh = wh + cm*(W(pp,j,k) + W(pm,j,k))

                    duc = duc + visc_stag(m)*(U(pp,j,k) - U(pm,j,k))
                    dvc = dvc + visc_stag(m)*(V(pp,j,k) - V(pm,j,k))
                    dwc = dwc + visc_stag(m)*(W(pp,j,k) - W(pm,j,k))
                    dtc = dtc + visc_stag(m)*(T(pp,j,k) - T(pm,j,k))

                    due = due + cm*(nd(ntp,1) + nd(ntm,1))
                    dve = dve + cm*(nd(ntp,2) + nd(ntm,2))
                    duz = duz + cm*(nd(ntp,3) + nd(ntm,3))
                    dwz = dwz + cm*(nd(ntp,4) + nd(ntm,4))
                  enddo

                  ux = csi_x*duc
                  uy = eta_y*due
                  uz = zeta_z*duz
                  vx = csi_x*dvc
                  vy = eta_y*dve
                  wx = csi_x*dwc
                  wz = zeta_z*dwz
                  tx = csi_x*dtc

                  div3 = two3*(ux + vy + wz)
                  t11  = muf*(2.0_rp*ux - div3)
                  t12  = muf*(uy + vx)
                  t13  = muf*(uz + wx)

                  flx1 = cof1*t11
                  flx2 = cof1*t12
                  flx3 = cof1*t13
                  qflx = lmf*(cof1*tx)
                  gflx = qflx*rti
                  pflx = gflx*dtc
                endif

                fc(ft,1) = flx1
                fc(ft,2) = flx2
                fc(ft,3) = flx3
                fc(ft,4) = gflx
                fc(ft,5) = uh
                fc(ft,6) = vh
                fc(ft,7) = wh
                fc(ft,8) = pflx
              enddo

              !$acc loop vector
              do i = cs0,cs1

                uu   = U(i,j,k)
                vv   = V(i,j,k)
                ww   = W(i,j,k)
                tnd  = T(i,j,k)
                dfl1 = 0.0_rp
                dfl2 = 0.0_rp
                dfl3 = 0.0_rp
                dhe  = 0.0_rp
                dpr  = 0.0_rp
                dph  = 0.0_rp

                !$acc loop seq
                do l = 1,visc_div_reach
                  bl = visc_divs(l)
                  cl = visc_intp(l)
                  !$acc loop seq
                  do s = 0,1
                    if(s == 1) then
                      p   =  i+l-1
                      fsg =  bl
                    else
                      p   =  i-l
                      fsg = -bl
                    endif
                    ft = p - cs0 + visc_div_reach

                    flx1 = fc(ft,1)
                    flx2 = fc(ft,2)
                    flx3 = fc(ft,3)
                    gflx = fc(ft,4)
                    pflx = fc(ft,8)
                    uh   = fc(ft,5)
                    vh   = fc(ft,6)
                    wh   = fc(ft,7)

                    dfl1 = dfl1 + fsg*flx1
                    dfl2 = dfl2 + fsg*flx2
                    dfl3 = dfl3 + fsg*flx3
                    dhe  = dhe  + fsg*gflx
                    dpr  = dpr  + cl*pflx
                    dph  = dph  + fsg*(flx1*(uh-uu) + flx2*(vh-vv) + flx3*(wh-ww))
                  enddo
                enddo

                Jac   = Jac2D(i,j)
                rhs_u = Jac*dfl1
                rhs_v = Jac*dfl2
                rhs_w = Jac*dfl3

                RHS(i,j,k,2) = RHS(i,j,k,2) + mu_inf*rhs_u
                RHS(i,j,k,3) = RHS(i,j,k,3) + mu_inf*rhs_v
                RHS(i,j,k,4) = RHS(i,j,k,4) + mu_inf*rhs_w
                RHS(i,j,k,5) = RHS(i,j,k,5) &
                             + mu_inf*(uu*rhs_u + vv*rhs_v + ww*rhs_w + Jac*(dph + tnd*dhe + dpr))

              enddo
            enddo
          enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine viscous_flux_3D_conservative_xi_cart

subroutine viscous_flux_3D_conservative_eta_curv
        implicit none

        real(rp), parameter :: two3 = 2.0_rp/3.0_rp

        integer , parameter :: ilanes    = 4
        integer , parameter :: svec      = 32
        integer , parameter :: vec_size  = ilanes*svec
        integer , parameter :: chunkSize = svec - 2*visc_div_reach + 1
        integer , parameter :: nnd_max   = svec + 2*visc_face_reach - 1

        real(rp), dimension(0:ilanes-1,0:svec-1,1:8)    :: fc
        real(rp), dimension(0:ilanes-1,0:nnd_max-1,1:8) :: nd

        real(rp) :: csi_x, csi_y, eta_x, eta_y, cof1, cof2, zeta_z
        real(rp) :: muf, lmf, Jac, bl, cm, cn, fsg
        real(rp) :: duc, dvc, dwc, dtc
        real(rp) :: due, dve, dwe, dte
        real(rp) :: duz, dvz, dwz
        real(rp) :: uh, vh, wh
        real(rp) :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tx, ty
        real(rp) :: t11, t12, t13, t22, t23, div3
        real(rp) :: flx1, flx2, flx3, qflx, gflx, pflx, rti
        real(rp) :: dfl1, dfl2, dfl3, dhe, dph, dpr, cl, tnd
        real(rp) :: rhs_u, rhs_v, rhs_w, uu, vv, ww
        integer  :: i, j, k, l, m, n, s, p, pm, pp
        integer  :: vfp_lo, vfp_hi
        integer  :: nt, q, q0, nnd, ntm, ntp
        integer  :: ft, ii, ib, cs0, cs1, nfa

        vfp_lo = sy - GN + visc_face_reach - 1
        vfp_hi = ey + GN - visc_face_reach

        call StartProfRange("viscous_flux_3D_conservative_eta_curv")

        !$acc parallel vector_length(vec_size) default(present) private(fc,nd)
        !$acc cache(fc,nd)
        !$acc loop gang collapse(2)
        do   k = sz,ez
          do  ib = sx,ex,ilanes

            zeta_z = ztastep_i(k)*iz_zta(k)

            do cs0 = sy,ey,chunkSize
              cs1 = min(cs0+chunkSize-1, ey)
              nfa = (cs1-cs0+1) + 2*visc_div_reach - 1
              nnd = nfa + 2*visc_face_reach - 1
              q0  = cs0 - visc_div_reach - visc_face_reach + 1

              !$acc loop vector collapse(2)
              do nt = 0,nnd-1
                do ii = 0,ilanes-1
                  i = ib + ii
                  q = q0 + nt
                  if(i <= ex .and. q >= lby .and. q <= uby) then
                    duc = 0.0_rp
                    dvc = 0.0_rp
                    dwc = 0.0_rp
                    dtc = 0.0_rp
                    duz = 0.0_rp
                    dvz = 0.0_rp
                    dwz = 0.0_rp
                    !$acc loop seq
                    do n = 1,3
                      cn  = central_1(n)
                      duc = duc + cn*(U(i+n,q,k) - U(i-n,q,k))
                      dvc = dvc + cn*(V(i+n,q,k) - V(i-n,q,k))
                      dwc = dwc + cn*(W(i+n,q,k) - W(i-n,q,k))
                      dtc = dtc + cn*(T(i+n,q,k) - T(i-n,q,k))
                      duz = duz + cn*(U(i,q,k+n) - U(i,q,k-n))
                      dvz = dvz + cn*(V(i,q,k+n) - V(i,q,k-n))
                      dwz = dwz + cn*(W(i,q,k+n) - W(i,q,k-n))
                    enddo
                    nd(ii,nt,1) = duc
                    nd(ii,nt,2) = dvc
                    nd(ii,nt,3) = dwc
                    nd(ii,nt,4) = dtc
                    nd(ii,nt,5) = duz
                    nd(ii,nt,6) = dvz
                    nd(ii,nt,7) = dwz
                    nd(ii,nt,8) = 1.0_rp/T(i,q,k)
                  endif
                enddo
              enddo

              !$acc loop vector collapse(2)
              do ft = 0,nfa-1
                do ii = 0,ilanes-1
                  i = ib + ii
                  if(i <= ex) then
                    p = cs0 - visc_div_reach + ft

                    if(p < vfp_lo .or. p > vfp_hi) then

                      if(p < vfp_lo) then
                        flx1 = vfclo_y(1,i,k)
                        flx2 = vfclo_y(2,i,k)
                        flx3 = vfclo_y(3,i,k)
                        gflx = vfclo_y(4,i,k)
                        pflx = vfclo_y(8,i,k)
                        uh   = vfclo_y(5,i,k)
                        vh   = vfclo_y(6,i,k)
                        wh   = vfclo_y(7,i,k)
                      else
                        flx1 = vfchi_y(1,i,k)
                        flx2 = vfchi_y(2,i,k)
                        flx3 = vfchi_y(3,i,k)
                        gflx = vfchi_y(4,i,k)
                        pflx = vfchi_y(8,i,k)
                        uh   = vfchi_y(5,i,k)
                        vh   = vfchi_y(6,i,k)
                        wh   = vfchi_y(7,i,k)
                      endif

                    else

                      csi_x = fmy_csix(i,p)
                      csi_y = fmy_csiy(i,p)
                      eta_x = fmy_etax(i,p)
                      eta_y = fmy_etay(i,p)
                      cof1  = fmy_cof1(i,p)
                      cof2  = fmy_cof2(i,p)
                      muf   = 0.0_rp
                      rti   = 0.0_rp
                      lmf   = 0.0_rp
                      uh    = 0.0_rp
                      vh    = 0.0_rp
                      wh    = 0.0_rp
                      duc   = 0.0_rp
                      dvc   = 0.0_rp
                      dwc   = 0.0_rp
                      dtc   = 0.0_rp
                      due   = 0.0_rp
                      dve   = 0.0_rp
                      dwe   = 0.0_rp
                      dte   = 0.0_rp
                      duz   = 0.0_rp
                      dvz   = 0.0_rp
                      dwz   = 0.0_rp

                      !$acc loop seq
                      do m = 1,visc_face_reach
                        cm  = visc_intp(m)
                        pm  = p+1-m
                        pp  = p+m
                        ntm = pm - q0
                        ntp = pp - q0

                        muf = muf + cm*(VIS(i,pp,k) + VIS(i,pm,k))
                        lmf = lmf + cm*(LMD(i,pp,k) + LMD(i,pm,k))
                        rti = rti + cm*(nd(ii,ntp,8) + nd(ii,ntm,8))

                        uh = uh + cm*(U(i,pp,k) + U(i,pm,k))
                        vh = vh + cm*(V(i,pp,k) + V(i,pm,k))
                        wh = wh + cm*(W(i,pp,k) + W(i,pm,k))

                        due = due + visc_stag(m)*(U(i,pp,k) - U(i,pm,k))
                        dve = dve + visc_stag(m)*(V(i,pp,k) - V(i,pm,k))
                        dwe = dwe + visc_stag(m)*(W(i,pp,k) - W(i,pm,k))
                        dte = dte + visc_stag(m)*(T(i,pp,k) - T(i,pm,k))

                        duc = duc + cm*(nd(ii,ntp,1) + nd(ii,ntm,1))
                        dvc = dvc + cm*(nd(ii,ntp,2) + nd(ii,ntm,2))
                        dwc = dwc + cm*(nd(ii,ntp,3) + nd(ii,ntm,3))
                        dtc = dtc + cm*(nd(ii,ntp,4) + nd(ii,ntm,4))
                        duz = duz + cm*(nd(ii,ntp,5) + nd(ii,ntm,5))
                        dvz = dvz + cm*(nd(ii,ntp,6) + nd(ii,ntm,6))
                        dwz = dwz + cm*(nd(ii,ntp,7) + nd(ii,ntm,7))
                      enddo

                      ux = csi_x*duc + eta_x*due
                      uy = csi_y*duc + eta_y*due
                      uz = zeta_z*duz
                      vx = csi_x*dvc + eta_x*dve
                      vy = csi_y*dvc + eta_y*dve
                      vz = zeta_z*dvz
                      wx = csi_x*dwc + eta_x*dwe
                      wy = csi_y*dwc + eta_y*dwe
                      wz = zeta_z*dwz
                      tx = csi_x*dtc + eta_x*dte
                      ty = csi_y*dtc + eta_y*dte

                      div3 = two3*(ux + vy + wz)
                      t11  = muf*(2.0_rp*ux - div3)
                      t12  = muf*(uy + vx)
                      t13  = muf*(uz + wx)
                      t22  = muf*(2.0_rp*vy - div3)
                      t23  = muf*(vz + wy)

                      flx1 = cof1*t11 + cof2*t12
                      flx2 = cof1*t12 + cof2*t22
                      flx3 = cof1*t13 + cof2*t23
                      qflx = lmf*(cof1*tx + cof2*ty)
                      gflx = qflx*rti
                      pflx = gflx*dte
                    endif

                    fc(ii,ft,1) = flx1
                    fc(ii,ft,2) = flx2
                    fc(ii,ft,3) = flx3
                    fc(ii,ft,4) = gflx
                    fc(ii,ft,5) = uh
                    fc(ii,ft,6) = vh
                    fc(ii,ft,7) = wh
                    fc(ii,ft,8) = pflx
                  endif
                enddo
              enddo

              !$acc loop vector collapse(2)
              do j = cs0,cs1
                do ii = 0,ilanes-1
                  i = ib + ii
                  if(i <= ex) then

                    uu   = U(i,j,k)
                    vv   = V(i,j,k)
                    ww   = W(i,j,k)
                    tnd  = T(i,j,k)
                    dfl1 = 0.0_rp
                    dfl2 = 0.0_rp
                    dfl3 = 0.0_rp
                    dhe  = 0.0_rp
                    dpr  = 0.0_rp
                    dph  = 0.0_rp

                    !$acc loop seq
                    do l = 1,visc_div_reach
                      bl = visc_divs(l)
                      cl = visc_intp(l)
                      !$acc loop seq
                      do s = 0,1
                        if(s == 1) then
                          p   =  j+l-1
                          fsg =  bl
                        else
                          p   =  j-l
                          fsg = -bl
                        endif
                        ft = p - cs0 + visc_div_reach

                        flx1 = fc(ii,ft,1)
                        flx2 = fc(ii,ft,2)
                        flx3 = fc(ii,ft,3)
                        gflx = fc(ii,ft,4)
                        pflx = fc(ii,ft,8)
                        uh   = fc(ii,ft,5)
                        vh   = fc(ii,ft,6)
                        wh   = fc(ii,ft,7)

                        dfl1 = dfl1 + fsg*flx1
                        dfl2 = dfl2 + fsg*flx2
                        dfl3 = dfl3 + fsg*flx3
                        dhe  = dhe  + fsg*gflx
                        dpr  = dpr  + cl*pflx
                        dph  = dph  + fsg*(flx1*(uh-uu) + flx2*(vh-vv) + flx3*(wh-ww))
                      enddo
                    enddo

                    Jac   = Jac2D(i,j)
                    rhs_u = Jac*dfl1
                    rhs_v = Jac*dfl2
                    rhs_w = Jac*dfl3

                    RHS(i,j,k,2) = RHS(i,j,k,2) + mu_inf*rhs_u
                    RHS(i,j,k,3) = RHS(i,j,k,3) + mu_inf*rhs_v
                    RHS(i,j,k,4) = RHS(i,j,k,4) + mu_inf*rhs_w
                    RHS(i,j,k,5) = RHS(i,j,k,5) &
                                 + mu_inf*(uu*rhs_u + vv*rhs_v + ww*rhs_w + Jac*(dph + tnd*dhe + dpr))

                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine viscous_flux_3D_conservative_eta_curv

subroutine viscous_flux_3D_conservative_eta_cart
        implicit none

        real(rp), parameter :: two3 = 2.0_rp/3.0_rp

        integer , parameter :: ilanes    = 4
        integer , parameter :: svec      = 32
        integer , parameter :: vec_size  = ilanes*svec
        integer , parameter :: chunkSize = svec - 2*visc_div_reach + 1
        integer , parameter :: nnd_max   = svec + 2*visc_face_reach - 1

        real(rp), dimension(0:ilanes-1,0:svec-1,1:8)    :: fc
        real(rp), dimension(0:ilanes-1,0:nnd_max-1,1:5) :: nd

        real(rp) :: csi_x, eta_y, cof2, zeta_z
        real(rp) :: muf, lmf, Jac, bl, cm, cn, fsg
        real(rp) :: duc, dvc, dwc, dtc
        real(rp) :: due, dve, dwe, dte
        real(rp) :: duz, dvz, dwz
        real(rp) :: uh, vh, wh
        real(rp) :: ux, uy, vx, vy, vz, wy, wz, ty
        real(rp) :: t12, t22, t23, div3
        real(rp) :: flx1, flx2, flx3, qflx, gflx, pflx, rti
        real(rp) :: dfl1, dfl2, dfl3, dhe, dph, dpr, cl, tnd
        real(rp) :: rhs_u, rhs_v, rhs_w, uu, vv, ww
        integer  :: i, j, k, l, m, n, s, p, pm, pp
        integer  :: vfp_lo, vfp_hi
        integer  :: nt, q, q0, nnd, ntm, ntp
        integer  :: ft, ii, ib, cs0, cs1, nfa

        vfp_lo = sy - GN + visc_face_reach - 1
        vfp_hi = ey + GN - visc_face_reach

        call StartProfRange("viscous_flux_3D_conservative_eta_cart")

        !$acc parallel vector_length(vec_size) default(present) private(fc,nd)
        !$acc cache(fc,nd)
        !$acc loop gang collapse(2)
        do   k = sz,ez
          do  ib = sx,ex,ilanes

            zeta_z = ztastep_i(k)*iz_zta(k)

            do cs0 = sy,ey,chunkSize
              cs1 = min(cs0+chunkSize-1, ey)
              nfa = (cs1-cs0+1) + 2*visc_div_reach - 1
              nnd = nfa + 2*visc_face_reach - 1
              q0  = cs0 - visc_div_reach - visc_face_reach + 1

              !$acc loop vector collapse(2)
              do nt = 0,nnd-1
                do ii = 0,ilanes-1
                  i = ib + ii
                  q = q0 + nt
                  if(i <= ex .and. q >= lby .and. q <= uby) then
                    duc = 0.0_rp
                    dvc = 0.0_rp
                    dvz = 0.0_rp
                    dwz = 0.0_rp
                    !$acc loop seq
                    do n = 1,3
                      cn  = central_1(n)
                      duc = duc + cn*(U(i+n,q,k) - U(i-n,q,k))
                      dvc = dvc + cn*(V(i+n,q,k) - V(i-n,q,k))
                      dvz = dvz + cn*(V(i,q,k+n) - V(i,q,k-n))
                      dwz = dwz + cn*(W(i,q,k+n) - W(i,q,k-n))
                    enddo
                    nd(ii,nt,1) = duc
                    nd(ii,nt,2) = dvc
                    nd(ii,nt,3) = dvz
                    nd(ii,nt,4) = dwz
                    nd(ii,nt,5) = 1.0_rp/T(i,q,k)
                  endif
                enddo
              enddo

              !$acc loop vector collapse(2)
              do ft = 0,nfa-1
                do ii = 0,ilanes-1
                  i = ib + ii
                  if(i <= ex) then
                    p = cs0 - visc_div_reach + ft

                    if(p < vfp_lo .or. p > vfp_hi) then

                      if(p < vfp_lo) then
                        flx1 = vfclo_y(1,i,k)
                        flx2 = vfclo_y(2,i,k)
                        flx3 = vfclo_y(3,i,k)
                        gflx = vfclo_y(4,i,k)
                        pflx = vfclo_y(8,i,k)
                        uh   = vfclo_y(5,i,k)
                        vh   = vfclo_y(6,i,k)
                        wh   = vfclo_y(7,i,k)
                      else
                        flx1 = vfchi_y(1,i,k)
                        flx2 = vfchi_y(2,i,k)
                        flx3 = vfchi_y(3,i,k)
                        gflx = vfchi_y(4,i,k)
                        pflx = vfchi_y(8,i,k)
                        uh   = vfchi_y(5,i,k)
                        vh   = vfchi_y(6,i,k)
                        wh   = vfchi_y(7,i,k)
                      endif

                    else

                      csi_x = fmy_csix(i,p)
                      eta_y = fmy_etay(i,p)
                      cof2  = fmy_cof2(i,p)
                      muf   = 0.0_rp
                      rti   = 0.0_rp
                      lmf   = 0.0_rp
                      uh    = 0.0_rp
                      vh    = 0.0_rp
                      wh    = 0.0_rp
                      duc   = 0.0_rp
                      dvc   = 0.0_rp
                      dwc   = 0.0_rp
                      dtc   = 0.0_rp
                      due   = 0.0_rp
                      dve   = 0.0_rp
                      dwe   = 0.0_rp
                      dte   = 0.0_rp
                      duz   = 0.0_rp
                      dvz   = 0.0_rp
                      dwz   = 0.0_rp

                      !$acc loop seq
                      do m = 1,visc_face_reach
                        cm  = visc_intp(m)
                        pm  = p+1-m
                        pp  = p+m
                        ntm = pm - q0
                        ntp = pp - q0

                        muf = muf + cm*(VIS(i,pp,k) + VIS(i,pm,k))
                        lmf = lmf + cm*(LMD(i,pp,k) + LMD(i,pm,k))
                        rti = rti + cm*(nd(ii,ntp,5) + nd(ii,ntm,5))

                        uh = uh + cm*(U(i,pp,k) + U(i,pm,k))
                        vh = vh + cm*(V(i,pp,k) + V(i,pm,k))
                        wh = wh + cm*(W(i,pp,k) + W(i,pm,k))

                        due = due + visc_stag(m)*(U(i,pp,k) - U(i,pm,k))
                        dve = dve + visc_stag(m)*(V(i,pp,k) - V(i,pm,k))
                        dwe = dwe + visc_stag(m)*(W(i,pp,k) - W(i,pm,k))
                        dte = dte + visc_stag(m)*(T(i,pp,k) - T(i,pm,k))

                        duc = duc + cm*(nd(ii,ntp,1) + nd(ii,ntm,1))
                        dvc = dvc + cm*(nd(ii,ntp,2) + nd(ii,ntm,2))
                        dvz = dvz + cm*(nd(ii,ntp,3) + nd(ii,ntm,3))
                        dwz = dwz + cm*(nd(ii,ntp,4) + nd(ii,ntm,4))
                      enddo

                      ux = csi_x*duc
                      uy = eta_y*due
                      vx = csi_x*dvc
                      vy = eta_y*dve
                      vz = zeta_z*dvz
                      wy = eta_y*dwe
                      wz = zeta_z*dwz
                      ty = eta_y*dte

                      div3 = two3*(ux + vy + wz)
                      t12  = muf*(uy + vx)
                      t22  = muf*(2.0_rp*vy - div3)
                      t23  = muf*(vz + wy)

                      flx1 = cof2*t12
                      flx2 = cof2*t22
                      flx3 = cof2*t23
                      qflx = lmf*(cof2*ty)
                      gflx = qflx*rti
                      pflx = gflx*dte
                    endif

                    fc(ii,ft,1) = flx1
                    fc(ii,ft,2) = flx2
                    fc(ii,ft,3) = flx3
                    fc(ii,ft,4) = gflx
                    fc(ii,ft,5) = uh
                    fc(ii,ft,6) = vh
                    fc(ii,ft,7) = wh
                    fc(ii,ft,8) = pflx
                  endif
                enddo
              enddo

              !$acc loop vector collapse(2)
              do j = cs0,cs1
                do ii = 0,ilanes-1
                  i = ib + ii
                  if(i <= ex) then

                    uu   = U(i,j,k)
                    vv   = V(i,j,k)
                    ww   = W(i,j,k)
                    tnd  = T(i,j,k)
                    dfl1 = 0.0_rp
                    dfl2 = 0.0_rp
                    dfl3 = 0.0_rp
                    dhe  = 0.0_rp
                    dpr  = 0.0_rp
                    dph  = 0.0_rp

                    !$acc loop seq
                    do l = 1,visc_div_reach
                      bl = visc_divs(l)
                      cl = visc_intp(l)
                      !$acc loop seq
                      do s = 0,1
                        if(s == 1) then
                          p   =  j+l-1
                          fsg =  bl
                        else
                          p   =  j-l
                          fsg = -bl
                        endif
                        ft = p - cs0 + visc_div_reach

                        flx1 = fc(ii,ft,1)
                        flx2 = fc(ii,ft,2)
                        flx3 = fc(ii,ft,3)
                        gflx = fc(ii,ft,4)
                        pflx = fc(ii,ft,8)
                        uh   = fc(ii,ft,5)
                        vh   = fc(ii,ft,6)
                        wh   = fc(ii,ft,7)

                        dfl1 = dfl1 + fsg*flx1
                        dfl2 = dfl2 + fsg*flx2
                        dfl3 = dfl3 + fsg*flx3
                        dhe  = dhe  + fsg*gflx
                        dpr  = dpr  + cl*pflx
                        dph  = dph  + fsg*(flx1*(uh-uu) + flx2*(vh-vv) + flx3*(wh-ww))
                      enddo
                    enddo

                    Jac   = Jac2D(i,j)
                    rhs_u = Jac*dfl1
                    rhs_v = Jac*dfl2
                    rhs_w = Jac*dfl3

                    RHS(i,j,k,2) = RHS(i,j,k,2) + mu_inf*rhs_u
                    RHS(i,j,k,3) = RHS(i,j,k,3) + mu_inf*rhs_v
                    RHS(i,j,k,4) = RHS(i,j,k,4) + mu_inf*rhs_w
                    RHS(i,j,k,5) = RHS(i,j,k,5) &
                                 + mu_inf*(uu*rhs_u + vv*rhs_v + ww*rhs_w + Jac*(dph + tnd*dhe + dpr))

                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine viscous_flux_3D_conservative_eta_cart

subroutine viscous_flux_3D_conservative_zta_curv
        implicit none

        real(rp), parameter :: two3 = 2.0_rp/3.0_rp

        integer , parameter :: ilanes    = 4
        integer , parameter :: svec      = 32
        integer , parameter :: vec_size  = ilanes*svec
        integer , parameter :: chunkSize = svec - 2*visc_div_reach + 1
        integer , parameter :: nnd_max   = svec + 2*visc_face_reach - 1

        real(rp), dimension(0:ilanes-1,0:svec-1,1:8)    :: fc
        real(rp), dimension(0:ilanes-1,0:nnd_max-1,1:7) :: nd

        real(rp) :: csi_x, csi_y, eta_x, eta_y, zeta_z, zetaf
        real(rp) :: muf, lmf, bl, cm, cn, fsg
        real(rp) :: duc, dvc, dwc
        real(rp) :: due, dve, dwe
        real(rp) :: duz, dvz, dwz, dtz
        real(rp) :: uh, vh, wh
        real(rp) :: ux, uz, vy, vz, wx, wy, wz, tz
        real(rp) :: t13, t23, t33, div3
        real(rp) :: flx1, flx2, flx3, qflx, gflx, pflx, rti
        real(rp) :: dfl1, dfl2, dfl3, dhe, dph, dpr, cl, tnd
        real(rp) :: rhs_u, rhs_v, rhs_w, uu, vv, ww
        integer  :: i, j, k, l, m, n, s, p, pm, pp
        integer  :: vfp_lo, vfp_hi
        integer  :: nt, q, q0, nnd, ntm, ntp
        integer  :: ft, ii, ib, cs0, cs1, nfa

        vfp_lo = sz - GN + visc_face_reach - 1
        vfp_hi = ez + GN - visc_face_reach

        call StartProfRange("viscous_flux_3D_conservative_zta_curv")

        !$acc parallel vector_length(vec_size) default(present) private(fc,nd)
        !$acc cache(fc,nd)
        !$acc loop gang collapse(2)
        do   j = sy,ey
          do  ib = sx,ex,ilanes

            do cs0 = sz,ez,chunkSize
              cs1 = min(cs0+chunkSize-1, ez)
              nfa = (cs1-cs0+1) + 2*visc_div_reach - 1
              nnd = nfa + 2*visc_face_reach - 1
              q0  = cs0 - visc_div_reach - visc_face_reach + 1

              !$acc loop vector collapse(2)
              do nt = 0,nnd-1
                do ii = 0,ilanes-1
                  i = ib + ii
                  q = q0 + nt
                  if(i <= ex .and. q >= lbz .and. q <= ubz) then
                    duc = 0.0_rp
                    dvc = 0.0_rp
                    dwc = 0.0_rp
                    due = 0.0_rp
                    dve = 0.0_rp
                    dwe = 0.0_rp
                    !$acc loop seq
                    do n = 1,3
                      cn  = central_1(n)
                      duc = duc + cn*(U(i+n,j,q) - U(i-n,j,q))
                      dvc = dvc + cn*(V(i+n,j,q) - V(i-n,j,q))
                      dwc = dwc + cn*(W(i+n,j,q) - W(i-n,j,q))
                      due = due + cn*(U(i,j+n,q) - U(i,j-n,q))
                      dve = dve + cn*(V(i,j+n,q) - V(i,j-n,q))
                      dwe = dwe + cn*(W(i,j+n,q) - W(i,j-n,q))
                    enddo
                    nd(ii,nt,1) = duc
                    nd(ii,nt,2) = dvc
                    nd(ii,nt,3) = dwc
                    nd(ii,nt,4) = due
                    nd(ii,nt,5) = dve
                    nd(ii,nt,6) = dwe
                    nd(ii,nt,7) = 1.0_rp/T(i,j,q)
                  endif
                enddo
              enddo

              !$acc loop vector collapse(2)
              do ft = 0,nfa-1
                do ii = 0,ilanes-1
                  i = ib + ii
                  if(i <= ex) then
                    p = cs0 - visc_div_reach + ft

                    csi_x = dxi_dx (i,j)
                    csi_y = dxi_dy (i,j)
                    eta_x = deta_dx(i,j)
                    eta_y = deta_dy(i,j)

                    if(p < vfp_lo .or. p > vfp_hi) then

                      if(p < vfp_lo) then
                        flx1 = vfclo_z(1,i,j)
                        flx2 = vfclo_z(2,i,j)
                        flx3 = vfclo_z(3,i,j)
                        gflx = vfclo_z(4,i,j)
                        pflx = vfclo_z(8,i,j)
                        uh   = vfclo_z(5,i,j)
                        vh   = vfclo_z(6,i,j)
                        wh   = vfclo_z(7,i,j)
                      else
                        flx1 = vfchi_z(1,i,j)
                        flx2 = vfchi_z(2,i,j)
                        flx3 = vfchi_z(3,i,j)
                        gflx = vfchi_z(4,i,j)
                        pflx = vfchi_z(8,i,j)
                        uh   = vfchi_z(5,i,j)
                        vh   = vfchi_z(6,i,j)
                        wh   = vfchi_z(7,i,j)
                      endif

                    else

                      zetaf = fmz_zeta(p)
                      muf   = 0.0_rp
                      rti   = 0.0_rp
                      lmf   = 0.0_rp
                      uh    = 0.0_rp
                      vh    = 0.0_rp
                      wh    = 0.0_rp
                      duz   = 0.0_rp
                      dvz   = 0.0_rp
                      dwz   = 0.0_rp
                      dtz   = 0.0_rp
                      duc   = 0.0_rp
                      dvc   = 0.0_rp
                      dwc   = 0.0_rp
                      due   = 0.0_rp
                      dve   = 0.0_rp
                      dwe   = 0.0_rp

                      !$acc loop seq
                      do m = 1,visc_face_reach
                        cm  = visc_intp(m)
                        pm  = p+1-m
                        pp  = p+m
                        ntm = pm - q0
                        ntp = pp - q0

                        muf = muf   + cm*(VIS(i,j,pp) + VIS(i,j,pm))
                        lmf = lmf   + cm*(LMD(i,j,pp) + LMD(i,j,pm))
                        rti = rti + cm*(nd(ii,ntp,7) + nd(ii,ntm,7))

                        uh = uh + cm*(U(i,j,pp) + U(i,j,pm))
                        vh = vh + cm*(V(i,j,pp) + V(i,j,pm))
                        wh = wh + cm*(W(i,j,pp) + W(i,j,pm))

                        duz = duz + visc_stag(m)*(U(i,j,pp) - U(i,j,pm))
                        dvz = dvz + visc_stag(m)*(V(i,j,pp) - V(i,j,pm))
                        dwz = dwz + visc_stag(m)*(W(i,j,pp) - W(i,j,pm))
                        dtz = dtz + visc_stag(m)*(T(i,j,pp) - T(i,j,pm))

                        duc = duc + cm*(nd(ii,ntp,1) + nd(ii,ntm,1))
                        dvc = dvc + cm*(nd(ii,ntp,2) + nd(ii,ntm,2))
                        dwc = dwc + cm*(nd(ii,ntp,3) + nd(ii,ntm,3))
                        due = due + cm*(nd(ii,ntp,4) + nd(ii,ntm,4))
                        dve = dve + cm*(nd(ii,ntp,5) + nd(ii,ntm,5))
                        dwe = dwe + cm*(nd(ii,ntp,6) + nd(ii,ntm,6))
                      enddo

                      ux = csi_x*duc + eta_x*due
                      uz = zetaf*duz
                      vy = csi_y*dvc + eta_y*dve
                      vz = zetaf*dvz
                      wx = csi_x*dwc + eta_x*dwe
                      wy = csi_y*dwc + eta_y*dwe
                      wz = zetaf*dwz
                      tz = zetaf*dtz

                      div3 = two3*(ux + vy + wz)
                      t13  = muf*(uz + wx)
                      t23  = muf*(vz + wy)
                      t33  = muf*(2.0_rp*wz - div3)

                      flx1 = t13
                      flx2 = t23
                      flx3 = t33
                      qflx = lmf*tz
                      gflx = qflx*rti
                      pflx = gflx*dtz
                    endif

                    fc(ii,ft,1) = flx1
                    fc(ii,ft,2) = flx2
                    fc(ii,ft,3) = flx3
                    fc(ii,ft,4) = gflx
                    fc(ii,ft,5) = uh
                    fc(ii,ft,6) = vh
                    fc(ii,ft,7) = wh
                    fc(ii,ft,8) = pflx
                  endif
                enddo
              enddo

              !$acc loop vector collapse(2)
              do k = cs0,cs1
                do ii = 0,ilanes-1
                  i = ib + ii
                  if(i <= ex) then

                    zeta_z = ztastep_i(k)*iz_zta(k)
                    uu     = U(i,j,k)
                    vv     = V(i,j,k)
                    ww     = W(i,j,k)
                    tnd    = T(i,j,k)
                    dfl1   = 0.0_rp
                    dfl2   = 0.0_rp
                    dfl3   = 0.0_rp
                    dhe    = 0.0_rp
                    dpr    = 0.0_rp
                    dph    = 0.0_rp

                    !$acc loop seq
                    do l = 1,visc_div_reach
                      bl = visc_divs(l)
                      cl = visc_intp(l)
                      !$acc loop seq
                      do s = 0,1
                        if(s == 1) then
                          p   =  k+l-1
                          fsg =  bl
                        else
                          p   =  k-l
                          fsg = -bl
                        endif
                        ft = p - cs0 + visc_div_reach

                        flx1 = fc(ii,ft,1)
                        flx2 = fc(ii,ft,2)
                        flx3 = fc(ii,ft,3)
                        gflx = fc(ii,ft,4)
                        pflx = fc(ii,ft,8)
                        uh   = fc(ii,ft,5)
                        vh   = fc(ii,ft,6)
                        wh   = fc(ii,ft,7)

                        dfl1 = dfl1 + fsg*flx1
                        dfl2 = dfl2 + fsg*flx2
                        dfl3 = dfl3 + fsg*flx3
                        dhe  = dhe  + fsg*gflx
                        dpr  = dpr  + cl*pflx
                        dph  = dph  + fsg*(flx1*(uh-uu) + flx2*(vh-vv) + flx3*(wh-ww))
                      enddo
                    enddo

                    rhs_u = zeta_z*dfl1
                    rhs_v = zeta_z*dfl2
                    rhs_w = zeta_z*dfl3

                    RHS(i,j,k,2) = RHS(i,j,k,2) + mu_inf*rhs_u
                    RHS(i,j,k,3) = RHS(i,j,k,3) + mu_inf*rhs_v
                    RHS(i,j,k,4) = RHS(i,j,k,4) + mu_inf*rhs_w
                    RHS(i,j,k,5) = RHS(i,j,k,5) &
                                 + mu_inf*(uu*rhs_u + vv*rhs_v + ww*rhs_w + zeta_z*(dph + tnd*dhe + dpr))

                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine viscous_flux_3D_conservative_zta_curv

subroutine viscous_flux_3D_conservative_zta_cart
        implicit none

        real(rp), parameter :: two3 = 2.0_rp/3.0_rp

        integer , parameter :: ilanes    = 4
        integer , parameter :: svec      = 32
        integer , parameter :: vec_size  = ilanes*svec
        integer , parameter :: chunkSize = svec - 2*visc_div_reach + 1
        integer , parameter :: nnd_max   = svec + 2*visc_face_reach - 1

        real(rp), dimension(0:ilanes-1,0:svec-1,1:8)    :: fc
        real(rp), dimension(0:ilanes-1,0:nnd_max-1,1:5) :: nd

        real(rp) :: csi_x, eta_y, zeta_z, zetaf
        real(rp) :: muf, lmf, bl, cm, cn, fsg
        real(rp) :: duc, dvc, dwc
        real(rp) :: due, dve, dwe
        real(rp) :: duz, dvz, dwz, dtz
        real(rp) :: uh, vh, wh
        real(rp) :: ux, uz, vy, vz, wx, wy, wz, tz
        real(rp) :: t13, t23, t33, div3
        real(rp) :: flx1, flx2, flx3, qflx, gflx, pflx, rti
        real(rp) :: dfl1, dfl2, dfl3, dhe, dph, dpr, cl, tnd
        real(rp) :: rhs_u, rhs_v, rhs_w, uu, vv, ww
        integer  :: i, j, k, l, m, n, s, p, pm, pp
        integer  :: vfp_lo, vfp_hi
        integer  :: nt, q, q0, nnd, ntm, ntp
        integer  :: ft, ii, ib, cs0, cs1, nfa

        vfp_lo = sz - GN + visc_face_reach - 1
        vfp_hi = ez + GN - visc_face_reach

        call StartProfRange("viscous_flux_3D_conservative_zta_cart")

        !$acc parallel vector_length(vec_size) default(present) private(fc,nd)
        !$acc cache(fc,nd)
        !$acc loop gang collapse(2)
        do   j = sy,ey
          do  ib = sx,ex,ilanes

            do cs0 = sz,ez,chunkSize
              cs1 = min(cs0+chunkSize-1, ez)
              nfa = (cs1-cs0+1) + 2*visc_div_reach - 1
              nnd = nfa + 2*visc_face_reach - 1
              q0  = cs0 - visc_div_reach - visc_face_reach + 1

              !$acc loop vector collapse(2)
              do nt = 0,nnd-1
                do ii = 0,ilanes-1
                  i = ib + ii
                  q = q0 + nt
                  if(i <= ex .and. q >= lbz .and. q <= ubz) then
                    duc = 0.0_rp
                    dwc = 0.0_rp
                    dve = 0.0_rp
                    dwe = 0.0_rp
                    !$acc loop seq
                    do n = 1,3
                      cn  = central_1(n)
                      duc = duc + cn*(U(i+n,j,q) - U(i-n,j,q))
                      dwc = dwc + cn*(W(i+n,j,q) - W(i-n,j,q))
                      dve = dve + cn*(V(i,j+n,q) - V(i,j-n,q))
                      dwe = dwe + cn*(W(i,j+n,q) - W(i,j-n,q))
                    enddo
                    nd(ii,nt,1) = duc
                    nd(ii,nt,2) = dwc
                    nd(ii,nt,3) = dve
                    nd(ii,nt,4) = dwe
                    nd(ii,nt,5) = 1.0_rp/T(i,j,q)
                  endif
                enddo
              enddo

              !$acc loop vector collapse(2)
              do ft = 0,nfa-1
                do ii = 0,ilanes-1
                  i = ib + ii
                  if(i <= ex) then
                    p = cs0 - visc_div_reach + ft

                    csi_x = dxi_dx (i,j)
                    eta_y = deta_dy(i,j)

                    if(p < vfp_lo .or. p > vfp_hi) then

                      if(p < vfp_lo) then
                        flx1 = vfclo_z(1,i,j)
                        flx2 = vfclo_z(2,i,j)
                        flx3 = vfclo_z(3,i,j)
                        gflx = vfclo_z(4,i,j)
                        pflx = vfclo_z(8,i,j)
                        uh   = vfclo_z(5,i,j)
                        vh   = vfclo_z(6,i,j)
                        wh   = vfclo_z(7,i,j)
                      else
                        flx1 = vfchi_z(1,i,j)
                        flx2 = vfchi_z(2,i,j)
                        flx3 = vfchi_z(3,i,j)
                        gflx = vfchi_z(4,i,j)
                        pflx = vfchi_z(8,i,j)
                        uh   = vfchi_z(5,i,j)
                        vh   = vfchi_z(6,i,j)
                        wh   = vfchi_z(7,i,j)
                      endif

                    else

                      zetaf = fmz_zeta(p)
                      muf   = 0.0_rp
                      rti   = 0.0_rp
                      lmf   = 0.0_rp
                      uh    = 0.0_rp
                      vh    = 0.0_rp
                      wh    = 0.0_rp
                      duz   = 0.0_rp
                      dvz   = 0.0_rp
                      dwz   = 0.0_rp
                      dtz   = 0.0_rp
                      duc   = 0.0_rp
                      dvc   = 0.0_rp
                      dwc   = 0.0_rp
                      due   = 0.0_rp
                      dve   = 0.0_rp
                      dwe   = 0.0_rp

                      !$acc loop seq
                      do m = 1,visc_face_reach
                        cm  = visc_intp(m)
                        pm  = p+1-m
                        pp  = p+m
                        ntm = pm - q0
                        ntp = pp - q0

                        muf = muf   + cm*(VIS(i,j,pp) + VIS(i,j,pm))
                        lmf = lmf   + cm*(LMD(i,j,pp) + LMD(i,j,pm))
                        rti = rti + cm*(nd(ii,ntp,5) + nd(ii,ntm,5))

                        uh = uh + cm*(U(i,j,pp) + U(i,j,pm))
                        vh = vh + cm*(V(i,j,pp) + V(i,j,pm))
                        wh = wh + cm*(W(i,j,pp) + W(i,j,pm))

                        duz = duz + visc_stag(m)*(U(i,j,pp) - U(i,j,pm))
                        dvz = dvz + visc_stag(m)*(V(i,j,pp) - V(i,j,pm))
                        dwz = dwz + visc_stag(m)*(W(i,j,pp) - W(i,j,pm))
                        dtz = dtz + visc_stag(m)*(T(i,j,pp) - T(i,j,pm))

                        duc = duc + cm*(nd(ii,ntp,1) + nd(ii,ntm,1))
                        dwc = dwc + cm*(nd(ii,ntp,2) + nd(ii,ntm,2))
                        dve = dve + cm*(nd(ii,ntp,3) + nd(ii,ntm,3))
                        dwe = dwe + cm*(nd(ii,ntp,4) + nd(ii,ntm,4))
                      enddo

                      ux = csi_x*duc
                      uz = zetaf*duz
                      vy = eta_y*dve
                      vz = zetaf*dvz
                      wx = csi_x*dwc
                      wy = eta_y*dwe
                      wz = zetaf*dwz
                      tz = zetaf*dtz

                      div3 = two3*(ux + vy + wz)
                      t13  = muf*(uz + wx)
                      t23  = muf*(vz + wy)
                      t33  = muf*(2.0_rp*wz - div3)

                      flx1 = t13
                      flx2 = t23
                      flx3 = t33
                      qflx = lmf*tz
                      gflx = qflx*rti
                      pflx = gflx*dtz
                    endif

                    fc(ii,ft,1) = flx1
                    fc(ii,ft,2) = flx2
                    fc(ii,ft,3) = flx3
                    fc(ii,ft,4) = gflx
                    fc(ii,ft,5) = uh
                    fc(ii,ft,6) = vh
                    fc(ii,ft,7) = wh
                    fc(ii,ft,8) = pflx
                  endif
                enddo
              enddo

              !$acc loop vector collapse(2)
              do k = cs0,cs1
                do ii = 0,ilanes-1
                  i = ib + ii
                  if(i <= ex) then

                    zeta_z = ztastep_i(k)*iz_zta(k)
                    uu     = U(i,j,k)
                    vv     = V(i,j,k)
                    ww     = W(i,j,k)
                    tnd    = T(i,j,k)
                    dfl1   = 0.0_rp
                    dfl2   = 0.0_rp
                    dfl3   = 0.0_rp
                    dhe    = 0.0_rp
                    dpr    = 0.0_rp
                    dph    = 0.0_rp

                    !$acc loop seq
                    do l = 1,visc_div_reach
                      bl = visc_divs(l)
                      cl = visc_intp(l)
                      !$acc loop seq
                      do s = 0,1
                        if(s == 1) then
                          p   =  k+l-1
                          fsg =  bl
                        else
                          p   =  k-l
                          fsg = -bl
                        endif
                        ft = p - cs0 + visc_div_reach

                        flx1 = fc(ii,ft,1)
                        flx2 = fc(ii,ft,2)
                        flx3 = fc(ii,ft,3)
                        gflx = fc(ii,ft,4)
                        pflx = fc(ii,ft,8)
                        uh   = fc(ii,ft,5)
                        vh   = fc(ii,ft,6)
                        wh   = fc(ii,ft,7)

                        dfl1 = dfl1 + fsg*flx1
                        dfl2 = dfl2 + fsg*flx2
                        dfl3 = dfl3 + fsg*flx3
                        dhe  = dhe  + fsg*gflx
                        dpr  = dpr  + cl*pflx
                        dph  = dph  + fsg*(flx1*(uh-uu) + flx2*(vh-vv) + flx3*(wh-ww))
                      enddo
                    enddo

                    rhs_u = zeta_z*dfl1
                    rhs_v = zeta_z*dfl2
                    rhs_w = zeta_z*dfl3

                    RHS(i,j,k,2) = RHS(i,j,k,2) + mu_inf*rhs_u
                    RHS(i,j,k,3) = RHS(i,j,k,3) + mu_inf*rhs_v
                    RHS(i,j,k,4) = RHS(i,j,k,4) + mu_inf*rhs_w
                    RHS(i,j,k,5) = RHS(i,j,k,5) &
                                 + mu_inf*(uu*rhs_u + vv*rhs_v + ww*rhs_w + zeta_z*(dph + tnd*dhe + dpr))

                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine viscous_flux_3D_conservative_zta_cart

subroutine viscous_face_plane_xi_curv(p, bnd, buf)
        implicit none

        integer , intent(in)                                :: p
        logical , intent(in)                                :: bnd
        real(rp), dimension(1:8,sy:ey,sz:ez), intent(inout) :: buf

        real(rp), parameter :: two3 = 2.0_rp/3.0_rp

        real(rp) :: csi_x, csi_y, eta_x, eta_y, cof1, cof2, zeta_z
        real(rp) :: duc, dvc, dwc, dtc
        real(rp) :: due, dve, dwe, dte
        real(rp) :: duz, dvz, dwz
        real(rp) :: uh, vh, wh
        real(rp) :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tx, ty
        real(rp) :: t11, t12, t13, t22, t23, div3
        real(rp) :: flx1, flx2, flx3, qflx, gflx, pflx, rti
        real(rp) :: nd1m, nd2m, nd3m, nd4m, nd5m, nd6m, nd7m
        real(rp) :: nd1p, nd2p, nd3p, nd4p, nd5p, nd6p, nd7p
        real(rp) :: muf, lmf, cm, cn, cs
        integer  :: k, j, m, n, pm, pp

        !$acc parallel default(present) present(buf)
        !$acc loop gang, vector collapse(2)
        do  k = sz,ez
          do j = sy,ey

            zeta_z = ztastep_i(k)*iz_zta(k)
            if(bnd) then
              csi_x = 0.0_rp
              csi_y = 0.0_rp
              eta_x = 0.0_rp
              eta_y = 0.0_rp
              cof1  = 0.0_rp
              cof2  = 0.0_rp
              !$acc loop seq
              do m = 1,visc_face_reach
                cm    = visc_intp_b(m)
                pm    = max(p+1-m, lbx)
                pp    = min(p+m  , ubx)
                csi_x = csi_x + cm*(dxi_dx (pp,j) + dxi_dx (pm,j))
                csi_y = csi_y + cm*(dxi_dy (pp,j) + dxi_dy (pm,j))
                eta_x = eta_x + cm*(deta_dx(pp,j) + deta_dx(pm,j))
                eta_y = eta_y + cm*(deta_dy(pp,j) + deta_dy(pm,j))
                cof1  = cof1  + cm*(dy_deta(pp,j) + dy_deta(pm,j))
                cof2  = cof2  - cm*(dx_deta(pp,j) + dx_deta(pm,j))
              enddo
            else
              csi_x = fmx_csix(p,j)
              csi_y = fmx_csiy(p,j)
              eta_x = fmx_etax(p,j)
              eta_y = fmx_etay(p,j)
              cof1  = fmx_cof1(p,j)
              cof2  = fmx_cof2(p,j)
            endif

            muf = 0.0_rp
            rti = 0.0_rp
            lmf = 0.0_rp
            uh  = 0.0_rp
            vh  = 0.0_rp
            wh  = 0.0_rp
            duc = 0.0_rp
            dvc = 0.0_rp
            dwc = 0.0_rp
            dtc = 0.0_rp
            due = 0.0_rp
            dve = 0.0_rp
            dwe = 0.0_rp
            dte = 0.0_rp
            duz = 0.0_rp
            dvz = 0.0_rp
            dwz = 0.0_rp

            !$acc loop seq
            do m = 1,visc_face_reach
              if(bnd) then
                cm = visc_intp_b(m)
                cs = visc_stag_b(m)
              else
                cm = visc_intp  (m)
                cs = visc_stag  (m)
              endif
              pm = max(p+1-m, lbx)
              pp = min(p+m  , ubx)

              muf = muf + cm*(VIS(pp,j,k) + VIS(pm,j,k))
              lmf = lmf + cm*(LMD(pp,j,k) + LMD(pm,j,k))
              rti = rti + cm*(1.0_rp/T(pp,j,k) + 1.0_rp/T(pm,j,k))

              uh = uh + cm*(U(pp,j,k) + U(pm,j,k))
              vh = vh + cm*(V(pp,j,k) + V(pm,j,k))
              wh = wh + cm*(W(pp,j,k) + W(pm,j,k))

              duc = duc + cs*(U(pp,j,k) - U(pm,j,k))
              dvc = dvc + cs*(V(pp,j,k) - V(pm,j,k))
              dwc = dwc + cs*(W(pp,j,k) - W(pm,j,k))
              dtc = dtc + cs*(T(pp,j,k) - T(pm,j,k))

              nd1m = 0.0_rp
              nd2m = 0.0_rp
              nd3m = 0.0_rp
              nd4m = 0.0_rp
              nd5m = 0.0_rp
              nd6m = 0.0_rp
              nd7m = 0.0_rp
              nd1p = 0.0_rp
              nd2p = 0.0_rp
              nd3p = 0.0_rp
              nd4p = 0.0_rp
              nd5p = 0.0_rp
              nd6p = 0.0_rp
              nd7p = 0.0_rp
              !$acc loop seq
              do n = 1,3
                cn   = central_1(n)
                nd1m = nd1m + cn*(U(pm,j+n,k) - U(pm,j-n,k))
                nd1p = nd1p + cn*(U(pp,j+n,k) - U(pp,j-n,k))
                nd2m = nd2m + cn*(V(pm,j+n,k) - V(pm,j-n,k))
                nd2p = nd2p + cn*(V(pp,j+n,k) - V(pp,j-n,k))
                nd3m = nd3m + cn*(W(pm,j+n,k) - W(pm,j-n,k))
                nd3p = nd3p + cn*(W(pp,j+n,k) - W(pp,j-n,k))
                nd4m = nd4m + cn*(T(pm,j+n,k) - T(pm,j-n,k))
                nd4p = nd4p + cn*(T(pp,j+n,k) - T(pp,j-n,k))
                nd5m = nd5m + cn*(U(pm,j,k+n) - U(pm,j,k-n))
                nd5p = nd5p + cn*(U(pp,j,k+n) - U(pp,j,k-n))
                nd6m = nd6m + cn*(V(pm,j,k+n) - V(pm,j,k-n))
                nd6p = nd6p + cn*(V(pp,j,k+n) - V(pp,j,k-n))
                nd7m = nd7m + cn*(W(pm,j,k+n) - W(pm,j,k-n))
                nd7p = nd7p + cn*(W(pp,j,k+n) - W(pp,j,k-n))
              enddo
              due = due + cm*(nd1p + nd1m)
              dve = dve + cm*(nd2p + nd2m)
              dwe = dwe + cm*(nd3p + nd3m)
              dte = dte + cm*(nd4p + nd4m)
              duz = duz + cm*(nd5p + nd5m)
              dvz = dvz + cm*(nd6p + nd6m)
              dwz = dwz + cm*(nd7p + nd7m)
            enddo

            ux = csi_x*duc + eta_x*due
            uy = csi_y*duc + eta_y*due
            uz = zeta_z*duz
            vx = csi_x*dvc + eta_x*dve
            vy = csi_y*dvc + eta_y*dve
            vz = zeta_z*dvz
            wx = csi_x*dwc + eta_x*dwe
            wy = csi_y*dwc + eta_y*dwe
            wz = zeta_z*dwz
            tx = csi_x*dtc + eta_x*dte
            ty = csi_y*dtc + eta_y*dte

            div3 = two3*(ux + vy + wz)
            t11  = muf*(2.0_rp*ux - div3)
            t12  = muf*(uy + vx)
            t13  = muf*(uz + wx)
            t22  = muf*(2.0_rp*vy - div3)
            t23  = muf*(vz + wy)

            flx1 = cof1*t11 + cof2*t12
            flx2 = cof1*t12 + cof2*t22
            flx3 = cof1*t13 + cof2*t23
            qflx = lmf*(cof1*tx + cof2*ty)
            gflx = qflx*rti
            pflx = gflx*dtc

            buf(1,j,k) = flx1
            buf(2,j,k) = flx2
            buf(3,j,k) = flx3
            buf(4,j,k) = gflx
            buf(5,j,k) = uh
            buf(6,j,k) = vh
            buf(7,j,k) = wh
            buf(8,j,k) = pflx

          enddo
        enddo
        !$acc end parallel

        return
end subroutine viscous_face_plane_xi_curv

subroutine viscous_face_plane_xi_cart(p, bnd, buf)
        implicit none

        integer , intent(in)                                :: p
        logical , intent(in)                                :: bnd
        real(rp), dimension(1:8,sy:ey,sz:ez), intent(inout) :: buf

        real(rp), parameter :: two3 = 2.0_rp/3.0_rp

        real(rp) :: csi_x, eta_y, cof1, zeta_z
        real(rp) :: duc, dvc, dwc, dtc
        real(rp) :: due, dve, dwe, dte
        real(rp) :: duz, dvz, dwz
        real(rp) :: uh, vh, wh
        real(rp) :: ux, uy, uz, vx, vy, wx, wz, tx
        real(rp) :: t11, t12, t13, div3
        real(rp) :: flx1, flx2, flx3, qflx, gflx, pflx, rti
        real(rp) :: nd1m, nd2m, nd3m, nd4m, nd1p, nd2p, nd3p, nd4p
        real(rp) :: muf, lmf, cm, cn, cs
        integer  :: k, j, m, n, pm, pp

        !$acc parallel default(present) present(buf)
        !$acc loop gang, vector collapse(2)
        do  k = sz,ez
          do j = sy,ey

            zeta_z = ztastep_i(k)*iz_zta(k)
            if(bnd) then
              csi_x = 0.0_rp
              eta_y = 0.0_rp
              cof1  = 0.0_rp
              !$acc loop seq
              do m = 1,visc_face_reach
                cm    = visc_intp_b(m)
                pm    = max(p+1-m, lbx)
                pp    = min(p+m  , ubx)
                csi_x = csi_x + cm*(dxi_dx (pp,j) + dxi_dx (pm,j))
                eta_y = eta_y + cm*(deta_dy(pp,j) + deta_dy(pm,j))
                cof1  = cof1  + cm*(dy_deta(pp,j) + dy_deta(pm,j))
              enddo
            else
              csi_x = fmx_csix(p,j)
              eta_y = fmx_etay(p,j)
              cof1  = fmx_cof1(p,j)
            endif

            muf = 0.0_rp
            rti = 0.0_rp
            lmf = 0.0_rp
            uh  = 0.0_rp
            vh  = 0.0_rp
            wh  = 0.0_rp
            duc = 0.0_rp
            dvc = 0.0_rp
            dwc = 0.0_rp
            dtc = 0.0_rp
            due = 0.0_rp
            dve = 0.0_rp
            dwe = 0.0_rp
            dte = 0.0_rp
            duz = 0.0_rp
            dvz = 0.0_rp
            dwz = 0.0_rp

            !$acc loop seq
            do m = 1,visc_face_reach
              if(bnd) then
                cm = visc_intp_b(m)
                cs = visc_stag_b(m)
              else
                cm = visc_intp  (m)
                cs = visc_stag  (m)
              endif
              pm = max(p+1-m, lbx)
              pp = min(p+m  , ubx)

              muf = muf + cm*(VIS(pp,j,k) + VIS(pm,j,k))
              lmf = lmf + cm*(LMD(pp,j,k) + LMD(pm,j,k))
              rti = rti + cm*(1.0_rp/T(pp,j,k) + 1.0_rp/T(pm,j,k))

              uh = uh + cm*(U(pp,j,k) + U(pm,j,k))
              vh = vh + cm*(V(pp,j,k) + V(pm,j,k))
              wh = wh + cm*(W(pp,j,k) + W(pm,j,k))

              duc = duc + cs*(U(pp,j,k) - U(pm,j,k))
              dvc = dvc + cs*(V(pp,j,k) - V(pm,j,k))
              dwc = dwc + cs*(W(pp,j,k) - W(pm,j,k))
              dtc = dtc + cs*(T(pp,j,k) - T(pm,j,k))

              nd1m = 0.0_rp
              nd2m = 0.0_rp
              nd3m = 0.0_rp
              nd4m = 0.0_rp
              nd1p = 0.0_rp
              nd2p = 0.0_rp
              nd3p = 0.0_rp
              nd4p = 0.0_rp
              !$acc loop seq
              do n = 1,3
                cn   = central_1(n)
                nd1m = nd1m + cn*(U(pm,j+n,k) - U(pm,j-n,k))
                nd1p = nd1p + cn*(U(pp,j+n,k) - U(pp,j-n,k))
                nd2m = nd2m + cn*(V(pm,j+n,k) - V(pm,j-n,k))
                nd2p = nd2p + cn*(V(pp,j+n,k) - V(pp,j-n,k))
                nd3m = nd3m + cn*(U(pm,j,k+n) - U(pm,j,k-n))
                nd3p = nd3p + cn*(U(pp,j,k+n) - U(pp,j,k-n))
                nd4m = nd4m + cn*(W(pm,j,k+n) - W(pm,j,k-n))
                nd4p = nd4p + cn*(W(pp,j,k+n) - W(pp,j,k-n))
              enddo
              due = due + cm*(nd1p + nd1m)
              dve = dve + cm*(nd2p + nd2m)
              duz = duz + cm*(nd3p + nd3m)
              dwz = dwz + cm*(nd4p + nd4m)
            enddo

            ux = csi_x*duc
            uy = eta_y*due
            uz = zeta_z*duz
            vx = csi_x*dvc
            vy = eta_y*dve
            wx = csi_x*dwc
            wz = zeta_z*dwz
            tx = csi_x*dtc

            div3 = two3*(ux + vy + wz)
            t11  = muf*(2.0_rp*ux - div3)
            t12  = muf*(uy + vx)
            t13  = muf*(uz + wx)

            flx1 = cof1*t11
            flx2 = cof1*t12
            flx3 = cof1*t13
            qflx = lmf*(cof1*tx)
            gflx = qflx*rti
            pflx = gflx*dtc

            buf(1,j,k) = flx1
            buf(2,j,k) = flx2
            buf(3,j,k) = flx3
            buf(4,j,k) = gflx
            buf(5,j,k) = uh
            buf(6,j,k) = vh
            buf(7,j,k) = wh
            buf(8,j,k) = pflx

          enddo
        enddo
        !$acc end parallel

        return
end subroutine viscous_face_plane_xi_cart

subroutine viscous_face_plane_eta_curv(p, bnd, buf)
        implicit none

        integer , intent(in)                                :: p
        logical , intent(in)                                :: bnd
        real(rp), dimension(1:8,sx:ex,sz:ez), intent(inout) :: buf

        real(rp), parameter :: two3 = 2.0_rp/3.0_rp

        real(rp) :: csi_x, csi_y, eta_x, eta_y, cof1, cof2, zeta_z
        real(rp) :: duc, dvc, dwc, dtc
        real(rp) :: due, dve, dwe, dte
        real(rp) :: duz, dvz, dwz
        real(rp) :: uh, vh, wh
        real(rp) :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tx, ty
        real(rp) :: t11, t12, t13, t22, t23, div3
        real(rp) :: flx1, flx2, flx3, qflx, gflx, pflx, rti
        real(rp) :: nd1m, nd2m, nd3m, nd4m, nd5m, nd6m, nd7m
        real(rp) :: nd1p, nd2p, nd3p, nd4p, nd5p, nd6p, nd7p
        real(rp) :: muf, lmf, cm, cn, cs
        integer  :: k, i, m, n, pm, pp

        !$acc parallel default(present) present(buf)
        !$acc loop gang, vector collapse(2)
        do  k = sz,ez
          do i = sx,ex

            zeta_z = ztastep_i(k)*iz_zta(k)
            if(bnd) then
              csi_x = 0.0_rp
              csi_y = 0.0_rp
              eta_x = 0.0_rp
              eta_y = 0.0_rp
              cof1  = 0.0_rp
              cof2  = 0.0_rp
              !$acc loop seq
              do m = 1,visc_face_reach
                cm    = visc_intp_b(m)
                pm    = max(p+1-m, lby)
                pp    = min(p+m  , uby)
                csi_x = csi_x + cm*(dxi_dx (i,pp) + dxi_dx (i,pm))
                csi_y = csi_y + cm*(dxi_dy (i,pp) + dxi_dy (i,pm))
                eta_x = eta_x + cm*(deta_dx(i,pp) + deta_dx(i,pm))
                eta_y = eta_y + cm*(deta_dy(i,pp) + deta_dy(i,pm))
                cof1  = cof1  - cm*(dy_dxi (i,pp) + dy_dxi (i,pm))
                cof2  = cof2  + cm*(dx_dxi (i,pp) + dx_dxi (i,pm))
              enddo
            else
              csi_x = fmy_csix(i,p)
              csi_y = fmy_csiy(i,p)
              eta_x = fmy_etax(i,p)
              eta_y = fmy_etay(i,p)
              cof1  = fmy_cof1(i,p)
              cof2  = fmy_cof2(i,p)
            endif

            muf = 0.0_rp
            rti = 0.0_rp
            lmf = 0.0_rp
            uh  = 0.0_rp
            vh  = 0.0_rp
            wh  = 0.0_rp
            duc = 0.0_rp
            dvc = 0.0_rp
            dwc = 0.0_rp
            dtc = 0.0_rp
            due = 0.0_rp
            dve = 0.0_rp
            dwe = 0.0_rp
            dte = 0.0_rp
            duz = 0.0_rp
            dvz = 0.0_rp
            dwz = 0.0_rp

            !$acc loop seq
            do m = 1,visc_face_reach
              if(bnd) then
                cm = visc_intp_b(m)
                cs = visc_stag_b(m)
              else
                cm = visc_intp  (m)
                cs = visc_stag  (m)
              endif
              pm = max(p+1-m, lby)
              pp = min(p+m  , uby)

              muf = muf + cm*(VIS(i,pp,k) + VIS(i,pm,k))
              lmf = lmf + cm*(LMD(i,pp,k) + LMD(i,pm,k))
              rti = rti + cm*(1.0_rp/T(i,pp,k) + 1.0_rp/T(i,pm,k))

              uh = uh + cm*(U(i,pp,k) + U(i,pm,k))
              vh = vh + cm*(V(i,pp,k) + V(i,pm,k))
              wh = wh + cm*(W(i,pp,k) + W(i,pm,k))

              due = due + cs*(U(i,pp,k) - U(i,pm,k))
              dve = dve + cs*(V(i,pp,k) - V(i,pm,k))
              dwe = dwe + cs*(W(i,pp,k) - W(i,pm,k))
              dte = dte + cs*(T(i,pp,k) - T(i,pm,k))

              nd1m = 0.0_rp
              nd2m = 0.0_rp
              nd3m = 0.0_rp
              nd4m = 0.0_rp
              nd5m = 0.0_rp
              nd6m = 0.0_rp
              nd7m = 0.0_rp
              nd1p = 0.0_rp
              nd2p = 0.0_rp
              nd3p = 0.0_rp
              nd4p = 0.0_rp
              nd5p = 0.0_rp
              nd6p = 0.0_rp
              nd7p = 0.0_rp
              !$acc loop seq
              do n = 1,3
                cn   = central_1(n)
                nd1m = nd1m + cn*(U(i+n,pm,k) - U(i-n,pm,k))
                nd1p = nd1p + cn*(U(i+n,pp,k) - U(i-n,pp,k))
                nd2m = nd2m + cn*(V(i+n,pm,k) - V(i-n,pm,k))
                nd2p = nd2p + cn*(V(i+n,pp,k) - V(i-n,pp,k))
                nd3m = nd3m + cn*(W(i+n,pm,k) - W(i-n,pm,k))
                nd3p = nd3p + cn*(W(i+n,pp,k) - W(i-n,pp,k))
                nd4m = nd4m + cn*(T(i+n,pm,k) - T(i-n,pm,k))
                nd4p = nd4p + cn*(T(i+n,pp,k) - T(i-n,pp,k))
                nd5m = nd5m + cn*(U(i,pm,k+n) - U(i,pm,k-n))
                nd5p = nd5p + cn*(U(i,pp,k+n) - U(i,pp,k-n))
                nd6m = nd6m + cn*(V(i,pm,k+n) - V(i,pm,k-n))
                nd6p = nd6p + cn*(V(i,pp,k+n) - V(i,pp,k-n))
                nd7m = nd7m + cn*(W(i,pm,k+n) - W(i,pm,k-n))
                nd7p = nd7p + cn*(W(i,pp,k+n) - W(i,pp,k-n))
              enddo
              duc = duc + cm*(nd1p + nd1m)
              dvc = dvc + cm*(nd2p + nd2m)
              dwc = dwc + cm*(nd3p + nd3m)
              dtc = dtc + cm*(nd4p + nd4m)
              duz = duz + cm*(nd5p + nd5m)
              dvz = dvz + cm*(nd6p + nd6m)
              dwz = dwz + cm*(nd7p + nd7m)
            enddo

            ux = csi_x*duc + eta_x*due
            uy = csi_y*duc + eta_y*due
            uz = zeta_z*duz
            vx = csi_x*dvc + eta_x*dve
            vy = csi_y*dvc + eta_y*dve
            vz = zeta_z*dvz
            wx = csi_x*dwc + eta_x*dwe
            wy = csi_y*dwc + eta_y*dwe
            wz = zeta_z*dwz
            tx = csi_x*dtc + eta_x*dte
            ty = csi_y*dtc + eta_y*dte

            div3 = two3*(ux + vy + wz)
            t11  = muf*(2.0_rp*ux - div3)
            t12  = muf*(uy + vx)
            t13  = muf*(uz + wx)
            t22  = muf*(2.0_rp*vy - div3)
            t23  = muf*(vz + wy)

            flx1 = cof1*t11 + cof2*t12
            flx2 = cof1*t12 + cof2*t22
            flx3 = cof1*t13 + cof2*t23
            qflx = lmf*(cof1*tx + cof2*ty)
            gflx = qflx*rti
            pflx = gflx*dte

            buf(1,i,k) = flx1
            buf(2,i,k) = flx2
            buf(3,i,k) = flx3
            buf(4,i,k) = gflx
            buf(5,i,k) = uh
            buf(6,i,k) = vh
            buf(7,i,k) = wh
            buf(8,i,k) = pflx

          enddo
        enddo
        !$acc end parallel

        return
end subroutine viscous_face_plane_eta_curv

subroutine viscous_face_plane_eta_cart(p, bnd, buf)
        implicit none

        integer , intent(in)                                :: p
        logical , intent(in)                                :: bnd
        real(rp), dimension(1:8,sx:ex,sz:ez), intent(inout) :: buf

        real(rp), parameter :: two3 = 2.0_rp/3.0_rp

        real(rp) :: csi_x, eta_y, cof2, zeta_z
        real(rp) :: duc, dvc, dwc, dtc
        real(rp) :: due, dve, dwe, dte
        real(rp) :: duz, dvz, dwz
        real(rp) :: uh, vh, wh
        real(rp) :: ux, uy, vx, vy, vz, wy, wz, ty
        real(rp) :: t12, t22, t23, div3
        real(rp) :: flx1, flx2, flx3, qflx, gflx, pflx, rti
        real(rp) :: nd1m, nd2m, nd3m, nd4m, nd1p, nd2p, nd3p, nd4p
        real(rp) :: muf, lmf, cm, cn, cs
        integer  :: k, i, m, n, pm, pp

        !$acc parallel default(present) present(buf)
        !$acc loop gang, vector collapse(2)
        do  k = sz,ez
          do i = sx,ex

            zeta_z = ztastep_i(k)*iz_zta(k)
            if(bnd) then
              csi_x = 0.0_rp
              eta_y = 0.0_rp
              cof2  = 0.0_rp
              !$acc loop seq
              do m = 1,visc_face_reach
                cm    = visc_intp_b(m)
                pm    = max(p+1-m, lby)
                pp    = min(p+m  , uby)
                csi_x = csi_x + cm*(dxi_dx (i,pp) + dxi_dx (i,pm))
                eta_y = eta_y + cm*(deta_dy(i,pp) + deta_dy(i,pm))
                cof2  = cof2  + cm*(dx_dxi (i,pp) + dx_dxi (i,pm))
              enddo
            else
              csi_x = fmy_csix(i,p)
              eta_y = fmy_etay(i,p)
              cof2  = fmy_cof2(i,p)
            endif

            muf = 0.0_rp
            rti = 0.0_rp
            lmf = 0.0_rp
            uh  = 0.0_rp
            vh  = 0.0_rp
            wh  = 0.0_rp
            duc = 0.0_rp
            dvc = 0.0_rp
            dwc = 0.0_rp
            dtc = 0.0_rp
            due = 0.0_rp
            dve = 0.0_rp
            dwe = 0.0_rp
            dte = 0.0_rp
            duz = 0.0_rp
            dvz = 0.0_rp
            dwz = 0.0_rp

            !$acc loop seq
            do m = 1,visc_face_reach
              if(bnd) then
                cm = visc_intp_b(m)
                cs = visc_stag_b(m)
              else
                cm = visc_intp  (m)
                cs = visc_stag  (m)
              endif
              pm = max(p+1-m, lby)
              pp = min(p+m  , uby)

              muf = muf + cm*(VIS(i,pp,k) + VIS(i,pm,k))
              lmf = lmf + cm*(LMD(i,pp,k) + LMD(i,pm,k))
              rti = rti + cm*(1.0_rp/T(i,pp,k) + 1.0_rp/T(i,pm,k))

              uh = uh + cm*(U(i,pp,k) + U(i,pm,k))
              vh = vh + cm*(V(i,pp,k) + V(i,pm,k))
              wh = wh + cm*(W(i,pp,k) + W(i,pm,k))

              due = due + cs*(U(i,pp,k) - U(i,pm,k))
              dve = dve + cs*(V(i,pp,k) - V(i,pm,k))
              dwe = dwe + cs*(W(i,pp,k) - W(i,pm,k))
              dte = dte + cs*(T(i,pp,k) - T(i,pm,k))

              nd1m = 0.0_rp
              nd2m = 0.0_rp
              nd3m = 0.0_rp
              nd4m = 0.0_rp
              nd1p = 0.0_rp
              nd2p = 0.0_rp
              nd3p = 0.0_rp
              nd4p = 0.0_rp
              !$acc loop seq
              do n = 1,3
                cn   = central_1(n)
                nd1m = nd1m + cn*(U(i+n,pm,k) - U(i-n,pm,k))
                nd1p = nd1p + cn*(U(i+n,pp,k) - U(i-n,pp,k))
                nd2m = nd2m + cn*(V(i+n,pm,k) - V(i-n,pm,k))
                nd2p = nd2p + cn*(V(i+n,pp,k) - V(i-n,pp,k))
                nd3m = nd3m + cn*(V(i,pm,k+n) - V(i,pm,k-n))
                nd3p = nd3p + cn*(V(i,pp,k+n) - V(i,pp,k-n))
                nd4m = nd4m + cn*(W(i,pm,k+n) - W(i,pm,k-n))
                nd4p = nd4p + cn*(W(i,pp,k+n) - W(i,pp,k-n))
              enddo
              duc = duc + cm*(nd1p + nd1m)
              dvc = dvc + cm*(nd2p + nd2m)
              dvz = dvz + cm*(nd3p + nd3m)
              dwz = dwz + cm*(nd4p + nd4m)
            enddo

            ux = csi_x*duc
            uy = eta_y*due
            vx = csi_x*dvc
            vy = eta_y*dve
            vz = zeta_z*dvz
            wy = eta_y*dwe
            wz = zeta_z*dwz
            ty = eta_y*dte

            div3 = two3*(ux + vy + wz)
            t12  = muf*(uy + vx)
            t22  = muf*(2.0_rp*vy - div3)
            t23  = muf*(vz + wy)

            flx1 = cof2*t12
            flx2 = cof2*t22
            flx3 = cof2*t23
            qflx = lmf*(cof2*ty)
            gflx = qflx*rti
            pflx = gflx*dte

            buf(1,i,k) = flx1
            buf(2,i,k) = flx2
            buf(3,i,k) = flx3
            buf(4,i,k) = gflx
            buf(5,i,k) = uh
            buf(6,i,k) = vh
            buf(7,i,k) = wh
            buf(8,i,k) = pflx

          enddo
        enddo
        !$acc end parallel

        return
end subroutine viscous_face_plane_eta_cart

subroutine viscous_face_plane_zta_curv(p, bnd, buf)
        implicit none

        integer , intent(in)                                :: p
        logical , intent(in)                                :: bnd
        real(rp), dimension(1:8,sx:ex,sy:ey), intent(inout) :: buf

        real(rp), parameter :: two3 = 2.0_rp/3.0_rp

        real(rp) :: csi_x, csi_y, eta_x, eta_y, zetaf
        real(rp) :: duc, dvc, dwc
        real(rp) :: due, dve, dwe
        real(rp) :: duz, dvz, dwz, dtz
        real(rp) :: uh, vh, wh
        real(rp) :: ux, uz, vy, vz, wx, wy, wz, tz
        real(rp) :: t13, t23, t33, div3
        real(rp) :: flx1, flx2, flx3, qflx, gflx, pflx, rti
        real(rp) :: nd1m, nd2m, nd3m, nd4m, nd5m, nd6m
        real(rp) :: nd1p, nd2p, nd3p, nd4p, nd5p, nd6p
        real(rp) :: muf, lmf, cm, cn, cs
        integer  :: j, i, m, n, pm, pp

        !$acc parallel default(present) present(buf)
        !$acc loop gang, vector collapse(2)
        do  j = sy,ey
          do i = sx,ex

            csi_x  = dxi_dx (i,j)
            csi_y  = dxi_dy (i,j)
            eta_x  = deta_dx(i,j)
            eta_y  = deta_dy(i,j)
            if(bnd) then
              zetaf = 0.0_rp
              !$acc loop seq
              do m = 1,visc_face_reach
                cm    = visc_intp_b(m)
                pm    = max(p+1-m, lbz)
                pp    = min(p+m  , ubz)
                zetaf = zetaf + cm*(ztastep_i(pp)*iz_zta(pp) + ztastep_i(pm)*iz_zta(pm))
              enddo
            else
              zetaf = fmz_zeta(p)
            endif

            muf = 0.0_rp
            rti = 0.0_rp
            lmf = 0.0_rp
            uh  = 0.0_rp
            vh  = 0.0_rp
            wh  = 0.0_rp
            duz = 0.0_rp
            dvz = 0.0_rp
            dwz = 0.0_rp
            dtz = 0.0_rp
            duc = 0.0_rp
            dvc = 0.0_rp
            dwc = 0.0_rp
            due = 0.0_rp
            dve = 0.0_rp
            dwe = 0.0_rp

            !$acc loop seq
            do m = 1,visc_face_reach
              if(bnd) then
                cm = visc_intp_b(m)
                cs = visc_stag_b(m)
              else
                cm = visc_intp  (m)
                cs = visc_stag  (m)
              endif
              pm = max(p+1-m, lbz)
              pp = min(p+m  , ubz)

              muf = muf   + cm*(VIS(i,j,pp) + VIS(i,j,pm))
              lmf = lmf   + cm*(LMD(i,j,pp) + LMD(i,j,pm))
              rti = rti + cm*(1.0_rp/T(i,j,pp) + 1.0_rp/T(i,j,pm))

              uh = uh + cm*(U(i,j,pp) + U(i,j,pm))
              vh = vh + cm*(V(i,j,pp) + V(i,j,pm))
              wh = wh + cm*(W(i,j,pp) + W(i,j,pm))

              duz = duz + cs*(U(i,j,pp) - U(i,j,pm))
              dvz = dvz + cs*(V(i,j,pp) - V(i,j,pm))
              dwz = dwz + cs*(W(i,j,pp) - W(i,j,pm))
              dtz = dtz + cs*(T(i,j,pp) - T(i,j,pm))

              nd1m = 0.0_rp
              nd2m = 0.0_rp
              nd3m = 0.0_rp
              nd4m = 0.0_rp
              nd5m = 0.0_rp
              nd6m = 0.0_rp
              nd1p = 0.0_rp
              nd2p = 0.0_rp
              nd3p = 0.0_rp
              nd4p = 0.0_rp
              nd5p = 0.0_rp
              nd6p = 0.0_rp
              !$acc loop seq
              do n = 1,3
                cn   = central_1(n)
                nd1m = nd1m + cn*(U(i+n,j,pm) - U(i-n,j,pm))
                nd1p = nd1p + cn*(U(i+n,j,pp) - U(i-n,j,pp))
                nd2m = nd2m + cn*(V(i+n,j,pm) - V(i-n,j,pm))
                nd2p = nd2p + cn*(V(i+n,j,pp) - V(i-n,j,pp))
                nd3m = nd3m + cn*(W(i+n,j,pm) - W(i-n,j,pm))
                nd3p = nd3p + cn*(W(i+n,j,pp) - W(i-n,j,pp))
                nd4m = nd4m + cn*(U(i,j+n,pm) - U(i,j-n,pm))
                nd4p = nd4p + cn*(U(i,j+n,pp) - U(i,j-n,pp))
                nd5m = nd5m + cn*(V(i,j+n,pm) - V(i,j-n,pm))
                nd5p = nd5p + cn*(V(i,j+n,pp) - V(i,j-n,pp))
                nd6m = nd6m + cn*(W(i,j+n,pm) - W(i,j-n,pm))
                nd6p = nd6p + cn*(W(i,j+n,pp) - W(i,j-n,pp))
              enddo
              duc = duc + cm*(nd1p + nd1m)
              dvc = dvc + cm*(nd2p + nd2m)
              dwc = dwc + cm*(nd3p + nd3m)
              due = due + cm*(nd4p + nd4m)
              dve = dve + cm*(nd5p + nd5m)
              dwe = dwe + cm*(nd6p + nd6m)
            enddo

            ux = csi_x*duc + eta_x*due
            uz = zetaf*duz
            vy = csi_y*dvc + eta_y*dve
            vz = zetaf*dvz
            wx = csi_x*dwc + eta_x*dwe
            wy = csi_y*dwc + eta_y*dwe
            wz = zetaf*dwz
            tz = zetaf*dtz

            div3 = two3*(ux + vy + wz)
            t13  = muf*(uz + wx)
            t23  = muf*(vz + wy)
            t33  = muf*(2.0_rp*wz - div3)

            flx1 = t13
            flx2 = t23
            flx3 = t33
            qflx = lmf*tz
            gflx = qflx*rti
            pflx = gflx*dtz

            buf(1,i,j) = flx1
            buf(2,i,j) = flx2
            buf(3,i,j) = flx3
            buf(4,i,j) = gflx
            buf(5,i,j) = uh
            buf(6,i,j) = vh
            buf(7,i,j) = wh
            buf(8,i,j) = pflx

          enddo
        enddo
        !$acc end parallel

        return
end subroutine viscous_face_plane_zta_curv

subroutine viscous_face_plane_zta_cart(p, bnd, buf)
        implicit none

        integer , intent(in)                                :: p
        logical , intent(in)                                :: bnd
        real(rp), dimension(1:8,sx:ex,sy:ey), intent(inout) :: buf

        real(rp), parameter :: two3 = 2.0_rp/3.0_rp

        real(rp) :: csi_x, eta_y, zetaf
        real(rp) :: duc, dvc, dwc
        real(rp) :: due, dve, dwe
        real(rp) :: duz, dvz, dwz, dtz
        real(rp) :: uh, vh, wh
        real(rp) :: ux, uz, vy, vz, wx, wy, wz, tz
        real(rp) :: t13, t23, t33, div3
        real(rp) :: flx1, flx2, flx3, qflx, gflx, pflx, rti
        real(rp) :: nd1m, nd2m, nd3m, nd4m, nd1p, nd2p, nd3p, nd4p
        real(rp) :: muf, lmf, cm, cn, cs
        integer  :: j, i, m, n, pm, pp

        !$acc parallel default(present) present(buf)
        !$acc loop gang, vector collapse(2)
        do  j = sy,ey
          do i = sx,ex

            csi_x  = dxi_dx (i,j)
            eta_y  = deta_dy(i,j)
            if(bnd) then
              zetaf = 0.0_rp
              !$acc loop seq
              do m = 1,visc_face_reach
                cm    = visc_intp_b(m)
                pm    = max(p+1-m, lbz)
                pp    = min(p+m  , ubz)
                zetaf = zetaf + cm*(ztastep_i(pp)*iz_zta(pp) + ztastep_i(pm)*iz_zta(pm))
              enddo
            else
              zetaf = fmz_zeta(p)
            endif

            muf = 0.0_rp
            rti = 0.0_rp
            lmf = 0.0_rp
            uh  = 0.0_rp
            vh  = 0.0_rp
            wh  = 0.0_rp
            duz = 0.0_rp
            dvz = 0.0_rp
            dwz = 0.0_rp
            dtz = 0.0_rp
            duc = 0.0_rp
            dvc = 0.0_rp
            dwc = 0.0_rp
            due = 0.0_rp
            dve = 0.0_rp
            dwe = 0.0_rp

            !$acc loop seq
            do m = 1,visc_face_reach
              if(bnd) then
                cm = visc_intp_b(m)
                cs = visc_stag_b(m)
              else
                cm = visc_intp  (m)
                cs = visc_stag  (m)
              endif
              pm = max(p+1-m, lbz)
              pp = min(p+m  , ubz)

              muf = muf   + cm*(VIS(i,j,pp) + VIS(i,j,pm))
              lmf = lmf   + cm*(LMD(i,j,pp) + LMD(i,j,pm))
              rti = rti + cm*(1.0_rp/T(i,j,pp) + 1.0_rp/T(i,j,pm))

              uh = uh + cm*(U(i,j,pp) + U(i,j,pm))
              vh = vh + cm*(V(i,j,pp) + V(i,j,pm))
              wh = wh + cm*(W(i,j,pp) + W(i,j,pm))

              duz = duz + cs*(U(i,j,pp) - U(i,j,pm))
              dvz = dvz + cs*(V(i,j,pp) - V(i,j,pm))
              dwz = dwz + cs*(W(i,j,pp) - W(i,j,pm))
              dtz = dtz + cs*(T(i,j,pp) - T(i,j,pm))

              nd1m = 0.0_rp
              nd2m = 0.0_rp
              nd3m = 0.0_rp
              nd4m = 0.0_rp
              nd1p = 0.0_rp
              nd2p = 0.0_rp
              nd3p = 0.0_rp
              nd4p = 0.0_rp
              !$acc loop seq
              do n = 1,3
                cn   = central_1(n)
                nd1m = nd1m + cn*(U(i+n,j,pm) - U(i-n,j,pm))
                nd1p = nd1p + cn*(U(i+n,j,pp) - U(i-n,j,pp))
                nd2m = nd2m + cn*(W(i+n,j,pm) - W(i-n,j,pm))
                nd2p = nd2p + cn*(W(i+n,j,pp) - W(i-n,j,pp))
                nd3m = nd3m + cn*(V(i,j+n,pm) - V(i,j-n,pm))
                nd3p = nd3p + cn*(V(i,j+n,pp) - V(i,j-n,pp))
                nd4m = nd4m + cn*(W(i,j+n,pm) - W(i,j-n,pm))
                nd4p = nd4p + cn*(W(i,j+n,pp) - W(i,j-n,pp))
              enddo
              duc = duc + cm*(nd1p + nd1m)
              dwc = dwc + cm*(nd2p + nd2m)
              dve = dve + cm*(nd3p + nd3m)
              dwe = dwe + cm*(nd4p + nd4m)
            enddo

            ux = csi_x*duc
            uz = zetaf*duz
            vy = eta_y*dve
            vz = zetaf*dvz
            wx = csi_x*dwc
            wy = eta_y*dwe
            wz = zetaf*dwz
            tz = zetaf*dtz

            div3 = two3*(ux + vy + wz)
            t13  = muf*(uz + wx)
            t23  = muf*(vz + wy)
            t33  = muf*(2.0_rp*wz - div3)

            flx1 = t13
            flx2 = t23
            flx3 = t33
            qflx = lmf*tz
            gflx = qflx*rti
            pflx = gflx*dtz

            buf(1,i,j) = flx1
            buf(2,i,j) = flx2
            buf(3,i,j) = flx3
            buf(4,i,j) = gflx
            buf(5,i,j) = uh
            buf(6,i,j) = vh
            buf(7,i,j) = wh
            buf(8,i,j) = pflx

          enddo
        enddo
        !$acc end parallel

        return
end subroutine viscous_face_plane_zta_cart

subroutine viscous_faces_conservative
        implicit none

        integer, parameter :: iW = 1, iE = 2, iS = 3, iN = 4, iB = 5, iF = 6

        integer :: plo, phi

        if(inplane_axis_aligned) then

            plo = sx - visc_div_reach
            phi = ex + visc_div_reach - 1
            call viscous_face_plane_xi_cart(sx+visc_div_reach-2, .false., vfcsl_x)
            call viscous_face_plane_xi_cart(ex-visc_div_reach+1, .false., vfcsh_x)
            if(my_neighbour(iW) < 0) call viscous_face_plane_xi_cart(plo, .true., vfclo_x)
            if(my_neighbour(iE) < 0) call viscous_face_plane_xi_cart(phi, .true., vfchi_x)

            plo = sy - visc_div_reach
            phi = ey + visc_div_reach - 1
            call viscous_face_plane_eta_cart(sy+visc_div_reach-2, .false., vfcsl_y)
            call viscous_face_plane_eta_cart(ey-visc_div_reach+1, .false., vfcsh_y)
            if(my_neighbour(iS) < 0) call viscous_face_plane_eta_cart(plo, .true., vfclo_y)
            if(my_neighbour(iN) < 0) call viscous_face_plane_eta_cart(phi, .true., vfchi_y)

            if(dims == 3) then
              plo = sz - visc_div_reach
              phi = ez + visc_div_reach - 1
              call viscous_face_plane_zta_cart(sz+visc_div_reach-2, .false., vfcsl_z)
              call viscous_face_plane_zta_cart(ez-visc_div_reach+1, .false., vfcsh_z)
              if(my_neighbour(iB) < 0) call viscous_face_plane_zta_cart(plo, .true., vfclo_z)
              if(my_neighbour(iF) < 0) call viscous_face_plane_zta_cart(phi, .true., vfchi_z)
            endif

        else

            plo = sx - visc_div_reach
            phi = ex + visc_div_reach - 1
            call viscous_face_plane_xi_curv(sx+visc_div_reach-2, .false., vfcsl_x)
            call viscous_face_plane_xi_curv(ex-visc_div_reach+1, .false., vfcsh_x)
            if(my_neighbour(iW) < 0) call viscous_face_plane_xi_curv(plo, .true., vfclo_x)
            if(my_neighbour(iE) < 0) call viscous_face_plane_xi_curv(phi, .true., vfchi_x)

            plo = sy - visc_div_reach
            phi = ey + visc_div_reach - 1
            call viscous_face_plane_eta_curv(sy+visc_div_reach-2, .false., vfcsl_y)
            call viscous_face_plane_eta_curv(ey-visc_div_reach+1, .false., vfcsh_y)
            if(my_neighbour(iS) < 0) call viscous_face_plane_eta_curv(plo, .true., vfclo_y)
            if(my_neighbour(iN) < 0) call viscous_face_plane_eta_curv(phi, .true., vfchi_y)

            if(dims == 3) then
              plo = sz - visc_div_reach
              phi = ez + visc_div_reach - 1
              call viscous_face_plane_zta_curv(sz+visc_div_reach-2, .false., vfcsl_z)
              call viscous_face_plane_zta_curv(ez-visc_div_reach+1, .false., vfcsh_z)
              if(my_neighbour(iB) < 0) call viscous_face_plane_zta_curv(plo, .true., vfclo_z)
              if(my_neighbour(iF) < 0) call viscous_face_plane_zta_curv(phi, .true., vfchi_z)
            endif

        endif

        call mpi_share_planes(mpi_comm_cart, my_neighbour, dims, &
                              vfcsl_x, vfcsh_x, vfclo_x, vfchi_x, &
                              vfcsl_y, vfcsh_y, vfclo_y, vfchi_y, &
                              vfcsl_z, vfcsh_z, vfclo_z, vfchi_z)

        return
end subroutine viscous_faces_conservative

subroutine divergence_3D_conservative
        implicit none

        real(rp) :: duc, due, dvc, dve, dwz, cn
        integer  :: i, j, k, n

        call StartProfRange("divergence_3D_conservative")

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(3)
        do   k = sz,ez
          do  j = sy,ey
            do i = sx,ex

              duc = 0.0_rp
              due = 0.0_rp
              dvc = 0.0_rp
              dve = 0.0_rp
              dwz = 0.0_rp

              !$acc loop seq
              do n = 1,3
                cn  = central_1(n)
                duc = duc + cn*(U(i+n,j,k) - U(i-n,j,k))
                due = due + cn*(U(i,j+n,k) - U(i,j-n,k))
                dvc = dvc + cn*(V(i+n,j,k) - V(i-n,j,k))
                dve = dve + cn*(V(i,j+n,k) - V(i,j-n,k))
                dwz = dwz + cn*(W(i,j,k+n) - W(i,j,k-n))
              enddo

              DIV(i,j,k) = dxi_dx(i,j)*duc + deta_dx(i,j)*due &
                + dxi_dy(i,j)*dvc + deta_dy(i,j)*dve &
                + ztastep_i(k)*iz_zta(k)*dwz

            enddo
          enddo
        enddo
        !$acc end parallel

        call EndProfRange

        return
end subroutine divergence_3D_conservative

end module viscous_conservative_module
