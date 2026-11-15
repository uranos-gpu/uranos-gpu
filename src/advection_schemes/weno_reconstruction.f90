module weno_reconstruction_module
use parameters_module
implicit none
private
public :: weno_reconstruction

contains

subroutine weno_reconstruction(weno_num,shock_recon,fp,fm,fhatp,fhatm)
        !$acc routine seq
        implicit none

        integer                                       , intent(in ) :: shock_recon, weno_num
        real(rp), dimension(5,-3:4)                   , intent(in ) :: fp, fm
        real(rp), dimension(5)                        , intent(out) :: fhatp, fhatm

        real(rp), parameter :: varepsilon = 1.0E-12_rp
        real(rp), parameter :: is1par = 13.0_rp/12.0_rp, is2par = 0.25_rp
        real(rp), parameter :: a1 = 10.5_rp, a2 = 3.5_rp, cr = 0.25_rp, csi = 1.0E-3_rp
        real(rp), parameter :: eps = 0.9_rp*cr/(1.0_rp-0.9_rp*cr)*csi**2
        real(rp), parameter :: ct_teno = 1.0E-7_rp

        integer  :: s, r
        real(rp) :: qp(0:3), qm(0:3), ISp(0:3), ISm(0:3), ap(0:3), am(0:3)
        real(rp) :: tp(0:1), tm(0:1), sp, sm, dsp, dsm, ctp, ctm, rmp, rmm
        real(rp) :: e0p, e1p, e2p, e0m, e1m, e2m
        integer  :: dp(0:3), dm(0:3)

        if(shock_recon == 1) then

        if    (weno_num == 2) then
        do s = 1, 5
           qp(0) = aweno3(0,0)*fp(s,-1) + aweno3(0,1)*fp(s,0)
           qm(0) = aweno3(0,0)*fm(s, 2) + aweno3(0,1)*fm(s,1)
           qp(1) = aweno3(1,0)*fp(s, 0) + aweno3(1,1)*fp(s,1)
           qm(1) = aweno3(1,0)*fm(s, 1) + aweno3(1,1)*fm(s,0)

           ISp(0) = (fp(s,-1)-fp(s,0))**2 ; ISm(0) = (fm(s,2)-fm(s,1))**2
           ISp(1) = (fp(s, 0)-fp(s,1))**2 ; ISm(1) = (fm(s,1)-fm(s,0))**2

           tp(0) = abs(ISp(0)-ISp(1))     ; tm(0) = abs(ISm(0)-ISm(1))
           do r = 0, 1
              ap(r) = cweno3(r)*(1.0_rp + tp(0)/(ISp(r) + varepsilon))
              am(r) = cweno3(r)*(1.0_rp + tm(0)/(ISm(r) + varepsilon))
           enddo

           sp = 1.0_rp/(ap(0) + ap(1)) ; sm = 1.0_rp/(am(0) + am(1))
           fhatp(s) = (ap(0)*sp)*qp(0) + (ap(1)*sp)*qp(1)
           fhatm(s) = (am(0)*sm)*qm(0) + (am(1)*sm)*qm(1)
        enddo

        elseif(weno_num == 3) then
        do s = 1, 5
           qp(0) = aweno5(0,0)*fp(s,-2) + aweno5(0,1)*fp(s,-1) + aweno5(0,2)*fp(s,0)
           qm(0) = aweno5(0,0)*fm(s, 3) + aweno5(0,1)*fm(s, 2) + aweno5(0,2)*fm(s,1)
           qp(1) = aweno5(1,0)*fp(s,-1) + aweno5(1,1)*fp(s, 0) + aweno5(1,2)*fp(s,1)
           qm(1) = aweno5(1,0)*fm(s, 2) + aweno5(1,1)*fm(s, 1) + aweno5(1,2)*fm(s,0)
           qp(2) = aweno5(2,0)*fp(s, 0) + aweno5(2,1)*fp(s, 1) + aweno5(2,2)*fp(s,2)
           qm(2) = aweno5(2,0)*fm(s, 1) + aweno5(2,1)*fm(s, 0) + aweno5(2,2)*fm(s,-1)

           ISp(0) = is1par*(fp(s,-2)-2.0_rp*fp(s,-1)+fp(s,0))**2 + is2par*(    fp(s,-2)-4.0_rp*fp(s,-1)+3.0_rp*fp(s,0))**2
           ISm(0) = is1par*(fm(s, 3)-2.0_rp*fm(s, 2)+fm(s,1))**2 + is2par*(    fm(s, 3)-4.0_rp*fm(s, 2)+3.0_rp*fm(s,1))**2
           ISp(1) = is1par*(fp(s,-1)-2.0_rp*fp(s, 0)+fp(s,1))**2 + is2par*(    fp(s,-1)                       -fp(s,1))**2
           ISm(1) = is1par*(fm(s, 2)-2.0_rp*fm(s, 1)+fm(s,0))**2 + is2par*(    fm(s, 2)                       -fm(s,0))**2
           ISp(2) = is1par*(fp(s, 0)-2.0_rp*fp(s, 1)+fp(s,2))**2 + is2par*(3.0_rp*fp(s, 0)-4.0_rp*fp(s, 1)    +fp(s,2))**2
           ISm(2) = is1par*(fm(s, 1)-2.0_rp*fm(s, 0)+fm(s,-1))**2 + is2par*(3.0_rp*fm(s, 1)-4.0_rp*fm(s, 0)    +fm(s,-1))**2

           tp(0) = abs(ISp(0)-ISp(2)) ; tm(0) = abs(ISm(0)-ISm(2))
           do r = 0, 2
              ap(r) = tp(0)/(ISp(r) + varepsilon) ; ap(r) = cweno5(r)*(1.0_rp + ap(r)*ap(r))
              am(r) = tm(0)/(ISm(r) + varepsilon) ; am(r) = cweno5(r)*(1.0_rp + am(r)*am(r))
           enddo

           sp = 1.0_rp/(ap(0) + ap(1) + ap(2)) ; sm = 1.0_rp/(am(0) + am(1) + am(2))
           fhatp(s) = (ap(0)*sp)*qp(0) + (ap(1)*sp)*qp(1) + (ap(2)*sp)*qp(2)
           fhatm(s) = (am(0)*sm)*qm(0) + (am(1)*sm)*qm(1) + (am(2)*sm)*qm(2)
        enddo

        elseif(weno_num == 4) then
        do s = 1, 5
           qp(0) = aweno7(0,0)*fp(s,-3) + aweno7(0,1)*fp(s,-2) + aweno7(0,2)*fp(s,-1) + aweno7(0,3)*fp(s,0)
           qm(0) = aweno7(0,0)*fm(s, 4) + aweno7(0,1)*fm(s, 3) + aweno7(0,2)*fm(s, 2) + aweno7(0,3)*fm(s, 1)
           qp(1) = aweno7(1,0)*fp(s,-2) + aweno7(1,1)*fp(s,-1) + aweno7(1,2)*fp(s, 0) + aweno7(1,3)*fp(s,1)
           qm(1) = aweno7(1,0)*fm(s, 3) + aweno7(1,1)*fm(s, 2) + aweno7(1,2)*fm(s, 1) + aweno7(1,3)*fm(s, 0)
           qp(2) = aweno7(2,0)*fp(s,-1) + aweno7(2,1)*fp(s, 0) + aweno7(2,2)*fp(s, 1) + aweno7(2,3)*fp(s,2)
           qm(2) = aweno7(2,0)*fm(s, 2) + aweno7(2,1)*fm(s, 1) + aweno7(2,2)*fm(s, 0) + aweno7(2,3)*fm(s,-1)
           qp(3) = aweno7(3,0)*fp(s, 0) + aweno7(3,1)*fp(s, 1) + aweno7(3,2)*fp(s, 2) + aweno7(3,3)*fp(s,3)
           qm(3) = aweno7(3,0)*fm(s, 1) + aweno7(3,1)*fm(s, 0) + aweno7(3,2)*fm(s,-1) + aweno7(3,3)*fm(s,-2)

           ISp(0) = fp(s,-3)*( 547.0_rp*fp(s,-3)- 3882.0_rp*fp(s,-2)+ 4642.0_rp*fp(s,-1)-1854.0_rp*fp(s,0)) + &
                    fp(s,-2)*(                    7043.0_rp*fp(s,-2)-17246.0_rp*fp(s,-1)+7042.0_rp*fp(s,0)) + &
                    fp(s,-1)*(                                       11003.0_rp*fp(s,-1)-9402.0_rp*fp(s,0)) + &
                    fp(s, 0)*(                                                           2107.0_rp*fp(s,0))
           ISm(0) = fm(s, 4)*( 547.0_rp*fm(s,4)- 3882.0_rp*fm(s,3)+ 4642.0_rp*fm(s, 2)-1854.0_rp*fm(s, 1)) + &
                    fm(s, 3)*(                   7043.0_rp*fm(s,3)-17246.0_rp*fm(s, 2)+7042.0_rp*fm(s, 1)) + &
                    fm(s, 2)*(                                     11003.0_rp*fm(s, 2)-9402.0_rp*fm(s, 1)) + &
                    fm(s, 1)*(                                                         2107.0_rp*fm(s, 1))
           ISp(1) = fp(s,-2)*( 267.0_rp*fp(s,-2)- 1642.0_rp*fp(s,-1)+ 1602.0_rp*fp(s, 0)- 494.0_rp*fp(s,1)) + &
                    fp(s,-1)*(                    2843.0_rp*fp(s,-1)- 5966.0_rp*fp(s, 0)+1922.0_rp*fp(s,1)) + &
                    fp(s, 0)*(                                        3443.0_rp*fp(s, 0)-2522.0_rp*fp(s,1)) + &
                    fp(s, 1)*(                                                            547.0_rp*fp(s,1))
           ISm(1) = fm(s, 3)*( 267.0_rp*fm(s,3)- 1642.0_rp*fm(s,2)+ 1602.0_rp*fm(s, 1)- 494.0_rp*fm(s, 0)) + &
                    fm(s, 2)*(                   2843.0_rp*fm(s,2)- 5966.0_rp*fm(s, 1)+1922.0_rp*fm(s, 0)) + &
                    fm(s, 1)*(                                      3443.0_rp*fm(s, 1)-2522.0_rp*fm(s, 0)) + &
                    fm(s, 0)*(                                                          547.0_rp*fm(s, 0))
           ISp(2) = fp(s,-1)*( 547.0_rp*fp(s,-1)- 2522.0_rp*fp(s, 0)+ 1922.0_rp*fp(s, 1)- 494.0_rp*fp(s,2)) + &
                    fp(s, 0)*(                    3443.0_rp*fp(s, 0)- 5966.0_rp*fp(s, 1)+1602.0_rp*fp(s,2)) + &
                    fp(s, 1)*(                                        2843.0_rp*fp(s, 1)-1642.0_rp*fp(s,2)) + &
                    fp(s, 2)*(                                                            267.0_rp*fp(s,2))
           ISm(2) = fm(s, 2)*( 547.0_rp*fm(s,2)- 2522.0_rp*fm(s,1)+ 1922.0_rp*fm(s, 0)- 494.0_rp*fm(s,-1)) + &
                    fm(s, 1)*(                   3443.0_rp*fm(s,1)- 5966.0_rp*fm(s, 0)+1602.0_rp*fm(s,-1)) + &
                    fm(s, 0)*(                                      2843.0_rp*fm(s, 0)-1642.0_rp*fm(s,-1)) + &
                    fm(s,-1)*(                                                          267.0_rp*fm(s,-1))
           ISp(3) = fp(s, 0)*(2107.0_rp*fp(s, 0)- 9402.0_rp*fp(s, 1)+ 7042.0_rp*fp(s, 2)-1854.0_rp*fp(s,3)) + &
                    fp(s, 1)*(                   11003.0_rp*fp(s, 1)-17246.0_rp*fp(s, 2)+4642.0_rp*fp(s,3)) + &
                    fp(s, 2)*(                                        7043.0_rp*fp(s, 2)-3882.0_rp*fp(s,3)) + &
                    fp(s, 3)*(                                                            547.0_rp*fp(s,3))
           ISm(3) = fm(s, 1)*(2107.0_rp*fm(s,1)- 9402.0_rp*fm(s,0)+ 7042.0_rp*fm(s,-1)-1854.0_rp*fm(s,-2)) + &
                    fm(s, 0)*(                  11003.0_rp*fm(s,0)-17246.0_rp*fm(s,-1)+4642.0_rp*fm(s,-2)) + &
                    fm(s,-1)*(                                      7043.0_rp*fm(s,-1)-3882.0_rp*fm(s,-2)) + &
                    fm(s,-2)*(                                                          547.0_rp*fm(s,-2))

           tp(0) = abs(ISp(0)             - ISp(3)) ; tm(0) = abs(ISm(0)             - ISm(3))
           tp(1) = abs(ISp(0)-ISp(1)-ISp(2)+ISp(3)) ; tm(1) = abs(ISm(0)-ISm(1)-ISm(2)+ISm(3))
           do r = 0, 3
              ap(r) = tp(mod(r,2))/(ISp(r) + varepsilon) ; ap(r) = cweno7(r)*(1.0_rp + ap(r)*ap(r))
              am(r) = tm(mod(r,2))/(ISm(r) + varepsilon) ; am(r) = cweno7(r)*(1.0_rp + am(r)*am(r))
           enddo

           sp = 1.0_rp/(ap(0) + ap(1) + ap(2) + ap(3)) ; sm = 1.0_rp/(am(0) + am(1) + am(2) + am(3))
           fhatp(s) = (ap(0)*sp)*qp(0) + (ap(1)*sp)*qp(1) + (ap(2)*sp)*qp(2) + (ap(3)*sp)*qp(3)
           fhatm(s) = (am(0)*sm)*qm(0) + (am(1)*sm)*qm(1) + (am(2)*sm)*qm(2) + (am(3)*sm)*qm(3)
        enddo

        endif

        elseif(shock_recon == 2) then

        if(weno_num == 3) then
        do s = 1, 5
           qp(0) = aweno5(0,0)*fp(s,-2) + aweno5(0,1)*fp(s,-1) + aweno5(0,2)*fp(s,0)
           qm(0) = aweno5(0,0)*fm(s, 3) + aweno5(0,1)*fm(s, 2) + aweno5(0,2)*fm(s,1)
           qp(1) = aweno5(1,0)*fp(s,-1) + aweno5(1,1)*fp(s, 0) + aweno5(1,2)*fp(s,1)
           qm(1) = aweno5(1,0)*fm(s, 2) + aweno5(1,1)*fm(s, 1) + aweno5(1,2)*fm(s,0)
           qp(2) = aweno5(2,0)*fp(s, 0) + aweno5(2,1)*fp(s, 1) + aweno5(2,2)*fp(s,2)
           qm(2) = aweno5(2,0)*fm(s, 1) + aweno5(2,1)*fm(s, 0) + aweno5(2,2)*fm(s,-1)

           ISp(0) = is1par*(fp(s,-2)-2.0_rp*fp(s,-1)+fp(s,0))**2 + is2par*(    fp(s,-2)-4.0_rp*fp(s,-1)+3.0_rp*fp(s,0))**2
           ISm(0) = is1par*(fm(s, 3)-2.0_rp*fm(s, 2)+fm(s,1))**2 + is2par*(    fm(s, 3)-4.0_rp*fm(s, 2)+3.0_rp*fm(s,1))**2
           ISp(1) = is1par*(fp(s,-1)-2.0_rp*fp(s, 0)+fp(s,1))**2 + is2par*(    fp(s,-1)                       -fp(s,1))**2
           ISm(1) = is1par*(fm(s, 2)-2.0_rp*fm(s, 1)+fm(s,0))**2 + is2par*(    fm(s, 2)                       -fm(s,0))**2
           ISp(2) = is1par*(fp(s, 0)-2.0_rp*fp(s, 1)+fp(s,2))**2 + is2par*(3.0_rp*fp(s, 0)-4.0_rp*fp(s, 1)    +fp(s,2))**2
           ISm(2) = is1par*(fm(s, 1)-2.0_rp*fm(s, 0)+fm(s,-1))**2 + is2par*(3.0_rp*fm(s, 1)-4.0_rp*fm(s, 0)    +fm(s,-1))**2

           tp(0) = abs(ISp(0)-ISp(2)) ; tm(0) = abs(ISm(0)-ISm(2))
           do r = 0, 2
              ap(r) = 1.0_rp + (tp(0)/(ISp(r) + varepsilon))**6
              am(r) = 1.0_rp + (tm(0)/(ISm(r) + varepsilon))**6
           enddo
           sp = 1.0_rp/(ap(0) + ap(1) + ap(2)) ; sm = 1.0_rp/(am(0) + am(1) + am(2))

           do r = 0, 2
              dp(r) = 0 ; if(ap(r)*sp >= ct_teno) dp(r) = 1
              dm(r) = 0 ; if(am(r)*sm >= ct_teno) dm(r) = 1
           enddo
           dsp = 1.0_rp/(dp(0)*cweno5(0) + dp(1)*cweno5(1) + dp(2)*cweno5(2))
           dsm = 1.0_rp/(dm(0)*cweno5(0) + dm(1)*cweno5(1) + dm(2)*cweno5(2))
           fhatp(s) = (dp(0)*cweno5(0)*dsp)*qp(0) + (dp(1)*cweno5(1)*dsp)*qp(1) + (dp(2)*cweno5(2)*dsp)*qp(2)
           fhatm(s) = (dm(0)*cweno5(0)*dsm)*qm(0) + (dm(1)*cweno5(1)*dsm)*qm(1) + (dm(2)*cweno5(2)*dsm)*qm(2)
        enddo

        elseif(weno_num == 4) then
        do s = 1, 5
           qp(0) = aweno7(0,0)*fp(s,-3) + aweno7(0,1)*fp(s,-2) + aweno7(0,2)*fp(s,-1) + aweno7(0,3)*fp(s,0)
           qm(0) = aweno7(0,0)*fm(s, 4) + aweno7(0,1)*fm(s, 3) + aweno7(0,2)*fm(s, 2) + aweno7(0,3)*fm(s, 1)
           qp(1) = aweno7(1,0)*fp(s,-2) + aweno7(1,1)*fp(s,-1) + aweno7(1,2)*fp(s, 0) + aweno7(1,3)*fp(s,1)
           qm(1) = aweno7(1,0)*fm(s, 3) + aweno7(1,1)*fm(s, 2) + aweno7(1,2)*fm(s, 1) + aweno7(1,3)*fm(s, 0)
           qp(2) = aweno7(2,0)*fp(s,-1) + aweno7(2,1)*fp(s, 0) + aweno7(2,2)*fp(s, 1) + aweno7(2,3)*fp(s,2)
           qm(2) = aweno7(2,0)*fm(s, 2) + aweno7(2,1)*fm(s, 1) + aweno7(2,2)*fm(s, 0) + aweno7(2,3)*fm(s,-1)
           qp(3) = aweno7(3,0)*fp(s, 0) + aweno7(3,1)*fp(s, 1) + aweno7(3,2)*fp(s, 2) + aweno7(3,3)*fp(s,3)
           qm(3) = aweno7(3,0)*fm(s, 1) + aweno7(3,1)*fm(s, 0) + aweno7(3,2)*fm(s,-1) + aweno7(3,3)*fm(s,-2)

           ISp(0) = fp(s,-3)*( 547.0_rp*fp(s,-3)- 3882.0_rp*fp(s,-2)+ 4642.0_rp*fp(s,-1)-1854.0_rp*fp(s,0)) + &
                    fp(s,-2)*(                    7043.0_rp*fp(s,-2)-17246.0_rp*fp(s,-1)+7042.0_rp*fp(s,0)) + &
                    fp(s,-1)*(                                       11003.0_rp*fp(s,-1)-9402.0_rp*fp(s,0)) + &
                    fp(s, 0)*(                                                           2107.0_rp*fp(s,0))
           ISm(0) = fm(s, 4)*( 547.0_rp*fm(s,4)- 3882.0_rp*fm(s,3)+ 4642.0_rp*fm(s, 2)-1854.0_rp*fm(s, 1)) + &
                    fm(s, 3)*(                   7043.0_rp*fm(s,3)-17246.0_rp*fm(s, 2)+7042.0_rp*fm(s, 1)) + &
                    fm(s, 2)*(                                     11003.0_rp*fm(s, 2)-9402.0_rp*fm(s, 1)) + &
                    fm(s, 1)*(                                                         2107.0_rp*fm(s, 1))
           ISp(1) = fp(s,-2)*( 267.0_rp*fp(s,-2)- 1642.0_rp*fp(s,-1)+ 1602.0_rp*fp(s, 0)- 494.0_rp*fp(s,1)) + &
                    fp(s,-1)*(                    2843.0_rp*fp(s,-1)- 5966.0_rp*fp(s, 0)+1922.0_rp*fp(s,1)) + &
                    fp(s, 0)*(                                        3443.0_rp*fp(s, 0)-2522.0_rp*fp(s,1)) + &
                    fp(s, 1)*(                                                            547.0_rp*fp(s,1))
           ISm(1) = fm(s, 3)*( 267.0_rp*fm(s,3)- 1642.0_rp*fm(s,2)+ 1602.0_rp*fm(s, 1)- 494.0_rp*fm(s, 0)) + &
                    fm(s, 2)*(                   2843.0_rp*fm(s,2)- 5966.0_rp*fm(s, 1)+1922.0_rp*fm(s, 0)) + &
                    fm(s, 1)*(                                      3443.0_rp*fm(s, 1)-2522.0_rp*fm(s, 0)) + &
                    fm(s, 0)*(                                                          547.0_rp*fm(s, 0))
           ISp(2) = fp(s,-1)*( 547.0_rp*fp(s,-1)- 2522.0_rp*fp(s, 0)+ 1922.0_rp*fp(s, 1)- 494.0_rp*fp(s,2)) + &
                    fp(s, 0)*(                    3443.0_rp*fp(s, 0)- 5966.0_rp*fp(s, 1)+1602.0_rp*fp(s,2)) + &
                    fp(s, 1)*(                                        2843.0_rp*fp(s, 1)-1642.0_rp*fp(s,2)) + &
                    fp(s, 2)*(                                                            267.0_rp*fp(s,2))
           ISm(2) = fm(s, 2)*( 547.0_rp*fm(s,2)- 2522.0_rp*fm(s,1)+ 1922.0_rp*fm(s, 0)- 494.0_rp*fm(s,-1)) + &
                    fm(s, 1)*(                   3443.0_rp*fm(s,1)- 5966.0_rp*fm(s, 0)+1602.0_rp*fm(s,-1)) + &
                    fm(s, 0)*(                                      2843.0_rp*fm(s, 0)-1642.0_rp*fm(s,-1)) + &
                    fm(s,-1)*(                                                          267.0_rp*fm(s,-1))
           ISp(3) = fp(s, 0)*(2107.0_rp*fp(s, 0)- 9402.0_rp*fp(s, 1)+ 7042.0_rp*fp(s, 2)-1854.0_rp*fp(s,3)) + &
                    fp(s, 1)*(                   11003.0_rp*fp(s, 1)-17246.0_rp*fp(s, 2)+4642.0_rp*fp(s,3)) + &
                    fp(s, 2)*(                                        7043.0_rp*fp(s, 2)-3882.0_rp*fp(s,3)) + &
                    fp(s, 3)*(                                                            547.0_rp*fp(s,3))
           ISm(3) = fm(s, 1)*(2107.0_rp*fm(s,1)- 9402.0_rp*fm(s,0)+ 7042.0_rp*fm(s,-1)-1854.0_rp*fm(s,-2)) + &
                    fm(s, 0)*(                  11003.0_rp*fm(s,0)-17246.0_rp*fm(s,-1)+4642.0_rp*fm(s,-2)) + &
                    fm(s,-1)*(                                      7043.0_rp*fm(s,-1)-3882.0_rp*fm(s,-2)) + &
                    fm(s,-2)*(                                                          547.0_rp*fm(s,-2))

           tp(0) = abs(ISp(0)             - ISp(3)) ; tm(0) = abs(ISm(0)             - ISm(3))
           tp(1) = abs(ISp(0)-ISp(1)-ISp(2)+ISp(3)) ; tm(1) = abs(ISm(0)-ISm(1)-ISm(2)+ISm(3))
           do r = 0, 3
              ap(r) = 1.0_rp + (tp(mod(r,2))/(ISp(r) + varepsilon))**6
              am(r) = 1.0_rp + (tm(mod(r,2))/(ISm(r) + varepsilon))**6
           enddo
           sp = 1.0_rp/(ap(0) + ap(1) + ap(2) + ap(3)) ; sm = 1.0_rp/(am(0) + am(1) + am(2) + am(3))

           do r = 0, 3
              dp(r) = 0 ; if(ap(r)*sp >= ct_teno) dp(r) = 1
              dm(r) = 0 ; if(am(r)*sm >= ct_teno) dm(r) = 1
           enddo
           dsp = 1.0_rp/(dp(0)*cweno7(0) + dp(1)*cweno7(1) + dp(2)*cweno7(2) + dp(3)*cweno7(3))
           dsm = 1.0_rp/(dm(0)*cweno7(0) + dm(1)*cweno7(1) + dm(2)*cweno7(2) + dm(3)*cweno7(3))
           fhatp(s) = (dp(0)*cweno7(0)*dsp)*qp(0) + (dp(1)*cweno7(1)*dsp)*qp(1) + (dp(2)*cweno7(2)*dsp)*qp(2) + (dp(3)*cweno7(3)*dsp)*qp(3)
           fhatm(s) = (dm(0)*cweno7(0)*dsm)*qm(0) + (dm(1)*cweno7(1)*dsm)*qm(1) + (dm(2)*cweno7(2)*dsm)*qm(2) + (dm(3)*cweno7(3)*dsm)*qm(3)
        enddo

        endif

        elseif(shock_recon == 3) then

        do s = 1, 5
           qp(0) = aweno5(0,0)*fp(s,-2) + aweno5(0,1)*fp(s,-1) + aweno5(0,2)*fp(s,0)
           qm(0) = aweno5(0,0)*fm(s, 3) + aweno5(0,1)*fm(s, 2) + aweno5(0,2)*fm(s,1)
           qp(1) = aweno5(1,0)*fp(s,-1) + aweno5(1,1)*fp(s, 0) + aweno5(1,2)*fp(s,1)
           qm(1) = aweno5(1,0)*fm(s, 2) + aweno5(1,1)*fm(s, 1) + aweno5(1,2)*fm(s,0)
           qp(2) = aweno5(2,0)*fp(s, 0) + aweno5(2,1)*fp(s, 1) + aweno5(2,2)*fp(s,2)
           qm(2) = aweno5(2,0)*fm(s, 1) + aweno5(2,1)*fm(s, 0) + aweno5(2,2)*fm(s,-1)

           ISp(0) = is1par*(fp(s,-2)-2.0_rp*fp(s,-1)+fp(s,0))**2 + is2par*(    fp(s,-2)-4.0_rp*fp(s,-1)+3.0_rp*fp(s,0))**2
           ISm(0) = is1par*(fm(s, 3)-2.0_rp*fm(s, 2)+fm(s,1))**2 + is2par*(    fm(s, 3)-4.0_rp*fm(s, 2)+3.0_rp*fm(s,1))**2
           ISp(1) = is1par*(fp(s,-1)-2.0_rp*fp(s, 0)+fp(s,1))**2 + is2par*(    fp(s,-1)                       -fp(s,1))**2
           ISm(1) = is1par*(fm(s, 2)-2.0_rp*fm(s, 1)+fm(s,0))**2 + is2par*(    fm(s, 2)                       -fm(s,0))**2
           ISp(2) = is1par*(fp(s, 0)-2.0_rp*fp(s, 1)+fp(s,2))**2 + is2par*(3.0_rp*fp(s, 0)-4.0_rp*fp(s, 1)    +fp(s,2))**2
           ISm(2) = is1par*(fm(s, 1)-2.0_rp*fm(s, 0)+fm(s,-1))**2 + is2par*(3.0_rp*fm(s, 1)-4.0_rp*fm(s, 0)    +fm(s,-1))**2

           tp(0) = abs(ISp(0)-ISp(2)) ; tm(0) = abs(ISm(0)-ISm(2))
           do r = 0, 2
              ap(r) = 1.0_rp + (tp(0)/(ISp(r) + varepsilon))**6
              am(r) = 1.0_rp + (tm(0)/(ISm(r) + varepsilon))**6
           enddo
           sp = 1.0_rp/(ap(0) + ap(1) + ap(2)) ; sm = 1.0_rp/(am(0) + am(1) + am(2))

           e0p = (2.0_rp*abs((fp(s, 0)-fp(s, -1))*(fp(s, -1)-fp(s, -2))) + eps)/ &
                        ((fp(s, 0)-fp(s, -1))**2+(fp(s, -1)-fp(s, -2))**2+ eps)
           e1p = (2.0_rp*abs((fp(s, 1)-fp(s,  0))*(fp(s,  0)-fp(s, -1))) + eps)/ &
                        ((fp(s, 1)-fp(s,  0))**2+(fp(s,  0)-fp(s, -1))**2+ eps)
           e2p = (2.0_rp*abs((fp(s, 2)-fp(s,  1))*(fp(s,  1)-fp(s,  0))) + eps)/ &
                        ((fp(s, 2)-fp(s,  1))**2+(fp(s,  1)-fp(s,  0))**2+ eps)
           rmp = 1.0_rp-min(1.0_rp,min(e0p,e1p,e2p)/cr)
           ctp = 1.0_rp*10.0_rp**(-floor(a1-a2*(1.0_rp-((1.0_rp-rmp)**4*(1.0_rp+4.0_rp*rmp)))))

           e0m = (2.0_rp*abs((fm(s,  1)-fm(s, 2))*(fm(s, 2)-fm(s, 3))) + eps)/ &
                        ((fm(s,  1)-fm(s, 2))**2+(fm(s, 2)-fm(s, 3))**2+ eps)
           e1m = (2.0_rp*abs((fm(s,  0)-fm(s, 1))*(fm(s, 1)-fm(s, 2))) + eps)/ &
                        ((fm(s,  0)-fm(s, 1))**2+(fm(s, 1)-fm(s, 2))**2+ eps)
           e2m = (2.0_rp*abs((fm(s, -1)-fm(s, 0))*(fm(s, 0)-fm(s, 1))) + eps)/ &
                        ((fm(s, -1)-fm(s, 0))**2+(fm(s, 0)-fm(s, 1))**2+ eps)
           rmm = 1.0_rp-min(1.0_rp,min(e0m,e1m,e2m)/cr)
           ctm = 1.0_rp*10.0_rp**(-floor(a1-a2*(1.0_rp-((1.0_rp-rmm)**4*(1.0_rp+4.0_rp*rmm)))))

           do r = 0, 2
              dp(r) = 0 ; if(ap(r)*sp >= ctp) dp(r) = 1
              dm(r) = 0 ; if(am(r)*sm >= ctm) dm(r) = 1
           enddo
           dsp = 1.0_rp/(dp(0)*cweno5(0) + dp(1)*cweno5(1) + dp(2)*cweno5(2))
           dsm = 1.0_rp/(dm(0)*cweno5(0) + dm(1)*cweno5(1) + dm(2)*cweno5(2))
           fhatp(s) = (dp(0)*cweno5(0)*dsp)*qp(0) + (dp(1)*cweno5(1)*dsp)*qp(1) + (dp(2)*cweno5(2)*dsp)*qp(2)
           fhatm(s) = (dm(0)*cweno5(0)*dsm)*qm(0) + (dm(1)*cweno5(1)*dsm)*qm(1) + (dm(2)*cweno5(2)*dsm)*qm(2)
        enddo

        endif

return
end subroutine weno_reconstruction

end module weno_reconstruction_module
