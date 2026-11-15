module thermodynamics_module
use parameters_module
implicit none

  real(rp), parameter :: R_0 = 1.0_rp
  real(rp), parameter :: T_0 = 1.0_rp
  real(rp), parameter :: rho_0 = 1.0_rp
  real(rp), parameter :: cp0 = 3.5_rp
  real(rp), parameter :: cv0 = 2.5_rp
  real(rp), parameter :: P_0 = rho_0 * R_0 * T_0
  real(rp)            :: cp_Tref, hL_Tref, hH_Tref, cv_Tref, gamma_Tref

  real(rp), parameter :: aL_air_m2 =  1.009950160E+04_rp
  real(rp), parameter :: aL_air_m1 = -1.968275610E+02_rp
  real(rp), parameter :: aL_air__0 =  5.009155110E+00_rp
  real(rp), parameter :: aL_air_p1 = -5.761013730E-03_rp
  real(rp), parameter :: aL_air_p2 =  1.066859930E-05_rp
  real(rp), parameter :: aL_air_p3 = -7.940297970E-09_rp
  real(rp), parameter :: aL_air_p4 =  2.185231910E-12_rp
  real(rp), parameter :: aH_air_m2 =  2.415214430E+05_rp
  real(rp), parameter :: aH_air_m1 = -1.257874600E+03_rp
  real(rp), parameter :: aH_air__0 =  5.144558670E+00_rp
  real(rp), parameter :: aH_air_p1 = -2.138541790E-04_rp
  real(rp), parameter :: aH_air_p2 =  7.065227840E-08_rp
  real(rp), parameter :: aH_air_p3 = -1.071483490E-11_rp
  real(rp), parameter :: aH_air_p4 =  6.577800150E-16_rp
  real(rp), parameter :: aL_N2_m2  =  2.210371497E+04_rp
  real(rp), parameter :: aL_N2_m1  = -3.818461820E+02_rp
  real(rp), parameter :: aL_N2__0  =  6.082738360E+00_rp
  real(rp), parameter :: aL_N2_p1  = -8.530914410E-03_rp
  real(rp), parameter :: aL_N2_p2  =  1.384646189E-05_rp
  real(rp), parameter :: aL_N2_p3  = -9.625793620E-09_rp
  real(rp), parameter :: aL_N2_p4  =  2.519705809E-12_rp
  real(rp), parameter :: aH_N2_m2  =  5.877124060E+05_rp
  real(rp), parameter :: aH_N2_m1  = -2.239249073E+03_rp
  real(rp), parameter :: aH_N2__0  =  6.066949220E+00_rp
  real(rp), parameter :: aH_N2_p1  = -6.139685500E-04_rp
  real(rp), parameter :: aH_N2_p2  =  1.491806679E-07_rp
  real(rp), parameter :: aH_N2_p3  = -1.923105485E-11_rp
  real(rp), parameter :: aH_N2_p4  =  1.061954386E-15_rp
  real(rp), parameter :: aL_O2_m2  = -3.425563420E+04_rp
  real(rp), parameter :: aL_O2_m1  =  4.847000970E+02_rp
  real(rp), parameter :: aL_O2__0  =  1.119010961E+00_rp
  real(rp), parameter :: aL_O2_p1  =  4.293889240E-03_rp
  real(rp), parameter :: aL_O2_p2  = -6.836300520E-07_rp
  real(rp), parameter :: aL_O2_p3  = -2.023372700E-09_rp
  real(rp), parameter :: aL_O2_p4  =  1.039040018E-12_rp
  real(rp), parameter :: aH_O2_m2  = -1.037939022E+06_rp
  real(rp), parameter :: aH_O2_m1  =  2.344830282E+03_rp
  real(rp), parameter :: aH_O2__0  =  1.819732036E+00_rp
  real(rp), parameter :: aH_O2_p1  =  1.267847582E-03_rp
  real(rp), parameter :: aH_O2_p2  = -2.188067988E-07_rp
  real(rp), parameter :: aH_O2_p3  =  2.053719572E-11_rp
  real(rp), parameter :: aH_O2_p4  = -8.193467050E-16_rp

  real(rp), parameter :: suth_air = 110.4_rp
  real(rp), parameter :: suth_N2  = 107.0_rp
  real(rp), parameter :: suth_O2  = 139.0_rp
  real(rp) :: suth_S = suth_air
  !$acc declare create(suth_S)

  real(rp) :: aL_adim_m2, aL_adim_m1, aL_adim__0, aL_adim_p1, aL_adim_p2, aL_adim_p3, aL_adim_p4
  real(rp) :: aH_adim_m2, aH_adim_m1, aH_adim__0, aH_adim_p1, aH_adim_p2, aH_adim_p3, aH_adim_p4
  real(rp) :: hL_adim_m2, hL_adim_p1, hL_adim_p2, hL_adim_p3, hL_adim_p4
  real(rp) :: hH_adim_m2, hH_adim_p1, hH_adim_p2, hH_adim_p3, hH_adim_p4
  real(rp) :: sL_adim_m2, sL_adim_m1, sL_adim_p1, sL_adim_p2, sL_adim_p3, sL_adim_p4
  real(rp) :: sH_adim_m2, sH_adim_m1, sH_adim_p1, sH_adim_p2, sH_adim_p3, sH_adim_p4
  real(rp) :: sL_Tref, sH_Tref
  !$acc declare create(aL_adim_m2, aL_adim_m1, aL_adim__0, aL_adim_p1, aL_adim_p2, aL_adim_p3, aL_adim_p4)
  !$acc declare create(aH_adim_m2, aH_adim_m1, aH_adim__0, aH_adim_p1, aH_adim_p2, aH_adim_p3, aH_adim_p4)
  !$acc declare create(hL_adim_m2, hL_adim_p1, hL_adim_p2, hL_adim_p3, hL_adim_p4)
  !$acc declare create(hH_adim_m2, hH_adim_p1, hH_adim_p2, hH_adim_p3, hH_adim_p4)
  !$acc declare create(sL_adim_m2, sL_adim_m1, sL_adim_p1, sL_adim_p2, sL_adim_p3, sL_adim_p4)
  !$acc declare create(sH_adim_m2, sH_adim_m1, sH_adim_p1, sH_adim_p2, sH_adim_p3, sH_adim_p4)
  !$acc declare create(hL_Tref, hH_Tref, sL_Tref, sH_Tref)
  !$acc declare create(cv_Tref,gamma_Tref,cp_Tref)
  real(rp) :: T_switch
  !$acc declare create(T_switch)
  real(rp) :: T_cold, aC_adim__0, aC_adim_p4, hC_Tref, sC_Tref
  !$acc declare create(T_cold, aC_adim__0, aC_adim_p4, hC_Tref, sC_Tref)
  real(rp) :: cp_coef(-2:4,0:2), h_coef(0:2)
  !$acc declare create(cp_coef, h_coef)

contains

subroutine init_non_calorically_coefficients
   use mpi_module, only: rank, root
   implicit none
   real(rp) :: invTref, Tref2, Tref3, Tref4, invTref2
   real(rp) :: h_at_1
   real(rp) :: aL(-2:4), aH(-2:4)

   select case (species_name)
      case ("air")
         suth_S = suth_air
         aL = [aL_air_m2, aL_air_m1, aL_air__0, aL_air_p1, aL_air_p2, aL_air_p3, aL_air_p4]
         aH = [aH_air_m2, aH_air_m1, aH_air__0, aH_air_p1, aH_air_p2, aH_air_p3, aH_air_p4]
      case ("N2")
         suth_S = suth_N2
         aL = [aL_N2_m2, aL_N2_m1, aL_N2__0, aL_N2_p1, aL_N2_p2, aL_N2_p3, aL_N2_p4]
         aH = [aH_N2_m2, aH_N2_m1, aH_N2__0, aH_N2_p1, aH_N2_p2, aH_N2_p3, aH_N2_p4]
      case ("O2")
         suth_S = suth_O2
         aL = [aL_O2_m2, aL_O2_m1, aL_O2__0, aL_O2_p1, aL_O2_p2, aL_O2_p3, aL_O2_p4]
         aH = [aH_O2_m2, aH_O2_m1, aH_O2__0, aH_O2_p1, aH_O2_p2, aH_O2_p3, aH_O2_p4]
      case default
         stop "Unknown species_name in init_non_calorically_coefficients"
   end select
   !$acc update device(suth_S)

   if(thermally_perfect_gas > 0) then
     invTref   = 1.0_rp  / Tref
     Tref2     = Tref    * Tref
     Tref3     = Tref2   * Tref
     Tref4     = Tref3   * Tref
     invTref2  = invTref * invTref

     aL_adim_m2 = aL(-2) * invTref2
     aL_adim_m1 = aL(-1) * invTref
     aL_adim__0 = aL(0)
     aL_adim_p1 = aL(1) * Tref
     aL_adim_p2 = aL(2) * Tref2
     aL_adim_p3 = aL(3) * Tref3
     aL_adim_p4 = aL(4) * Tref4

     aH_adim_m2 = aH(-2) * invTref2
     aH_adim_m1 = aH(-1) * invTref
     aH_adim__0 = aH(0)
     aH_adim_p1 = aH(1) * Tref
     aH_adim_p2 = aH(2) * Tref2
     aH_adim_p3 = aH(3) * Tref3
     aH_adim_p4 = aH(4) * Tref4
     hL_adim_m2 =             - aL(-2) * invTref2
     hL_adim_p1 =   0.5_rp *  ( aL(1) * Tref )
     hL_adim_p2 =             ( aL(2) * Tref2) / 3.0_rp
     hL_adim_p3 =   0.25_rp * ( aL(3) * Tref3)
     hL_adim_p4 =   0.2_rp  * ( aL(4) * Tref4)

     hH_adim_m2 =             - aH(-2) * invTref2
     hH_adim_p1 =   0.5_rp *  ( aH(1) * Tref )
     hH_adim_p2 =             ( aH(2) * Tref2) / 3.0_rp
     hH_adim_p3 =   0.25_rp * ( aH(3) * Tref3)
     hH_adim_p4 =   0.2_rp *  ( aH(4) * Tref4)
     sL_adim_m2 = -0.5_rp  *  aL_adim_m2
     sL_adim_m1 =          -  aL_adim_m1
     sL_adim_p1 =             aL_adim_p1
     sL_adim_p2 =  0.5_rp  *  aL_adim_p2
     sL_adim_p3 =             aL_adim_p3 / 3.0_rp
     sL_adim_p4 =  0.25_rp *  aL_adim_p4

     sH_adim_m2 = -0.5_rp  *  aH_adim_m2
     sH_adim_m1 =          -  aH_adim_m1
     sH_adim_p1 =             aH_adim_p1
     sH_adim_p2 =  0.5_rp  *  aH_adim_p2
     sH_adim_p3 =             aH_adim_p3 / 3.0_rp
     sH_adim_p4 =  0.25_rp *  aH_adim_p4

     T_switch   = 1000.0_rp / Tref

     T_cold = tail_anchor_temperature()
     if (T_cold > 0.0_rp) then
        aC_adim_p4 = 0.25_rp*dcp_low(T_cold)/T_cold**3
        aC_adim__0 = cp_low(T_cold) - aC_adim_p4*T_cold**4
        if (rank == root) write(*,'(A,F8.2,A,ES10.3)') ' init_non_calorically_coefficients: '// &
             'cold tail below ', T_cold*Tref, ' K, cp drop to 0 K ', aC_adim_p4*T_cold**4
     else
        aC_adim__0 = 0.0_rp
        aC_adim_p4 = 0.0_rp
        if (rank == root) write(*,'(A)') ' init_non_calorically_coefficients: no anchor for '// &
             'the cold tail in 200-1000 K, the LOW cp fit is extrapolated below 200 K'
     endif

     if (T_cold >= 1.0_rp) then
        cp_Tref = aC_adim__0 + aC_adim_p4
        hC_Tref = cp_Tref - h_cold_kernel(1.0_rp)
        sC_Tref = -s_cold_kernel(1.0_rp)
        hL_Tref = hC_Tref + h_cold_kernel(T_cold) - h_low_kernel(T_cold)
        sL_Tref = sC_Tref + s_cold_kernel(T_cold) - s_low_kernel(T_cold)
        call enforce_branch_continuity(T_switch,                   &
             hL_Tref, sL_Tref,                                     &
             aL_adim_m1, aL_adim__0,                               &
             hL_adim_m2, hL_adim_p1, hL_adim_p2,                   &
             hL_adim_p3, hL_adim_p4,                               &
             sL_adim_m2, sL_adim_m1,                               &
             sL_adim_p1, sL_adim_p2, sL_adim_p3, sL_adim_p4,       &
             aH_adim_m1, aH_adim__0,                               &
             hH_adim_m2, hH_adim_p1, hH_adim_p2,                   &
             hH_adim_p3, hH_adim_p4,                               &
             sH_adim_m2, sH_adim_m1,                               &
             sH_adim_p1, sH_adim_p2, sH_adim_p3, sH_adim_p4,       &
             hH_Tref, sH_Tref)
     else if (Tref <= 1000.0_rp) then
        cp_Tref = aL_adim_m2 + aL_adim_m1 + aL_adim__0 + &
                 aL_adim_p1 + aL_adim_p2 + aL_adim_p3 + aL_adim_p4
        h_at_1 = hL_adim_m2 + aL_adim__0 + hL_adim_p1 + &
                 hL_adim_p2 + hL_adim_p3 + hL_adim_p4
        hL_Tref = cp_Tref - h_at_1
        sL_Tref = -(sL_adim_m2 + sL_adim_m1 + &
                 sL_adim_p1 + sL_adim_p2 + sL_adim_p3 + sL_adim_p4)
        call enforce_branch_continuity(T_switch,                   &
             hL_Tref, sL_Tref,                                     &
             aL_adim_m1, aL_adim__0,                               &
             hL_adim_m2, hL_adim_p1, hL_adim_p2,                   &
             hL_adim_p3, hL_adim_p4,                               &
             sL_adim_m2, sL_adim_m1,                               &
             sL_adim_p1, sL_adim_p2, sL_adim_p3, sL_adim_p4,       &
             aH_adim_m1, aH_adim__0,                               &
             hH_adim_m2, hH_adim_p1, hH_adim_p2,                   &
             hH_adim_p3, hH_adim_p4,                               &
             sH_adim_m2, sH_adim_m1,                               &
             sH_adim_p1, sH_adim_p2, sH_adim_p3, sH_adim_p4,       &
             hH_Tref, sH_Tref)
     else
        cp_Tref = aH_adim_m2 + aH_adim_m1 + aH_adim__0 + &
                 aH_adim_p1 + aH_adim_p2 + aH_adim_p3 + aH_adim_p4
        h_at_1 = hH_adim_m2 + aH_adim__0 + hH_adim_p1 + &
                 hH_adim_p2 + hH_adim_p3 + hH_adim_p4
        hH_Tref = cp_Tref - h_at_1
        sH_Tref = -(sH_adim_m2 + sH_adim_m1 + &
                 sH_adim_p1 + sH_adim_p2 + sH_adim_p3 + sH_adim_p4)
        call enforce_branch_continuity(T_switch,                   &
             hH_Tref, sH_Tref,                                     &
             aH_adim_m1, aH_adim__0,                               &
             hH_adim_m2, hH_adim_p1, hH_adim_p2,                   &
             hH_adim_p3, hH_adim_p4,                               &
             sH_adim_m2, sH_adim_m1,                               &
             sH_adim_p1, sH_adim_p2, sH_adim_p3, sH_adim_p4,       &
             aL_adim_m1, aL_adim__0,                               &
             hL_adim_m2, hL_adim_p1, hL_adim_p2,                   &
             hL_adim_p3, hL_adim_p4,                               &
             sL_adim_m2, sL_adim_m1,                               &
             sL_adim_p1, sL_adim_p2, sL_adim_p3, sL_adim_p4,       &
             hL_Tref, sL_Tref)
     endif
     if (T_cold > 0.0_rp .and. T_cold < 1.0_rp) then
        hC_Tref = hL_Tref + h_low_kernel(T_cold) - h_cold_kernel(T_cold)
        sC_Tref = sL_Tref + s_low_kernel(T_cold) - s_cold_kernel(T_cold)
     else if (T_cold <= 0.0_rp) then
        hC_Tref = 0.0_rp
        sC_Tref = 0.0_rp
     endif
     cv_Tref   = cp_Tref - 1.0_rp
     gamma_Tref = cp_Tref / cv_Tref

   else
     gamma_Tref = gamma0
     cv_Tref    = cv0
     cp_Tref    = cp0
     hL_Tref    = cp0
     hH_Tref    = cp0
     aL_adim_m2 = 0.0_rp; aL_adim_m1 = 0.0_rp; aL_adim__0 = 0.0_rp
     aL_adim_p1 = 0.0_rp; aL_adim_p2 = 0.0_rp; aL_adim_p3 = 0.0_rp; aL_adim_p4 = 0.0_rp
     aH_adim_m2 = 0.0_rp; aH_adim_m1 = 0.0_rp; aH_adim__0 = 0.0_rp
     aH_adim_p1 = 0.0_rp; aH_adim_p2 = 0.0_rp; aH_adim_p3 = 0.0_rp; aH_adim_p4 = 0.0_rp
     hL_adim_m2 = 0.0_rp; hL_adim_p1 = 0.0_rp; hL_adim_p2 = 0.0_rp
     hL_adim_p3 = 0.0_rp; hL_adim_p4 = 0.0_rp
     hH_adim_m2 = 0.0_rp; hH_adim_p1 = 0.0_rp; hH_adim_p2 = 0.0_rp
     hH_adim_p3 = 0.0_rp; hH_adim_p4 = 0.0_rp
     sL_adim_m2 = 0.0_rp; sL_adim_m1 = 0.0_rp; sL_adim_p1 = 0.0_rp
     sL_adim_p2 = 0.0_rp; sL_adim_p3 = 0.0_rp; sL_adim_p4 = 0.0_rp
     sH_adim_m2 = 0.0_rp; sH_adim_m1 = 0.0_rp; sH_adim_p1 = 0.0_rp
     sH_adim_p2 = 0.0_rp; sH_adim_p3 = 0.0_rp; sH_adim_p4 = 0.0_rp
     sL_Tref    = 0.0_rp
     sH_Tref    = 0.0_rp
     T_switch   = huge(1.0_rp)
     T_cold     = 0.0_rp
     aC_adim_p4 = 0.0_rp
     aC_adim__0 = 0.0_rp
     hC_Tref    = cp0
     sC_Tref    = 0.0_rp
   endif

   cp_coef(:,0) = [0.0_rp, 0.0_rp, aC_adim__0, 0.0_rp, 0.0_rp, 0.0_rp, aC_adim_p4]
   cp_coef(:,1) = [aL_adim_m2, aL_adim_m1, aL_adim__0, aL_adim_p1, aL_adim_p2, aL_adim_p3, aL_adim_p4]
   cp_coef(:,2) = [aH_adim_m2, aH_adim_m1, aH_adim__0, aH_adim_p1, aH_adim_p2, aH_adim_p3, aH_adim_p4]
   h_coef       = [hC_Tref, hL_Tref, hH_Tref]

   !$acc update device(aL_adim_m2, aL_adim_m1, aL_adim__0, aL_adim_p1, aL_adim_p2, aL_adim_p3, aL_adim_p4)
   !$acc update device(aH_adim_m2, aH_adim_m1, aH_adim__0, aH_adim_p1, aH_adim_p2, aH_adim_p3, aH_adim_p4)
   !$acc update device(hL_adim_m2, hL_adim_p1, hL_adim_p2, hL_adim_p3, hL_adim_p4)
   !$acc update device(hH_adim_m2, hH_adim_p1, hH_adim_p2, hH_adim_p3, hH_adim_p4)
   !$acc update device(sL_adim_m2, sL_adim_m1, sL_adim_p1, sL_adim_p2, sL_adim_p3, sL_adim_p4)
   !$acc update device(sH_adim_m2, sH_adim_m1, sH_adim_p1, sH_adim_p2, sH_adim_p3, sH_adim_p4)
   !$acc update device(sL_Tref, sH_Tref)
   !$acc update device(hL_Tref, hH_Tref, cv_Tref, cp_Tref, gamma_Tref)
   !$acc update device(T_switch)
   !$acc update device(T_cold, aC_adim__0, aC_adim_p4, hC_Tref, sC_Tref)
   !$acc update device(cp_coef, h_coef)

   return
end subroutine init_non_calorically_coefficients

function tail_anchor_temperature() result(Ta)
   implicit none
   real(rp) :: Ta
   real(rp), parameter :: Tfloor_K = 200.0_rp, Thi_K = 1000.0_rp
   real(rp) :: Tf, Th, T1, Tz

   Ta = 0.0_rp
   Tf = Tfloor_K/Tref
   Th = Thi_K/Tref
   T1 = Tf
   if (dcp_low(Tf) <= 0.0_rp) then
      if (dcp_low(Th) <= 0.0_rp .or. cp_low(Th) <= cp_low(Tf)) return
      Tz = rising_crossing(1, 0.0_rp, Tf, Th)
      T1 = rising_crossing(2, cp_low(Tf), Tz, Th)
   endif
   Ta = T1

contains

   function rising_crossing(kind, target, a, b) result(Tc)
      integer , intent(in) :: kind
      real(rp), intent(in) :: target, a, b
      real(rp) :: Tc, lo, hi, f
      integer  :: it

      lo = a
      hi = b
      Tc = 0.5_rp*(lo + hi)
      do it = 1, 200
         Tc = 0.5_rp*(lo + hi)
         if (Tc <= lo .or. Tc >= hi) exit
         if (kind == 1) then
            f = dcp_low(Tc)
         else
            f = cp_low(Tc)
         endif
         if (f < target) then
            lo = Tc
         else
            hi = Tc
         endif
      end do
   end function rising_crossing

end function tail_anchor_temperature

pure function cp_low(T) result(c)
   real(rp), intent(in) :: T
   real(rp) :: c
   c = (aL_adim_m1 + aL_adim_m2/T)/T + aL_adim__0 + &
       T*(aL_adim_p1 + T*(aL_adim_p2 + T*(aL_adim_p3 + T*aL_adim_p4)))
end function cp_low

pure function dcp_low(T) result(d)
   real(rp), intent(in) :: T
   real(rp) :: d
   d = -(2.0_rp*aL_adim_m2/T + aL_adim_m1)/(T*T) + &
       aL_adim_p1 + T*(2.0_rp*aL_adim_p2 + T*(3.0_rp*aL_adim_p3 + T*4.0_rp*aL_adim_p4))
end function dcp_low

pure function h_low_kernel(T) result(h)
   real(rp), intent(in) :: T
   real(rp) :: h
   h = hL_adim_m2/T + aL_adim_m1*log(T) + &
       T*(aL_adim__0 + T*(hL_adim_p1 + T*(hL_adim_p2 + T*(hL_adim_p3 + T*hL_adim_p4))))
end function h_low_kernel

pure function s_low_kernel(T) result(sk)
   real(rp), intent(in) :: T
   real(rp) :: sk
   sk = sL_adim_m2/(T*T) + sL_adim_m1/T + aL_adim__0*log(T) + &
        T*(sL_adim_p1 + T*(sL_adim_p2 + T*(sL_adim_p3 + T*sL_adim_p4)))
end function s_low_kernel

pure function h_cold_kernel(T) result(h)
   real(rp), intent(in) :: T
   real(rp) :: h
   h = T*(aC_adim__0 + 0.2_rp*aC_adim_p4*T**4)
end function h_cold_kernel

pure function s_cold_kernel(T) result(sk)
   real(rp), intent(in) :: T
   real(rp) :: sk
   sk = aC_adim__0*log(T) + 0.25_rp*aC_adim_p4*T**4
end function s_cold_kernel

pure elemental function heat_capacity_t(T) result(cp_t)
!$acc routine seq
  implicit none
  real(rp), intent(in) :: T
  real(rp) :: cp_t, invT

  invT = 1.0_rp / T
  if (T <= T_cold) then
     cp_t = aC_adim__0 + aC_adim_p4*T**4
  else if (T <= T_switch) then
     cp_t = invT*(aL_adim_m1 + invT*aL_adim_m2) + aL_adim__0 + &
            T*(aL_adim_p1 + T*(aL_adim_p2 + T*(aL_adim_p3 + T*aL_adim_p4)))
  else
     cp_t = invT*(aH_adim_m1 + invT*aH_adim_m2) + aH_adim__0 + &
            T*(aH_adim_p1 + T*(aH_adim_p2 + T*(aH_adim_p3 + T*aH_adim_p4)))
  end if

end function heat_capacity_t

pure elemental function enthalpy_t(T) result(enthalpy)
!$acc routine seq
  implicit none
  real(rp), intent(in) :: T
  real(rp) :: enthalpy, invT

  invT = 1.0_rp / T
  if (T <= T_cold) then
     enthalpy = hC_Tref + T*(aC_adim__0 + 0.2_rp*aC_adim_p4*T**4)
  else if (T <= T_switch) then
     enthalpy = hL_Tref + invT*hL_adim_m2 + aL_adim_m1*log(T) + &
                T*(aL_adim__0 + T*(hL_adim_p1 + T*(hL_adim_p2 + T*(hL_adim_p3 + T*hL_adim_p4))))
  else
     enthalpy = hH_Tref + invT*hH_adim_m2 + aH_adim_m1*log(T) + &
                T*(aH_adim__0 + T*(hH_adim_p1 + T*(hH_adim_p2 + T*(hH_adim_p3 + T*hH_adim_p4))))
  end if

end function enthalpy_t

pure function cp_local(T) result(cp)
  implicit none
  real(rp), intent(in) :: T
  real(rp) :: cp

  if(thermally_perfect_gas > 0) then
     cp = heat_capacity_t(T)
  else
     cp = cp_Tref
  endif

end function cp_local

pure function sound_speed(T) result(a)
  implicit none
  real(rp), intent(in) :: T
  real(rp) :: a, cp_loc, gamma_loc

  if(thermally_perfect_gas > 0) then
     cp_loc    = heat_capacity_t(T)
     gamma_loc = cp_loc / (cp_loc - 1.0_rp)
  else
     gamma_loc = gamma_Tref
  endif
  a = sqrt(gamma_loc * T)

end function sound_speed

pure elemental function entropy_t(T) result(s_t)
!$acc routine seq
  implicit none
  real(rp), intent(in) :: T
  real(rp) :: s_t, invT, invT2

  invT  = 1.0_rp / T
  invT2 = invT * invT
  if (T <= T_cold) then
     s_t = sC_Tref + aC_adim__0*log(T) + 0.25_rp*aC_adim_p4*T**4
  else if (T <= T_switch) then
     s_t = sL_Tref + invT2*sL_adim_m2 + invT*sL_adim_m1 + aL_adim__0*log(T) + &
           T*(sL_adim_p1 + T*(sL_adim_p2 + T*(sL_adim_p3 + T*sL_adim_p4)))
  else
     s_t = sH_Tref + invT2*sH_adim_m2 + invT*sH_adim_m1 + aH_adim__0*log(T) + &
           T*(sH_adim_p1 + T*(sH_adim_p2 + T*(sH_adim_p3 + T*sH_adim_p4)))
  end if

end function entropy_t

elemental function inverse_temperature(intern_energy) result(Tnew)
!$acc routine seq
   implicit none
   real(rp), intent(in) :: intern_energy
   real(rp)             :: Tnew

   real(rp) :: t, dt, denom, cp_t, h_t, invT, r2
   integer,  parameter :: MAXIT = 50
   real(rp), parameter :: toll2 = 1.0e-12_rp
   real(rp), parameter :: T_min = 1.0e-2_rp
   integer :: iter

   t = intern_energy / cv_Tref

   do iter = 1, MAXIT

      invT = 1.0_rp / t
      if (t <= T_cold) then
         r2   = aC_adim_p4*t**4
         cp_t = aC_adim__0 + r2
         h_t  = hC_Tref + t*(aC_adim__0 + 0.2_rp*r2)
      else if (t <= T_switch) then
         r2   = aL_adim_m2 * invT
         cp_t = invT*(aL_adim_m1 + r2) + aL_adim__0 + &
                t*(aL_adim_p1 + t*(aL_adim_p2 + t*(aL_adim_p3 + t*aL_adim_p4)))
         h_t  = hL_Tref - r2 + aL_adim_m1*log(t) + &
                t*(aL_adim__0 + t*(hL_adim_p1 + t*(hL_adim_p2 + t*(hL_adim_p3 + t*hL_adim_p4))))
      else
         r2   = aH_adim_m2 * invT
         cp_t = invT*(aH_adim_m1 + r2) + aH_adim__0 + &
                t*(aH_adim_p1 + t*(aH_adim_p2 + t*(aH_adim_p3 + t*aH_adim_p4)))
         h_t  = hH_Tref - r2 + aH_adim_m1*log(t) + &
                t*(aH_adim__0 + t*(hH_adim_p1 + t*(hH_adim_p2 + t*(hH_adim_p3 + t*hH_adim_p4))))
      end if

      denom = 1.0_rp - cp_t
      dt    = -(intern_energy - h_t + t) / denom

      t = max(t + dt, T_min)

      if (dt*dt <= toll2) exit
   end do

   Tnew = t
end function inverse_temperature

function recovery_temperature() result(Trec)
   use mpi_module, only: rank, root
   implicit none
   real(rp) :: Trec

   real(rp), parameter :: toll_h = 1.0E-12_rp
   integer , parameter :: itmx_h = 50
   real(rp) :: rec_fac, h_inf, h_rec, h_T, cp_T, dT
   integer  :: itr
   logical  :: converged

   rec_fac = Prandtl**(1.0_rp/3.0_rp)

   Trec = 1.0_rp + 0.5_rp*rec_fac*(gamma_Tref-1.0_rp)*Mach*Mach

   if(thermally_perfect_gas > 0) then
      h_inf     = enthalpy_t(1.0_rp)
      h_rec     = h_inf + 0.5_rp*rec_fac*gamma_Tref*Mach*Mach
      converged = .false.
      dT        = huge(1.0_rp)
      do itr = 1,itmx_h
         h_T  = enthalpy_t(Trec)
         cp_T = heat_capacity_t(Trec)
         dT   = (h_rec - h_T)/cp_T
         Trec = Trec + dT
         if(abs(dT) < toll_h) then
            converged = .true.
            exit
         endif
      enddo
      if(.not. converged .and. rank == root) then
         write(*,'(A)')        ' ========== recovery_temperature: WARNING =========='
         write(*,'(A,I0,A)')   '   Newton did not converge in ', itr, ' iter'
         write(*,'(A,ES12.4)') '   residual |dT|   = ', abs(dT)
         write(*,'(A,ES12.4)') '   tolerance       = ', toll_h
         write(*,'(A)')        ' ==================================================='
      endif
   endif

end function recovery_temperature

subroutine enforce_branch_continuity(T_sw,                        &
   hP_Tref, sP_Tref,                                           &
   aP_m1,   aP_0,                                              &
   hP_m2,   hP_p1, hP_p2, hP_p3, hP_p4,                        &
   sP_m2,   sP_m1,                                             &
   sP_p1,   sP_p2, sP_p3, sP_p4,                               &
   aS_m1,   aS_0,                                              &
   hS_m2,   hS_p1, hS_p2, hS_p3, hS_p4,                        &
   sS_m2,   sS_m1,                                             &
   sS_p1,   sS_p2, sS_p3, sS_p4,                               &
   hS_Tref, sS_Tref)
   implicit none
   real(rp), intent(in)  :: T_sw, hP_Tref, sP_Tref
   real(rp), intent(in)  :: aP_m1, aP_0
   real(rp), intent(in)  :: hP_m2, hP_p1, hP_p2, hP_p3, hP_p4
   real(rp), intent(in)  :: sP_m2, sP_m1
   real(rp), intent(in)  :: sP_p1, sP_p2, sP_p3, sP_p4
   real(rp), intent(in)  :: aS_m1, aS_0
   real(rp), intent(in)  :: hS_m2, hS_p1, hS_p2, hS_p3, hS_p4
   real(rp), intent(in)  :: sS_m2, sS_m1
   real(rp), intent(in)  :: sS_p1, sS_p2, sS_p3, sS_p4
   real(rp), intent(out) :: hS_Tref, sS_Tref

   real(rp) :: invT, invT2, logT
   real(rp) :: Kh_P, Kh_S, Ks_P, Ks_S

   invT  = 1.0_rp / T_sw
   invT2 = invT * invT
   logT  = log(T_sw)

   Kh_P = invT*hP_m2 + aP_m1*logT + &
          T_sw*(aP_0 + T_sw*(hP_p1 + T_sw*(hP_p2 + T_sw*(hP_p3 + T_sw*hP_p4))))
   Kh_S = invT*hS_m2 + aS_m1*logT + &
          T_sw*(aS_0 + T_sw*(hS_p1 + T_sw*(hS_p2 + T_sw*(hS_p3 + T_sw*hS_p4))))

   Ks_P = invT2*sP_m2 + invT*sP_m1 + aP_0*logT + &
          T_sw*(sP_p1 + T_sw*(sP_p2 + T_sw*(sP_p3 + T_sw*sP_p4)))
   Ks_S = invT2*sS_m2 + invT*sS_m1 + aS_0*logT + &
          T_sw*(sS_p1 + T_sw*(sS_p2 + T_sw*(sS_p3 + T_sw*sS_p4)))

   hS_Tref = hP_Tref + Kh_P - Kh_S
   sS_Tref = sP_Tref + Ks_P - Ks_S
end subroutine enforce_branch_continuity

subroutine write_thermodynamics_profiles
   use mpi_module, only: rank, root
   implicit none
   real(rp),         parameter :: Tmin_K = 50.0_rp
   real(rp),         parameter :: Tmax_K = 6000.0_rp
   real(rp),         parameter :: dT_K   = 1.0_rp
   character(len=:), allocatable :: fname
   integer  :: unit, i, nT
   real(rp) :: T_K, T_nd, cp_nd, h_nd, s_nd

   if (thermally_perfect_gas <= 0) return
   if (rank /= root)              return

   fname = 'DATA/'//trim(data_dir)//'/thermo_profiles.txt'
   nT = nint((Tmax_K - Tmin_K)/dT_K) + 1
   open(newunit=unit, file=fname, status='replace', action='write')
   write(unit,'(A)') '# T[K]            T_nondim         cp/R             h/(R*Tref)       s/R'
   do i = 1, nT
      T_K   = Tmin_K + real(i-1, rp)*dT_K
      T_nd  = T_K / Tref
      cp_nd = heat_capacity_t(T_nd)
      h_nd  = enthalpy_t(T_nd)
      s_nd  = entropy_t(T_nd)
      write(unit,'(5(ES16.8,1X))') T_K, T_nd, cp_nd, h_nd, s_nd
   end do
   close(unit)

   write(*,'(A,A)') ' write_thermodynamics_profiles: wrote ', fname
end subroutine write_thermodynamics_profiles

end module thermodynamics_module
