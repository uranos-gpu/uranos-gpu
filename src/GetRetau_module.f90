module GetRetau_module

use parameters_module
use mpi_module
implicit none

private
public :: init_Reynolds

contains

subroutine init_Reynolds

        implicit none
        real(rp), parameter :: reTauMin = 5.0e1_rp
        real(rp), parameter :: reTauMax = 5.0e4_rp

        select case(ic)
        case('turbulent_channel')
          call check_ReTau_input(reTauMin, reTauMax)
          call GetReynolds_TCH

        case('turbulent_BL')
          call check_ReTau_input(reTauMin, reTauMax)
          call GetReynolds_TBL

        case default
          ReTau = 0.0_rp

        end select

        return
end subroutine init_Reynolds

subroutine check_ReTau_input(reTauMin, reTauMax)

        implicit none
        real(rp), intent(in) :: reTauMin, reTauMax
        logical              :: bad

        bad = (Reynolds /= Reynolds)                    &
         .or. (abs(Reynolds) > huge(1.0_rp))            &
         .or. (Reynolds <= 0.0_rp)                      &
         .or. (Reynolds < reTauMin)                     &
         .or. (Reynolds > reTauMax)

        if(.not. bad) return

        if(rank == root) then
           write(*,'(A)')        ' ============== init_Reynolds: INPUT CHECK FAILED =============='
           write(*,'(2A)')       '   ic             = ', trim(ic)
           write(*,'(A,ES12.4)') '   Reynolds (dat) = ', Reynolds
           write(*,'(A,ES10.2,A,ES10.2,A)') &
                                 '   expected ReTau in [', reTauMin,' , ',reTauMax,']'
           write(*,'(A)')        '   Please verify that a bulk or free-stream Reynolds was'
           write(*,'(A)')        '   not provided instead of ReTau.'
           write(*,'(A)')        ' ==============================================================='
        endif

        call secure_stop

        return
end subroutine check_ReTau_input

subroutine GetReynolds_TCH

        use fluid_functions_module, only: ReichardtLaw
        use thermodynamics_module , only: gamma_Tref

        implicit none

        real(rp), dimension(:), allocatable :: yTmp
        real(rp), dimension(:), allocatable :: uPlus, ucPlus, ucPlusOld
        real(rp), dimension(:), allocatable :: rhoY, tmpY
        real(rp), parameter                 :: toll  = 1.0E-14_rp
        integer , parameter                 :: itMax = 100
        real(rp), parameter                 :: sRey  = 1.1_rp
        real(rp)                            :: yP, tWall, rWall, tRec, alf
        real(rp)                            :: uu, fuu, du, uci, tEdge
        real(rp)                            :: errLoc, reOut
        real(rp)                            :: rhoSum, rhouSum, dy, uBulk
        integer                             :: n, j, iter, err = 0

        ReTau = Reynolds

        n = 1000000

        allocate(yTmp(0:n), uPlus(0:n), ucPlus(0:n),                      &
                 ucPlusOld(0:n), rhoY(0:n), tmpY(0:n), stat=err)
        if(err /= 0) stop ' Allocation error in GetReynolds_TCH'

        do j = 0,n
           yTmp(j) = j*2/real(2*n,rp)
        enddo

        uPlus = 0.0_rp
        do j = 1,n
           yP       = yTmp(j)*ReTau
           uPlus(j) = ReichardtLaw(yP)
        enddo
        ucPlus = uPlus

        tWall  = 1.0_rp
        rWall  = 1.0_rp
        alf    = sRey*Prandtl

        tRec  = T_recovery
        tEdge = tRec

        iter   = 0
        errLoc = 2*toll
        do while(errLoc > toll .and. iter < itMax)

           ucPlusOld = ucPlus
           tmpY(0)   = tWall
           rhoY(0)   = rWall
           do j = 1,n
              uu      = ucPlus(j)/ucPlus(n)
              fuu     = alf*uu + (1.0_rp - alf)*uu**2
              tmpY(j) = tWall + (tRec - tWall)*fuu + (tEdge - tRec)*uu**2
              rhoY(j) = 1.0_rp/tmpY(j)
           enddo

           do j = 1,n
              du        = uPlus(j) - uPlus(j-1)
              uci       = 0.5_rp*(sqrt(rWall/rhoY(j)) + sqrt(rWall/rhoY(j-1)))
              ucPlus(j) = ucPlus(j-1) + uci*du
           enddo

           errLoc = maxval(abs(ucPlus - ucPlusOld))
           iter   = iter + 1

        enddo

        rhoSum  = 0.0_rp
        rhouSum = 0.0_rp
        do j = 1,n
           dy      = yTmp(j) - yTmp(j-1)
           rhoSum  = rhoSum  + rhoY(j)*dy
           rhouSum = rhouSum + rhoY(j)*ucPlus(j)*dy
        enddo
        uBulk = rhouSum/rhoSum

        reOut    = ReTau*uBulk*rhoSum
        Reynolds = reOut
        mu_inf   = sqrt(gamma_Tref)*Mach/Reynolds

        deallocate(yTmp, uPlus, ucPlus, ucPlusOld, rhoY, tmpY)

        return
end subroutine GetReynolds_TCH

subroutine GetReynolds_TBL

        use fluid_functions_module, only: MuskerProfile, CompressibleCorrection, &
                                          laminar_viscosity
        use thermodynamics_module , only: gamma_Tref

        implicit none

        integer , parameter                 :: nAux        = 100000
        real(rp), parameter                 :: yMaxAux     = 10.0_rp
        real(rp), parameter                 :: yPlusTarget = 0.001_rp

        real(rp), dimension(:), allocatable :: yTmp, yMusker
        real(rp), dimension(:), allocatable :: uPlusInc, uPlusCmp
        real(rp), dimension(:), allocatable :: rhoY, tmpY
        real(rp)                            :: tWall, rWall, mWall, uEdge
        integer                             :: err = 0

        ReTau = Reynolds

        allocate(yTmp(0:nAux), yMusker(0:nAux),                              &
                 uPlusInc(0:nAux), uPlusCmp(0:nAux),                         &
                 rhoY(0:nAux), tmpY(0:nAux), stat=err)
        if(err /= 0) stop ' Allocation error in GetReynolds_TBL'

        call BuildAuxYGrid(nAux, yMaxAux, ReTau, yPlusTarget, yTmp)

        yMusker = yTmp + ymin
        call MuskerProfile(nAux, yMusker, ReTau, uPlusInc)

        call CompressibleCorrection(nAux, uPlusInc, uPlusCmp, rhoY, tmpY, tWall, rWall)

        uEdge = uPlusCmp(nAux)
        mWall = laminar_viscosity(tWall, Tref, vis_flag)

        Reynolds = ReTau*uEdge/rWall*mWall
        mu_inf   = sqrt(gamma_Tref)*Mach/Reynolds

        deallocate(yTmp, yMusker, uPlusInc, uPlusCmp, rhoY, tmpY)

        return
end subroutine GetReynolds_TBL

subroutine BuildAuxYGrid(n, Ly, reTauIn, yPlusTarget, y)

        implicit none
        integer , intent(in)  :: n
        real(rp), intent(in)  :: Ly, reTauIn, yPlusTarget
        real(rp), intent(out) :: y(0:n)

        real(rp), parameter :: toll  = 1.0E-15_rp
        integer , parameter :: itMax = 1000
        real(rp)            :: alpha, alphaOld, xi1, errLoc, kRatio
        integer             :: iter, j

        xi1 = 1.0_rp/real(n,rp)

        kRatio = reTauIn*Ly*xi1/yPlusTarget

        if(kRatio <= 1.0_rp) then
           alpha = 0.0_rp
        else
           alpha  = 1.0_rp
           errLoc = 2.0_rp*toll
           iter   = 0
           do while(errLoc > toll .and. iter < itMax)
              alphaOld = alpha
              alpha    = asinh(reTauIn*Ly*sinh(alphaOld*xi1)/yPlusTarget)
              errLoc   = abs(alpha - alphaOld)
              iter     = iter + 1
           enddo
           if(alpha <= 0.0_rp) stop ' BuildAuxYGrid: stretching search failed'
        endif

        y(0) = 0.0_rp
        if(alpha == 0.0_rp) then
           do j = 1,n
              y(j) = Ly*real(j,rp)/real(n,rp)
           enddo
        else
           do j = 1,n
              y(j) = Ly*sinh(alpha*real(j,rp)/real(n,rp))/sinh(alpha)
           enddo
        endif

        return
end subroutine BuildAuxYGrid

end module GetRetau_module
