module ic_module
use parameters_module
use storage_module
use mpi_module
use fluid_functions_module
use thermodynamics_module, only: gamma_Tref, heat_capacity_t, sound_speed, cp_local
implicit none

contains

subroutine ambient(status)
        implicit none
        integer, intent(out) :: status

        if(rank == root) write(*,'(A)', advance = 'no') ' Initializing ambient initial condition: '

        do k = lbz,ubz
           do j = lby,uby
              do i = lbx,ubx

                 PHI(i,j,k,1) = 1.0_rp
                 U  (i,j,k)   = 0.0_rp
                 V  (i,j,k)   = 0.0_rp
                 W  (i,j,k)   = 0.0_rp
                 P  (i,j,k)   = 1.0_rp

              enddo
           enddo
        enddo

        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return

end subroutine ambient

subroutine TGV(status)
        use mesh_module, only: xx, yy
        implicit none
        integer, intent(out) :: status
        real(rp)             :: x_, y_, z_, rho_z, u_z, p_z

        if(rank == root) write(*,'(A,A,A)', advance = 'no') 'Initializing TGV: '

        p_z   = 100.0_rp
        rho_z = 1.0_rp
        u_z   = Mach*sqrt(gamma_Tref*p_z/rho_z)

        do k = lbz,ubz
           do j = lby,uby
              do i = lbx,ubx

                 x_ = xx(i,j)
                 y_ = yy(i,j)
                 z_ = z(k)

                 P  (i,j,k)   =   p_z + u_z*u_z*(cos(2*x_)+cos(2*y_))*(2.0_rp+cos(2*z_))/16.0_rp
                 PHI(i,j,k,1) =   P(i,j,k)*rho_z/p_z
                 U  (i,j,k)   =   u_z*sin(x_)*cos(y_)*cos(z_)
                 V  (i,j,k)   = - u_z*cos(x_)*sin(y_)*cos(z_)
                 W  (i,j,k)   =   0.0_rp

              enddo
           enddo
        enddo

        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return

end subroutine TGV

subroutine init_periodic_euler(dir,status)
        use mesh_module, only: xx, yy
        implicit none
        character(len=1), intent(in)  :: dir
        integer         , intent(out) :: status

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing periodic euler along ', dir, ': '

        select case(dir)
        case('x')
           do k = lbz,ubz
              do j = lby,uby
                 do i = lbx,ubx

                    PHI(i,j,k,1) = 1.0_rp + 0.2_rp*sin(pi*xx(i,j))
                    U  (i,j,k)   = Mach
                    V  (i,j,k)   = 0.0_rp
                    W  (i,j,k)   = 0.0_rp
                    P  (i,j,k)   = 1.0_rp

                 enddo
              enddo
           enddo

        case('y')
           do k = lbz,ubz
              do j = lby,uby
                 do i = lbx,ubx

                    PHI(i,j,k,1) = 1.0_rp + 0.2_rp*sin(pi*yy(i,j))
                    U  (i,j,k)   = 0.0_rp
                    V  (i,j,k)   = Mach
                    W  (i,j,k)   = 0.0_rp
                    P  (i,j,k)   = 1.0_rp

                 enddo
              enddo
           enddo

        case('z')
           do k = lbz,ubz
              do j = lby,uby
                 do i = lbx,ubx

                    PHI(i,j,k,1) = 1.0_rp + 0.2_rp*sin(pi*z(k))
                    U  (i,j,k)   = 0.0_rp
                    V  (i,j,k)   = 0.0_rp
                    W  (i,j,k)   = Mach
                    P  (i,j,k)   = 1.0_rp

                 enddo
              enddo
           enddo

        endselect

        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return

end subroutine init_periodic_euler

subroutine init_shock_tube(dir, status)
        implicit none
        character(len=1), intent(in)  :: dir
        integer         , intent(out) :: status

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing shock tube along ', dir, ': '

        select case(dir)
        case('x')
           do k = lbz,ubz
              do j = lby,uby
                 do i = lbx,ubx

                    if(x(i) <= 0.0_rp) then
                       PHI(i,j,k,1) = 1.0_rp
                       U  (i,j,k)   = 0.0_rp
                       V  (i,j,k)   = 0.0_rp
                       W  (i,j,k)   = 0.0_rp
                       P  (i,j,k)   = 1.0_rp
                    else
                       PHI(i,j,k,1) = 0.125_rp
                       U  (i,j,k)   = 0.0_rp
                       V  (i,j,k)   = 0.0_rp
                       W  (i,j,k)   = 0.0_rp
                       P  (i,j,k)   = 0.1_rp
                    endif

                 enddo
              enddo
           enddo

        case('y')
           do k = lbz,ubz
              do j = lby,uby
                 do i = lbx,ubx

                    if(y(j) <= 0.0_rp) then
                       PHI(i,j,k,1) = 1.0_rp
                       U  (i,j,k)   = 0.0_rp
                       V  (i,j,k)   = 0.0_rp
                       W  (i,j,k)   = 0.0_rp
                       P  (i,j,k)   = 1.0_rp
                    else
                       PHI(i,j,k,1) = 0.125_rp
                       U  (i,j,k)   = 0.0_rp
                       V  (i,j,k)   = 0.0_rp
                       W  (i,j,k)   = 0.0_rp
                       P  (i,j,k)   = 0.1_rp
                    endif

                 enddo
              enddo
           enddo

        case('z')
           do k = lbz,ubz
              do j = lby,uby
                 do i = lbx,ubx

                    if(z(k) <= 0.0_rp) then
                       PHI(i,j,k,1) = 1.0_rp
                       U  (i,j,k)   = 0.0_rp
                       V  (i,j,k)   = 0.0_rp
                       W  (i,j,k)   = 0.0_rp
                       P  (i,j,k)   = 1.0_rp
                    else
                       PHI(i,j,k,1) = 0.125_rp
                       U  (i,j,k)   = 0.0_rp
                       V  (i,j,k)   = 0.0_rp
                       W  (i,j,k)   = 0.0_rp
                       P  (i,j,k)   = 0.1_rp
                    endif

                 enddo
              enddo
           enddo

        endselect

        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return

end subroutine init_shock_tube

subroutine init_isentropic_vortex(dir,status)
        use mesh_module, only: xx, yy
        implicit none
        character(len=1), intent(in ) :: dir
        integer         , intent(out) :: status
        real(rp)                      :: x0, y0, alpha, beta, r2
        real(rp)                      :: exp_, x_, y_
        real(rp)                      :: gm1, exp_rho, exp_p

        x0   = 0.0_rp
        y0   = 0.0_rp
        beta = 5.0_rp

        gm1     = gamma_Tref - 1.0_rp
        exp_rho = 1.0_rp     / gm1
        exp_p   = gamma_Tref / gm1
        alpha   = beta*beta * gm1 / (8.0_rp*gamma_Tref*pi*pi)

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing isentropic vortex along ', dir, ': '

        select case(dir)
        case('x')
           do k = lbz,ubz
              do j = lby,uby
                 do i = lbx,ubx

                    x_   = xx(i,j) - x0
                    y_   = yy(i,j) - y0
                    r2   = x_*x_ + y_*y_
                    exp_ = exp(1.0_rp - r2)

                    PHI(i,j,k,1) = (1.0_rp - alpha*exp_)**exp_rho
                    U  (i,j,k)   = u_inf - beta/(2.0_rp*pi) * y_ * sqrt(exp_)
                    V  (i,j,k)   =         beta/(2.0_rp*pi) * x_ * sqrt(exp_)
                    W  (i,j,k)   = 0.0_rp
                    P  (i,j,k)   = (1.0_rp - alpha*exp_)**exp_p

                 enddo
              enddo
           enddo

        case('y')
           do k = lbz,ubz
              do j = lby,uby
                 do i = lbx,ubx

                    x_   = xx(i,j) - x0
                    y_   = yy(i,j) - y0
                    r2   = x_*x_ + y_*y_
                    exp_ = exp(1.0_rp - r2)

                    PHI(i,j,k,1) = (1.0_rp - alpha*exp_)**exp_rho
                    U  (i,j,k)   =         beta/(2.0_rp*pi) * y_ * sqrt(exp_)
                    V  (i,j,k)   = u_inf - beta/(2.0_rp*pi) * x_ * sqrt(exp_)
                    W  (i,j,k)   = 0.0_rp
                    P  (i,j,k)   = (1.0_rp - alpha*exp_)**exp_p

                 enddo
              enddo
           enddo

        endselect

        if(rank == root) write(*, '(A)') 'done!'
        status = 0
        return

end subroutine init_isentropic_vortex

subroutine init_turbulent_channel(status)
        implicit none
        integer, intent(out) :: status

        integer , parameter :: nAux   =  100000
        integer , parameter :: nBand  =       4
        real(rp), parameter :: ypStrk =    15.0_rp
        real(rp), parameter :: rollAR =     5.0_rp
        real(rp), parameter :: lamZp  =   100.0_rp
        real(rp), parameter :: lamXp  =  1000.0_rp
        real(rp), parameter :: lamZot =     3.0_rp
        real(rp), parameter :: oNorm  = 3.4938562148_rp
        real(rp), parameter :: aStrk  =     2.7_rp
        real(rp), parameter :: aRoll  =     1.0_rp
        real(rp), parameter :: wiggle =     0.6_rp
        real(rp), parameter :: sRey   =     1.1_rp
        real(rp), parameter :: tollVD = 1.0E-14_rp
        integer , parameter :: itMxVD =     100

        real(rp), dimension(nBand), parameter :: zBand = [0.7_rp, 1.0_rp, 1.4_rp, 2.0_rp]
        real(rp), dimension(nBand), parameter :: xBand = [1.0_rp, 2.0_rp, 1.0_rp, 2.0_rp]

        integer , parameter :: nMode = 2*nBand + 1
        real(rp), dimension(nMode) :: betaM, alphaM, phiM, chiM, yprM
        integer , dimension(nMode) :: wallM
        integer , dimension(nBand) :: nzB, nxB

        real(rp), dimension(:), allocatable :: yAux, upAux, ucAux, ucOld, tAux, rAux
        real(rp), dimension(:), allocatable :: uMean, rMean
        real(rp), dimension(:), allocatable :: cosBz, sinBz, cosAx, sinAx
        real(rp), dimension(:), allocatable :: fRoll, dRoll, gStrk

        real(rp) :: nuWall, uTau, tWall, rWall, tRec, alfVD, errVD
        real(rp) :: uu, fuu, du, uci, gold, silv, wDist, tPos, wgt
        real(rp) :: ypl, ypr, dypdy, ang, csT, snT, amp, bandNorm, lamZm
        real(rp) :: lclRU, lclVol, gblRU, gblVol, uBulk, scal, dv
        integer  :: imode, iband, iw, jtab, mcount, iterVD, nzM, nxM, ierrA

        if(rank == root) write(*,'(A)', advance = 'no') &
        ' Initializing turbulent channel: '

        status = 0

        if(dims /= 3) then
           if(rank == root) print*, ' turbulent_channel requires dims = 3 (Lz unset in 2D)'
           call secure_stop
        endif

        nuWall = mu_inf
        uTau   = ReTau*nuWall

        allocate(yAux (0:nAux), upAux(0:nAux), ucAux(0:nAux), &
                 ucOld(0:nAux), tAux (0:nAux), rAux (0:nAux), stat = ierrA)
        if(ierrA /= 0) then
           print*, ' Allocation error in init_turbulent_channel, rank ', rank
           call secure_stop
        endif

        do jtab = 0,nAux
           yAux (jtab) = real(jtab,rp)/real(nAux,rp)
           upAux(jtab) = ReichardtLaw(yAux(jtab)*ReTau)
        enddo
        ucAux = upAux

        tWall = 1.0_rp
        rWall = 1.0_rp
        tRec  = T_recovery
        alfVD = sRey*Prandtl

        iterVD = 0
        errVD  = 2.0_rp*tollVD
        do while(errVD > tollVD .and. iterVD < itMxVD)

           ucOld     = ucAux
           tAux(0)   = tWall
           rAux(0)   = rWall
           do jtab = 1,nAux
              uu         = ucAux(jtab)/ucAux(nAux)
              fuu        = alfVD*uu + (1.0_rp - alfVD)*uu**2
              tAux(jtab) = tWall + (tRec - tWall)*fuu
              rAux(jtab) = 1.0_rp/tAux(jtab)
           enddo

           do jtab = 1,nAux
              du          = upAux(jtab) - upAux(jtab-1)
              uci         = 0.5_rp*(sqrt(rWall/rAux(jtab)) + sqrt(rWall/rAux(jtab-1)))
              ucAux(jtab) = ucAux(jtab-1) + uci*du
           enddo

           errVD  = maxval(abs(ucAux - ucOld))
           iterVD = iterVD + 1

        enddo

        allocate(uMean(lby:uby), rMean(lby:uby), stat = ierrA)
        if(ierrA /= 0) then
           print*, ' Allocation error in init_turbulent_channel, rank ', rank
           call secure_stop
        endif

        do j = lby,uby
           uMean(j) = 0.0_rp
           rMean(j) = 1.0_rp
           if(abs(y(j)) < 1.0_rp) then
              wDist    = 1.0_rp - abs(y(j))
              tPos     = wDist*real(nAux,rp)
              jtab     = min(nAux-1, max(0, int(tPos)))
              wgt      = tPos - real(jtab,rp)
              uMean(j) = uTau*((1.0_rp - wgt)*ucAux(jtab) + wgt*ucAux(jtab+1))
              rMean(j) =       (1.0_rp - wgt)*rAux (jtab) + wgt*rAux (jtab+1)
           endif
        enddo

        uu = 0.0_rp
        do jtab = 1,nAux
           uu = uu + rAux(jtab)*ucAux(jtab)
        enddo
        uu = uu/real(nAux,rp)
        du = Reynolds/ReTau

        deallocate(yAux, upAux, ucAux, ucOld, tAux, rAux)

        do k       = lbz,ubz
           do j    = lby,uby
              do i = lbx,ubx
                 phi(i,j,k,1) = rMean(j)
                 U  (i,j,k)   = uMean(j)
                 V  (i,j,k)   = 0.0_rp
                 W  (i,j,k)   = 0.0_rp
                 P  (i,j,k)   = 1.0_rp
              enddo
           enddo
        enddo

        gold = 0.5_rp*(sqrt(5.0_rp) - 1.0_rp)
        silv =         sqrt(2.0_rp) - 1.0_rp

        do iband = 1,nBand
           nzB(iband) = max(1, min(nz/4, nint(Lz*ReTau/(lamZp*zBand(iband)))))
           nxB(iband) = max(1, min(nx/4, nint(Lx*ReTau/(lamXp*xBand(iband)))))
           if(iband > 1) nzB(iband) = max(1, min(nzB(iband), nzB(iband-1) - 1))
        enddo

        if(rank == root .and. nzB(2) /= max(1, nint(Lz*ReTau/lamZp))) then
           write(*,'(/,A)') '  WARNING: spanwise seed clamped by the mesh -- the near-wall'
           write(*,'(A)')   '           cycle is NOT resolved on this grid; the lambda_z+'
           write(*,'(A)')   '           reported below is what the grid carries, not what'
           write(*,'(A)')   '           was asked for.'
        endif

        mcount = 0
        do iw = 1,2
           do iband = 1,nBand
              mcount = mcount + 1
              nzM    = nzB(iband)
              nxM    = nxB(iband)
              lamZm  = Lz*ReTau/real(nzM,rp)
              betaM (mcount) = 2.0_rp*pi*real(nzM,rp)/Lz
              alphaM(mcount) = 2.0_rp*pi*real(nxM,rp)/Lx
              yprM  (mcount) = lamZm/rollAR
              wallM (mcount) = iw
           enddo
        enddo

        mcount = mcount + 1
        nzM    = max(1, nint(Lz/lamZot))
        nxM    = max(1, nint(Lx*ReTau/lamXp))
        betaM (mcount) = 2.0_rp*pi*real(nzM,rp)/Lz
        alphaM(mcount) = 2.0_rp*pi*real(nxM,rp)/Lx
        yprM  (mcount) = 0.0_rp
        wallM (mcount) = 0

        bandNorm = sqrt(2.0_rp/real(nBand,rp))

        do imode = 1,nMode
           phiM(imode) = 2.0_rp*pi*(real(imode,rp)*gold - floor(real(imode,rp)*gold))
           chiM(imode) = 2.0_rp*pi*(real(imode,rp)*silv - floor(real(imode,rp)*silv))
        enddo

        allocate(cosBz(lbz:ubz), sinBz(lbz:ubz), cosAx(lbx:ubx), sinAx(lbx:ubx), &
                 fRoll(lby:uby), dRoll(lby:uby), gStrk(lby:uby), stat = ierrA)
        if(ierrA /= 0) then
           print*, ' Allocation error in init_turbulent_channel, rank ', rank
           call secure_stop
        endif

        do imode = 1,nMode

           do k = lbz,ubz
              ang      = betaM(imode)*z(k) + phiM(imode)
              cosBz(k) = cos(ang)
              sinBz(k) = sin(ang)
           enddo

           do i = lbx,ubx
              ang      = wiggle*sin(alphaM(imode)*x(i) + chiM(imode))
              cosAx(i) = cos(ang)
              sinAx(i) = sin(ang)
           enddo

           do j = lby,uby

              fRoll(j) = 0.0_rp
              dRoll(j) = 0.0_rp
              gStrk(j) = 0.0_rp
              if(abs(y(j)) >= 1.0_rp) cycle

              if(wallM(imode) == 0) then

                 fRoll(j) = -oNorm*y(j)*(1.0_rp - y(j)**2)**2
                 dRoll(j) =  oNorm*(1.0_rp - y(j)**2)*(5.0_rp*y(j)**2 - 1.0_rp)
                 gStrk(j) =   1.0_rp - y(j)**2

              else

                 if(wallM(imode) == 1) then
                    ypl   = (1.0_rp + y(j))*ReTau
                    dypdy =  ReTau
                 else
                    ypl   = (1.0_rp - y(j))*ReTau
                    dypdy = -ReTau
                 endif

                 ypr      = yprM(imode)
                 fRoll(j) = (ypl/ypr)**2*exp(2.0_rp - 2.0_rp*ypl/ypr)
                 dRoll(j) = dypdy*(2.0_rp/ypr)*(ypl/ypr)*(1.0_rp - ypl/ypr) &
                          * exp(2.0_rp - 2.0_rp*ypl/ypr)
                 gStrk(j) = (ypl/ypStrk)*exp(1.0_rp - ypl/ypStrk)

              endif

           enddo

           amp = aRoll*bandNorm*uTau
           if(wallM(imode) == 2) amp = -amp

           do k       = lbz,ubz
              do j    = lby,uby
                 do i = lbx,ubx

                    csT = cosBz(k)*cosAx(i) - sinBz(k)*sinAx(i)
                    snT = sinBz(k)*cosAx(i) + cosBz(k)*sinAx(i)

                    U(i,j,k) = U(i,j,k) + aStrk*bandNorm*uTau*gStrk(j)*sinBz(k)

                    V(i,j,k) = V(i,j,k) - amp*fRoll(j)*snT
                    W(i,j,k) = W(i,j,k) - amp/betaM(imode)*dRoll(j)*csT

                 enddo
              enddo
           enddo

        enddo

        deallocate(cosBz, sinBz, cosAx, sinAx, fRoll, dRoll, gStrk)
        deallocate(uMean, rMean)

        lclRU  = 0.0_rp
        lclVol = 0.0_rp
        do k       = sz,ez
           do j    = sy,ey
              do i = sx,ex
                 if(abs(y(j)) < 1.0_rp) then
                    dv     = xstep(i)*ystep(j)*zstep(k)
                    lclRU  = lclRU  + phi(i,j,k,1)*U(i,j,k)*dv
                    lclVol = lclVol + dv
                 endif
              enddo
           enddo
        enddo

        call MPI_allreduce(lclRU ,gblRU ,1,MPI_RP,mpi_sum,mpi_comm_cart,mpi_err)
        call MPI_allreduce(lclVol,gblVol,1,MPI_RP,mpi_sum,mpi_comm_cart,mpi_err)

        uBulk = gblRU/gblVol
        scal  = u_inf/uBulk

        do k       = lbz,ubz
           do j    = lby,uby
              do i = lbx,ubx
                 U(i,j,k) = scal*U(i,j,k)
              enddo
           enddo
        enddo

        if(rank == root) then
           write(*,'(A)') 'done!'
           write(*,'(A,ES13.6)') '    u_tau                   : ', uTau
           write(*,'(A,F10.3,A,F9.3)') &
                                 '    bulk u+  table/GetRetau : ', uu, '  /', du
           write(*,'(A,F10.6)')  '    mass-flux rescaling     : ', scal
           write(*,'(A,F10.2,A,F7.2,A)') '    streak lambda_z+        : ', &
                                 Lz*ReTau/real(nzB(2),rp), '   (', real(nz,rp)/real(nzB(2),rp), ' pts/wave)'
           write(*,'(A,F10.2,A,F7.2,A)') '    sinuous lambda_x+       : ', &
                                 Lx*ReTau/real(nxB(1),rp), '   (', real(nx,rp)/real(nxB(1),rp), ' pts/wave)'
           if(abs(uu - du) > 1.0E-3_rp*du) write(*,'(A)') &
           '  WARNING: local mean-profile table disagrees with GetReynolds_TCH.'
        endif

        return
end subroutine init_turbulent_channel

subroutine init_TurbulentBoundaryLayerMean(status)

        use fluid_functions_module, only: MuskerProfile, CompressibleCorrection, &
                                          BoundaryLayerQuantities, laminar_viscosity
        use mesh_module           , only: west_wall_normal_distance, &
                                          bottom_wall_streamwise_distance, &
                                          yy, dxi_dx, deta_dx, nx_xi, ny_xi
        use real_to_integer_module, only: locate
        use interpolation_module  , only: polint
        use FileModule

        implicit none
        integer, intent(out) :: status

        real(rp), dimension(:,:,:), allocatable :: tempPrim
        real(rp), dimension(:)    , allocatable :: delta_array, deltaVArray
        real(rp), dimension(:)    , allocatable :: theta_array, frict_array
        real(rp), dimension(:)    , allocatable :: y_tmp

        real(rp), dimension(0:ny) :: UPlusInc, UPlusCmp, RhoY, TmpY, velY
        real(rp), dimension(0:ny) :: wdist
        real(rp), dimension(0:ny) :: dru_dx

        real(rp), dimension(-GN:nx+GN) :: sdist

        real(rp), parameter :: toll   = 1.0E-14_rp
        integer , parameter :: itrmax = 30
        real(rp) :: tWall, rWall
        real(rp) :: delta, deltaV, thrat, deltaStar, theta, cf
        real(rp) :: dy, ReTau_new, ReTau_old, err
        real(rp) :: ue, ReTauI, yl, errR, errU
        real(rp) :: ri, ui, drxi, dreta
        real(rp) :: r_, ir, u_, v_
        real(rp) :: u_y, mu_wall
        integer  :: i, j, jj, jjj, m, jg, ji
        integer  :: itr, lbx_loc, k, isign, nbr

        type(FileType) :: blData, inflow, frictionFile

        if(rank == root) write(*,'(A)', advance='no') ' Initializing Turbulent boundary layer: '

        call AllocateReal(delta_array, -GN,nx+GN)
        call AllocateReal(deltaVArray, -GN,nx+GN)
        call AllocateReal(theta_array, -GN,nx+GN)
        call AllocateReal(frict_array, -GN,nx+GN)

        lbx_loc = lbx - 1
        call AllocateReal(tempPrim, lbx_loc,ubx, lby,uby, 1,3)

        call AllocateReal(y_tmp, 0,ny)
        call west_wall_normal_distance(wdist)
        y_tmp(0:ny) = ymin + wdist(0:ny)

        call bottom_wall_streamwise_distance(sdist)

        call MuskerProfile         (ny, y_tmp, ReTau, UPlusInc)
        call CompressibleCorrection(ny, UPlusInc, UPlusCmp, RhoY, TmpY, tWall, rWall)
        call BoundaryLayerQuantities(y_tmp, RhoY, UPlusCmp, rWall, thrat, deltaStar, cf)

        delta_array(0) = 1.0_rp
        deltaVArray(0) = 1.0_rp/ReTau
        frict_array(0) = cf
        theta_array(0) = thrat

        ReTau_new = ReTau
        ReTau_old = ReTau - 1.0_rp

        sweep_k: do k = 1, max(GN, nx+GN)
           sweep_dir: do isign = 1, -1, -2

              i = isign*k
              if (i > nx+GN .or. i < -GN) cycle sweep_dir
              nbr = i - isign

              itr = 0
              err = huge(1.0_rp)
              picard: do while(err > toll .and. itr < itrmax)

                 call MuskerProfile         (ny, y_tmp, ReTau_new, UPlusInc)
                 call CompressibleCorrection(ny, UPlusInc, UPlusCmp, RhoY, TmpY, tWall, rWall)
                 call BoundaryLayerQuantities(y_tmp, RhoY, UPlusCmp, rWall, thrat, deltaStar, cf)

                 theta     = theta_array(nbr) + 0.25_rp*(sdist(i)-sdist(nbr))*(cf + frict_array(nbr))
                 delta     = theta/thrat
                 deltaV    = deltaVArray(nbr)*sqrt(frict_array(nbr)/cf)
                 ReTau_new = delta/deltaV

                 err       = abs(ReTau_new - ReTau_old)
                 ReTau_old = ReTau_new
                 itr       = itr + 1

              enddo picard

              theta_array(i) = theta
              delta_array(i) = delta
              deltaVArray(i) = deltaV
              frict_array(i) = cf

           enddo sweep_dir
        enddo sweep_k

        do i = lbx_loc, ubx

           delta  = delta_array(i)
           deltaV = deltaVArray(i)
           ReTauI = delta/deltaV

           call MuskerProfile         (ny, y_tmp, ReTauI, UPlusInc)
           call CompressibleCorrection(ny, UPlusInc, UPlusCmp, RhoY, TmpY, tWall, rWall)

           ue   = maxval(UPlusCmp)
           velY = UPlusCmp/ue

           do j = 0, ny
              yl  = wdist(j)/delta
              call locate(wdist, 0, ny, yl, jj)

              m   = 4
              jjj = min(max(jj-(m-1)/2, 1), ny+1-m)
              call polint(wdist(jjj), RhoY(jjj), m, yl, ri, errR)
              call polint(wdist(jjj), velY(jjj), m, yl, ui, errU)

              tempPrim(i,j,1) = ri
              tempPrim(i,j,2) = ri*ui
              tempPrim(i,j,3) = 0.0_rp
           enddo

        enddo

        do i = lbx, ubx

           do j = 0, ny
              drxi  = tempPrim(i,j,2) - tempPrim(i-1,j,2)
              if    (j == 0 ) then
                 dreta = tempPrim(i,1 ,2) - tempPrim(i,0   ,2)
              elseif(j == ny) then
                 dreta = tempPrim(i,ny,2) - tempPrim(i,ny-1,2)
              else
                 dreta = 0.5_rp*(tempPrim(i,j+1,2) - tempPrim(i,j-1,2))
              endif
              dru_dx(j) = dxi_dx(i,j)*drxi + deta_dx(i,j)*dreta
           enddo

           do j = 1, ny
              dy = yy(i,j) - yy(i,j-1)
              tempPrim(i,j,3) = tempPrim(i,j-1,3) - 0.5_rp*(dru_dx(j)+dru_dx(j-1))*dy
           enddo

        enddo

        do    j = 1,   GN
           do i = lbx, ubx
              jg = ny + j
              ji = ny - (j-1)
              tempPrim(i,jg,1) =  tempPrim(i,ji,1)
              tempPrim(i,jg,2) =  tempPrim(i,ji,2)
              tempPrim(i,jg,3) =  tempPrim(i,ji,3)
           enddo
        enddo

        do    j = 1,   GN
           do i = lbx, ubx
              jg = 1 - j
              ji = 1 + (j-1)
              tempPrim(i,jg,1) =  tempPrim(i,ji,1)
              tempPrim(i,jg,2) = -tempPrim(i,ji,2)
              tempPrim(i,jg,3) = -tempPrim(i,ji,3)
           enddo
        enddo

        do       k = lbz, ubz
           do    j = lby, uby
              do i = lbx, ubx
                 r_ = tempPrim(i,j,1)
                 ir = 1.0_rp/r_
                 u_ = tempPrim(i,j,2)*ir
                 v_ = tempPrim(i,j,3)*ir

                 phi(i,j,k,1) = r_
                 U  (i,j,k)   = u_inf*(nx_xi(i,j)*u_ - ny_xi(i,j)*v_)
                 V  (i,j,k)   = u_inf*(ny_xi(i,j)*u_ + nx_xi(i,j)*v_)
                 W  (i,j,k)   = 0.0_rp
                 P  (i,j,k)   = 1.0_rp
              enddo
           enddo
        enddo

        if(rank == root) then
           call MuskerProfile         (ny, y_tmp, ReTau, UPlusInc)
           call CompressibleCorrection(ny, UPlusInc, UPlusCmp, RhoY, TmpY, tWall, rWall)

           inflow%name = 'inflowProfile'
           inflow%dir  = trim(data_dir)//'/INITIAL_CONDITION'
           call OpenNewFile(inflow, it)
           do j = 0, ny
              write(inflow%unit,*) y_tmp(j), RhoY(j), UPlusCmp(j)/maxval(UPlusCmp), TmpY(j)
           enddo
           call CloseFile(inflow)

           blData%name = 'blData'
           blData%dir  = trim(data_dir)//'/INITIAL_CONDITION'
           call OpenNewFile(blData, it)
           write(blData%unit,*) '# s_wall    theta   delta   delta_v   cf   ReTau'
           do i = -GN, nx+GN
              write(blData%unit,*) sdist(i), theta_array(i), delta_array(i), &
                                   deltaVArray(i), frict_array(i), delta_array(i)/deltaVArray(i)
           enddo
           call CloseFile(blData)
        endif

        frictionFile%name = 'InitFriction'//trim(str(rank))
        frictionFile%dir  = trim(data_dir)//'/INITIAL_CONDITION'
        call OpenNewFile(frictionFile, it)
        mu_wall = laminar_viscosity(tWall, Tref, vis_flag)
        do i = lbx, ubx
           u_y = U(i,sy,sz)/y_tmp(1)
           write(frictionFile%unit,*) x(i), mu_wall*mu_inf*u_y/q_inf
        enddo
        call CloseFile(frictionFile)

        call DeallocateReal(delta_array)
        call DeallocateReal(deltaVArray)
        call DeallocateReal(theta_array)
        call DeallocateReal(frict_array)
        call DeallocateReal(tempPrim)
        call DeallocateReal(y_tmp)

        if(rank == root) write(*,'(A)') 'done!'
        status = 0

        return
end subroutine init_TurbulentBoundaryLayerMean

function tbl_noise_temperature(Tm, us_m, us_p, j) result(T_)
        use df_module, only: DF_LundMatrix

        implicit none
        real(rp), intent(in) :: Tm, us_m, us_p
        integer , intent(in) :: j
        real(rp)             :: T_, dTdu, th, sg

        dTdu = GRATemperatureSlope(us_m/u_inf, Tm)/u_inf
        th   = dTdu*us_p/Tm
        sg   = dTdu*u_inf*DF_LundMatrix(1,1,max(0,min(ny,j)))/Tm
        T_   = Tm*exp(th - 0.5_rp*sg*sg)

        return
end function tbl_noise_temperature

subroutine init_TBLNoise_fast(status)
        use df_module
        use random_module
        use mesh_module, only: nx_xi, ny_xi

        implicit none

        integer, intent(out) :: status

        real(rp), dimension(:,:,:), allocatable :: rnd_tmp, fy, vf, uf, uf_old

        real(rp), dimension(3) :: mean, rmsq, rms
        real(rp), dimension(3) :: exppar1, exppar2, par1, par2
        real(rp) :: r, rm, um, vm, wm, Tm, up, vp, wp, T_, dx
        real(rp) :: us_p, vn_p, us_m
        integer :: j,k,m, jj, kk, l
        integer :: req_recv, plane_count

        if(rank == root) write(*,'(A,A,A)', advance = 'no') ' Initializing a turbulent-consistent field: '

        call AllocateReal(rnd_tmp, 1,3, 1-DF_N,ny+DF_N, 1-DF_N,nz+DF_N)
        call AllocateReal(fy     , 1,3, 1-DF_N,ny+DF_N, 1-DF_N,nz+DF_N)
        call AllocateReal(vf     , 1,3, 1,ny, 1,nz)
        call AllocateReal(uf     , 1,3, 1-GN,ny+GN, 1-GN,nz+GN)
        call AllocateReal(uf_old , 1,3, 1-GN,ny+GN, 1-GN,nz+GN)

        uf_old = 0.0_rp

        plane_count = 3*(ny+2*GN)*(nz+2*GN)
        req_recv = MPI_REQUEST_NULL
        if (my_neighbour(1) /= MPI_PROC_NULL) then
           call MPI_IRECV(uf_old, plane_count, MPI_RP, my_neighbour(1), &
                          700, mpi_comm_cart, req_recv, mpi_err)
           call MPI_WAIT(req_recv, MPI_STATUS_IGNORE, mpi_err)
        endif

        do i = sx,ex

           do k = 1,nz
              do j = 1-DF_N,ny+DF_N
                 do m = 1,3

                    call random_normal(r)
                    rnd_tmp(m,j,k) = r

                 enddo
              enddo
           enddo

           do j = 1-DF_N,ny+DF_N

              mean = 0.0_rp
              rmsq = 0.0_rp
              do k = 1,nz
                 do m = 1,3
                    mean(m) = mean(m) + rnd_tmp(m,j,k)
                    rmsq(m) = rmsq(m) + rnd_tmp(m,j,k)**2
                 enddo
              enddo
              mean = mean/real(nz,rp)
              rmsq = rmsq/real(nz,rp)
              rms  = sqrt(rmsq - mean**2)

              do k = 1,nz
                 do m = 1,3
                    rnd_tmp(m,j,k) = (rnd_tmp(m,j,k) - mean(m))/rms(m)
                 enddo
              enddo
           enddo

           do k = 1,DF_N
              do j = 1-DF_N,ny+DF_N
                 do m = 1,3
                    rnd_tmp(m,j,nz+k) = rnd_tmp(m,j,k)
                    rnd_tmp(m,j,1 -k) = rnd_tmp(m,j,nz+1-k)
                 enddo
              enddo
           enddo

           fy = 0.0_rp
           do k = 1-DF_N,nz+DF_N
              do j = 1,ny
                 do jj = -DF_N,DF_N
                    do m = 1,3
                       fy(m,j,k) = fy(m,j,k) + DF_By(m,jj,j)*rnd_tmp(m,j+jj,k)
                    enddo
                 enddo
              enddo
           enddo

           vf = 0.0_rp
           do k = 1,nz
              do j = 1,ny
                 do kk = -DF_N,DF_N
                    do m = 1,3
                       vf(m,j,k) = vf(m,j,k) + DF_Bz(m,kk,j)*fy(m,j,k+kk)
                    enddo
                 enddo
              enddo
           enddo

           uf = 0.0_rp
           do k = 1,nz
              do j = 1,ny
                 do m = 1,3
                    do l = 1,3
                       uf(m,j,k) = uf(m,j,k) + DF_LundMatrix(m,l,j)*vf(l,j,k)
                    enddo
                 enddo
              enddo
           enddo

           dx      = gbl_min_step(1)
           exppar1(1) = -0.5_rp*pi*dx/DF_xlen1
           exppar1(2) = -0.5_rp*pi*dx/DF_xlen2
           exppar1(3) = -0.5_rp*pi*dx/DF_xlen3
           exppar2(1) = -pi*dx/DF_xlen1
           exppar2(2) = -pi*dx/DF_xlen2
           exppar2(3) = -pi*dx/DF_xlen3
           par1    = exp(exppar1)
           par2    = sqrt(1.0_rp-exp(exppar2))

           do m = 1,3
              uf(m,:,:) = uf_old(m,:,:)*par1(m) + uf(m,:,:)*par2(m)
           enddo
           uf_old = uf

           do j = 1,GN
              uf(:,ny + j,:) = uf(:,ny - (j-1),:)
           enddo
           do j = 1,GN
              uf(:,1 - j,:) = - uf(:,1 + (j-1),:)
           enddo
           do k = 1,GN
              uf(:,:,nz+k) = uf(:,:,k)
              uf(:,:,1-k)  = uf(:,:,nz+1-k)
           enddo

           do    k = lbz,ubz
              do j = lby,uby

                 rm = phi(i,j,k,1)
                 um = U(i,j,k)
                 vm = V(i,j,k)
                 wm = W(i,j,k)
                 Tm = 1.0_rp/rm

                 us_p = u_inf*uf(1,j,k)
                 vn_p = u_inf*uf(2,j,k)
                 wp   = u_inf*uf(3,j,k)

                 us_m = nx_xi(i,j)*um + ny_xi(i,j)*vm
                 T_   = tbl_noise_temperature(Tm, us_m, us_p, j)

                 up = nx_xi(i,j)*us_p - ny_xi(i,j)*vn_p
                 vp = ny_xi(i,j)*us_p + nx_xi(i,j)*vn_p

                 phi(i,j,k,1) = 1.0_rp/T_
                 U(i,j,k)     = um + up
                 V(i,j,k)     = vm + vp
                 W(i,j,k)     = wm + wp

              enddo
           enddo

        enddo

        if (my_neighbour(2) /= MPI_PROC_NULL) then
           call MPI_SEND(uf_old, plane_count, MPI_RP, my_neighbour(2), &
                         700, mpi_comm_cart, mpi_err)
        endif

        call DeallocateReal(rnd_tmp)
        call DeallocateReal(fy)
        call DeallocateReal(vf)
        call DeallocateReal(uf)
        call DeallocateReal(uf_old)

        if(rank == root) write(*, '(A)') 'done!'

        status = 0

        return
end subroutine init_TBLNoise_fast

end module ic_module
