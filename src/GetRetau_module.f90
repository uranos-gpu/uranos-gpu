module GetRetau_module
use parameters_module
use mpi_module
implicit none

private
public init_reynolds


contains
subroutine init_Reynolds


        implicit none

        selectcase(ic)
        case('turbulent_channel')
          call GetReynolds_TCH

        case('turbulent_BL')
          call GetReynolds_TBL

        case('smooth_body', 'supersonic_ramp')
          call GetReynolds_SmoothBody

        case default
          ReTau = 0.0_rp

        endselect

        return
end subroutine init_Reynolds



subroutine GetReynolds_TCH

        use fluid_functions_module, only: ReichardtLaw
        use mesh_module           , only: compute_grid_point

        implicit none

        real(rp), dimension(:), allocatable :: y_tmp
        real(rp), dimension(:), allocatable :: uPlus, ucPlus, ucPlusOld
        real(rp), dimension(:), allocatable :: RhoY, TmpY
        real(rp), parameter                 :: toll = 1.0E-14_rp
        integer , parameter                 :: itmx = 100
        real(rp)                            :: yp, t_wl, r_wl,  refc, trec, s, alf
        real(rp)                            :: uu, fuu, du, uci, te
        real(rp)                            :: rerr, ReOut
        real(rp)                            :: rhosum,  rhousum, dy, ubulk
        integer                             :: n, j, iter, err = 0

        ReTau = Reynolds

        n = 1000000

        allocate(y_tmp(0:n),uPlus(0:n),ucPlus(0:n), &
                 ucPlusOld(0:n),RhoY(0:n),TmpY(0:n), stat=err)
        if(err .ne. 0) stop ' Allocation error get Reynolds Channel'

        ! build a temporary and very fine grid
        do j = 0,n
           y_tmp(j) = j*2/real(2*n,rp)
        enddo
       
        ! init uPlus
        uPlus = 0.0_rp
        do j = 1,n
           yP = (y_tmp(j))*ReTau
           uPlus(j) = ReichardtLaw(yP)
        enddo
        ucPlus = uPlus
        
        ! get wall properties
        t_Wl = 1.0_rp
        r_wl = 1.0_rp
        refc = Prandtl**(1.0_rp/3.0_rp)
        Trec = 1.0_rp + refc*0.5_rp*(gamma0-1.0_rp)*Mach**2
        s    = 1.1_rp
        te   = Trec
        alf  = s*Prandtl

        iter = 0
        rerr = 2*toll
        do while(rerr > toll .and. iter < itmx)
        
           ucPlusOld = ucPlus
           TmpY(0)   = t_wl
           RhoY(0)   = r_wl
           do j = 1,n
              uu      = ucplus(j)/ucplus(n)
              fuu     = alf*uu+(1.0_rp-alf)*uu**2
              TmpY(j) = t_wl+(trec-t_wl)*fuu+(te-trec)*uu**2
              RhoY(j) = 1.0_rp/TmpY(j)
           enddo

           do j=1,n
            du        = uplus(j)-uplus(j-1)
            uci       = 0.5_rp*(sqrt(r_wl/RhoY(j))+sqrt(r_wl/RhoY(j-1)))
            ucplus(j) = ucplus(j-1)+uci*du
           enddo

           rerr  = maxval(abs(ucPlus - ucplusOld))
           iter = iter + 1

        enddo

        rhosum   = 0.0_rp
        rhousum  = 0.0_rp
        do j=1,n
         dy      = y_tmp(j)-y_tmp(j-1)
         rhosum  = rhosum+RhoY(j)*dy
         rhousum = rhousum+RhoY(j)*ucplus(j)*dy
        enddo
        ubulk =  rhousum/rhosum

        ReOut = ReTau*ubulk*rhosum

        Reynolds = ReOut
        mu_inf   = sqrt(gamma0)*Mach/Reynolds

        deallocate(y_tmp,uPlus,ucPlus,ucPlusOld,RhoY,TmpY)

        return
end subroutine GetReynolds_TCH


subroutine GetReynolds_TBL
        
        use mesh_module           , only: compute_grid_point
        use fluid_functions_module, only: MuskerProfile, CompressibleCorrection, &
                                          Sutherland

        implicit none

        ! local declarations
        real(rp), dimension(:), allocatable :: y_tmp
        real(rp), dimension(0:ny) :: UPlusInc, UplusCmp, RhoY, TmpY
        real(rp)                  :: mWall,rWall,tWall,ue, ReOut
        integer                   :: err = 0

        ReTau = Reynolds
        
        allocate(y_tmp(0:ny),stat=err)
        if(err .ne. 0) stop ' Allocation error in GetRetau_TBL'
        call compute_grid_point(y_tmp,ny,0,ny,ymin,ymax,stretching_par,gridpoint_y)
        y_tmp(0) = 0.0_rp

        ! init incompressible profile
        call MuskerProfile(ny,y_tmp,ReTau,UplusInc)
        uPlusCmp = uPlusInc
        
        ! iterate to get the compressible correction
        call CompressibleCorrection(ny,.true.,Mach,Prandtl,Trat,tWall,rWall,&
                                     UPlusInc,UPlusCmp,RhoY,TmpY)
        
        ! obtain reference Reynolds number
        ue       = maxval(UPlusCmp)
        mWall    = Sutherland(tWall)

        ReOut    = ReTau*ue/rWall*mWall
        Reynolds = ReOut
        mu_inf   = sqrt(gamma0)*Mach/Reynolds

        deallocate(y_tmp)

        return
end subroutine GetReynolds_TBL



subroutine GetReynolds_SmoothBody

        use mesh_module, only: compute_grid_point
        use fluid_functions_module, only: MuskerProfile, CompressibleCorrection, &
                                          sutherland

        implicit none
        integer, parameter :: n = 10000
        real(rp), dimension(:), allocatable :: y_tmp
        real(rp), dimension(0:n)            :: uPlusInc, uPlusCmp, RhoY, TmpY

        real(rp) :: xi, a, mWall, tWall, rWall, ue, ReOut
        integer  :: j, err = 0

        ReTau = Reynolds

        allocate(y_tmp(0:n),stat=err)
        if(err .ne. 0) stop ' Allocation error in GetRetau_SmoothBody'
        
        a = 10.0_rp
        do j = 0,n
           xi = (j)/real(n,rp)
           y_tmp(j) = (ymax-ymin)*(1.0_rp+tanh((xi-1.0_rp)*a)/tanh(a))
        enddo

        call MuskerProfile(n,y_tmp,ReTau,UplusInc)
        UPlusCmp = uPlusInc

        call CompressibleCorrection(n,.true.,Mach,Prandtl,Trat,tWall,rWall,&
                                     UPlusInc,UPlusCmp,RhoY,TmpY)

        ue = maxval(UPlusCmp)
        mWall    = Sutherland(tWall)

        ReOut    = ReTau*ue/rWall*mWall!*(1.0_rp/0.032_rp)
        Reynolds = ReOut
        mu_inf   = sqrt(gamma0)*Mach/Reynolds

        deallocate(y_tmp)

        return
end subroutine GetReynolds_SmoothBody
















end module GetRetau_module
