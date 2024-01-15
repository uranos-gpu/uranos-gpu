module inflow_module
use parameters_module
use storage_module
use profiling_module

implicit none
private

real(rp), allocatable, dimension(:,:,:), public :: iflow_mean               !< mean      inflow field
real(rp), allocatable, dimension(:,:,:), public :: iflow_turb               !< turbulent inflow field
real(rp), allocatable, dimension(:,:,:), public :: iflow_vf_old             !< old velocity fluctuations
real(rp), allocatable, dimension(:,:,:), public :: iflow_vf_new             !< new velocity fluctuations
logical                                , public :: iflow_TurbFlag = .false.
logical                                , public :: iflow_smthFlag = .false.

public init_inflow_profile, init_turbulent_inflow, turbulent_inflow


contains

subroutine init_inflow_profile

        implicit none
        integer :: err = 0
        real(rp) :: stepH 

        !
        ! === set the Inflow flags
        !
        if(smooth_inflow) iflow_SmthFlag = .true.
        if(turb_inflow  ) iflow_turbflag = .true.
        !
        ! === allocate the mean and turb inflow profiles
        !
        allocate(iflow_mean(-GN:0  , lby:uby,1:5), &
                 iflow_turb(1:3,sy:ey, sz:ez), stat=err)
        if(err.ne.0) stop ' Allocation error in init_inflow_profile!' 

        iflow_mean = 0.0_rp
        iflow_turb = 0.0_rp

        selectcase(inflow_profile)
          case('constant')
          call inflow_constant()

          case('parabolic')
          call inflow_parabolic()

          case('blasius')
          call inflow_blasius()

          case('pohlhausen')
          call inflow_pohlhausen()

          case('TurbulentBoundaryLayerInflow')
          call inflow_TurbulentBoundaryLayer()

          case('shock_inflow')
          call shock_inflow()

          case('smooth_body_inflow')
          stepH = 0.0_rp
          call inflow_step(stepH)

          case('supersonic_ramp')
          if(     viscous) call inflow_step(0.0_rp)
          if(.not.viscous) call inflow_constant()

          case default

          if(rank == root) print*, 'Inflow profile ', trim(inflow_profile), ' is not implemented.'
          call secure_stop

        end select

        if(iflow_turbFlag) call init_turbulent_inflow

        return
end subroutine init_inflow_profile




subroutine inflow_constant()

        implicit none
        real(rp) :: Ttot, Trat, T_

        Ttot = total_temperature_inlet
        Trat = (1.0_rp + 0.5_rp*(gamma0-1.0_rp)*Mach**2)
        T_   = Ttot/Trat

        iflow_mean(:,:,1) = 1.0_rp
        iflow_mean(:,:,2) = u_inf
        iflow_mean(:,:,3) = 0.0_rp
        iflow_mean(:,:,4) = 0.0_rp
        iflow_mean(:,:,5) = T_

        return
end subroutine inflow_constant

subroutine inflow_parabolic()
        implicit none

        ! local declaration
        real(rp) :: u_, T_, a, b, c
        integer  :: i,j

        ! parameter
        a = -u_inf/(ymin**2)
        b = 0.0_rp
        c = u_inf

        do    j = lby,uby
           do i = -GN,0
              
              u_ = a*y(j)**2 + b*y(j) + c
              T_ = 1.0_rp - 0.5_rp* (gamma0-1.0_rp)/gamma0 * Prandtl * u_**2

              iflow_mean(i,j,1) = 1.0_rp
              iflow_mean(i,j,2) = u_
              iflow_mean(i,j,3) = 0.0_rp
              iflow_mean(i,j,4) = 0.0_rp
              iflow_mean(i,j,5) = T_

           enddo
        enddo

        return
end subroutine inflow_parabolic


subroutine inflow_blasius()
        
        use fluid_functions_module, only: BlasiusProfile
        implicit none
        ! local declaration
        real(rp), dimension(1-GN:ny+GN) :: VelY
        real(rp), dimension(1-GN:ny+GN) :: RhoY
        real(rp), dimension(1-GN:ny+GN) :: TmpY
        integer                         :: i,j

        ! === get the nominal uncompressible inflow profile
        call BlasiusProfile(VelY,TmpY,RhoY)
        
        do j    = lby,uby
           do i = -GN,0

              iflow_mean(i,j,1) = RhoY(j)
              iflow_mean(i,j,2) = VelY(j)
              iflow_mean(i,j,3) = 0.0_rp
              iflow_mean(i,j,4) = 0.0_rp
              iflow_mean(i,j,5) = TmpY(j)

           enddo
        enddo
       
        return
end subroutine inflow_blasius


subroutine inflow_pohlhausen()
        
        use fluid_functions_module, only: pohlhausenProfile, CompressibleCorrection

        implicit none

        ! local declaration
        real(rp), dimension(1-GN:ny+GN) :: u_inc
        real(rp), dimension(1-GN:ny+GN) :: u_cmp, u_cmp_old
        real(rp), dimension(1-GN:ny+GN) :: RhoY
        real(rp), dimension(1-GN:ny+GN) :: TmpY
        real(rp)                        :: Tw, Rw,RecFac
        real(rp)                        :: alf, uu, du, fuu, Tr, uci, err
        real(rp), parameter             :: s      = 1.1_rp 
        integer , parameter             :: itmx   = 100
        real(rp), parameter             :: toll   = 1.0E-14_rp
        integer                         :: i,j, itr


        call PohlhausenProfile(u_inc)
        !call CompressibleCorrection(ny,.false.,Mach,Prandtl,Trat,Tw,Rw,VelY_inc,VelY_cmp,RhoY,TmpY)

        RecFac = sqrt(Prandtl)
        Tr  = 1.0_rp + 0.5_rp*RecFac*(gamma0-1.0_rp)*Mach**2
        Tw  = TRat*Tr
        Rw  = 1.0_rp/Tw
        alf = s*Prandtl
        
        RhoY = 1.0_rp
        tmpY = 1.0_rp

        itr   = 0
        err   = 2*toll
        u_cmp = u_inc
        do while(err > toll .and. itr < itmx)

           u_cmp_old = u_cmp
                
           ! compute Temperature and Density profiles
           do j = 0,ny+GN
              uu      = u_cmp(j)/u_cmp(ny)
              fuu     = alf*uu+(1.0_rp-alf)*uu**2
              TmpY(j) = tw+(tr-tw)*fuu+(1.0_rp-tr)*uu**2           !< Zhang
              RhoY(j) = 1.0_rp/TmpY(j)
           enddo
        
           ! apply the correction to uPlus_i according to VAN DRIEST
           do j = 1,ny+GN
              du       = u_inc(j)-u_inc(j-1)
              uci      = 0.5_rp*(sqrt(Rw/RhoY(j))+sqrt(Rw/RhoY(j-1))) 
              u_cmp(j) = u_cmp(j-1)+uci*du
           enddo

           err = sum(abs(u_cmp - u_cmp_old))
           itr = itr + 1

        enddo

        do j    = lby,uby
           do i = -GN,0

              iflow_mean(i,j,1) = RhoY(j)
              iflow_mean(i,j,2) = u_inc(j)
              iflow_mean(i,j,3) = 0.0_rp
              iflow_mean(i,j,4) = 0.0_rp
              iflow_mean(i,j,5) = TmpY(j)

           enddo
        enddo
       
        return
end subroutine inflow_pohlhausen



subroutine inflow_step(stepH)

        use fluid_functions_module, only: CompressibleCorrection

        implicit none
        real(rp) , intent(in)    :: stepH

        real(rp), parameter                 :: pi_wake = 0.434_rp
        real(rp), parameter                 :: toll    = 1.0E-14_rp
        real(rp), dimension(:), allocatable :: y_tmp, uPlus
        real(rp), dimension(1-GN:ny+GN)     :: VelY_inc, VelY_cmp, RhoY, TmpY
        real(rp), dimension(1-GN:ny+GN)     :: VelY_cmp_old
        real(rp) :: yj, eta,yp,up, ue, iReTau, Tr, Tw, Rw, alf, err
        real(rp) :: uu, fuu, du, uci
        integer  :: j, jw, itr, itmx = 100

        allocate(y_tmp(1-GN:ny+GN))
        allocate(UPlus(1-GN:ny+GN))
        call compute_grid_point(y_tmp,ny,lbound(y_tmp,1),ubound(y_tmp,1),&
                                ymin,ymax,stretching_par,gridpoint_y)
        
        !
        ! === find wall location
        !
        jw = 1-GN
        do j = 1-GN, ny+GN
           if(y_tmp(j) > stepH) then
             jw = j-1
             exit
           endif
        enddo
        !
        ! === get musker profile as a function of ReTau
        !
        iReTau = 1.0_rp/ReTau
        uPlus  = 0.0_rp
        do j = jw, ny+GN
        
             yj  = y_tmp(j) - y_tmp(jw)
             yp  = yj * ReTau

             eta = yp*iReTau
             eta = min(1.0_rp,eta)
             yp  = eta*retau

             up = 5.424_rp*atan((2*yp-8.15_rp)/16.7_rp)&
                     +log10((yp+10.6_rp)**9.6_rp/(yp**2-8.15_rp*yp+86)**2)-3.51132976630723_rp+&
                      2.44_rp*(pi_wake*(6*eta**2-4*eta**3)+(eta**2*(1-eta)))

             UPlus(j) = up

        enddo
        !
        ! === rescale Musker profile to get physical speed
        !
        ue = maxval(UPlus)
        do j = 1-GN,ny+GN
           VelY_inc(j) = uPlus(j)/ue*u_inf
        enddo
        !
        ! === get temperature and density profiles via SRA
        !
        Tr  = 1.0_rp + 0.5_rp*(Prandtl)**(1.0_rp/3.0_rp)*(gamma0-1.0_rp)*Mach**2
        Tw  = TRat*Tr
        Rw  = 1.0_rp/Tw
        alf = 1.1_rp*Prandtl

        itr = 0
        err = 2.0_rp*toll
        VelY_cmp = VelY_inc
        TmpY     = 1.0_rp
        RhoY     = 1.0_rp
        do while(err > toll .and. itr < itmx)
           
           VelY_cmp_old = VelY_cmp
           ! compute Temperature and Density profiles
           do j = jw,ny+GN
                uu      = VelY_cmp(j)/VelY_cmp(ny)
                fuu     = alf*uu+(1.0_rp-alf)*uu**2
                TmpY(j) = tw+(tr-tw)*fuu+(1.0_rp-tr)*uu**2           !< Zhang
                RhoY(j) = 1.0_rp/TmpY(j)
           enddo

           do j = jw+1,ny+GN
              du       = VelY_inc(j)-VelY_inc(j-1)
              uci      = 0.5_rp*(sqrt(Rw/RhoY(j))+sqrt(Rw/RhoY(j-1))) 
              VelY_cmp(j) = VelY_cmp(j-1)+uci*du
           enddo
        
           err = sum(abs(VelY_cmp - VelY_cmp_old))
           itr = itr + 1
        enddo

        !
        ! === assign to inflow
        !
        do j    = lby,uby
           do i = -GN,0

              iflow_mean(i,j,1) = RhoY(j)
              iflow_mean(i,j,2) = VelY_inc(j)
              iflow_mean(i,j,3) = 0.0_rp
              iflow_mean(i,j,4) = 0.0_rp
              iflow_mean(i,j,5) = TmpY(j)

           enddo
        enddo
        
        deallocate(y_tmp,uPlus)
        return
end subroutine inflow_step


subroutine shock_inflow()
! ----------------------------------------------------------------------------------
!       
!       This subroutine provides the exact inflow condition of a inflowing shock wave
!
! ----------------------------------------------------------------------------------
        
        use parameters_module     , only: rp, Mach
        use fluid_functions_module, only: inverse_rankine_hugoniot, rankine_hugoniot
        use math_tools_module     , only: newton_raphson

        implicit none

        ! === local declarations
        ! ambient conditions
        real(rp), parameter :: r1     =  1.0_rp 
        real(rp), parameter :: p1     =  1.0_rp
        real(rp), parameter :: T1     =  p1/r1
        ! other parameters
        real(rp), parameter :: gamma0 =  1.4_rp        
        real(rp), parameter :: toll   =  1.0E-14_rp
        integer , parameter :: itmax  =  10000
        ! post-shocked flow variables
        real(rp)            :: r2, p2, T2, u2, v2, w2, M2, M1

        ! compute post-shocked field based on the post shock Mach number
        call newton_raphson(inverse_rankine_hugoniot,itmax,toll, Mach+1.0_rp, M1)
        call Rankine_Hugoniot(r1,p1,T1,M1,r2,p2,T2,M2)

        u2 = Mach*sqrt(gamma0*T2)
        v2 = 0.0_rp
        w2 = 0.0_rp

        ! set variable inflow profile
        iflow_mean(:,:,1) = r2
        iflow_mean(:,:,2) = u2
        iflow_mean(:,:,3) = v2
        iflow_mean(:,:,4) = w2
        iflow_mean(:,:,5) = p2/r2

        if(rank == root) then
          write(*,'(A,e18.6)') ' The flow  is inflowing  at Mach: ', Mach
          write(*,'(A,e18.6)') ' The shock is travelling at Mach: ', M1
          write(*,'(A)')       ' ------------------------------------------------------ '
        endif

        return
end subroutine shock_inflow




subroutine inflow_TurbulentBoundaryLayer()
        
        use fluid_functions_module, only: MuskerProfile, CompressibleCorrection, laminar_viscosity, BoundaryLayerQuantities
        use real_to_integer_module, only: locate
        use interpolation_module  , only: polint
        use parameters_module     , only: ReTau
        use FileModule          

        implicit none

        ! local declarations
        real(rp), dimension(:,:,:), allocatable :: tempPrim
        real(rp), dimension(:)    , allocatable :: delta_array
        real(rp), dimension(:)    , allocatable :: deltaVArray
        real(rp), dimension(:)    , allocatable :: frict_array
        real(rp), dimension(:)    , allocatable :: theta_array
        real(rp), dimension(:)    , allocatable :: y_tmp
        real(rp), dimension(0:ny)               :: UPlusInc, UplusCmp, RhoY, TmpY, velY
        real(rp), parameter                     :: toll = 1.0E-14_rp
        real(rp)                                :: tWall, rWall
        real(rp)                                :: delta, deltaV, thrat, deltaStar, theta, cf
        real(rp)                                :: dx, dy, ReTau_new, ReTau_old, err
        real(rp)                                :: ue, ReTauI, yl, errR, errU
        real(rp)                                :: ri, ui, vi_j00, vi_jm1, vi_j
        integer                                 :: i,j,jj,jjj, m, jg, ji
        integer                                 :: itr, itrmax = 30, ierr = 0

        real(rp) :: delta_REF, theta_REF, deltaVREF, frict_REF
        
        character(7)   :: crank
        type(FileType) :: inflow

        !
        ! === allocate temporary variables
        !
        allocate(delta_array(-GN:0), &
                 deltaVArray(-GN:0), &
                 frict_Array(-GN:0), &
                 theta_array(-GN:0), stat = ierr)
        if(ierr.ne.0) stop ' Allocation error'
        
        allocate(tempPrim(-GN:0, 1-GN:ny+GN, 3), y_tmp(0:ny), stat = ierr)
        if(ierr.ne.0) stop ' Allocation error'
        
        call compute_grid_point(y_tmp,ny,0,ny,ymin,ymax,stretching_par,gridpoint_y)
        y_tmp(0) = 0.0_rp
        
        !
        ! === getting the REFERENCE inflow profile
        !
        call MuskerProfile(ny,y_tmp,ReTau,UplusInc)
        call CompressibleCorrection(ny,.true.,Mach,Prandtl,Trat,tWall,rWall,UPlusInc,UPlusCmp,RhoY,TmpY)

        !
        ! === computing the growth of the boundary layer
        !
        call BoundaryLayerQuantities(y_tmp,RhoY,UPlusCmp,rWall,thrat,deltaStar,cf)
        delta_REF = 1.0_rp
        deltaVREF = delta_REF/ReTau
        frict_REF = cf
        theta_REF = thrat*delta_REF

        ReTau_old = ReTau - 1.0_rp
        ReTau_new = ReTau
        dx        = Lx/real(nx,rp)

        ! get the node-0 values integrating from REF in half cell space
        itr = 0
        err = huge(1.0_rp)
        do while(err > toll .and. itr < itrmax)

           call MuskerProfile(ny,y_tmp,ReTau_new,UplusInc)
           call CompressibleCorrection(ny,.true.,Mach,Prandtl,Trat,tWall,rWall,UPlusInc,UPlusCmp,RhoY,TmpY)
           call BoundaryLayerQuantities(y_tmp,RhoY,UPlusCmp,rWall,thrat,deltaStar,cf)

           theta     = theta_REF - 0.25_rp*0.5_rp*dx*(cf+frict_REF)
           delta     = theta/thrat
           deltaV    = deltaVREF*sqrt(frict_REF/cf)
           reTau_new = delta/deltaV

           err = abs(ReTau_new - ReTau_old) 

           ReTau_old = ReTau_new
           itr = itr + 1

        enddo
        theta_array(0) = theta
        frict_array(0) = cf
        delta_array(0) = delta
        deltaVArray(0) = deltaV

        !
        ! integrating from -1 to -GN
        !
        do i = -1,-GN,-1

           itr = 0
           err = huge(1.0_rp)
           do while(err > toll .and. itr < itrmax)

              call MuskerProfile(ny,y_tmp,ReTau_new,UplusInc)
              call CompressibleCorrection(ny,.true.,Mach,Prandtl,Trat,tWall,rWall,UPlusInc,UPlusCmp,RhoY,TmpY)
              call BoundaryLayerQuantities(y_tmp,RhoY,UPlusCmp,rWall,thrat,deltaStar,cf)

              theta     = theta_array(i+1) - 0.25_rp*dx*(cf+frict_array(i+1))
              delta     = theta/thrat
              deltaV    = deltaVArray(i+1)*sqrt(frict_array(i+1)/cf)
              reTau_new = delta/deltaV

              err = abs(ReTau_new - ReTau_old) 

              ReTau_old = ReTau_new
              itr = itr + 1

           enddo
           theta_array(i) = theta
           frict_array(i) = cf
           delta_array(i) = delta
           deltaVArray(i) = deltaV

        enddo

        !
        ! === interpolation on the mesh
        !
        do i = -GN,0
           delta  = delta_array(i)
           deltaV = deltaVArray(i)
           ReTauI = delta/deltaV

           call MuskerProfile(ny,y_tmp,ReTauI,UplusInc)
           call CompressibleCorrection(ny,.true.,Mach,Prandtl,Trat,tWall,rWall,UPlusInc,UPlusCmp,RhoY,TmpY)

           ue = maxval(UPlusCmp)
           VelY = UPlusCmp/ue

           do j = 0, ny

              yl = y_tmp(j)/delta
              call locate(y_tmp,0,ny,yl,jj)

              m = 4
              jjj = min(max(jj-(m-1)/2,1),ny+1-m)
              call polint(y_tmp(jjj),rhoY(jjj),m,yl,ri,errR)
              call polint(y_tmp(jjj),velY(jjj),m,yl,ui,errU)

              tempPrim(i,j,1) = ri
              tempPrim(i,j,2) = ri*ui
              tempPrim(i,j,3) = 0.0_rp

           enddo
        enddo

        ! === integrate rhov
        do j    = 1   , ny
           do i = 1-GN, 0

              dy = y_tmp(j) - y_tmp(j-1)

              vi_j00 = (tempPrim(i,j  ,2) - tempPrim(i-1,j  ,2))/dx
              vi_jm1 = (tempPrim(i,j-1,2) - tempPrim(i-1,j-1,2))/dx
              vi_j   = 0.5_rp*(vi_j00 + vi_jm1)

              tempPrim(i,j,3) = tempPrim(i,j-1,3) - vi_j * dy

           enddo
        enddo

        ! extapolation on north face ghosts
        do    j = 1,GN
           do i = -GN,0

              jg = ny + j      !< ghost node
              ji = ny - (j-1)  !< inner node

              tempPrim(i,jg,1) =  tempPrim(i,ji,1)
              tempPrim(i,jg,2) =  tempPrim(i,ji,2)
              tempPrim(i,jg,3) =  tempPrim(i,ji,3)

           enddo
        enddo

        ! extapolation on south face ghosts
        do    j = 1,GN
           do i = -GN,0

              jg = 1 - j      !< ghost node
              ji = 1 + (j-1)  !< inner node

              tempPrim(i,jg,1) =  tempPrim(i,ji,1)
              tempPrim(i,jg,2) = -tempPrim(i,ji,2)
              tempPrim(i,jg,3) = -tempPrim(i,ji,3)

           enddo
        enddo


        do j = lby,uby
           do i = 1-GN,0

              iflow_mean(i,j,1) = tempPrim(i,j,1)
              iflow_mean(i,j,2) = u_inf*tempPrim(i,j,2)/tempPrim(i,j,1)
              iflow_mean(i,j,3) = u_inf*tempPrim(i,j,3)/tempPrim(i,j,1)
              iflow_mean(i,j,4) = 0.0_rp
              iflow_mean(i,j,5) = 1.0_rp/tempPrim(i,j,1)

           enddo
        enddo

        write(crank, '(i7.7)') rank
        inflow%name = 'InflowData'//trim(crank)
        inflow%dir  = trim(data_dir)//'/INITIAL_CONDITION'
        call OpenNewFile(inflow,it)
        do j = lby,uby
           write(inflow%unit,*) y(j), iflow_mean(0,j,1), iflow_mean(0,j,2)/u_inf, &
                                      iflow_mean(0,j,3)/u_inf, iflow_mean(0,j,5)
        enddo
        call CloseFile(inflow)


        deallocate(delta_array,deltaVArray,frict_Array,theta_array)
        deallocate(y_tmp,tempPrim)


        return
end subroutine inflow_TurbulentBoundaryLayer


subroutine init_turbulent_inflow

        use df_module
        use fileModule

        implicit none
        real(rp), dimension(:), allocatable :: y_tmp
        integer                             :: err = 0
#ifdef DEBUG
        type(FileType)         :: reyStressFile
        real(rp), dimension(3) :: rms, rmean, rmeansq
        real(rp)               :: rey12, rey13, rey23
        integer                :: m
#endif
        
        ! allocations
        allocate(y_tmp(0:ny), stat = err)
        if(err .ne. 0) stop ' Allocation error in init_turbulent_inflow'

        allocate(DF_ylen(3,0:ny), DF_zlen(3,0:ny), stat = err)
        if(err .ne. 0) stop ' Allocation error in init_turbulent_inflow'

        allocate(DF_By(3,-DF_N:DF_N, 0:ny), DF_Bz(3,-DF_N:DF_N, 0:ny), stat = err)
        if(err .ne. 0) stop ' Allocation error in init_turbulent_inflow'

        allocate(DF_Rnd2D(3, 1-DF_N:ny+DF_N, 1-DF_N:nz+DF_N), stat = err)
#ifdef AMD
        rnd2Dnumshape = shape(DF_Rnd2D)
        rnd2Dnumsize = rnd2Dnumshape(1) * rnd2Dnumshape(2) * rnd2Dnumshape(3)
        allocate(rnd2Dnum(rnd2Dnumsize), stat = err)
#endif
        if(err .ne. 0) stop ' Allocation error in init_turbulent_inflow'

        allocate(iflow_vf_old(3, sy:ey, sz:ez), stat = err)
        if(err .ne. 0) stop ' Allocation error in init_turbulent_inflow'

        allocate(iflow_vf_new(3, sy:ey, sz:ez), stat = err)
        if(err .ne. 0) stop ' Allocation error in init_turbulent_inflow'

        ! step 1) create a temporary grid 
        call compute_grid_point(y_tmp,ny,0,ny,ymin,ymax,stretching_par,gridpoint_y)
        y_tmp(0) = 0.0_rp
        
        selectcase(inflow_profile)

          case('TurbulentBoundaryLayerInflow', 'smooth_body_inflow', 'supersonic_ramp')
          ! step 2) Read input data from database
          call DFInputData_TBL(ReTau,Mach,Prandtl,y_tmp,ny)

          ! step 3) Get integral lenght scale
          call DFIntegralLenght_TBL(y_tmp,ny)
        
          case('constant')
          ! step 2) Read input data from database
          call DFInputData_HTURB(y_tmp,ny)

          ! step 3) Get integral lenght scale
          call DFIntegralLenght_HTURB(Lz)

        case default
          print*, 'Turbulent inflow is not implemented for inflow "', trim(inflow_profile), '"'
          stop
        endselect


        ! step 4) compute DF coefficients
        call DFCoefficients(y_tmp,ny,gbl_min_step)
#ifdef AMD
        !$acc data copy(rnd2Dnum, rnd2Dnumshape,rnd2Dnumptr) &
        !$acc copy(DF_Rnd2D,DF_N,DF_ylen,DF_zlen,DF_By,DF_Bz,DF_Fy,DF_LundMatrix) &
        !$acc copy(iflow_vf_old, iflow_vf_new, iflow_turb)
#endif
#ifdef NVIDIA
        !$acc data copy(DF_Rnd2D,DF_ylen,DF_zlen,DF_By,DF_Bz,DF_Fy,DF_LundMatrix) &
        !$acc copy(iflow_vf_old, iflow_vf_new, iflow_turb)
#endif


        ! step 5) get a random slice
        call DFRandomField2D(ny,nz)

        ! step 6) compute velocity flunctuations based on the extracted random slice
        call DFConvolution2D(sy,ey,sz,ez,iflow_vf_old)

        ! step 7) get a new random slice
        call DFRandomField2D(ny,nz)
        
        ! step 8) compute a new velocity fluctuation based on the extracted random slice
        call DFConvolution2D(sy,ey,sz,ez,iflow_vf_new)

        ! step 9) provide time-correlation via Castro's formula
        call DFCastroTimeCorrelation(ik,c_rk,u_inf,dt,sy,ey,sz,ez,iflow_vf_old,iflow_turb)

        ! step 6) Enforce Reynolds stress and obtain turb inflow profile
        call DFEnforceReynoldsStresses2D(sy,ey,sz,ez,iflow_vf_new,iflow_turb)

        !$acc end data

#ifdef DEBUG
        if(nprocs == 1) then
        reyStressFile%name = 'rey_computed_inflow' 
        reyStressFile%dir  = trim(data_dir)//'/INFLOW_REYNOLDS_STRESS'
        call OpenNewFile(ReyStressFile,it)
        do j=1,ny
         rmean   = 0._rp
         rmeansq = 0._rp
         rey12   = 0._rp
         rey13   = 0._rp
         rey23   = 0._rp
         do k=1,nz
           do m=1,3
            rmean(m)   = rmean(m)   + iflow_turb(m,j,k)
            rmeansq(m) = rmeansq(m) + iflow_turb(m,j,k)**2
           enddo
           rey12 = rey12+iflow_turb(1,j,k)*iflow_turb(2,j,k)
           rey13 = rey13+iflow_turb(1,j,k)*iflow_turb(3,j,k)
           rey23 = rey23+iflow_turb(2,j,k)*iflow_turb(3,j,k)
         enddo
         rmean   = rmean/(nz)
         rmeansq = rmeansq/(nz)
         rms     = sqrt(rmeansq-rmean*rmean)
         rey12   = rey12/(nz)
         rey13   = rey13/(nz)
         rey23   = rey23/(nz)
         rey12   = rey12-rmean(1)*rmean(2)
         rey13   = rey13-rmean(1)*rmean(3)
         rey23   = rey23-rmean(2)*rmean(3)
         write(reyStressFile%unit,*) &
         y(j),(rms(m)**2,m=1,3),rey12,rey13,rey23,(rmean(m),m=1,3)
        enddo
        call CloseFile(ReyStressFile)
        endif
#endif
        deallocate(y_tmp)

        return
end subroutine init_turbulent_inflow





subroutine turbulent_inflow()
        
        use df_module
        use FileModule

        implicit none
        integer :: j,k,l

#ifdef DEBUG
        type(fileType)         :: turbInflow
        real(rp)               :: lcl_ekt, gbl_ekt
        integer                :: lcl_pts, gbl_pts
        real(rp)               :: lcl_ekt_max, gbl_ekt_max
        real(rp)               :: ek
        integer                :: err = 0
#endif

        call StartProfRange('turbulent_inflow')

        if(iflow_TurbFlag) then

          ! step 1) extract a new random slice
          call DFRandomField2D(ny,nz)

          ! step 2) get velocity flunctuations based on the new random slice
          call DFConvolution2D(sy,ey,sz,ez,iflow_vf_new)

          ! step 3) provide time correlation of the inflow slices
          call DFCastroTimeCorrelation(ik,c_rk,u_inf,dt,sy,ey,sz,ez,iflow_vf_old,iflow_vf_new)

          ! step 4) enforce Reynolds Stress
          call DFEnforceReynoldsStresses2D(sy,ey,sz,ez,iflow_vf_new,iflow_turb)

        else

          !$acc parallel default(present)
          !$acc loop gang, vector collapse(3)
          do       k = sz,ez
             do    j = sy,ey
                do l = 1,3
                   iflow_turb(l,j,k) = 0.0_rp
                enddo
             enddo
          enddo
          !$acc end parallel

        endif

#ifdef DEBUG
        turbInflow%name = 'inflowFluctuation'//trim(str(rank))
        turbInflow%dir  = trim(data_dir)//'/TURB_INFLOW/'
        call OpenNewFile(turbInflow,it)

        lcl_ekt     = 0.0_rp
        lcl_pts     = 0
        lcl_ekt_max = 0.0_rp
        do k    = sz,ez
           do j = sy,ey 
              ek = (u_inf*iflow_turb(1,j,k))**2 + (u_inf*iflow_turb(2,j,k))**2 + (u_inf*iflow_turb(3,j,k))**2

              lcl_ekt     = lcl_ekt + ek
              lcl_ekt_max = max(ek,lcl_ekt_max)
              lcl_pts     = lcl_pts + 1
           enddo
        enddo
        if(mpi_flag) then
          call MPI_allreduce(lcl_ekt    ,gbl_ekt    ,1,MPI_RP, MPI_SUM, mpi_comm_cart, err)
          call MPI_allreduce(lcl_pts    ,gbl_pts    ,1,MPI_INTEGER         , MPI_SUM, mpi_comm_cart, err)
          call MPI_allreduce(lcl_ekt_max,gbl_ekt_max,1,MPI_RP, MPI_MAX, mpi_comm_cart, err)
        else
          gbl_ekt     = lcl_ekt
          gbl_pts     = lcl_pts
          gbl_ekt_max = lcl_ekt_max
        endif

        write(turbInflow%unit,*) '# turbulent mach nmber ', sqrt(gbl_ekt/real(gbl_pts,rp))
        write(turbInflow%unit,*) '# max turbulent energy ', gbl_ekt_max
        do k    = sz,ez
           do j = sy,ey
              write(turbInflow%unit,*) y(j), z(k), u_inf*iflow_turb(:,j,k)
           enddo
        enddo
        call CloseFile(turbInflow)
#endif

        call EndProfRange

        return
end subroutine turbulent_inflow






end module inflow_module
