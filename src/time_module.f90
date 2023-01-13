module time_module
use parameters_module
use storage_module
use mpi_module
use shock_detection_module
use profiling_module

implicit none
private
public init_runge_kutta, runge_kutta, compute_dt, last_iteration, explicit_filter_wmles

contains

subroutine init_runge_kutta
! ------------------------------------------------------------------
!
!       This subroutine init the solution for Runge-Kutta substeps
!       
! ------------------------------------------------------------------
        implicit none
        integer :: i, j, k, l
#ifdef TIME
        call mpi_stime(s_irk_time)
#endif
        call StartProfRange("init_runge_kutta") 

        !$acc parallel default(present)
        !$acc loop gang, vector collapse(4)
        do          l = 1,eqs
           do       k = lbz,ubz
              do    j = lby,uby
                 do i = lbx,ubx
                    phi_n(i,j,k,l) = phi(i,j,k,l)
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        ! init ducros sensor if the scheme is hybrid
        if(hybrid_weno) call compute_hybrid_weno_flag

        call EndProfRange
#ifdef TIME
        call mpi_etime(s_irk_time,t_irk_calls,t_irk_time)
#endif
        
return
end subroutine init_runge_kutta


subroutine runge_kutta
! ----------------------------------------------------------------------------------------
!
!       This subroutine evolves conservative variable via a a Runge-Kutta scheme. 
!
!       REF: SIGAL GOTTLIEB AND CHI-WANG SHU, 
!            "TOTAL VARIATION DIMINISHING RUNGE-KUTTA SCHEMES"
!
! ----------------------------------------------------------------------------------------
        implicit none
        real(rp) :: aik, bik
        integer  :: i,j,k,n,l

#ifdef TIME
        call mpi_stime(s_crk_time)
#endif
        call StartProfRange("runge_kutta")

        ! Runge - Kutta time
        time = time + Dt * c_rk(ik)

        ! Runge - Kutta coefficient
        aik = a_rk(ik)
        bik = b_rk(ik)

        ! Runge - Kutta scheme
        !$acc parallel default(present)
        !$acc loop gang, vector collapse(4)
        do          l = 1,eqs
           do       k = sz,ez
              do    j = sy,ey
                 do i = sx,ex
                    phi(i,j,k,l) = aik * phi_n(i,j,k,l) &
                                 + bik * (phi(i,j,k,l) + dt * RHS(i,j,k,l))
                 enddo
              enddo
           enddo
        enddo
        !$acc end parallel

        call check_for_NAN_values

        call EndProfRange
#ifdef TIME
        call mpi_etime(s_crk_time,t_crk_calls,t_crk_time)
#endif
        
return
end subroutine runge_kutta


subroutine compute_dt
! -------------------------------------------------------------------------------
!
!       This subroutine update time and iteration and calculate the time-step dt
!       in order to ensure CFL condition. 
!
!       REFERENCE: 2002, Pirozzoli, "Conservative Hybrid Compact-WENO Schemes 
!                                    for Schock-Turbolence Interaction",
!                                    pag. 105
!
! -------------------------------------------------------------------------------
        implicit none
        real(rp) :: my_dt          !< proc's time step
        real(rp) :: my_dt_CFL      !< proc's dt courant condition
        real(rp) :: my_dt_DIF1     !< proc's dt diffusion
        real(rp) :: my_dt_DIF2     !< proc's dt diffusion
       
        real(rp) :: max_speed1, max_speed2, max_speed3
        real(rp) :: dt_i1, dt_i2, dt_i3
        real(rp) :: max_vis, max_lmd, min_st2
        real(rp) :: ru, rv, rw, ir, T_, c_
        integer  :: i,j,k

#ifdef TIME
        call mpi_stime(s_cfl_time)
#endif
        call StartProfRange("compute_dt") 

        if(time+dt < tmax) then

          max_speed1 = 0.0_rp
          max_speed2 = 0.0_rp
          max_speed3 = 0.0_rp
          !$acc parallel default(present) copy(max_speed1,max_speed2,max_speed3)
          !$acc loop gang, vector collapse(3) &
          !$acc reduction(max:max_speed1,max_speed2,max_speed3)
          do       k = sz,ez
             do    j = sy,ey
                do i = sx,ex
                
                   ru = phi(i,j,k,2)
                   rv = phi(i,j,k,3)
                   rw = phi(i,j,k,4)

                   ir = 1.0_rp/phi(i,j,k,1)
                   T_ = ir*cv_i*(phi(i,j,k,5) - 0.5_rp*ir*(ru*ru+rv*rv+rw*rw))
                   c_ = sqrt(gamma0*T_)

                   max_speed1 = max(abs(ru*ir) + c_, max_speed1)
                   max_speed2 = max(abs(rv*ir) + c_, max_speed2)
                   max_speed3 = max(abs(rw*ir) + c_, max_speed3)

                enddo
             enddo
          enddo
          !$acc end parallel
        
          ! critical dt for ith direction
          dt_i1 = gbl_min_step(1)/max_speed1
          dt_i2 = gbl_min_step(2)/max_speed2
          dt_i3 = gbl_min_step(3)/max_speed3
        
          ! convective criteria
          my_dt_CFL = CFL * min(dt_i1, dt_i2, dt_i3)

          my_dt = my_dt_CFL

          if(viscous) then

            max_vis = 0.0_rp
            max_lmd = 0.0_rp
            !$acc parallel default(present) copy(max_vis,max_lmd) 
            !$acc loop gang, vector collapse(3) &
            !$acc reduction(max:max_vis,max_lmd) 
            do       k = sz,ez
               do    j = sy,ey
                  do i = sx,ex
                     max_vis = max(VIS(i,j,k),max_vis)
                     max_lmd = max(LMD(i,j,k),max_lmd)
                  enddo
               enddo
            enddo
            !$acc end parallel

            min_st2 = min(gbl_min_step(1),gbl_min_step(2),gbl_min_step(3))**2
        
            my_dt_DIF1 = 0.1_rp * min_st2 / (mu_inf*max_vis)
            my_dt_DIF2 = 0.1_rp * min_st2 / (mu_inf*max_lmd)

            my_dt = min(my_dt_CFL, my_dt_DIF1, my_dt_DIF2)

          endif

        else ! get exactly to tmax

          my_dt = (tmax - time)
          istop = 1

        endif

        ! compute GLOBAL time step
        call MPI_allreduce(my_dt, dt, 1, MPI_RP, MPI_MIN, mpi_comm_cart, mpi_err)

        call EndProfRange
#ifdef TIME
        call mpi_etime(s_cfl_time,t_cfl_calls,t_cfl_time)
#endif

        return
end subroutine compute_dt



subroutine explicit_filter_wmles
        implicit none
        
        real(rp), allocatable, dimension(:,:,:,:) :: tmp
        real(rp), dimension(5)                    :: phif
        real(rp), parameter                       :: alpha = 0.999_rp
        real(rp), dimension(-2:2, -2:2, -2:2)     :: a
        logical                                   :: wmflag
        character(3)                              :: bchar
        integer                                   :: f,i,j,k,m,n,l

        wmflag = .false.
        do f = 1,6
           bchar = trim(bc(f))
           if(bchar == 'aws') WMflag = .true.
        enddo

        if(wmflag) then

          a = 0.0_rp
          a(0,0,0) = 1.0_rp/ 2.0_rp

          a(1,0,0) = 1.0_rp/16.0_rp; a(-1,0,0) = a(1,0,0)
          a(0,1,0) = a(1,0,0); a(0,-1,0) = a(1,0,0)
          a(0,0,1) = a(1,0,0); a(0,0,-1) = a(1,0,0)

          a(2,0,0) = -1.0_rp/16.0_rp; a(-2,0,0) = a(2,0,0)
          a(0,2,0) =  a(2,0,0); a(0,-2,0) = a(2,0,0)
          a(0,0,2) =  a(2,0,0); a(0,0,-2) = a(2,0,0)

          a(1,1,0) = 1.0_rp/32.0_rp; a(-1,-1,0) = a(1,1,0)
          a(0,1,1) = a(1,1,0); a(0,-1,-1) = a(1,1,0)
          a(1,0,1) = a(1,1,0); a(-1,0,-1) = a(1,1,0)

          a(1,-1,0) = a(1,1,0); a(-1,1,0) = a(1,1,0)
          a(0,-1,1) = a(1,1,0); a(0,1,-1) = a(1,1,0)
          a(1,0,-1) = a(1,1,0); a(-1,0,1) = a(1,1,0)

          a( 1 ,1,1) = 1.0_rp/64.0_rp
          a(-1 ,1,1) = a(1,1,1)
          a( 1,-1,1) = a(1,1,1)
          a( 1,1,-1) = a(1,1,1)
          a(1,-1,-1) = a(1,1,1)
          a(-1,-1,1) = a(1,1,1)
          a(-1,1,-1) = a(1,1,1)
          a(-1,-1,-1) = a(1,1,1)

          allocate(tmp(lbx:ubx, lby:uby, lbz:ubz, 5))

          tmp = phi
          do k = sz,ez
             do j = sy,ey
                do i = sx,ex

                   phif = 0.0_rp
                   do       l = -2,2
                      do    m = -2,2
                         do n = -2,2
                            phif = phif + a(n,m,l)*tmp(i+n,j+m,k+l,:)
                         enddo
                      enddo
                   enddo

                   phi(i,j,k,:) = alpha*phi(i,j,k,:) + (1.0_rp-alpha)*phif

                enddo
             enddo
          enddo

          deallocate(tmp)

        endif

        return
end subroutine explicit_filter_wmles



subroutine check_for_NAN_values
        use, intrinsic ::ieee_arithmetic
        implicit none
        !implicit none
        !character(3)   :: check_nan
        !character(3), parameter   :: nan_val = 'NaN' 

        !do iv = 1,5
        !   write(check_nan,'(F3.3)') phi(sx,sy,sz,iv)
        !   if (check_nan.eq.nan_val) stop 'Divergency has been detected' 
        !   write(check_nan,'(F3.3)') phi(ex,ey,ez,iv)
        !   if (check_nan.eq.nan_val) stop 'Divergency has been detected' 
        !enddo

        do iv = 1,5
           if (ieee_is_nan(phi(sx,sy,sz,iv))) stop 'Divergency has been detected'
           if (ieee_is_nan(phi(ex,ey,ez,iv))) stop 'Divergency has been detected'
        enddo

        return
end subroutine check_for_NAN_values


subroutine last_iteration(sTime,istop)
        implicit none
        real(rp), intent(in)  :: sTime
        integer , intent(out) :: iStop
        
        character(dl), parameter :: marconi100 = '/m100/home/userexternal/fdevanna/uranos_mpi'
        character(dl), parameter :: sapienza   = '/scratch/mbernard/uranos_mpi'
        character(dl), parameter :: irene      = '/ccc/scratch/cont005/ra0111/devannaf/uranos_mpi'
        character(dl), parameter :: galileo100 = '/g100_work/IscrC_iWLES/uranosToGPU/fdevanna/uranos_mpi'
        character(dl)            :: wrkdir
        real(rp)                 :: ActualTime, MachinTmax

        call StartProfRange("last_iteration") 

        call getcwd(wrkdir)
        
        selectcase(trim(wrkdir))

          case(trim(marconi100), trim(irene))

            ActualTime = MPI_WTIME()
            MachinTmax = 23.00_rp*3600.0_rp

            if(ActualTime - sTime > MachinTmax) then
              if(rank == root) print*, ' Last iteration!'
              istop = 1
            endif

          case(trim(galileo100))

            ActualTime = MPI_WTIME()
            MachinTmax = 7.50_rp*3600.0_rp

            if(ActualTime - sTime > MachinTmax) then
              if(rank == root) print*, ' Last iteration!'
              istop = 1
            endif

          case(trim(sapienza))

            ActualTime = MPI_WTIME()
            MachinTmax = 71.0_rp*3600.0_rp

            if(ActualTime - sTime > MachinTmax) then
              if(rank == root) print*, ' Last iteration!'
              istop = 1
            endif

          case default
            continue

        endselect

        call EndProfRange

        return
end subroutine last_iteration




end module time_module
